library(tidyverse)

create_globals <- function () {
  
  # globals were previously defined in plasma_omics.wdl which is bad practice
  # from a security point of view. Now, we'll pull these same variables out
  # of environmental variables set in the .Renviron
  
  required_global_mapping <- tibble::tribble(
    ~ global_var, ~ env_var_name,
    "project_path", "doomics_zfs_root",
    "repo_path", "doomics_repo_path",
    "do_drive_token", "doomics_drive_working_uri",
    "shock_token", "doomics_shock_token",
    "out_base_dir", "doomics_out_base_dir"
  )
  
  global_values <- required_global_mapping %>%
    mutate(env_var = purrr::map_chr(env_var_name, Sys.getenv))
  
  stopifnot(all(global_values$env_var != ""))
  
  globals_list <- as.list(global_values$env_var)
  names(globals_list) <-  required_global_mapping$global_var
  
  return(globals_list)
}

test_config_globals <- function(do_config_json) {
  
  stopifnot("globals" %in% names(do_config_json))
  defined_globals <- do_config_json[["globals"]]
  
  required_globals <- c("project_path", "repo_path", "do_drive_token", "shock_token", "out_base_dir")
  missing_required_globals <- setdiff(required_globals, names(defined_globals))
  if (length(missing_required_globals) != 0) {
    stop ("missing required global params in json: ", paste(missing_required_globals, collapse = ", "))
  }
  
  # check path validity
  global_paths <- required_globals[required_globals %>% stringr::str_detect('path$')]
  existing_paths <- defined_globals[global_paths] %>%
    purrr::map_lgl(file.exists)
  if (sum(!existing_paths) != 0) {
    stop("Some paths are missing: ", paste(names(existing_paths)[!existing_paths], collapse = ", "))
  }
  
  # can n_slurm_jobs be coerced to integer  
  n_slurm_jobs = defined_globals[["n_slurm_jobs"]]
  stopifnot(!is.na(as.integer(n_slurm_jobs)))
}

upload_knitted_html <- function (out_dir, doc_name, connect_id = NULL) {
  
  out_htmls <- list.files(path = out_dir, pattern = ".html$")
  
  if (length(out_htmls) == 1) {
    full_html_path <- file.path(out_dir, out_htmls)
    print(glue::glue("Uploading {full_html_path} to Connect"))
    
    calibase::connect_deploy_doc(
      full_html_path,
      doc_name,
      server = "public",
      appId = connect_id,
      forceUpdate = TRUE
      )
  } else {
    cli::cli_abort("{length(out_htmls)} .html files exist in {.file {out_dir}}")
  }
}

featurization_scripts = tibble::tribble(
  ~ data_type, ~ script,
  "proteomics", "featurization_proteomics.Rmd",
  "small-molecules", "featurization_small_molecules.Rmd",
  "o-link", "featurization_olink.Rmd"
  )

# format inputs

if (interactive()) {
  # interactive
  formatted_flags <- c(
    process_feature_id = "small-molecules.lipids-neg.lipids-neg-3",
    json_path_table = "~/cromwell-executions/biom/f1b7548a-f8d7-4064-9d56-c3a48d5507e5/call-featurization/shard-0/inputs/149462423/featurization_table.tsv",
    path_do_config = file.path(Sys.getenv("doomics_repo_path"), "bioinformatics", "1_wdl_processes", "do_config.json"),
    overwrite = "true"
    )
} else {
  # command line
  input_args <- commandArgs()
  formatted_flags <- quahog::format_flagged_arguments(input_args)
}

if ("overwrite" %in% names(formatted_flags)) {
  overwrite <- as.logical(formatted_flags["overwrite"])
  if (is.na(overwrite)) {
    stop ("overwrite was mis-specified:", formatted_flags["overwrite"], "it must be false or true if provided")
    }
} else {
  overwrite <- FALSE
}

print(formatted_flags)

# test required args

required_args = c("process_feature_id", "json_path_table", "path_do_config")
missing_required_args <- setdiff(required_args, names(formatted_flags))
if (length(missing_required_args) != 0) {
  stop ("missing required arguments: ", paste(missing_required_args, collapse = ", "))
}

# read json config and test globals

do_config_json <- jsonlite::fromJSON(formatted_flags["path_do_config"])
do_config_json$globals <- create_globals()

test_config_globals(do_config_json)

# What run mode are we in?
json_path_table_basename <- basename(formatted_flags["json_path_table"])
script_mode <- switch(json_path_table_basename, 
                      "featurization_table.tsv" = "featurization",
                      "integration_table.tsv" = "integration",
                      "analysis_table.tsv" = "analysis")

if (is.null(script_mode)) {
  stop ("Inappropriate json_path_table provided; ", json_path_table_basename, "; valid types are
        \"featurization_table.tsv\", \"integration_table.tsv\", \"analysis_table.tsv\"")
}

# select DO json path
json_path_table_tsv <- readr::read_tsv(formatted_flags["json_path_table"], show_col_types = FALSE)
if (!(formatted_flags["process_feature_id"] %in% json_path_table_tsv$run_id)) {
  stop ("No information for run_id: ", formatted_flags["process_feature_id"], " in featurization_table; valid IDs are: ", paste(json_path_table_tsv$run_id, collapse = ", ")) 
}
selected_json_path <- json_path_table_tsv %>%
  dplyr::filter(run_id == formatted_flags["process_feature_id"]) %>%
  unlist()

if (script_mode == "featurization") {
  stopifnot(colnames(json_path_table_tsv) == c("run_id", "featurization", "datatype", "runspecs"))
  run_id = unname(selected_json_path["run_id"])
  data_type = unname(selected_json_path["featurization"])
  call_script <- featurization_scripts$script[featurization_scripts$data_type == data_type]
  call_script_path <- file.path(
    do_config_json[["globals"]][["repo_path"]],
    "bioinformatics",
    "1_wdl_processes",
    call_script)
  
  # read run specs
  run_specs <- try(do_config_json[["featurization"]][[data_type]][[selected_json_path["datatype"]]][[selected_json_path["runspecs"]]], silent = TRUE)
  if (class(run_specs) %in% c("NULL", "try-error")) {
    stop ("json path specs for run_id are malformed: ", paste(selected_json_path, collapse = ", "))
  }
  
  # setup output and test for existence and whether to eval or overwrite
  run_outdir = file.path(do_config_json[["globals"]][["project_path"]], do_config_json[["globals"]][["out_base_dir"]], run_id)
  
  existing_analysis <- file.exists(run_outdir) && overwrite == FALSE
  ignore_analysis <- "eval" %in% names(run_specs) && run_specs$eval == "FALSE"
  skip_analysis <- (existing_analysis | ignore_analysis)
  
  # run update
  if (skip_analysis) {
    message("Not run due to non-overwrite or non-eval")
  } else {
    
    if (file.exists(run_outdir)) {
      unlink(run_outdir, recursive = TRUE) 
    }
    dir.create(run_outdir)
    
    # set up parameters to pass into Rmd as a YAML header
    run_params <- append(run_specs, do_config_json[["globals"]])
    run_params <- run_params[!(names(run_params) %in% c("eval", "out_base_dir"))]
    run_params$run_outdir = run_outdir
    
    rmd_paramers <- run_params[!(names(run_params) %in% c("connect_id"))]
    rmarkdown::render(input = call_script_path,
                      output_dir = run_outdir,
                      params = rmd_paramers,
                      envir = new.env())
    
    # upload rendered .html
    doc_name = glue::glue("DO Mice : {run_id}")
    if ("connect_id" %in% names(run_specs)) {
      # provide a connect ID via the config to update an existing connect report
      connect_id <- as.integer(unname(run_specs["connect_id"]))
    } else {
      connect_id <- NULL
    }
    
    print("Uploading knitted document to Connect")
    upload_knitted_html(run_outdir, doc_name, connect_id = connect_id) 
  }}

# write mtcars if successful :)
write.table(mtcars, file = "outputA.ext")
