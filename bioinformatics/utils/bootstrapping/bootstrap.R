if (interactive()) {
  job_id_value <- 1
  nbs <- 1000
} else {
  # read in command line arguements
  args <- commandArgs()
  job_id_value = as.numeric(unlist(strsplit(args[grep("job_id", args)], "="))[2])
  nbs = as.numeric(unlist(strsplit(args[grep("nbs", args)], "="))[2])
}

print(paste("JOB ID: ", job_id_value, ", # OF BOOTSTRAPS: ", nbs, sep = ""))
print(sessionInfo())

# env setup
library(tidyverse)
future::plan("sequential")

repo_path <- Sys.getenv("doomics_repo_path")
source(file.path(repo_path, "bioinformatics", "utils", "linear_model_fxns.R"))

# path to data containing feature-level data and models-to-fit
diffex_root_path <- file.path(
  Sys.getenv("doomics_zfs_root"),
  "WDL_20220919",
  "differential_abundances",
  "out_files"
)
lm_setup_data <- file.path(diffex_root_path, "lm_setup_data.Rds")

out_file <- glue::glue("fit_lms_job{job_id_value}_bs{nbs}.Rds")
output_path <- file.path(diffex_root_path, "bootstrap_results", out_file)

overwrite = TRUE
if (!overwrite && file.exists(output_path)) {
  stop (glue::glue("Output exists and overwrite is FALSE: {output_path}"))
}

lm_setup_data_job <- readRDS(lm_setup_data) %>%
  dplyr::filter(job_id == job_id_value)

fit_lms <- lm_setup_data_job %>%
  dplyr::mutate(lms = furrr::future_map2(
    model_data,
    model_params,
    lm_fitter,
    n_bootstraps = nbs,
    .progress = FALSE,
    .options = furrr::furrr_options(seed = TRUE)
    )) %>%
  dplyr::select(-model_data, -model_params) %>%
  tidyr::unnest(lms)

saveRDS(fit_lms, output_path)

