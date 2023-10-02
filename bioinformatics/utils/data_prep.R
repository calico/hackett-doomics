load_sample_labels <- function(path_to_do_github) {
  
  supporting_files_path <- file.path(path_to_do_github, "supporting_files")
  
  # bin by time of experiment (date of blood draw)
  goal_age <- tibble::tibble(
    blood_draw_age_order = 1:3,
    blood_draw_age_coarse = c("8 months", "14 months", "20 months")
    ) %>%
    dplyr::mutate(blood_draw_age_coarse = ordered(blood_draw_age_coarse, levels = blood_draw_age_coarse))
  
  # set and sample uniquely identify a sample. Batches are taken directly from the data.
  
  set12_design <- suppressMessages(readr::read_csv(
    file.path(supporting_files_path, "design_matrix_set1_and_set2.csv")
    )) %>%
    dplyr::mutate(sample = paste0("s", `sample ID`)) %>%
    dplyr::select(set, sample, Mouse.ID, blood_draw_age_actual = Age.at.Exp) %>%
    dplyr::group_by(Mouse.ID) %>%
    dplyr::arrange(blood_draw_age_actual) %>%
    dplyr::mutate(blood_draw_age_order = 1:n()) %>%
    dplyr::left_join(goal_age, by = "blood_draw_age_order") %>%
    dplyr::select(-blood_draw_age_order) %>%
    dplyr::ungroup()
  
  # use with_edition so that that newlines are recognized
  # https://community.rstudio.com/t/readr-2-0-0-now-vroom-doesnt-seem-to-recognise-files-with-lines-delimited-by-cr-instead-of-cr-lf/113462
  set3_design <- suppressMessages(readr::with_edition(1, readr::read_tsv(
    file.path(supporting_files_path, "DO set3 anonymous IDs and batches.txt")
    ))) %>%
    dplyr::mutate(
      set = 3,
      batch = reassigned_batch,
      sample = paste0("s", anon.num),
      blood_draw_age_actual = lubridate::time_length(
        lubridate::mdy(Harvest.Date) - lubridate::mdy(DOB),
        unit = "days"
        )) %>%
    dplyr::group_by(Mouse.ID.x) %>%
    dplyr::arrange(blood_draw_age_actual) %>%
    dplyr::mutate(blood_draw_age_order = 1:n()) %>%
    dplyr::left_join(goal_age, by = "blood_draw_age_order") %>%
    dplyr::select(set, sample, Mouse.ID = Mouse.ID.x, blood_draw_age_coarse, blood_draw_age_actual) %>%
    dplyr::ungroup()
  
  all_designs <- dplyr::bind_rows(set12_design,
                                  set3_design) %>%
    dplyr:::rename(Age = blood_draw_age_coarse)
  
  sample_count <- all_designs %>%
    dplyr::count(Mouse.ID) %>%
    dplyr::count(n, name = "nn")
  
  if (nrow(sample_count) != 1) {
    stop("mice have different numbers of samples") 
  }
  
  #ggplot(all_designs, aes(x = Age, y = blood_draw_age_actual)) + geom_boxplot()
  
  return(all_designs)
}

load_mouse_phenotypes <- function(shock_token, reduced_phenotypes = TRUE) {
  
  # Setup additional phenotypes
  
  local_shock_file = "/tmp/longitudinal.Rdata"
  
  # download file
  
  if (!file.exists(local_shock_file)) {
    googledrive::as_id(shock_token) %>%
      googledrive::drive_ls() %>%
      dplyr::filter(name == "Data") %>%
      googledrive::drive_ls() %>%
      dplyr::filter(name == "longitudinal.Rdata") %>%
      googledrive::drive_download(path = local_shock_file)
  }
  
  # load shock R tables
  
  load(local_shock_file)
  
  mouse_level_phenotypes <- if (reduced_phenotypes) {
    data.life %>%
      dplyr::select(Mouse.ID, Generation, Sex, Birth_date, lifespan)
  } else {
    
    death_categories <- data.life %>%
      dplyr::count(Death.type) %>%
      dplyr::arrange(desc(n)) %>%
      # found dead
      dplyr::mutate(Death.type.coarse = ifelse(grepl('^F\\.?D\\.?', Death.type), "F.D.", Death.type),
                    # euthenized
                    Death.type.coarse = ifelse(grepl('^E\\.?S\\.?', Death.type.coarse), "E.S.", Death.type.coarse),
                    # flag other mice since most are missing or NAs
                    Death.type.coarse = ifelse(Death.type.coarse %in% c("F.D.", "E.S."), Death.type.coarse, NA)) %>%
      dplyr::mutate(death_FD = ifelse(!is.na(Death.type.coarse) & Death.type.coarse == "F.D.", 1, 0),
                    death_ES = ifelse(!is.na(Death.type.coarse) & Death.type.coarse == "E.S.", 1, 0)) %>%
      dplyr::select(-n)
    
    data.life %>%
      dplyr::left_join(death_categories, by = "Death.type") %>%
      # only keep lifespan for F.D. and E.S. mice
      dplyr::mutate(lifespan = ifelse(is.na(Death.type.coarse), NA, lifespan),
                    Phenotype.Category = "mouse") %>%
      dplyr::select(Mouse.ID, Sex, Generation, Phenotype.Category, death_FD, death_ES, lifespan, Tail_Length, Birth_date) %>%
      # lifepsans specific to a mode of death
      dplyr::mutate(lifespan_FD = ifelse(death_FD == 1, lifespan, NA),
                    lifespan_ES = ifelse(death_ES == 1, lifespan, NA)) %>%
      tidyr::gather(Phenotype, Value, -Mouse.ID, -Sex, -Generation, -Phenotype.Category, -Birth_date)
  }
  
  # time specific covariates
  
  age_based_phenotypes <- if (reduced_phenotypes) {
    
    data.bw.long %>%
      dplyr::select(Mouse.ID, Age, bw, gluc) %>%
      dplyr::mutate(blood_draw_age_coarse = ordered(paste0(as.character(as.numeric(as.character(Age)) + 2), " months"),
                                                    c("8 months", "14 months", "20 months"))) %>%
      dplyr::select(-Age)
    
  } else {
    
    dplyr::bind_rows(
      gather_phenotypes(data.bw.long, "bw"),
      gather_phenotypes(data.advia.long, "advia"),
      gather_phenotypes(data.chem.long, "chem"),
      gather_phenotypes(data.facs.long, "facs"),
      gather_phenotypes(data.urine.long, "urine")) %>%
      dplyr::mutate(blood_draw_age_coarse = ordered(paste0(as.character(as.numeric(as.character(Age)) + 2), " months"),
                                                    c("8 months", "14 months", "20 months"))) %>%
      dplyr::select(-Age)
    
  }
  
  output <- list()
  output$mouse_level_phenotypes <- mouse_level_phenotypes
  output$age_based_phenotypes <- age_based_phenotypes
  output
}

gather_phenotypes <- function(phenotype.category.long, label) {
  phenotype.category.long %>%
    dplyr::mutate(Phenotype.Category = label) %>%
    tidyr::gather(Phenotype, Value, -Phenotype.Category, -Mouse.ID, -Sex, -Generation, -Age)
}

load_wide_age_based_phenotypes <- function (shock_token) {
  
  full_mouse_phenotypes <- load_mouse_phenotypes(shock_token = shock_token,
                                                 reduced_phenotypes = FALSE)
  
  wide_mouse_level <- full_mouse_phenotypes$mouse_level_phenotypes %>%
    tidyr::spread(Phenotype, Value) %>%
    dplyr::select(Mouse.ID, death_ES, death_FD)
  
  wide_age_based <- full_mouse_phenotypes$age_based_phenotypes %>%
    dplyr::filter(Phenotype.Category %in% c("bw", "chem", "urine"),
                  Phenotype != "mg") %>%
    dplyr::select(-Phenotype.Category, -Sex) %>%
    tidyr::spread(Phenotype, Value)
  
  return (
    wide_age_based %>%
      dplyr::left_join(wide_mouse_level, by = "Mouse.ID")
  )
}

read_design_matrix <- function(design_matrix_path) {
  
  readr::read_csv(design_matrix_path) %>%
    dplyr::mutate(sample_id = glue::glue('{set}_{`sample ID`}')) %>%
    dplyr::select(set, batch, `batch order`, sample_id) %>%
    dplyr::arrange(set, batch) %>%
    dplyr::mutate(unique_batch = paste0("set", set, "_b", batch),
                  unique_batch = factor(unique_batch, levels = unique(unique_batch))) 
  
}

read_design_matrix_protein_quant <- function(design_matrix_path, set = 'DOBATCH1') {
  
  if(set == 'DOBATCH1'){
    design_matrix <- readr::read_csv(design_matrix_path) %>%
      dplyr::mutate(sample_id = glue::glue('{set}_{`sample ID`}')) %>%
      dplyr::select(set, batch, `batch order`, sample_id) %>%
      dplyr::arrange(set, batch) %>%
      dplyr::mutate(unique_batch = paste0("set", set, "_b", batch),
                    unique_batch = factor(unique_batch, levels = unique(unique_batch))) 
  }else if(set == 'DOBATCH2'){
    # TODO (MM): clean this up, a bit messy because I was confused about sample_id and anon.num 
    design_matrix <- readr::read_tsv(design_matrix_path) %>%
      select(Mouse.ID.x, reassigned_batch, anon.num) %>%
      rename(c('sample_id' = 'anon.num', 'batch' = 'reassigned_batch')) %>%
      mutate(set = 3,
             `batch order` = sample_id,
             sample_id = as.character(sample_id),
             unique_batch = paste0('set', set, '_b', batch),
             unique_batch = factor(unique_batch, levels = unique(unique_batch))) %>%
      select(set, batch, `batch order`, sample_id, unique_batch) %>%
      dplyr::arrange(set, batch, `batch order`)
    
    # We have more measurements than we do samples for some reason.
    # looking at the design matrix in masspike conditions 90 and 91 are missing
    # which correspond to plex 17, channels TMT9 and TMT10.
    # measurements from these channels are labeled as dummy variables
    # to be removed later.
    design_matrix <- bind_rows(design_matrix,
                               list('set' = 3, 'batch' = 17, 'batch order' = 9, 'sample_id' = 'dummy_1', 'unique_batch' = 'set3_b17'),
                               list('set' = 3, 'batch' = 17, 'batch order' = 10, 'sample_id' = 'dummy_2', 'unique_batch' = 'set3_b17')) %>%
      mutate(`batch order` = rep(1:10, times = length(unique(batch)))) %>%
      mutate(sample_id = stringr::str_c(set, '_', sample_id))
    
  }
  
  design_matrix
}


read_protein_quant <- function(path, design_matrix, ssn_threshold = 200, sn_sample_cutoff = 10, iso_spec_threshold = 0.7) {
  # correspondence between tags and common labels
  
  TMT_tags <- tibble::tribble(~ tag, ~ label,
                              "rq_126_sn", "TMT1",
                              "rq_127n_sn", "TMT2",
                              "rq_127c_sn", "TMT3",
                              "rq_128n_sn", "TMT4",
                              "rq_128c_sn", "TMT5",
                              "rq_129n_sn", "TMT6",
                              "rq_129c_sn", "TMT7",
                              "rq_130n_sn", "TMT8",
                              "rq_130c_sn", "TMT9",
                              "rq_131_sn", "TMT10",
                              # Batch 2 (Set 3) was measured using 11 channels
                              "rq_131n_sn", "TMT10", 
                              "rq_131c_sn", "TMT11")
  
  
  # correspondence between TMT samples and set / batch
  
  class_lookup_table <- tibble::tribble(~ Class, ~ set, ~ batch,
                                        "DO1TMT5", 1, 1,
                                        "DO1TMT6", 1, 2,
                                        "DO1TMT7", 1, 3,
                                        "DO1TMT8", 1, 4,
                                        "DO1TMT9", 1, 5,
                                        "DO1TMT10", 1, 6,
                                        "DO1TMT11", 2, 1,
                                        "DO1TMT12", 2, 2,
                                        "DO1TMT13", 2, 3,
                                        "DO1TMT14", 2, 4,
                                        "DO1TMT15", 2, 5,
                                        "DO1TMT16", 2, 6,
                                        "DO1TMT17", 2, 7,
                                        "DO1TMT18", 2, 8,
                                        "DO1TMT19", 2, 9,
                                        "DO1TMT20", 2, 10,
                                        "DO1TMT21", 2, 11,
                                        "DO1TMT22", 2, 12,
                                        "DO2TMT1", 3, 1,
                                        "DO2TMT2", 3, 2,
                                        "DO2TMT3", 3, 3,
                                        "DO2TMT4", 3, 4,
                                        "DO2TMT5", 3, 5,
                                        "DO2TMT6", 3, 6,
                                        "DO2TMT7", 3, 7,
                                        "DO2TMT8", 3, 8,
                                        "DO2TMT9", 3, 9,
                                        "DO2TMT10", 3, 10,
                                        "DO2TMT11", 3, 11,
                                        "DO2TMT12", 3, 12,
                                        "DO2TMT13", 3, 13,
                                        "DO2TMT14", 3, 14,
                                        "DO2TMT15", 3, 15,
                                        "DO2TMT16", 3, 16,
                                        "DO2TMT17", 3, 17)
  
  # match TMT channels to correct batch order (it varies by set)
  
  tag_order_index <- tibble::tribble(~ set, ~ label, ~ `batch order`,
                                     1, "TMT1", "1",
                                     1, "TMT2", "2",
                                     1, "TMT3", "3",
                                     1, "TMT4", "4",
                                     1, "TMT5", "5",
                                     1, "TMT6", "6",
                                     1, "TMT7", "7",
                                     1, "TMT8", "8",
                                     1, "TMT9", "9",
                                     1, "TMT10", "bridge",
                                     2, "TMT1", "bridge",
                                     2, "TMT2", "1",
                                     2, "TMT3", "2",
                                     2, "TMT4", "3",
                                     2, "TMT5", "4",
                                     2, "TMT6", "5",
                                     2, "TMT7", "6",
                                     2, "TMT8", "7",
                                     2, "TMT9", "8",
                                     2, "TMT10", "9",
                                     3, "TMT1", "1",
                                     3, "TMT2", "2",
                                     3, "TMT3", "3",
                                     3, "TMT4", "4",
                                     3, "TMT5", "5",
                                     3, "TMT6", "6",
                                     3, "TMT7", "7",
                                     3, "TMT8", "8",
                                     3, "TMT9", "9",
                                     3, "TMT10", "10",
                                     3, "TMT11", "bridge")
  
  mp_dataset_raw <- readr::read_tsv(path)
  
  # filter TMT tags only found in other batch i.e in Set1/2 or Set3
  TMT_tags <- TMT_tags %>% filter(tag %in% colnames(mp_dataset_raw))
  TMT_tags_vec <- TMT_tags$tag; names(TMT_tags_vec) <- TMT_tags$label
  
  mp_dataset_raw <- mp_dataset_raw %>%
    # rename TMT tags
    dplyr::rename(!!TMT_tags_vec) %>%
    # remove reverse hits and low isolation specificity
    dplyr::filter(!stringr::str_detect(ProteinId, '^##'),
                  !stringr::str_detect(ProteinId, 'contaminant')) %>%
    dplyr::filter(IsolationSpecificity > iso_spec_threshold) %>%
    dplyr::mutate(.row_number = 1:n())
  
  # lookup table for ProteinId + GeneSymbol -> GeneName
  protein_identifiers <- mp_dataset_raw %>%
    # deal with missing gene symbols and inconsistent description
    # by defining a unique GeneName
    dplyr::distinct(ProteinId, GeneSymbol, Description) %>%
    dplyr::mutate(GeneName = dplyr::case_when(!is.na(GeneSymbol) ~ GeneSymbol,
                                              TRUE ~ ProteinId))
  
  unique_proteins <- protein_identifiers %>%
    # define proteins with a unique gene-name and choose one description
    dplyr::group_by(GeneName) %>%
    dplyr::summarize(Description = Description[1]) %>%
    dplyr::ungroup()
  
  unique_peptides <- mp_dataset_raw %>%
    dplyr::left_join(protein_identifiers %>% dplyr::select(-Description), by = c("ProteinId", "GeneSymbol")) %>%
    dplyr::distinct(GeneName, PeptideSequence) %>%
    dplyr::group_by(GeneName) %>%
    dplyr::mutate(pepindx = 1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Peptide_ID = glue::glue('{GeneName}_{pepindx}'))
  
  summed_SN <-   mp_dataset_raw %>%
    dplyr::select(!!TMT_tags$label, .row_number) %>%
    tidyr::gather(key = "label", value = "SN", -.row_number) %>%
    dplyr::group_by(.row_number) %>%
    dplyr::summarize(SSN = sum(SN)) %>%
    dplyr::filter(SSN >= ssn_threshold)
  
  unique_measurements <- mp_dataset_raw %>%
    # add row-level summed signal to noise
    dplyr::inner_join(summed_SN, by = ".row_number") %>%
    # add unique Peptide_ID
    dplyr::left_join(protein_identifiers %>% dplyr::select(-Description), by = c("ProteinId", "GeneSymbol")) %>%
    dplyr::left_join(unique_peptides %>% dplyr::select(GeneName, PeptideSequence, Peptide_ID), by = c("GeneName", "PeptideSequence")) %>%
    dplyr::group_by(Class, Peptide_ID) %>%
    dplyr::arrange(desc(SSN)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  tall_measurements <- unique_measurements %>%
    dplyr::select(Class, GeneName, Peptide_ID, any_of(paste0('TMT', 1:11))) %>%
    tidyr::gather("label", "SN", -Class, -Peptide_ID, -GeneName) %>%
    # convert from class to set + batch
    dplyr::left_join(class_lookup_table, by = "Class") %>%
    # add batch order / bridge
    dplyr::left_join(tag_order_index, by = c("set", "label")) %>%
    # add sample id (anonymized) based on set, batch and batch order    
    dplyr::left_join(design_matrix %>%
                       dplyr::mutate(`batch order` = as.character(`batch order`)),
                     by = c("set", "batch", "batch order")) %>%
    # missing samples are bridge
    dplyr::mutate(sample_type = dplyr::case_when(is.na(sample_id) ~ "bridge",
                                                 TRUE ~ "experimental"),
                  sample_id = dplyr::case_when(is.na(sample_id) ~ glue::glue("{set}_bridge"),
                                               TRUE ~ sample_id)) %>%
    dplyr::select(GeneName, Peptide_ID, set, batch, unique_batch, sample_id, sample_type, SN) %>%
    # remove low SN observations or references to low SN bridges
    dplyr::filter(SN > sn_sample_cutoff) %>%
    dplyr::group_by(Peptide_ID, set, batch) %>%
    dplyr::filter("bridge" %in% sample_type) %>%
    dplyr::ungroup() %>%
    # compare each observation to bridge
    dplyr::group_by(Peptide_ID, set, batch) %>%
    dplyr::mutate(log2.SN = log2(SN),
                  log2.SN.bridge = log2.SN[sample_type == "bridge"],
                  log2.RSN = log2.SN - log2.SN.bridge) %>%
    dplyr::ungroup() %>%
    dplyr::select(-log2.SN.bridge)
  
  reduced_peptides <- unique_peptides %>%
    dplyr::semi_join(tall_measurements, by = "Peptide_ID")
  
  reduced_proteins <- unique_proteins %>%
    dplyr::semi_join(reduced_peptides, by = "GeneName")
  
  # remove observations from 'empty' channels (channels which don't have a corresponding sample)
  tall_measurements <- tall_measurements %>%
    filter(!stringr::str_detect(sample_id, 'dummy'))
  
  list(unique_proteins = reduced_proteins,
       unique_peptides = reduced_peptides,
       tall_measurements = tall_measurements)
}

read_small_molecules <- function(mzroll_db_path, full_inclusion_file_path) {
  
  stopifnot(file.exists(mzroll_db_path))
  stopifnot(class(full_inclusion_file_path) == "character")
  if (!is.na(full_inclusion_file_path)) {
    stopifnot(file.exists(full_inclusion_file_path))
    
    inclusion_tbl <- readr::read_tsv(full_inclusion_file_path, show_col_types = FALSE) %>%
      dplyr::select(dplyr::one_of("filename", "keep?", "date", "time", "rebatched")) %>%
      dplyr::rename(keep = `keep?`) %>%
      dplyr::mutate(keep = ifelse(keep %in% c("yes", "y", "Y"), TRUE, FALSE))
    
    if (all(c("time", "date") %in% colnames(inclusion_tbl))) {
      inclusion_tbl <- inclusion_tbl %>%
        tidyr::unite(date, time, col = "datetime", sep = " ") %>%
        dplyr::mutate(datetime = lubridate::mdy_hms(datetime))
    }
  }
  
  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  peakgroups <- dplyr::tbl(mzroll_db_con, "peakgroups") %>%
    dplyr::collect() %>%
    dplyr::mutate(label = "")
  DBI::dbWriteTable(mzroll_db_con, "peakgroups", peakgroups, overwrite = TRUE)
  DBI::dbDisconnect(mzroll_db_con)
  
  small_molecules_mzroll <- claman::process_mzroll(
    mzroll_db_path,
    only_identified = FALSE,
    validate = TRUE
    )
  
  mode_label <- ifelse("[M-H]-" %in% small_molecules_mzroll$features$adductName, "(-) ", "(+) ")
  
  small_molecules_mzroll$features <- small_molecules_mzroll$features %>%
    # clean-up compound names
    dplyr::mutate(compoundName = paste0(mode_label, stringr::str_trim(compoundName)))
  
  # update samples
  
  parsed_samples <- small_molecules_mzroll$samples %>%
    dplyr::mutate(.set = stringr::str_extract(name, 'set[0-9]')) %>%
    # nest and iterate over sets to allow for different naming conventions in different sets
    tidyr::nest(setData = -.set) %>%
    dplyr::mutate(parsedSetData = purrr::map2(setData, .set, parse_sm_do_filestring)) %>%
    # reform parsed data across sets
    dplyr::select(-setData) %>%
    tidyr::unnest(parsedSetData) %>%
    dplyr::select(-.set) %>%
    # add sample_type
    dplyr::mutate(
      sample_type = dplyr::case_when(
        grepl('^s[0-9]+', sample) ~ "sample",
        grepl('blank', sample) ~ "blank",
        grepl('^TrueBlank', sample) ~ "blank",
        grepl('^pos', sample) ~ "posctl",
        sample == 'negctl' ~ "negctl",
        grepl('^std', sample) ~ "standards",
        grepl("LipidSplash", sample) ~ "standards",
        grepl('^pool', note) ~ "targeted_runs"
      ),
      sample_type = factor(sample_type, levels = c("blank", "negctl", "targeted_runs", "standards", "posctl", "sample")),
      set = as.numeric(sub('^set', '', set)),
      batch = as.numeric(gsub('(^b)|([A-Z]$)', '', batch)),
      injection = as.numeric(sub('^inj', '', injection)),
      sample_id = dplyr::case_when(
        sample_type == "sample" ~ glue::glue("{set}_{stringr::str_replace(sample, '^s', '')}"),
        sample_type == "posctl" ~ glue::glue("{set}_p{1:n()}"),
        TRUE ~ sample
        )
      ) %>%
    dplyr::arrange(set, batch) %>%
    dplyr::mutate(unique_batch = paste0("set", set, "_b", batch),
                  unique_batch = factor(unique_batch, levels = unique(unique_batch)))
  
  # remove samples which should be ignored
  if (!is.na(full_inclusion_file_path)) {
    
    parsed_samples <- parsed_samples %>%
      # reduce to included filenames
      dplyr::left_join(inclusion_tbl, by = c("name" = "filename")) %>%
      dplyr::filter(keep)
    
    if ("rebatched" %in% colnames(parsed_samples)) {
      
      parsed_samples <- parsed_samples %>%
        dplyr::arrange(set, rebatched) %>%
        dplyr::mutate(unique_batch = paste0("set", set, "_b", rebatched),
                      unique_batch = factor(unique_batch, levels = unique(unique_batch)))
      
    }
    
    # remove dropped peaks and peakgroups
    small_molecules_mzroll <- romic::filter_tomic(
      small_molecules_mzroll,
      filter_type = "quo",
      filter_table = "samples",
      filter_value = rlang::quo(sampleId %in% parsed_samples$sampleId)
      )
    
  } else {
    message ("No \"full_inclusion_file_path\" provided, all samples will be retained")
  }
  
  undefined_sample_types <- parsed_samples %>%
    dplyr::filter(is.na(sample_type))
  if (nrow(undefined_sample_types) != 0) {
    stop (nrow(undefined_sample_types), " samples with undefined \"sample_type\"")
  }
  
  small_molecules_mzroll <- romic::update_tomic(small_molecules_mzroll, parsed_samples)

  small_molecules_mzroll
}

parse_sm_do_filestring <- function(samples, set) {
  
  if (set == "set3") {
    samples %>%
      tidyr::separate(col = name, into = c("experiment", "method", "set", "batch", "sample", "injection", "column", "note"), sep = '[._-]', remove = FALSE)
  } else if (set %in% c("set1", "set2")) {
    samples %>%
      dplyr::filter(!grepl('pool', filename)) %>%
      tidyr::separate(col = name, into = c("experiment", "set", "batch", "sample", "injection", "column", "mode", "note"), sep = '[._-]', remove = FALSE) %>%
      dplyr::mutate(note = ifelse(note == "mzXML", NA, note)) %>%
      dplyr::bind_rows(
        samples %>%
          dplyr::filter(grepl('pool', filename)) %>%
          tidyr::separate(col = filename, into = c("experiment", "set", "note", "column", "mode", "sample"), sep = '[._-]', remove = FALSE))
  } else {
    stop ("Undefined set")
  }
}

read_small_molecule_annotations <- function (annotation_file, annotation_data_type) {
  
  annotations <- readr::read_tsv(annotation_file) %>%
    dplyr::filter(data_type == annotation_data_type)
  
  if (nrow(annotations) == 0) {
    stop ("No annotations available for ", annotation_data_type)
  }
  
  annotations <- annotations %>%
    dplyr::group_by(groupId) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(compoundName = dplyr::case_when(!is.na(systematicName) ~ systematicName,
                                                  TRUE ~ compoundName),
                  mode_label = dplyr::case_when(annotation_data_type == "metabolites-neg" ~ "(-) ",
                                                TRUE ~ "(+) "),
                  compoundName = paste0(mode_label, compoundName)) %>%
    dplyr::select(-mode_label)
  
  return (annotations)
}

read_olink <- function (full_dataset_path) {
  
  olink_header <- readxl::read_excel(full_dataset_path, col_names = FALSE, n_max = 6, skip = 2) %>%
    tidyr::gather(feature_index, value, -`...1`) %>%
    tidyr::spread(`...1`, value) %>%
    dplyr::mutate(feature_index_numeric = as.numeric(stringr::str_replace(feature_index, '^...', '')))
  
  olink_spacer_rows <- readxl::read_excel(full_dataset_path, col_names = FALSE) %>%
    dplyr::select(`...1`) %>%
    dplyr::mutate(NA_pos = is.na(`...1`)) %>%
    {which(.$NA_pos)}
  
  olink_footer <- readxl::read_excel(full_dataset_path, col_names = FALSE, skip = last(olink_spacer_rows)) %>%
    tidyr::gather(feature_index, value, -`...1`) %>%
    tidyr::spread(`...1`, value)
  
  olink_feature_attributes <- olink_header %>%
    dplyr::left_join(olink_footer, by = "feature_index", multiple = "all") %>%
    dplyr::arrange(feature_index_numeric) %>%
    dplyr::select(-feature_index) 
  
  olink_data <- readxl::read_excel(full_dataset_path, skip = olink_spacer_rows[2], n_max = olink_spacer_rows[3] - olink_spacer_rows[2], 
                                   col_names = c("sample", olink_feature_attributes$Assay)) %>%
    dplyr::select(-`Plate ID`, -`QC Warning`) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(replicate = dplyr::case_when(n() == 1 ~ NA_integer_,
                                               TRUE ~ 1:n())) %>%
    tidyr::gather(Assay, abundance, -sample, -replicate) %>%
    dplyr::ungroup() %>%
    dplyr::select(Assay, sample, normalized_log2_abundance = abundance) %>%
    dplyr::filter(!stringr::str_detect(sample, "CONTROL")) %>%
    dplyr::mutate(sampleId = sample,
                  sample = stringr::str_replace(sample, '^', 's'),
                  unique_batch = NA_character_,
                  dataset = 'olink',
                  data_type = 'olink',
                  set = 3)
  
  olink_data$groupId <- olink_data %>% 
    dplyr::group_by(Assay) %>% 
    dplyr::group_indices()
  
  olink_feature_attributes <- olink_feature_attributes %>%
    dplyr::filter(!(Assay %in% c('Plate ID', 'QC Warning'))) %>%
    dplyr::left_join(olink_data %>% dplyr::distinct(Assay, groupId), by = 'Assay') %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      data_type = "olink",
      NameUnique = Assay
    )
  
  olink_data <- olink_data %>% 
    dplyr::select(-Assay)
  
  list(feature_attributes = olink_feature_attributes,
       tall_measurements = olink_data)
}

featurization_parameterization <- function(path_do_config, path_out_base_dir) {
  
  do_config_json <- jsonlite::fromJSON(path_do_config)
  
  featurization_parameters <- tibble::tibble(data_type = names(do_config_json[["featurization"]]),
                                             data_val = unname(do_config_json[["featurization"]])) %>%
    dplyr::mutate(dataset = purrr::map(data_val, tibble_list)) %>%
    tidyr::unnest_legacy(dataset) %>%
    dplyr::mutate(dataset = purrr::map(values, tibble_list)) %>%
    tidyr::unnest_legacy(dataset) %>%
    dplyr::mutate(dataset = purrr::map(values, tibble_list)) %>%
    tidyr::unnest_legacy(dataset) %>%
    dplyr::mutate(value = purrr::map_chr(values, identity)) %>%
    dplyr::select(data_type, method = names...2, dataset = names...3, parameter = names...4, value) %>%
    tidyr::spread(parameter, value) %>%
    tidyr::unite(data_type, method, dataset, sep = ".", col = dataset_dir, remove = FALSE) %>%
    dplyr::mutate(dataset_path = file.path(path_out_base_dir, dataset_dir)) %>%
    tidyr::separate(dataset_dir, into = c("featurization", "data_type", "data_partition"), sep = "\\.")
  
  return (featurization_parameters)
}

tibble_list <- function (list) {
  tibble::tibble(names = names(list),
                 values = unname(list))
}

small_molecules_tall_do_abundances <- function(mzroll_list) {
  
  mzroll_list$measurements %>%
    dplyr::select(groupId, sampleId, normalized_log2_abundance) %>%
    dplyr::mutate(
      groupId = as.integer(as.character(groupId)),
      sampleId = as.character(sampleId)
      ) %>%
    dplyr::inner_join(
      mzroll_list$samples %>%
        dplyr::filter(sample_type == "sample") %>%
        dplyr::select(sampleId, set, sample, unique_batch) %>%
        dplyr::mutate(sampleId = as.character(sampleId)),
      by = "sampleId"
      ) %>%
    dplyr::left_join(
      mzroll_list$features %>%
        dplyr::select(groupId, peak_label, adductName) %>%
        dplyr::mutate(
          peak_label = stringr::str_trim(peak_label),
          groupId = as.integer(as.character(groupId))),
      by = "groupId"
      )
}

construct_features_with_design <- function (params, met_input = "valid") {
  
  # populate feature sets to look at
  path_do_config <- file.path(params$repo_path, "do_config.json")
  path_out_base_dir <- file.path(params$project_path, params$out_base_dir)
  featurization_parameters <- suppressMessages(featurization_parameterization(path_do_config, path_out_base_dir))
  
  # setup small molecule data
  small_molecules <- load_small_molecules(featurization_parameters, met_input)
  
  # load CompMS protein relative abundance measurements
  proteomics_dir <- file.path(params$project_path, 'data/proteomics/out_files')
  proteomics <- load_proteomics(proteomics_dir)
  
  # load peptide level data from Masspike Protein Quant
  peptidomics_dir <- file.path(params$project_path, 'data/peptidomics/out_files')
  peptidomics <- load_peptidomics(peptidomics_dir)
  
  # read Olink file from their excel speadsheet
  # TO DO - see if results differ w/ and w/o values that fail QC
  olink_file = file.path(params$project_path, 'data', 'Olink', '20170803_McAlister_NPX_LOD.xlsx')  
  olink <- suppressMessages(read_olink(olink_file))
  
  # combine data types
  all_normalized_features <- dplyr::bind_rows(
    proteomics$proteins,
    peptidomics$peptides,
    small_molecules$all_small_molecules,
    olink$tall_measurements
    )
  
  # load phenotypes
  sample_phenotypes <- load_sample_phenotypes(params)
  
  # combine features with phenotypes
  
  features_with_design <- all_normalized_features %>%
    dplyr::mutate(set = as.factor(set)) %>%
    # add phenotypes
    dplyr::left_join(sample_phenotypes, by = c("set", "sample"), multiple = "all") %>%
    # find fold-changes w.r.t. the earliest age possible (keep as NA if the 8 month draw is missing)
    dplyr::group_by(dataset, groupId, Mouse.ID) %>%
    dplyr::mutate(initial_abundance = ifelse("8 months" %in% Age, normalized_log2_abundance[Age == "8 months"], NA_real_)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log2_fold_change = normalized_log2_abundance - initial_abundance)
  
  # test whether the join missed anything
  
  sample_mismatches <- features_with_design %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.na(Sex) | is.na(Age) | is.na(Mouse.ID) | is.na(lifespan)) %>%
    dplyr::distinct(dataset, sampleId)
  
  if (nrow(sample_mismatches) != 0) {
    stop(nrow(sample_mismatches), " samples could not be matches within a dataset")
  }
  
  output <- list(sample_phenotypes = sample_phenotypes,
                 featurization_parameters = featurization_parameters,
                 features_with_design = features_with_design,
                 small_molecules_features = small_molecules$small_molecules_features,
                 proteomics_features = proteomics$proteins_features,
                 peptidomics_features = peptidomics$peptidomics_features,
                 olink_features = olink$feature_attributes)
  
  return(output)
}

load_sample_phenotypes <- function (params) {
  
  # phenotypes
  mouse_phenotypes <- load_mouse_phenotypes(
    shock_token = params$shock_token,
    reduced_phenotypes = TRUE
    )
  
  # de-anonymize samples
  sample_labels <- load_sample_labels(params$repo_path)
  
  sample_phenotypes <- sample_labels %>%
    dplyr::left_join(mouse_phenotypes$mouse_level_phenotypes, by = "Mouse.ID") %>%
    dplyr::mutate(lifespan_remaining = lifespan - blood_draw_age_actual,
                  fraction_of_life_lived = blood_draw_age_actual / lifespan,
                  is_ddm_sample = ifelse(lifespan_remaining <= 20, TRUE, FALSE),
                  Draw_date = Birth_date + blood_draw_age_actual,
                  is_lbd_sample = dplyr::case_when(Generation == "G11" & Draw_date >= "2013-10-09" ~ TRUE,
                                                   Generation %in% c("G9", "G10") & Age == "20 months" ~ TRUE,
                                                   TRUE ~ FALSE),
                  is_eg8_sample = dplyr::case_when(Generation == "G8" & Draw_date < "2012-12-10" ~ TRUE,
                                                   TRUE ~ FALSE),
                  set = as.factor(set),
                  unique_sample_id = glue::glue("set{set}_{sample}")) %>%
    # calculate a covariate difference for eg8 (since fold-changes of two
    # is_eg8_samples effects would cancel each other out) eg8 is the only
    # affected categorical covariate, since is_ddm_sample and is_lbd_sample
    # only affect 14/20 month samples
    dplyr::group_by(Mouse.ID) %>%
    dplyr::mutate(eg8_ref_sample = is_eg8_sample[Age == "8 months"]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(is_eg8_fc_sample = as.logical(eg8_ref_sample - is_eg8_sample)) %>%
    dplyr::select(-eg8_ref_sample) %>%
    # combine EG8 and LBD effects into a single batch covariate
    dplyr::mutate(xs_batch = dplyr::case_when(is_lbd_sample ~ "LBD",
                                              is_eg8_sample ~ "EG8",
                                              TRUE ~ "Ref"),
                  xs_batch = factor(xs_batch, levels = c("Ref", "LBD", "EG8")),
                  fc_batch = dplyr::case_when(is_lbd_sample ~ "LBD",
                                              is_eg8_fc_sample ~ "EG8",
                                              TRUE ~ "Ref"),
                  fc_batch = factor(fc_batch, levels = c("Ref", "LBD", "EG8")))
  
  return (sample_phenotypes)
}

load_peptidomics <- function(peptidomics_dir){
  
  peptides <- readRDS(file.path(peptidomics_dir, 'peptidomics_features.RDS'))
  
  peptidomics_features <- peptides %>%
    dplyr::select(PeptideSequence, GeneName, Description, groupId) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      data_type = "peptides",
      NameUnique = PeptideSequence
      )
  
  peptides <- peptides %>%
    dplyr::rename(sampleId = sample_id, normalized_log2_abundance = normalized.log2.RSN) %>%
    dplyr::mutate(dataset = 'peptides', 
                  data_type = 'peptides',
                  set = as.numeric(set),
                  sample = stringr::str_c('s', stringr::str_sub(sampleId, 3))) %>%
    dplyr::select(dataset, data_type, groupId, sampleId, normalized_log2_abundance, set, sample, unique_batch)
  
  output <- list(peptides = peptides,
                 peptidomics_features = peptidomics_features)
  
  return(output)
}

load_proteomics <- function(proteomics_dir){
  
  compms <- readRDS(file.path(proteomics_dir, 'compms_proteomics_features.RDS'))
  compms$groupId <- compms %>% group_by(NameUnique) %>% group_indices()
  
  protein_features <- compms %>% 
    select(NameUnique, GeneName, ProteinId, groupId) %>%
    distinct() %>%
    dplyr::mutate(data_type = "tryptic-compms") %>%
    dplyr::select(data_type, groupId, NameUnique, GeneName, ProteinId)
  
  proteomics_batches <- compms %>%
    dplyr::distinct(set, unique_batch) %>%
    dplyr::mutate(batch = as.numeric(stringr::str_extract(unique_batch, "[0-9]+$"))) %>%
    dplyr::arrange(set, batch) %>%
    dplyr::mutate(unique_batch_new = glue::glue('set{set}_b{batch}'),
                  unique_batch_new = factor(unique_batch_new, levels = unique_batch_new)) %>%
    dplyr::select(unique_batch, unique_batch_new)
  
  compms_formatted <- compms %>%
    mutate(dataset = 'tryptic-compms',
           data_type = 'tryptic-compms',
           set = as.numeric(set),
           sample = stringr::str_c('s', stringr::str_sub(sample_id, 3)),
           weight = 1/variance) %>%
    dplyr::left_join(proteomics_batches, by = "unique_batch") %>%
    dplyr::select(-unique_batch) %>%
    select(dataset, data_type, groupId, sampleId = sample_id,
           normalized_log2_abundance = normalized.log2.RA, weight,
           set, sample, unique_batch = unique_batch_new)
  
  output <- list(proteins = compms_formatted,
                 proteins_features = protein_features)
  
  return (output)
}

load_small_molecules <- function (featurization_parameters, met_input) {
  
  checkmate::assertDataFrame(featurization_parameters)
  checkmate::assertChoice(met_input, c("valid", "normalized"))
  
  met_input_file <- if (met_input == "valid") {
    "valid_small_molecules.Rds"
  } else if (met_input == "normalized") {
    "normalized_small_molecules.Rds"
  } else {
    stop (met_input, " is not a defined value for \"met_input\", valid values are ", paste(valid_met_inputs, collapse = ", "))
  }
  
  # read datasets
  
  featurized_small_molecules <- featurization_parameters %>%
    dplyr::filter(featurization == "small-molecules",
                  eval == "TRUE") %>%
    dplyr::mutate(valid_peakgroups = file.path(dataset_path, "out_files", met_input_file),
                  peakgroup_file_found = file.exists(valid_peakgroups))
  
  # check for missing input
  
  missing_small_molecules <- featurized_small_molecules %>%
    dplyr::filter(!peakgroup_file_found)
  
  if (nrow(missing_small_molecules) != 0) {
    stop (missing_small_molecules %>% glue::glue_data("{met_input_file} missing for {data_partition}.  "))
  }
  
  # load results
  featurized_small_molecules <- featurized_small_molecules %>%
    dplyr::mutate(dat = purrr::map(valid_peakgroups, readRDS))
  
  # validate results
  purrr::walk(featurized_small_molecules$dat, romic::check_tomic)
  
  # summary of abundances
  all_small_molecules <- featurized_small_molecules %>%
    dplyr::mutate(tall_dataset = purrr::map(dat, small_molecules_tall_do_abundances)) %>%
    dplyr::select(dataset = data_partition, data_type, tall_dataset) %>%
    tidyr::unnest_legacy(tall_dataset)
  
  # join metabolites by groupId and lipids by name
  
  cross_set_lipids <- all_small_molecules %>%
    dplyr::filter(data_type %in% c("lipids-neg", "lipids-pos")) %>%
    dplyr::distinct(dataset, peak_label, adductName) %>%
    dplyr::group_by(peak_label, adductName) %>%
    dplyr::filter(n() == 2) %>%
    dplyr::distinct(peak_label, adductName) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(groupId = 1:dplyr::n())
    
  # only retain peakgroups which exist in all datasets
  cross_set_metabolites <- all_small_molecules %>%
    dplyr::filter(data_type %in% c("metabolites-pos", "metabolites-neg")) %>%
    dplyr::distinct(data_type, dataset, groupId) %>%
    dplyr::count(data_type, groupId) %>%
    dplyr::filter(n == max(n))
  
  # filter a couple outlier samples based on downstream EDA
  outlier_samples <- all_small_molecules %>%
    dplyr::distinct(dataset, data_type, set, sample) %>%
    dplyr::mutate(is_outlier = dplyr::case_when(
      dataset == "metabolites-neg-3" & sample %in% c("s29", "s107") ~ TRUE,
      dataset == "lipids-neg-3" & sample %in% "s77" ~ TRUE,
      dataset == "lipids-pos-12" & sample %in% "s36" ~ TRUE,
      TRUE ~ FALSE)) %>%
    dplyr::filter(is_outlier)
  
  if (nrow(outlier_samples) != 0) {
    message(outlier_samples %>% glue::glue_data("filtered {sample} from {dataset}.  "))
  }
  
  # infer limit of detection
  log2_floor_value <- min(all_small_molecules$normalized_log2_abundance)
  
  # create set12 + set3 metabolites and lipids
  normalized_small_molecules <- dplyr::bind_rows(
    all_small_molecules %>%
      dplyr::semi_join(cross_set_metabolites, by = c("data_type", "groupId")),
    all_small_molecules %>%
      dplyr::select(-groupId) %>%
      dplyr::filter(data_type %in% c("lipids-pos", "lipids-neg")) %>%
      dplyr::inner_join(cross_set_lipids, by = c("peak_label", "adductName"))
  ) %>%
    dplyr::select(-peak_label, -adductName) %>%
    dplyr::anti_join(outlier_samples, by = c("dataset", "data_type", "set", "sample")) %>%
    # normalize across datasets for the same data_type
    dplyr::group_by(dataset, data_type, groupId) %>%
    dplyr::mutate(.dataset_center = median(normalized_log2_abundance)) %>%
    dplyr::group_by(data_type, groupId) %>%
    dplyr::mutate(.data_type_center = median(normalized_log2_abundance)) %>%
    dplyr::ungroup() %>%
    # center each batch (dataset) then add back the original abundance scale
    # retain previously floored values at limit of detection
    dplyr::mutate(normalized_log2_abundance = dplyr::case_when(
      normalized_log2_abundance < log2_floor_value + 0.001 ~ log2_floor_value,
      TRUE ~ normalized_log2_abundance - .dataset_center + .data_type_center)
      ) %>%
    dplyr::select(-.dataset_center, -.data_type_center) %>%
    # floor to limit of detection
    dplyr::mutate(
      normalized_log2_abundance = ifelse(
        normalized_log2_abundance < log2_floor_value,
        log2_floor_value,
        normalized_log2_abundance)
      )
  
  # summarize feature mz, rt and average ic
  
  featurized_small_molecules_abund <- normalized_small_molecules %>%
    dplyr::group_by(data_type, groupId) %>%
    dplyr::summarize(mean_log2_abundance = mean(normalized_log2_abundance), .groups = "drop")
  
  small_molecule_features <- featurized_small_molecules %>%
    dplyr::mutate(features = purrr::map(dat, function(x){
      x$features
    })) %>%
    dplyr::select(data_type, features) %>%
    tidyr::unnest(features) %>%
    dplyr::mutate(
      peak_label = stringr::str_trim(peak_label),
      groupId = as.integer(as.character(groupId))
      )
  
  small_molecules_features <- dplyr::bind_rows(
    small_molecule_features %>%
      dplyr::filter(data_type %in% c("metabolites-neg", "metabolites-pos")) %>%
      dplyr::semi_join(cross_set_metabolites, by = c("data_type", "groupId")) %>%
      dplyr::mutate(groupId = as.integer(as.character(groupId))),
    small_molecule_features %>%
      dplyr::filter(data_type %in% c("lipids-neg", "lipids-pos")) %>%
      dplyr::select(-groupId) %>%
      dplyr::inner_join(cross_set_lipids, by = c("peak_label", "adductName"))
  ) %>%
    # only retain one description of each peakgroup
    dplyr::group_by(data_type, groupId) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    # add extra summaries
    dplyr::left_join(featurized_small_molecules_abund, by = c("data_type", "groupId"))
  
  missing_mets <- normalized_small_molecules %>%
    distinct(data_type, groupId) %>%
    dplyr::anti_join(small_molecules_features, by = c("data_type", "groupId"))
  
  if (nrow(missing_mets) != 0) {
    stop(
      missing_mets %>%
        glue::glue_data(
        "no feature annotations for groupId {groupId} of {data_type}.
        - ")
      )
  }

  output <- list(all_small_molecules = normalized_small_molecules,
                 small_molecules_features = small_molecules_features)
  
  return(output)
}

check_regression_summaries <- function (
  features_with_design,
  lm_models,
  lm_params,
  model_signif,
  model_summaries
) {
  
  lm_models
  model_summaries
  
  # check that features agree across summary tables
  distinct_features_with_design <- features_with_design %>%
    dplyr::distinct(data_type, groupId) %>%
    dplyr::mutate(is_features = TRUE)
  
  distinct_lm_models <- lm_models %>%
    dplyr::distinct(data_type, groupId) %>%
    dplyr::mutate(is_models = TRUE)
  
  distinct_lm_params <- lm_params %>%
    dplyr::distinct(data_type, groupId) %>%
    dplyr::mutate(is_params = TRUE)
  
  distinct_model_signif <- model_signif %>%
    dplyr::distinct(data_type, groupId) %>%
    dplyr::mutate(is_signif = TRUE)
  
  distinct_model_summaries <- model_summaries %>%
    dplyr::distinct(data_type, groupId) %>%
    dplyr::mutate(is_comparison = TRUE)
  
  misaligned_stats <-  distinct_lm_models %>%
    dplyr::full_join(distinct_lm_params, by = c("data_type", "groupId")) %>%
    dplyr::full_join(distinct_model_signif, by = c("data_type", "groupId")) %>%
    dplyr::full_join(distinct_model_summaries, by = c("data_type", "groupId")) %>%
    dplyr::mutate(
      is_models =  tidyr::replace_na(is_models, FALSE),
      is_params =  tidyr::replace_na(is_params, FALSE),
      is_signif =  tidyr::replace_na(is_signif, FALSE),
      is_comparison = tidyr::replace_na(is_comparison, FALSE)
    ) %>%
    dplyr::filter(!is_models | !is_params | !is_signif, !is_comparison)
  
  if (nrow(misaligned_stats) != 0) {
    stop("Different features were found between stat parameters, significance and fits")
  }
  
  misaligned_features_signif <- distinct_features_with_design %>%
    dplyr::full_join(distinct_model_signif, by = c("data_type", "groupId")) %>%
    dplyr::mutate(
      is_features =  tidyr::replace_na(is_features, FALSE),
      is_signif =  tidyr::replace_na(is_signif, FALSE)
    ) %>%
    dplyr::filter(!is_features | !is_signif) %>%
    tidyr::gather(problem, value, is_features, is_signif) %>%
    dplyr::filter(!value) %>%
    dplyr::count(data_type, problem) %>%
    dplyr::arrange(problem)
  
  if (nrow(misaligned_features_signif) != 0) {
    print(knitr::kable(misaligned_features_signif))
  }
  
  if ("is_features" %in% misaligned_features_signif$problem) {
    stop("Some features with signifiance lacked data")
  }
}
