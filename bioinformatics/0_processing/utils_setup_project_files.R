#' Detect Biological Enriched Peakgroups
#'
#' Identify peakgroups where the signal in biological samples is not higher than in negative controls.
detect_bioenriched_peakgroups <- function (
  mzroll_db_path,
  biosample_regex = "s[0-9]+",
  negsample_regex = "negctl",
  above_limit_tresh = 1
) {
  
  mzroll_list <- calicomics::process_mzroll(
    mzroll_db_path,
    only_identified = FALSE
    )
  
  mzroll_list$samples <- mzroll_list$samples %>%
    dplyr::mutate(sample_type = dplyr::case_when(
      stringr::str_detect(name, biosample_regex) ~ "bio",
      stringr::str_detect(name, negsample_regex) ~ "neg",
      TRUE ~ NA_character_
      )) %>%
    dplyr::filter(!is.na(sample_type))
  
  if (length(unique(mzroll_list$samples$sample_type)) != 2) {
    stop ("bio and negative samples were not detected")
  }
  
  mzroll_list <- calicomics:::reconcile_mzroll_list(mzroll_list)
  mzroll_list <- calicomics::floor_peaks(mzroll_list, 12)
  
  group_averages = mzroll_list$peaks %>%
    dplyr::left_join(mzroll_list$samples %>%
                       dplyr::select(sampleId, sample_type),
                     by = "sampleId") %>%
    dplyr::group_by(groupId, sample_type) %>%
    dplyr::summarize(mean_log2_abundance = mean(log2_abundance), .groups = "drop")
  
  bio_enrichments <-  group_averages %>%
    tidyr::spread(sample_type, mean_log2_abundance) %>%
    dplyr::mutate(bio_signal_enrichment = bio - neg)
  
  #library(ggplot2)
  #ggplot(bio_enrichments, aes(x = bio_signal_enrichment)) + geom_histogram(binwidth = 0.1) +
  #  geom_vline(xintercept = above_limit_tresh, color = "BLUE")
  
  bioenriched_peaks <- bio_enrichments %>%
    dplyr::filter(bio_signal_enrichment >= above_limit_tresh) %>%
    {.$groupId}
  
  return (bioenriched_peaks)
}

#' Reduce MzRoll Peakgroups
#'
#' Reduce an mzrollDB to a desired set of peakgroups"
reduce_mzroll_peakgroups <- function(mzroll_db_con, retained_groupIds) {
  
  stopifnot(class(retained_groupIds) == "integer")
  
  present_tables <- DBI::dbListTables(mzroll_db_con)
  updated_tables <- c("peakgroups", "peaks", "matches")
  
  for (a_table in updated_tables) {
    
    if (!(a_table %in% present_tables)) {
      warning (a_table, " is missing")
      next
    }
    
    sql_query <- paste0("SELECT * FROM ", a_table)
    
    updated_table <- dplyr::tbl(mzroll_db_con, dbplyr::sql(sql_query)) %>%
      dplyr::collect() %>%
      dplyr::filter(groupId %in% retained_groupIds)
    
    DBI::dbWriteTable(mzroll_db_con, a_table, updated_table, overwrite = TRUE)
  }
  
  # shrink database to filtered data
  DBI::dbSendQuery(mzroll_db_con, "VACUUM")
  
  invisible(0)
}

#' DO Filter Non-Biological Peakgroups
#'
#' Remove non-biological peakgroups from an mzroll
do_filter_nonbio_peakgroups_mzroll <- function (mzroll_db_path) {
  
  # find peakgroups which are enriched in biological samples relative
  #   negative controls
  bioenriched_peaks <- detect_bioenriched_peakgroups(mzroll_db_path)
  
  print(glue::glue("{length(bioenriched_peaks)} bio-enriched peaks"))
  
  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  reduce_mzroll_peakgroups(mzroll_db_con, bioenriched_peaks) 
  DBI::dbDisconnect(mzroll_db_con)
  
}

#' Reduce MzRoll Samples via Regular Expressions
#'
#' Identify samples whose name matches a regular expression and retain them.
reduce_mzroll_samples_regex <- function(mzroll_db_path, pattern) {
  
  #mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  mzroll_db_con <- clamshell::create_sqlite_con(mzroll_db_path)
  
  samples <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT * FROM samples")) %>%
    dplyr::collect()
  
  retained_sampleIds <- samples %>%
    dplyr::filter(stringr::str_detect(name, pattern)) %>%
    {.$sampleId}
  
  message("\nretain ", length(retained_sampleIds), " samples")
  
  clamr::reduce_mzroll_samples(mzroll_db_con, retained_sampleIds)
  DBI::dbDisconnect(mzroll_db_con)
}

#' Reduce DO
#'
#" Find a small number of samples from each batch and remove all other samples from an mzroll
reduce_do <- function(mzroll_db_path) {
  
  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  #mzroll_db_con <- clamshell::create_sqlite_con(mzroll_db_path)
  
  suppressWarnings(
    samples <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT * FROM samples")) %>%
      dplyr::collect() %>%
      dplyr::mutate(name = stringr::str_replace(name, 'DOplasma-', ''),
                    name = stringr::str_replace(name, 'X[0-9]{4}-M[0-9]{3}[A-Za-z]{1}-', '')) %>%
      tidyr::separate(name, into = c("set", "batch", "sample"))
  )
  
  reduced_sampleIds <- reduce_do_samples(samples)
  
  message("\nretain ", length(reduced_sampleIds), " samples")
  
  clamr::reduce_mzroll_samples(mzroll_db_con, reduced_sampleIds)
  DBI::dbDisconnect(mzroll_db_con)
}

#' Reduce DO samples
#'
#' Select samples to be retained in reduce_do()
reduce_do_samples <- function(samples, seed = 1234) {
  
  set.seed(seed)
  
  reduced_sampleIds <- dplyr::bind_rows(
    # set12
    samples %>%
      dplyr::filter(set %in% c("set1", "set2"),
                    stringr::str_detect(sample, '^s[0-9]')) %>%
      dplyr::group_by(set, batch) %>%
      dplyr::sample_n(1) %>%
      dplyr::ungroup(),
    
    # set3
    samples %>%
      dplyr::filter(set %in% c("set3"),
                    stringr::str_detect(sample, '^pos')) %>%
      dplyr::group_by(set, batch) %>%
      dplyr::sample_n(2) %>%
      dplyr::ungroup()
    
  ) %>%
    {.$sampleId}
  
  return(reduced_sampleIds)
}

#' Hoist Samples
#'
#' Find samples files associated with an mzrollDB
hoist_samples <- function(mzroll_db_path, raw_db) {
  
  sample_mzxml_dir <- dirname(dirname(raw_db))
  all_sample_paths <- tibble::tibble(name = list.files(sample_mzxml_dir, pattern = '.mzXML$'),
                                     path = file.path(sample_mzxml_dir, name))
  
  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  samples <- dplyr::tbl(mzroll_db_con, dbplyr::sql("SELECT name FROM samples")) %>%
    dplyr::collect() %>%
    dplyr::left_join(all_sample_paths, by = "name")
  
  if (any(is.na(samples$path))) {
    stop (sum(is.na(samples$path)), " samples could not be found")
  }
  
  return (samples)
}

#' Update Reduced MzXML Directory
#'
#' Remove and recreate reduced_mxml_dir
update_reduced_mzxml_dir <- function(reduced_mxml_dir) {
  
  if (file.exists(reduced_mxml_dir)) {
    unlink(reduced_mxml_dir, recursive = TRUE)
  }
  
  dir.create(reduced_mxml_dir)
}

#' Update Reduced MzXML Directory on Drive
#'
#' Upload a directory of .mzXMLs to Google Drive.
update_reduced_mzxml_dir_drive <- function(method_dir, DO_drive_path) {
  
  method_dir_drive_path = googledrive::drive_ls(DO_drive_path) %>%
    dplyr::filter(name == method_dir)
  
  if (nrow(method_dir_drive_path) != 0) {
    # delete existing folder if found
    googledrive::drive_rm(method_dir_drive_path)
  }
  
  googledrive::drive_mkdir(method_dir, DO_drive_path)
  
  method_dir_drive_path <- googledrive::drive_ls(DO_drive_path) %>%
    dplyr::filter(name == method_dir) %>%
    {googledrive::as_id(.$id[1])}
  
  return(method_dir_drive_path)
}

#' DO Copy Split
#'
#' Move an mzrollDB to a new location and then retain samples whose name matches a regular expression.
do_copy_split <- function(from, to, regex, overwrite) {
  
  # use return code to determine whether reduce should be performed.
  copy_ind <- clamshell::sqlite_copy(from, to, overwrite)
  
  if (copy_ind == 0) {
    reduce_mzroll_samples_regex(to, regex)
  }
}

patch_malformed_mzrolls <- function(labelled_mzroll, scan_mzroll) {
  labelled_mzroll_db_con <- clamshell::create_sqlite_con(labelled_mzroll)
  scan_mzroll_db_con <- clamshell::create_sqlite_con(scan_mzroll)
  
  if (!("scans" %in% DBI::dbListTables(scan_mzroll_db_con))) {
    stop(glue::glue("{scan_mzroll_db_con} does not have scans"))
  }
  
  peakgroups <- tbl(labelled_mzroll_db_con, sql("SELECT * FROM peakgroups")) %>% collect()
  peaks <- tbl(labelled_mzroll_db_con, sql("SELECT * FROM peaks")) %>% collect()
  samples <- tbl(labelled_mzroll_db_con, sql("SELECT * FROM samples")) %>% collect()
  
  DBI::dbWriteTable(scan_mzroll_db_con, "peakgroups", peakgroups, overwrite = TRUE)
  DBI::dbWriteTable(scan_mzroll_db_con, "peaks", peaks, overwrite = TRUE)
  DBI::dbWriteTable(scan_mzroll_db_con, "samples", samples, overwrite = TRUE)
  
  if ("compounds" %in% DBI::dbListTables(labelled_mzroll_db_con)) {
    compounds <- tbl(labelled_mzroll_db_con, sql("SELECT * FROM compounds")) %>% collect()
    DBI::dbWriteTable(scan_mzroll_db_con, "compounds", compounds, overwrite = TRUE)
  }
}

find_mzroll_systematic_compounds <- function (
  peakgroups,
  standard_databases,
  additional_smiles
  ) {
  
  systematic_compounds <- clamshell::query_systematic_compounds(
    systematic_compounds_con = standard_databases$systematic_compounds_con,
  )
  
  augmented_peakgroups <- suppressWarnings(
    peakgroups %>%
      tidyr::separate(compoundName, into = c("name", "smiles"), sep = ": ", remove = FALSE)
    )
  
  # match peakgroups to systematic compounds by name and structure
  
  name_direct_matches <- augmented_peakgroups %>%
    dplyr::select(groupId, name) %>%
    clamdb::match_structures(
      systematic_compounds$distinct_compounds,
      systematic_compounds$compound_aliases,
      match_var_order = c("name", "alias")
    ) %>%
    dplyr::select(groupId, systematicCompoundId)
  
  # extra manually provided smiles
  
  if (length(additional_smiles) != 0) {
    
    smiles_manual_matches <- tibble::tibble(
      name = names(additional_smiles),
      manual_smiles = unname(additional_smiles)
    )
    
    canonical_structures <- clamdb::canonicalize_structures(
      input_vector = smiles_manual_matches$manual_smiles,
      input_type = "smiles",
      output_types = c("smiles", "inchikey"),
      conda_env = "rdkit37"
    ) %>%
      tidyr::spread(type, canonical_ID)
    
    smiles_manual_matches <- smiles_manual_matches %>%
      dplyr::left_join(canonical_structures, by = c("manual_smiles" = "input"))  %>%
      dplyr::inner_join(augmented_peakgroups %>%
                          dplyr::select(groupId, name),
                        by = "name")
    
    if (nrow(smiles_manual_matches) == 0) {
      smiles_manual_matches <- smiles_manual_matches <- tibble::tibble()
    } else{
      smiles_manual_matches <- smiles_manual_matches %>%
        clamdb::match_structures(
          systematic_compounds$distinct_compounds,
          match_var_order = c("smiles", "inchikey", "inchikey_noproto", "inchikey_connectivity")
        ) %>%
        dplyr::select(groupId, systematicCompoundId)
      
      missed_manual_matches <- smiles_manual_matches %>%
        dplyr::filter(is.na(systematicCompoundId))
      
      if (nrow(missed_manual_matches) > 0) {
        warning(glue::glue(
          "{nrow(missed_manual_matches)} manually provided groupIds did not
          match the systematic compounds database:
          {paste(missed_manual_matches$groupId, collapse = ', ')}
         -- "))
      }
    }
  } else {
    smiles_manual_matches <- tibble::tibble()
  }
  
  # manual name -> smiles matches
  
  id_matches <- do.call(dplyr::bind_rows, 
                        list(name_direct_matches,
                             #smiles_compounds_matches,
                             #smiles_direct_matches,
                             smiles_manual_matches
                        )
  ) %>%
    dplyr::distinct() 
  
  degenerate_matches <- id_matches %>%
    dplyr::filter(!is.na(systematicCompoundId)) %>%
    dplyr::count(groupId) %>%
    dplyr::filter(n > 1)
  
  if (nrow(degenerate_matches) != 0) {
    warning (paste(degenerate_matches$groupId, collapse = ", "), " - matched multiple distinct IDs")
  }
  
  missing_matches <- id_matches %>%
    dplyr::group_by(groupId) %>%
    dplyr::filter(all(is.na(systematicCompoundId))) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(peakgroups %>%
                       dplyr::select(groupId, compoundName),
                     by = "groupId")
  
  if (nrow(missing_matches) != 0) {
    missing_matches %>%
      dplyr::slice(1:10) %>%
      glue::glue_data(
        "groupId {groupId}, compoundName: {compoundName} - matched zero IDs
       -- ") %>%
      warning()
  }
  
  matches <- id_matches %>%
    dplyr::filter(!is.na(systematicCompoundId)) %>%
    dplyr::left_join(systematic_compounds$distinct_compounds %>%
                       dplyr::select(systematicCompoundId,
                                     systematicName = name),
                     by = "systematicCompoundId")
  
  return (matches)
}

process_annotation_mzroll <- function (
  annotation_file,
  standard_databases,
  additional_smiles = c()
  ) {
  
  # reduce to just bookmarks and combine the bookmarks of
  # peakgroups with the same annotation
  
  mzroll_db_con <- clamr::mzroll_db_sqlite(annotation_file)
  
  # only look at bookmarks and good labels
  
  annotated_peakgroups <- dplyr::tbl(mzroll_db_con, "peakgroups") %>%
    dplyr::collect() %>%
    # retain just bookmarks
    dplyr::filter(searchTableName == "Bookmarks") %>%
    # drop isotopoologues
    dplyr::filter(!stringr::str_detect(tagString, "^(C13)|(N15)")) %>%
    # remove invalid peakgroups
    dplyr::filter(!(is.na(compoundName) & is.na(displayName) & is.na(compoundId)))
  
  # match peakgroup annotations to the systematic compounds database
  
  systematic_ids <- find_mzroll_systematic_compounds(annotated_peakgroups, standard_databases, additional_smiles)
  
  # reduce peaks based on retained peakgroups
  
  annotated_peakgroups <- annotated_peakgroups %>%
    dplyr::left_join(systematic_ids, by = "groupId") %>%
    # rename compoundName based on systematic and display names
    dplyr::mutate(compoundName = dplyr::case_when(
      !is.na(systematicCompoundId) ~ systematicName,
      !is.na(displayName) ~ displayName,
      !is.na(compoundName) ~ compoundName)
    ) %>%
    # remove adduct name from compound if one is present
    dplyr::mutate(
      compoundName = stringr::str_replace(
        compoundName,
        "\\[M[ +A-Za-z0-9-]+\\][-+]$",
        ""
        ),
      compoundName = stringr::str_replace(
        compoundName,
        " \\(?[0-9]+\\)?$",
        ""
        ),
      compoundName = stringr::str_trim(compoundName)
      ) %>%
    # number entries of the same compound
    dplyr::group_by(adductName, compoundName) %>%
    # extra a number label for compounds which are present multiple times
    dplyr::mutate(
      suffix = stringr::str_match(displayName, " \\(?([0-9])\\)?$")[,2],
      suffix = ifelse(is.na(suffix), 1, as.numeric(suffix))
    ) %>%
    dplyr::ungroup() %>%
    # multiple annotations of the same groupId can occur
    dplyr::group_by(groupId) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    # define compounds with multiple annotations
    dplyr::group_by(adductName, compoundName) %>%
    dplyr::mutate(multiannot = ifelse(any(suffix > 1), TRUE, FALSE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(compoundName = dplyr::case_when(
      multiannot ~ glue::glue("{compoundName} ({suffix})"),
      TRUE ~ compoundName
    )) %>%
    dplyr::select(-multiannot)
      
  annotated_peaks <- dplyr::tbl(mzroll_db_con, "peaks") %>%
    dplyr::collect() %>%
    dplyr::semi_join(annotated_peakgroups, by = "groupId")
  
  clamshell::sqlite_write(mzroll_db_con, "peaks", annotated_peaks, overwrite = TRUE)
  clamshell::sqlite_write(mzroll_db_con, "peakgroups", annotated_peakgroups %>% dplyr::select(-systematicCompoundId, -systematicName), overwrite = TRUE)
  clamshell::sqlite_write(mzroll_db_con, "systematic_ids", systematic_ids, overwrite = TRUE)
  
  # Combine peakgroups with identical annotations
  
  compound_updates <- annotated_peakgroups %>%
    dplyr::select(groupId_old = groupId, adductName, compoundName) %>%
    dplyr::group_by(compoundName, adductName) %>%
    dplyr::mutate(groupId_new = groupId_old[1]) %>%
    dplyr::ungroup() %>%
    dplyr::select(groupId_old, groupId_new)
  
  annotated_peakgroups %>%
    dplyr::inner_join(
      compound_updates %>%
        dplyr::group_by(groupId_new) %>%
        dplyr::filter(n() > 1),
      by = c("groupId" = "groupId_old")) %>%
    dplyr::arrange(adductName, compoundName)
  
  clamr::merge_peakgroups(mzroll_db_con, compound_updates)
  
  # shrink database to filtered data
  DBI::dbSendQuery(mzroll_db_con, "VACUUM")
  DBI::dbDisconnect(mzroll_db_con)
  
  invisible(0)
}

create_preliminary_peakgroup_annotations_tbl <- function (nfs_root, patching_paths, overwrite = FALSE) {
  
  annotations_dir <- file.path(nfs_root, "annotations")
  
  # one file to store all groupId annotations
  peakgroup_annotation_path <- file.path(annotations_dir, "peakgroup_annotations.tsv")
  
  if (overwrite | !file.exists(peakgroup_annotation_path)) {
    
    annotation_files <- tibble::tibble(annotation_path = patching_paths$out_file) %>%
      dplyr::mutate(
        method = stringr::str_extract(annotation_path, 'M[0-9]{3}'),
        data_type = dplyr::case_when(
          method == "M001" ~ "metabolites-neg",
          method == "M002" ~ "metabolites-pos",
          method == "M004" ~ "lipids-neg",
          method == "M005" ~ "lipids-pos"
          ))
    
    if (any(is.na(annotation_files$method))) {
      stop ("undefined method, add logic")  
    }
    
    peakgroup_annotations <- patching_paths %>%
      dplyr::left_join(annotation_files, by = c("out_file" = "annotation_path")) %>%
      dplyr::mutate(peakgroup_annotations = purrr::map2(scan_file, out_file, annotate_peakgroups)) %>%
      tidyr::unnest(peakgroup_annotations)
    
    readr::write_tsv(peakgroup_annotations, file = file.path(annotations_dir, "peakgroup_annotations.tsv"))
  } else {
    peakgroup_annotations <- readr::read_tsv(file = file.path(annotations_dir, "peakgroup_annotations.tsv"))
  }
  
  return (peakgroup_annotations)
}

#' Patch metaGroupId
#' 
#' Some file have metaGroupId as a numeric instead of an integer which
#'   causes schema errors
patch_metaGroupId <- function(mzroll_db_path) {
  
  mzroll_db_con <- clamshell::create_sqlite_con(mzroll_db_path)
  
  peakgroups <- dplyr::tbl(mzroll_db_con, "peakgroups") %>%
    dplyr::collect()
  
  if (class(peakgroups$metaGroupId) != "integer") {
    print(glue::glue("patching metaGroupId in {mzroll_db_path}"))
    patched_peakgroups <- peakgroups %>%
      dplyr::mutate(metaGroupId = as.integer(metaGroupId))
    
    DBI::dbWriteTable(mzroll_db_con, "peakgroups", patched_peakgroups, overwrite = TRUE)
  }
  
  DBI::dbDisconnect(mzroll_db_con)
}

annotate_peakgroups <- function (dataset_mzroll_db, annotation_mzroll_db) {
  
  patch_metaGroupId(dataset_mzroll_db)
  dataset_mzroll_db_con <- clamr::mzroll_db_sqlite(dataset_mzroll_db)
  dataset_peakset <- clamr::extract_peakset(dataset_mzroll_db_con)
  DBI::dbDisconnect(dataset_mzroll_db_con)

  patch_metaGroupId(annotation_mzroll_db)
  annotation_mzroll_db_con <- clamr::mzroll_db_sqlite(annotation_mzroll_db)
  annotation_peakset <- clamr::extract_peakset(annotation_mzroll_db_con)
  annotation_systematic_ids <- dplyr::tbl(annotation_mzroll_db_con, "systematic_ids") %>%
    dplyr::collect()
  
  annotation_peakgroups <- annotation_peakset$peakgroups %>%
    dplyr::select(groupId, label, compoundName, adductName) %>%
    dplyr::left_join(annotation_systematic_ids, by = "groupId")
  
  DBI::dbDisconnect(annotation_mzroll_db_con)
  
  # search for peakgroups which are identical between dataset and annotation
  annotation_matches <- clamr::match_peaksets(dataset_peakset, annotation_peakset, cutoff = 2) %>% 
    dplyr::left_join(
      dataset_peakset$peaks %>%
        dplyr::count(groupId),
      by = c("groupId_query" = "groupId")) %>%
    # enough samples overlap for correlation to be reliably assessed
    # peakgroups contain all samples' peaks
    dplyr::mutate(status = ifelse(n_samples >= 10 & n == max(n), "keep", "discard"))
  
  group_annotations <- annotation_matches %>%
    dplyr::select(
      groupId = groupId_query,
      annotation_groupId = groupId_target,
      n_samples = n,
      status
    ) %>%
    dplyr::left_join(annotation_peakgroups, by = c("annotation_groupId" = "groupId")) %>%
    dplyr::arrange(systematicName, compoundName)
  
  unusable_peakgroups <- annotation_peakgroups %>%
    dplyr::anti_join(group_annotations, by = c("groupId" = "annotation_groupId"))
  
  if (nrow(unusable_peakgroups) != 0) {
    warning(nrow(unusable_peakgroups), " annotations couldn't be used: ", paste(unusable_peakgroups$compoundName, collapse = ", "))
  }
  
  return (group_annotations)
}

#' Create Peakset
#"
#" Create a peakset from a provided table of peakgroups and peaks (which may be missing some features)
create_peakset <- function (peakgroups, peaks) {
  
  checkmate::assertDataFrame(peakgroups)
  checkmate::assertDataFrame(peaks)
  
  # create a table containing default fields required by the mzroll format.
  # required variables are not present since they must be user-specified.
  peakgroup_null_fields <- tibble::tibble(parentGroupId = 0L,
                                          tagString = "",
                                          metaGroupId = 0L,
                                          expectedRtDiff = -1,
                                          groupRank = 0,
                                          label = "",
                                          type = NA_integer_,
                                          srmId = "",
                                          ms2EventCount = NA_integer_,
                                          ms2Score = 0,
                                          adductName = NA_character_,
                                          compoundId = NA_character_,
                                          compoundName = NA_character_,
                                          compoundDB = NA_character_,
                                          searchTableName = NA_character_)
  
  missing_peakgroup_null_fields <- dplyr::select(peakgroup_null_fields, !!!rlang::syms(setdiff(colnames(peakgroup_null_fields), colnames(peakgroups))))
  
  peaks_null_fields <- tibble::tibble(pos = NA_integer_,
                                      minpos = NA_integer_,
                                      maxpos = NA_integer_,
                                      rt = NA_real_,
                                      rtmin = NA_real_,
                                      rtmax = NA_real_,
                                      mzmin = NA_real_,
                                      mzmax = NA_real_,
                                      scan = NA_integer_,
                                      minscan = NA_integer_,
                                      maxscan = NA_integer_,
                                      peakArea = NA_real_,
                                      peakAreaCorrected = NA_real_,
                                      peakAreaTop = NA_real_,
                                      peakAreaFractional = NA_real_,
                                      peakRank = NA_real_,
                                      peakIntensity = NA_real_,
                                      peakBaseLineLevel = NA_real_,
                                      peakMz = NA_real_,
                                      medianMz = NA_real_,
                                      baseMz = NA_real_,
                                      quality = NA_real_,
                                      width = NA_integer_,
                                      gaussFitSigma = NA_real_,
                                      gaussFitR2 = NA_real_,
                                      noNoiseObs = NA_integer_,
                                      noNoiseFraction = NA_real_,
                                      symmetry = NA_real_,
                                      signalBaselineRatio = NA_real_,
                                      groupOverlap = NA_real_,
                                      groupOverlapFrac = NA_real_,
                                      localMaxFlag = NA_real_,
                                      fromBlankSample = NA_integer_,
                                      label = NA_integer_)
  
  missing_peaks_null_fields <- dplyr::select(
    peaks_null_fields,
    !!!rlang::syms(setdiff(colnames(peaks_null_fields), colnames(peaks)))
  )
  
  mzroll_peakgroups <- tidyr::crossing(peakgroups, missing_peakgroup_null_fields)
  mzroll_peaks <- tidyr::crossing(peaks, missing_peaks_null_fields)
  
  return (clamr:::peakset_format_for_mzroll(list(peakgroups = mzroll_peakgroups,
                                                peaks = mzroll_peaks)))
}

add_eic_peaks <- function(annot_mzroll_db, method_eics) {
  
  ### Summarize the status of the TO BE updated mzrollDB file.
  
  annot_mzroll_con <- clamr::mzroll_db_sqlite(annot_mzroll_db)
  
  mzroll_samples <- tbl(annot_mzroll_con, "samples") %>%
    dplyr::collect() %>%
    dplyr::select(sampleId, filename) %>%
    dplyr::mutate(sample = basename(filename))
  
  max_groupId <- DBI::dbGetQuery(annot_mzroll_con, "SELECT MAX(groupId) FROM peakgroups")[[1]]
  max_peakId <- DBI::dbGetQuery(annot_mzroll_con, "SELECT COUNT(*) FROM peaks")[[1]]
  
  eic_peakgroups <- method_eics %>%
    dplyr::distinct(old_groupId = groupId, compoundName) %>%
    dplyr::mutate(
      adductName = stringr::str_extract(compoundName, "\\[M[+0-9A-Za-z-]*\\][+-]"),
      compoundName = stringr::str_replace(compoundName, "\\[M[+0-9A-Za-z-]*\\][+-]", ""),
      # update peakgroups
      groupId = max_groupId + 1:n(),
      searchTableName = "EICs"
      )
  
  # are there samples present which don't have any extracted EICs
  mismatched_samples <- mzroll_samples %>%
    dplyr::anti_join(method_eics, by = "sample")
  
  if (nrow(mismatched_samples) != 0) {
    warning(glue::glue(
      "reextracted peaks missing for {nrow(mismatched_samples)} samples 
      - in {annot_mzroll_db}
      - {paste(mismatched_samples$sample, collapse = ', ')}"
      ))
  }
  
  eic_peaks <- method_eics %>%
    dplyr::rename(old_groupId = groupId) %>%
    dplyr::left_join(eic_peakgroups %>%
                       dplyr::select(old_groupId, groupId),
                     by = "old_groupId") %>%
    dplyr::select(groupId, sampleId, peakAreaTop = intensity, mzmin, mzmax, rtmin = aligned_rtmin, rtmax = aligned_rtmax) %>%
    dplyr::mutate(peakId = max_peakId + 1:n())
  
  eic_peakset <- create_peakset(eic_peakgroups, eic_peaks)
  
  ### Update Annotation DB
  
  clamshell::sqlite_write(annot_mzroll_con, "peaks", eic_peakset$peaks, append = TRUE)
  clamshell::sqlite_write(annot_mzroll_con, "peakgroups", eic_peakset$peakgroups, append = TRUE)
  DBI::dbDisconnect(annot_mzroll_con)
  
  return (invisible(0))
}

#' Update Existing Peakgroups
#'
#' Add systematic annotations to confidently matched peakgroups and 
#'   remove borderline ones as they'll be re-quantified as EICs
update_existing_peakgroups <- function (mzroll_db_path, annotations) {
  
  mzroll_db_con <- clamr::mzroll_db_sqlite(mzroll_db_path)
  
  peakgroups <- dplyr::tbl(mzroll_db_con, "peakgroups") %>%
    dplyr::collect()
  
  updated_peakgroups <- peakgroups %>%
    dplyr::left_join(
      annotations %>%
        dplyr::select(
          groupId,
          status,
          annot_compoundName = compoundName,
          annot_adductName = adductName
        ) %>%
        dplyr::mutate(groupId = as.integer(groupId)),
      by = "groupId") %>%
    dplyr::mutate(
      compoundName = ifelse(
        !is.na(annot_compoundName),
        annot_compoundName,
        compoundName
      ),
      adductName = 
        dplyr::case_when(
          !is.na(annot_compoundName) ~ annot_adductName, # if there is an annotation but no adduct then keep that NA
          TRUE ~ adductName
        ),
      searchTableName = ifelse(
        !is.na(annot_compoundName),
        "Bookmarks",
        searchTableName
      )) %>%
    # drop peakgroups which match an annotation but are incomplete
    # these will be replaced by the reextracted intensities
    dplyr::filter(is.na(status) | status != "discard") %>%
    dplyr::select(groupId:displayName) %>%
    # drop all peakgroups which are duplicated
    # these represent cases where 2+ annotations are associatd with the
    # same peakgroup
    dplyr::group_by(groupId) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup()
  
  if (nrow(updated_peakgroups) == 0) {
    message("No update needed; continuing")
  } else {
    
    message(glue::glue(
      "{nrow(peakgroups) - nrow(updated_peakgroups)} peakgroups dropped due to
      - matching existing annotations. These annotated features will be
      - added using directly extracted EICs."
    ))
    
    updated_peaks <- dplyr::tbl(mzroll_db_con, "peaks") %>%
      dplyr::collect() %>%
      dplyr::semi_join(updated_peakgroups, by = "groupId")
    
    clamshell::sqlite_write(mzroll_db_con, "peaks", updated_peaks, overwrite = TRUE)
    clamshell::sqlite_write(mzroll_db_con, "peakgroups", updated_peakgroups, overwrite = TRUE)
    
  }
  
  DBI::dbDisconnect(mzroll_db_con)
  return (invisible(0))
}
