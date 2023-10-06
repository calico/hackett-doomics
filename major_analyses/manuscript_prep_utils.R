setup_figure_params <- function (config_path = "config.json") {
  # read config
  
  if (!file.exists(config_path) && !file.exists(config_path)) {
    stop (glue::glue("config.json not present at {config_path}"))
  }
  config <- jsonlite::read_json(config_path)
  
  required_json_vars <- c("cache_dir", "update_figures")
  missing_required_json_vars <- setdiff(required_json_vars, names(config))
  if (length(missing_required_json_vars) != 0) {
    stop (glue::glue("{length(missing_required_json_vars)} required variables are not specified in config.json: {paste(missing_required_json_vars, collapse = ', ')}"))
  }
  
  for (x in names(config)) {
    assign(x, config[[x]])
  }
  
  # Setup global parameters 
  
  if (!(dir.exists(cache_dir))) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  params <- list(
    # thks will be an empty string if you have not set this variable
    manuscript_drive_root = Sys.getenv("doomics_manuscript_uri"),
    repo_path = dirname(getwd())
  )
  
  # update options to avoid timeout when downloading large files
  options("timeout" = 10000)
  
  if (basename(params$repo_path) != "hackett-doomics") {
    cli::cli_abort(
      "Something went wrong, the directory above the working
      directory should be the repo root named {.var hackett-doomics} but
      it was {.path {params$repo_path}}
      ")
  }
  
  params$cache_dir <- cache_dir
  checkmate::assertDirectoryExists(params$cache_dir, access = "w")
  
  params$figures_dir <- file.path(params$cache_dir, "figures")
  if (!(dir.exists(params$figures_dir))) {
    dir.create(params$figures_dir)
  }
  checkmate::assertDirectoryExists(params$figures_dir, access = "w")
  
  checkmate::assertLogical(update_figures, len = 1)
  
  # additional setup for uploading results to google drive
  if (update_figures) {
    
    params$update_figures <- TRUE
    
    # drive ID for manuscript figures and other assets
    if (params$manuscript_drive_root == "") {
      cli::cli_abort(
        "{.var manuscript_drive_root} must be set in the .Renviron
        to export results to gDrive. This process is only supported
        for the study's authors. If you are an end-user then please set
        {.var update_figures} to FALSE in {.path config.json}"
        )
    }
    
    calibase::configure_google_access()
    # suppress drive printing
    options(googledrive_quiet = TRUE)
  } else {
    params$update_figures <- FALSE
  }
  
  return (params)
}

load_doomics <- function (
  an_asset = NULL,
  outdir = "/tmp/doomics/data",
  overwrite = FALSE
  ) {
  
  # Download data and modeling results
  
  checkmate::assertCharacter(outdir, len = 1)
  checkmate::assertLogical(overwrite, len = 1)
  
  if (!(dir.exists(outdir))) {
    dir.create(outdir, recursive = TRUE)
  }
  checkmate::assertDirectoryExists(outdir, access = "w")
  
  doomics_available_assets <- tibble::tribble(
    ~ asset, ~ object, ~ public_url, ~ description,
    "all_model_signif", "all_model_signif.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/all_model_signif.Rds", "FDR-controlled terms for all models",
    "normalized_peaks", "all_normalized_peaks.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/all_normalized_peaks.Rds", "peak level summaries of alternative normalization approaches for metabolomics and lipidomics",
    "gene_sets", "gene_sets", list(
      "m2.cp.reactome.v2022.1.Mm.entrez.gmt" = "https://storage.googleapis.com/calico-hackett-doomics-public/genesets/m2.cp.reactome.v2022.1.Mm.entrez.gmt",
      "m5.go.bp.v2022.1.Mm.entrez.gmt" = "https://storage.googleapis.com/calico-hackett-doomics-public/genesets/m5.go.bp.v2022.1.Mm.entrez.gmt"
    ), "Mouse genesets from MSigDB; these are bundled here to facilitate reproducible research. The primary source is http://www.gsea-msigdb.org/gsea/downloads.jsp.",
    "qtl_data", "all_qtl_results.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/all_qtl_results.Rds", "QTL lod profiles and null distributions for ~2200 molecules' abundances",   
    "feature_design_list", "feature_design_list.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/feature_design_list.Rds", "Extra metadata for individual data modalities",
    "features_extended", "features_with_design.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/features_with_design.Rds", "Measurements of all features in all samples (including data modalities not discussed in manuscript - Olink, and peptide-level analysis)",
    "features", "features_with_design_core.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/features_with_design_core.Rds", "Measurements of all features in all samples",
    "mouse_phenotypes", "mouse_phenotypes.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/mouse_phenotypes.Rds", "Mouse- and blood-draw-level phenotypes in both a reduced and comprehensive format.",
    "genetic_map", "genetic_map.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/genetic_map.Rds", "Genetic map mapping markers onto genomic coordinates",
    "protein_annotations", "proteomics_features_annotated.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/proteomics_features_annotated.Rds", "Protein annotatoins based on curating proteins associated with aging",
    "signif_all_models", "lm_models.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/lm_models.Rds", "Model-level summaries of all regressions",
    "signif", "model_signif.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/model_signif.Rds", "FDR-controlled terms for selected coefficients",
    "select_informatics_plots", "select_informatics_plots", list(
      "batch_post_correction-1.png" = "https://storage.googleapis.com/calico-hackett-doomics-public/select_informatics_plots/batch_post_correction-1.png",
      "metabolites-neg_injection_comparison.png" = "https://storage.googleapis.com/calico-hackett-doomics-public/select_informatics_plots/metabolites-neg_injection_comparison.png",
      "metabolites-pos_injection_comparison.png" = "https://storage.googleapis.com/calico-hackett-doomics-public/select_informatics_plots/metabolites-pos_injection_comparison.png",
      "tmt_ra_correlations-1.png" = "https://storage.googleapis.com/calico-hackett-doomics-public/select_informatics_plots/tmt_ra_correlations-1.png",
      "tmt_bridge_correlations-1.png" = "https://storage.googleapis.com/calico-hackett-doomics-public/select_informatics_plots/tmt_bridge_correlations-1.png"
    ), "Figures generated upstream which can be downloaded from GCS",
    "tomic_extended", "tidy_features.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/tidy_features.Rds", "features, samples, and measurements list (including data modalities not discussed in manuscript - Olink, and peptide-level analysis)",
    "tomic", "tidy_features_core.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/tidy_features_core.Rds", "features, samples, and measurements list",
    "tomic_peptides", "tidy_features_peptides.Rds", "https://storage.googleapis.com/calico-hackett-doomics-public/tidy_features_peptides.Rds", "features, samples, and measurements list (for peptide-level measurements only)"
  )
  
  available_assets_table <- function(doomics_available_assets) {
    # internal function
    doomics_available_assets %>%
      dplyr::select(asset, description) %>%
      mutate(description = stringr::str_wrap(description, 50)) %>%
      knitr::kable() %>%
      kableExtra::kable_styling() %>%
      print()
  }
  
  if (class(an_asset) == "NULL") {
    available_assets_table(doomics_available_assets)
    cli::cli_abort("You did not provide an input for {.var an_asset}. Please select a valid asset from {.val {paste(doomics_available_assets$asset, collapse = ', ')}}")
  }
  
  checkmate::assertString(an_asset)
  if (!(an_asset %in% doomics_available_assets$asset)) {
    available_assets_table(doomics_available_assets)
    cli::cli_abort("{.val {an_asset}} is not a valid value for {.var asset}. Please select a valid asset from {.val {paste(doomics_available_assets$asset, collapse = ', ')}}")
  }
  
  # does the asset already exist?
  asset_metadata <- doomics_available_assets %>% dplyr::filter(asset == an_asset)
  stopifnot(nrow(asset_metadata) == 1)
  
  asset_local_path <- file.path(outdir, asset_metadata$object[1])
  
  return_functional <- function (asset_local_path) {
    if (stringr::str_detect(asset_local_path, "\\.Rds$")) {
      return(readRDS(asset_local_path))
    } else {
      return(asset_local_path)
    }
  }
  
  if (file.exists(asset_local_path) && !overwrite) {
    return(return_functional(asset_local_path))
  }
  
  # if the object is a file we can directly load it, otherwise download a
  # set of files
  
  payload <- asset_metadata$public_url[[1]]
  if (class(payload) == "character") {
    utils::download.file(payload, destfile = asset_local_path)
  } else if (class(payload) == "list") {
    
    if (!(dir.exists(asset_local_path))) {
      dir.create(asset_local_path, recursive = TRUE)
    }
    
    purrr::walk2(
      names(payload),
      unname(payload),
      ~ utils::download.file(.y, file.path(asset_local_path, .x))
      )
  } else {
    cli::cli_abort("public_url is a {.val {class(payload)[1]}} and must be a character or list")
  }
  
  return(return_functional(asset_local_path))
}

create_and_upload_figure <- function (
  grob = NULL,
  params,
  name,
  drive_path,
  width,
  height,
  extensions = ".pdf"
  ) {
  
  # ggsave to a local file and add local plots to gDrive if update_figures is TRUE
  
  if (is.null(grob)) {
    grob = last_plot()
  }
  checkmate::assertClass(grob, "ggplot")
  
  checkmate::assertLogical(params$update_figures, len = 1)
  if (!params$update_figures) {
    # no upload needed
    return (invisible(0))
  }
  
  checkmate::assertCharacter(extensions)
  purrr::walk(extensions, checkmate::assertChoice, choices = c(".pdf", ".png", ".eps"))
  
  checkmate::assertString(params$figures_dir)
  checkmate::assertPathForOutput(file.path(params$figures_dir, "foo.txt"))
  checkmate::assertString(name)
  checkmate::assertNumber(width)
  checkmate::assertNumber(height)
  
  manuscript_drive_root_con <- googledrive::as_id(params$manuscript_drive_root)
  
  plots_to_save <- tibble::tibble(
    extension = extensions
  ) %>%
    mutate(
      filename = glue::glue("{name}{extension}"),
      localpath = file.path(params$figures_dir, filename)
    )
  
  if (".pdf" %in% extensions) {
    ggsave(
      plot = grob,
      filename = plots_to_save$localpath[plots_to_save$extension == ".pdf"],
      height = height,
      width = width,
      useDingbats = FALSE
    )
  }
  
  vanilla_extensions <- setdiff(extensions, ".pdf")
  if (length(vanilla_extensions) != 0) {
    plots_to_save %>%
      dplyr::filter(extension %in% vanilla_extensions) %>%
      {purrr::walk(.$localpath, ggsave, height = height, width = width, plot = grob)}
  }
  
  purrr::walk(
    plots_to_save$localpath,
    calibase::upload_to_drive,
    drive_path = drive_path,
    drive_root = manuscript_drive_root_con
  )
}

add_feature_label <- function(feature_df, feature_design_list) {
  
  # Rename features based on feature_name and proteomics_features metadata
  
  checkmate::assertDataFrame(feature_df)
  checkmate::assertChoice("feature_name", colnames(feature_df))
  
  update_protein_labels <- feature_design_list$proteomics_features %>%
    separate(ProteinId, into = c("src", "uniprot", "common", "species")) %>%
    mutate(feature_label = case_when(
      !is.na(GeneName) & stringr::str_length(GeneName) > 2 ~ stringr::str_to_title(GeneName),
      TRUE ~ uniprot
    ))
  
  stopifnot(update_protein_labels %>% filter(is.na(feature_label)) %>% nrow() == 0)
  
  feature_df_w_labels <- feature_df %>%
    left_join(
      update_protein_labels %>%
        select(feature_name = NameUnique, feature_label),
      by = "feature_name"
    ) %>%
    mutate(feature_label = ifelse(is.na(feature_label), feature_name, feature_label))
  
  return(feature_df_w_labels)
}

add_term_labels <- function(terms_df) {
  terms_df %>%
    mutate(
      term_label = stringr::str_wrap(stringr::str_to_title(term), 15),
      term_label = ifelse(term == "ddm", "DDM", term_label),
      term_label = factor(term_label, levels = term_label),
      term_label_simple = stringr::str_replace(term_label, " ?[\\(\\|][A-Za-z \\)]+$", "")
    )
}

get_subset_correlation <- function (
  feature_names,
  all_phenotype_abundances,
  cor_method = "spearman"
  ) {
  
  checkmate::assertDataFrame(all_phenotype_abundances)
  purrr::walk(feature_names, ~ checkmate::assertChoice(., unique(all_phenotype_abundances$feature_name)))
  feature_names <- unique(feature_names)
  stopifnot(length(feature_names) > 1)  
  
  X <- all_phenotype_abundances %>%
    filter(feature_name %in% feature_names) %>%
    reshape2::acast(
      feature_name ~ unique_sample_id,
      value.var = "corrected_log2_abundance"
    )
  
  all_feature_correlations <- cor(t(X), use = "pairwise.complete.obs", method = cor_method) %>%
    as.data.frame() %>%
    mutate(feature_2 = rownames(.)) %>%
    tidyr::gather(feature_1, corr, -feature_2) %>%
    tibble::as_tibble() %>%
    filter(feature_1 != feature_2)
  
  return(all_feature_correlations)
}

# Figure 3

calculate_category_enrichment <- function (
  category_members,
  all_analytes,
  term_hits
  ) {
  
  # generate 2x2 table for a functional enrichment and apply a Fisher Exact
  # test of enrichment for category membership among differential
  # expression/abundance features
  
  # significant changes among analytes being considered
  valid_term_hits <- term_hits %>%
    semi_join(all_analytes, c("data_type", "groupId"))
  
  # df of category members which are also diffex
  category_hits_df <- category_members %>%
    semi_join(valid_term_hits, by = c("data_type", "groupId"))
  
  two_by_two_vals <- c(
    category_hits_df %>% nrow(), # category hit
    category_members %>% anti_join(valid_term_hits, by = c("data_type", "groupId")) %>% nrow(), # category non-hit
    valid_term_hits %>% anti_join(category_members, by = c("data_type", "groupId")) %>% nrow() # non-category hit
  )
  
  contingency_table <- c(two_by_two_vals, nrow(all_analytes) - sum(two_by_two_vals)) %>%
    matrix(byrow = TRUE, ncol = 2, nrow = 2)
  
  fisher.test(contingency_table, alternative = "greater") %>%
    broom::tidy() %>%
    select(estimate, p.value, conf.low) %>%
    mutate(
      category_hits = nrow(category_hits_df),
      category_size = nrow(category_members)
    ) %>%
    # add a table with category hits
    crossing(
      tidyr::nest(
        category_hits_df %>%
          select(data_type, groupId),
        category_hits_df = c(data_type, groupId)
        ))
}

drop_degenerate_categories <- function (
  term_discoveries,
  term_signif_hits,
  jaccard_cutoff = 0.7
  ) {
  
  # flag degenerate functional categories which are less significant
  # than another category and have highly overlapping membership
  
  # generate all pairs of significant categories
  go_pairs <- crossing(
    term_discoveries %>%
      select(category1 = category, value1 = value),
    term_discoveries %>%
      select(category2 = category, value2 = value)
  ) %>%
    # add members
    left_join(protein_category_members %>% dplyr::rename(category1 = category, value1 = value, members1 = members), by = c("category1", "value1"), multiple = "all") %>%
    left_join(protein_category_members %>% dplyr::rename(category2 = category, value2 = value, members2 = members), by = c("category2", "value2"), multiple = "all") %>%
    # ignore self-self
    filter(value1 != value2) %>%
    mutate(jaccard_index = purrr::map2_dbl(members1, members2, find_jaccard_index, term_signif_hits = term_signif_hits))
  
  # filter to degenerate pathways
  degenerate_go_terms <- go_pairs %>%
    filter(jaccard_index >= jaccard_cutoff)
  
  # order by p/qvalue and starting with the strongest enrichments drops degenerate categories
  running_terms <- term_discoveries %>%
    arrange(qvalue, p.value)
  nondegenerate_terms <- NULL
  
  while(nrow(running_terms) > 0) {
    
    # take the most signifiant term
    new_term <- running_terms %>% dplyr::slice(1)
    nondegenerate_terms <- dplyr::bind_rows(nondegenerate_terms, new_term %>% dplyr::select(category, value))
    running_terms <- running_terms %>% dplyr::slice(-1)
    
    # drop all lower significance degenerate term
    terms_to_drop <- new_term %>%
      inner_join(
        degenerate_go_terms,
        by = c("category" = "category1", "value" = "value1"),
        multiple = "all"
        )
    
    if (nrow(terms_to_drop) > 0) {
      running_terms <- running_terms %>%
        anti_join(terms_to_drop, by = c("category" = "category2", "value" = "value2"))
    }
  }
  
  term_discoveries %>%
    left_join(
      nondegenerate_terms %>%
        select(category, value) %>%
        mutate(distinct_term = TRUE),
      by = c("category", "value")) %>%
    mutate(distinct_term = ifelse(is.na(distinct_term), FALSE, TRUE))
}

find_jaccard_index <- function(members1, members2, term_signif_hits) {
  genes1 <- members1$gene_id[members1$groupId %in% term_signif_hits$groupId]
  genes2 <- members2$gene_id[members2$groupId %in% term_signif_hits$groupId]
  # intersection / union
  length(intersect(genes1, genes2)) / length(union(genes1, genes2))
}

plot_fxnl_enrichment_beeswarm <- function (
    gsea_category_plot_df_subset,
    scale_for_per_day_effect_sizes = FALSE,
    facet_formula = as.formula(term_label_simple ~ .)
) {
  
  checkmate::assertLogical(scale_for_per_day_effect_sizes, len = 1)
  if (scale_for_per_day_effect_sizes) {
    estimate_var <- "scaled_effect_size"
    x_axis_label <- expression("Regression" ~ log[2] ~ "fold-change Per 100 Days")
  } else {
    estimate_var <- "estimate"
    x_axis_label <- expression("Regression" ~ log[2] ~ "fold-change")
  }
  
  beeswarm_theme <- theme_bw() +
    theme(
      text = element_text(size = 10),
      axis.text.y = ggtext::element_markdown(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  gsea_category_plot_df_subset %>%
    filter(!is.na(!!rlang::sym(estimate_var))) %>%
    ggplot(aes(x = term_go_pair, y = !!rlang::sym(estimate_var), color = is_signif)) +
    ggbeeswarm::geom_beeswarm(size = 1, cex = 1) +
    geom_hline(yintercept = 0, color = "gray40", linetype = 3) +
    facet_grid(facet_formula, scales = "free", space = "free") +
    scale_x_discrete(NULL, breaks = gsea_category_plot_df_subset$term_go_pair, labels = gsea_category_plot_df_subset$category_label_w_star_md, limits=rev, expand = c(0,1.5)) +
    scale_y_continuous(x_axis_label) +
    scale_color_manual("Significance", values = c("q < 0.1" = "#5BB867", "n.s." = "gray40")) +
    guides(color=guide_legend(override.aes = list(size=5))) +
    coord_flip() +
    beeswarm_theme
}

calculate_signed_nlog10p = function(estimate, pvalue) {
  ifelse(estimate < 0, -1, 1) * -1*log10(pvalue)
}

# figure 4

bootstrap_relative_prob <- function (category_data, nbs = 1000) {
  
  purrr::map(
    seq(nbs),
    ~ category_data %>%
      dplyr::sample_n(size = dplyr::n(), replace = TRUE) %>%
      mutate(bs_num = .x) %>%
      unnest(feature_model_results)
  ) %>%
    # combine all bootstraps
    dplyr::bind_rows() %>%
    # sum up each model's AIC across features
    group_by(bs_num, model_name) %>%
    summarize(
      aic_sum = sum(AICc),
      .groups = "drop"
    ) %>%
    # summarize the relative fits of models for each bootstrap
    group_by(bs_num) %>%
    mutate(
      relative_lik = exp((min(aic_sum) - aic_sum)/2),
      relative_prob = relative_lik / sum(relative_lik)
    )
}

data_modality_colorscheme <- tibble::tribble(
  ~ data_modality, ~ color,
  "proteomics", "royalblue1",
  "metabolomics", "red3",
  "lipidomics", "limegreen",
  "non-omics", "gray50"
)

light_label_categories <- c(
  "Metabolomics Knowns",
  "Metabolism",
  "Signaling",
  "Stress",
  "Other Proteomics",
  "Chronological Aging Lipids",
  "Lipids with 24:1 (FLL-associated)"
)
# Figure 5


min_dmeasure <- function (
    sample_dat,
    dmeasure_method = "lm",
    dmeasure_vals = seq(-10, 10, by = 0.001)
    ) {
  
  checkmate::assertDataFrame(sample_dat)
  checkmate::assertChoice(dmeasure_method, c("mle", "lm"))
  checkmate::assertNumeric(dmeasure_vals)
  
  if (dmeasure_method == "lm") {
    
    min_dmeasure <- lm(formula = .resid ~ estimate + 0, data = sample_dat) %>%
      broom::tidy() %>%
      {.$estimate}
    
    # the Maximum Likelihood Estimator (mle) method confirms that lm solves
    # the appropriate problem and generates
    # equivalent results. There should be no reason to use the mle method.
    # the only reason this method is it took me a little while to realize
    # that this was just a regression problem.
  } else if (dmeasure_method == "mle") {
    
    # scaled dmeasure is a range of possible values of delta-age
    scaled_dmeasure = dmeasure_vals / sqrt(mean(sample_dat$estimate^2))
    
    dmeasure_df <- tibble::tibble(dmeasure = scaled_dmeasure) %>%
      tidyr::crossing(sample_dat) %>%
      mutate(adjusted_resid = .resid - dmeasure * estimate) %>%
      summarize(
        RSS = sum(adjusted_resid^2),
        .by = dmeasure
      )
    
    min_dmeasure <- dmeasure_df$dmeasure[which.min(dmeasure_df$RSS)]
    if (min_dmeasure %in% c(min(dmeasure_df$dmeasure), max(dmeasure_df$dmeasure))) {
      stop(glue::glue(
        "Adjusted residuals are minimized at the boundary of delta measures = {min_dmeasure}
      Please expand the optimization window
      "
      ))
    } 
  } else {
    stop("unexpected behavior")
  }
  
  return(min_dmeasure)
}

infer_age_measure_given_go <- function(
    members,
    an_aging_archetype,
    signif_all_models
    ) {
  
  aging_archetype_attrs <- age_driver_model_definitions %>%
    dplyr::filter(aging_archetype == an_aging_archetype) %>%
    as.list()
  
  stopifnot(length(aging_archetype_attrs) == 4)
  
  # attributes to multiple to recover original entry in design matrix
  design_measure <- unlist(aging_archetype_attrs$measure_list)
  
  member_model_fits <- signif_all_models %>%
    # filter to correct model
    dplyr::filter(
      model_name == aging_archetype_attrs[["model_name"]]
    ) %>%
    # reduce to relevant members
    semi_join(members, by = c("data_type", "groupId"))
  
  member_fits <- member_model_fits %>%
    dplyr::select(data_type, groupId, params, fits) %>%
    unnest(params) %>%
    filter(term == aging_archetype_attrs[["relevent_term"]]) %>%
    unnest(fits)
  
  go_category_deltaage <- member_fits %>%
    # sample IDs are inconsistent in this fits
    # since they occurred in the analysis before data types were harmonized
    dplyr::select(-sampleId) %>%
    tidyr::nest(sample_data = -c(set:is_ddm_sample)) %>%
    dplyr::mutate(
      dmeasure = purrr::map_dbl(sample_data, min_dmeasure)
    )
  
  stopifnot(nrow(go_category_deltaage) <= 320) # less than or equal to number of non-DDM samples
  
  # add dmeasure to original aging measure
  
  normalized_aging_measure <- go_category_deltaage %>%
    rowwise() %>%
    mutate(age_measure = prod(!!!rlang::syms(design_measure))) %>%
    ungroup() %>%
    mutate(age_measure_w_dmeasure = age_measure + dmeasure) %>%
    # standardize age measures across individuals
    # rescale since measures are in original units (e.g., days {of lifespan} x days of {age})
    mutate(across(c(dmeasure, age_measure, age_measure_w_dmeasure), ~ c(scale(.x, center = FALSE))))
  
  return(normalized_aging_measure)
}

# figure 5
create_aging_measures_heatmap <- function (go_age_measures, samples) {
  
  stratify_by_measures <- "aging_archetype" %in% colnames(go_age_measures)  
  if (!stratify_by_measures) {
    aging_archetype <- "fraction of life lived"
  }
  
  go_age_measures_tomic <- go_age_measures %>%
    select(data_modality:category_color, any_of("aging_archetype"), go_age_measure) %>%
    unnest(go_age_measure) %>%
    select(-sample_data, -unique_batch) %>%
    # add additional sample metadata (a common primary key)
    left_join(
      samples %>% select(unique_sample_id, set, sample),
      by = c("set", "sample")
    ) %>%
    dplyr::mutate(
      consensus_measure = dplyr::case_when(
        aging_archetype == "age x lifespan" ~ -1 * dmeasure,
        aging_archetype == "lifespan remaining" ~ -1 * age_measure_w_dmeasure,
        TRUE ~ age_measure_w_dmeasure
      ),
      category_color_hex = gplots::col2hex(category_color),
      category_label_md = glue::glue("<span style='color:{category_color_hex}'><b>{category_label}</b></span>")
    ) %>%
    romic::create_tidy_omic(
      feature_pk = "category_label_md",
      feature_vars = c("category_label", "data_modality", "category", "value", "category_general_label", "category_color") %>%
        {c(., if(stratify_by_measures){"aging_archetype"}else{NULL})}, 
      sample_pk = "unique_sample_id",
      sample_vars = c("sample", "set", "Mouse.ID", "Age", "Sex", "blood_draw_age_actual", "lifespan", "lifespan_remaining", "fraction_of_life_lived", "is_ddm_sample")
    ) %>%
    romic::center_tomic()
  
  aging_measures_sample_heatmap <- romic::plot_heatmap(
    go_age_measures_tomic,
    value_var = "consensus_measure",
    sample_var = "fraction_of_life_lived",
    cluster_dim = "rows",
    change_threshold = 3,
    x_label = "Samples ordered by fraction of life lived",
    y_label = "",
    colorbar_label = "Relative age given\npathway abundances"
  ) + 
    theme(
      text = element_text(size = 10),
      legend.title = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.text.y.right = ggtext::element_markdown(),
      strip.text = element_text(size = 12)
    )
  
  if (stratify_by_measures) {
    aging_measures_sample_heatmap <- aging_measures_sample_heatmap +
      facet_grid(aging_archetype ~ ., scales = "free", space = "free", switch = "y")
  } else{
    aging_measures_sample_heatmap <- aging_measures_sample_heatmap +
      facet_grid(. ~ Age, scales = "free", space = "free", switch = "y")
  }
  
  return(aging_measures_sample_heatmap)
}

calculate_cross_group_distance <- function(
  category_data_1,
  category_data_2,
  feature_feature_similarity,
  n_permutations = 100
  ) {
  
  # add intersection to the smaller category so that we can compare
  # a nested gene set to its parent category (i.e., protein catabolism to proteolysis)
  
  if (nrow(category_data_1) > nrow(category_data_2)) {
    major_category <- category_data_1
    minor_category <- category_data_2
  } else {
    major_category <- category_data_2
    minor_category <- category_data_1
  }
  
  # drop entries from major category which are present in minor category
  major_category <- anti_join(
    major_category,
    minor_category,
    by = "feature_name"
  )
  
  # create major and minor clusters permuting initial category labels
  compared_categories <- bind_rows(
    major_category %>%
      mutate(category = "major"),
    minor_category %>%
      mutate(category = "minor")
  ) %>%
    mutate(null_i = -1)
  
  permuted_categories <- tibble::tibble(
    null_i = seq_len(n_permutations)
  ) %>%
    mutate(null_categories = purrr::map(null_i, ~ compared_categories %>% select(feature_name) %>% mutate(category = sample(compared_categories$category, replace = FALSE)))) %>%
    unnest(null_categories)
  
  all_categories <- bind_rows(
    compared_categories,
    permuted_categories
  ) %>%
    calculate_category_homogeneity(feature_feature_similarity = feature_feature_similarity)
  
  # for the real categories determine whether the major or minor
  # category is more coherent for each measure
  
  category_coherence <- all_categories %>%
    filter(null_i == -1) %>%
    select(-null_i) %>%
    gather(measure, value, -category)
  
  top_coherent_category <- bind_rows(
    category_coherence %>%
      filter(measure == "median_distance") %>%
      arrange(value) %>%
      slice(1),
    category_coherence %>%
      filter(measure != "median_distance") %>%
      arrange(desc(value)) %>%
      group_by(measure) %>%
      slice(1)
  )
  
  summary_versus_permutation <- all_categories %>%
    gather(measure, value, -category, -null_i) %>%
    # average major and minor summary for each null_i
    semi_join(top_coherent_category, by = c("category", "measure")) %>%
    arrange(value) %>%
    group_by(measure) %>%
    mutate(i_rank = 1:n()) %>%
    filter(null_i == -1) %>%
    select(-value, -category, -null_i) %>%
    spread(measure, i_rank) %>%
    # flip the order
    mutate(
      mean_abs_corr = (n_permutations + 1) - mean_abs_corr + 1,
      median_abs_corr = (n_permutations + 1) - median_abs_corr + 1
    ) %>%
    mutate(across(mean_abs_corr:median_distance, ~ . / n_permutations))
  
  coherence_summary <- top_coherent_category %>% select(measure, coherence = value) %>%
    left_join(
      summary_versus_permutation %>%
        gather(measure, coherence_null_quantile),
      by = "measure")
  
  return(coherence_summary)
}

calculate_category_homogeneity <- function (categories_w_members, feature_feature_similarity) {
  
  homogeneity <- categories_w_members %>%
    rename(feature_1 = feature_name) %>%
    # compare all members - all members
    left_join(
      categories_w_members %>% rename(feature_2 = feature_name),
      by = c("category", "null_i"),
      relationship = "many-to-many"
    ) %>%
    filter(feature_1 < feature_2) %>%
    inner_join(feature_feature_similarity, by = c("feature_1", "feature_2")) %>%
    group_by(category, null_i) %>%
    summarize(
      median_distance = median(distance),
      mean_abs_corr = mean(abs(corr)),
      median_abs_corr = median(abs(corr)),
      .groups = "drop"
    )
  
  return(homogeneity)
}

# figure 6
visualize_feature_subsets <- function (feature_subset, features_with_design, top_partial_corrs, facet_by = NULL) {
  
  # validate inputs
  stopifnot(all(c("feature_name", "feature_label", "feature_aging_archetype", "change_w_age") %in% colnames(feature_subset)))
  stopifnot(length(feature_subset$feature_name) == length(unique(feature_subset$feature_name)))
  
  if (!all(feature_subset$feature_name %in% features_with_design$feature_name)) {
    missing_features <- setdiff(feature_subset$feature_name, features_with_design$feature_name)
    stop(glue::glue(
      "{length(missing_features)} features were present in feature_subset but not features_with_design including:\n",
      paste(missing_features[1:pmin(5, length(missing_features))], collapse = ", ")
    ))
  }
  
  # cluster features
  centered_subset_abundances <- features_with_design %>%
    semi_join(feature_subset, by = "feature_name") %>%
    group_by(data_type, groupId) %>%
    mutate(corrected_log2_abundance = corrected_log2_abundance - median(corrected_log2_abundance, na.rm = TRUE)) %>%
    ungroup()
  
  hclust_feature_orders <- centered_subset_abundances %>%
    romic::hclust_order(
      feature_pk = "feature_name",
      sample_pk = "unique_sample_id",
      value_var = "corrected_log2_abundance",
      cluster_dim = "rows",
      distance_measure = "corr"
    )
  
  # order features and add other metadata
  feature_subset_annot <- feature_subset %>%
    mutate(
      change_w_age_symbol = ifelse(change_w_age == "up", "+", "-"),
      feature_name = factor(feature_name, levels = hclust_feature_orders$rows)
    )
  
  if (class(facet_by) != "NULL") {
    checkmate::assertChoice(facet_by, colnames(feature_subset))
    
    # decide order of categories based on general ordering of features
    ordered_facet_categories <- feature_subset_annot %>%
      mutate(feature_order_int = as.integer(feature_name)) %>%
      group_by(!!rlang::sym(facet_by)) %>%
      summarize(category_order = mean(feature_order_int), .groups = "drop") %>%
      arrange(category_order) %>%
      select(-category_order) %>%
      mutate(
        ordered_category = !!rlang::sym(facet_by),
        # replace double backslashes which were escaped when loading
        # categories from googlesheets. We want the \n to be treated as a
        # newline
        ordered_category = stringr::str_replace(ordered_category, "\\\\n", "\n"),
        ordered_category = factor(ordered_category, levels = ordered_category)
      )
    
    # reaorder features based on category orders; they should be grouped together
    # first by category and then by hclust ordering so the upper-diagonal -> corr;
    # lower-diagonal -> partial-corr visualization works. This is setup
    # by distinguishing corrs and partials in corr_or_partial based on feature order
    
    feature_subset_annot <- feature_subset_annot %>%
      left_join(ordered_facet_categories, by = facet_by) %>%
      mutate(
        category1 = ordered_category,
        category2 = forcats::fct_rev(ordered_category)
      ) %>%
      # reorder by ordered_category -> feature_name
      arrange(ordered_category, feature_name) %>%
      mutate(feature_name = factor(feature_name, levels = feature_name))
  }
  
  feature_corrs <- widyr::pairwise_cor(
    centered_subset_abundances,
    feature_name,
    unique_sample_id,
    corrected_log2_abundance,
    method = "pearson",
    use = "pairwise.complete.obs"
  ) %>%
    # add partial correlations
    left_join(top_partial_corrs, by = c("item1" = "feature_1", "item2" = "feature_2"))
  
  if (class(facet_by) != "NULL") {
    
    feature_corrs <- feature_corrs %>%
      left_join(
        feature_subset_annot %>%
          select(item1 = feature_name, category1),
        by = "item1"
      ) %>%
      left_join(
        feature_subset_annot %>%
          select(item2 = feature_name, category2),
        by = "item2"
      )
  }
  
  ordered_feature_corrs <- feature_corrs %>%
    # order by hclust
    mutate(
      item1 = ordered(item1, levels = feature_subset_annot$feature_name),
      item2 = ordered(item2, levels = feature_subset_annot$feature_name),
      corr_or_partial = ifelse(item1 < item2, correlation, partial_corr)
    ) 
  
  corr_heatmap_grob <- ggplot(ordered_feature_corrs, aes(x = item1, y = item2)) +
    geom_tile(aes(fill = corr_or_partial)) +
    scale_fill_gradientn(
      "Pearson Correlation",
      colours = c("steelblue1", "steelblue2", "steelblue3", "steelblue4", "black", "yellow4", "yellow3", "yellow2", "yellow"),
      limits = c(-1,1)
    ) +
    #scale_fill_gradient2(
    #  "Pearson Correlation",
    #  low = "steelblue1", 
    #  mid = "black",
    #  high = "yellow",
    #  midpoint = 0
    #) +
    geom_point(data = feature_subset_annot, aes(x = feature_name, y = feature_name, shape = feature_aging_archetype, color = change_w_age), size = 4) +
    scale_shape_manual(
      "Aging archetype", 
      values = c(
        "lifespan" = 15,
        "chronological age" = 16,
        "fraction of life lived" = 17,
        "age x lifespan" = 13,
        "lifespan-remaining" = 18,
        "non-omics" = 8)
    ) +
    scale_color_manual(
      "Change with age",
      values = c(
        "up" = "magenta",
        "down" = "orchid4"
      )
    ) +
    # use pretty labels
    scale_x_discrete(breaks = feature_subset_annot$feature_name, labels = feature_subset_annot$feature_label) +
    scale_y_discrete(breaks = rev(feature_subset_annot$feature_name), labels = rev(feature_subset_annot$feature_label)) +
    #coord_equal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_blank(),
      strip.background = element_rect(fill = "gray30"),
      strip.text = element_text(color = "gray90")
    )
  
  if (class(facet_by) != "NULL") {
    corr_heatmap_grob <- corr_heatmap_grob + facet_grid(category2 ~ category1, scales = "free", space = "free")
  }
  
  corr_heatmap_grob
}
