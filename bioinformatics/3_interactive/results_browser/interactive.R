#' @param outdir Directory where output should be saved
#' @param overwrite Overwrite existing files?
read_gcp_modeling_results <- function (outdir = "/tmp", overwrite = FALSE) {
  
  checkmate::assertCharacter(outdir, len = 1)
  checkmate::assertDirectoryExists(outdir, access = "w")
  checkmate::assertLogical(overwrite, len = 1)
  
  do_dir <- file.path(outdir, "domics")
  if (!(dir.exists(do_dir))) {
    dir.create(do_dir)
  }
  
  domic_files <- tibble::tribble(
    ~ type, ~ object,
    "tomic", "tidy_features.Rds",
    "tomic_nopeptides", "tidy_features_no_peptides.Rds",
    "tomic_peptides", "tidy_features_peptides.Rds",
    "signif", "model_signif.Rds",
    "features", "features_with_design.Rds"
  ) %>%
    dplyr::mutate(
      local_path = file.path(do_dir, object),
      exists = file.exists(local_path)
    )
  
  if (all(domic_files$exists) && !overwrite) {
    return (domic_files %>% dplyr::select(type, object, local_path))
  }
  
  # GCP service token is saved with the package
  # <<TO DO - update auth>>
  available_objects <- googleCloudStorageR::gcs_list_objects(prefix = "domics")
  
  domic_files <- domic_files %>%
    dplyr::mutate(gcs_object_path = file.path("domics", object),
                  gcs_exists = gcs_object_path %in% available_objects$name)
  
  unsatisfiable_requirments <- domic_files %>%
    dplyr::filter(!exists & !gcs_exists)
  if (nrow(unsatisfiable_requirments) != 0) {
    stop (glue::glue(
      "{nrow(unsatisfiable_requirments)} required objects were missing from GCS:
      {paste(unsatisfiable_requirments$gcs_object_path, collapse = ', ')}"
    ))
  }
  
  if (overwrite) {
    assets_to_get <- domic_files
  } else {
    assets_to_get <- domic_files %>%
      dplyr::filter(!exists)
  }
  
  for (i in 1:nrow(assets_to_get)) {
    googleCloudStorageR::gcs_get_object(
      object = assets_to_get$gcs_object_path[i],
      saveToDisk = assets_to_get$local_path[i],
      overwrite = overwrite
    )
  }
  
  return (domic_files %>% dplyr::select(type, object, local_path))
}

stratify_by <- tibble::tribble(
  ~ stratify_var, ~ stratify_label,
  "fc_batch", "known biological batch effects",
  "Generation", "mouse generation",
  "set", "sample analysis sets"
)

stratify_grob <- function (grob, stratify_var, stratify_label, grouping_var_formula = NULL) {
  
  facet_rhs <- ifelse(class(grouping_var_formula) == "NULL", ".", all.vars(grouping_var_formula))
  updated_formula <- as.formula(paste(stratify_var, " ~ ", facet_rhs))
  
  updated_title <- glue::glue("{grob$labels$title} stratified by {stratify_label}")
  
  stratified_grob <- grob +
    facet_grid(updated_formula, scale = "free_y") +
    ggtitle(updated_title)
  
  return(stratified_grob)
}

plot_ddm <- function (
    features_with_design,
    model_signif,
    feature_names,
    facet_ncol = 3
    ) {
  
  # feature_names <- sample(model_signif$groupName, 6)
  
  ddm_signif <- model_signif %>%
    dplyr::filter(
      term == "ddm",
      groupName %in% feature_names
      )
  
  if (nrow(ddm_signif) == 0) {
    # o-link doesn't have ddm so it is possible not to have ddm results
    # also, all of the DDM samples could have missing values
    return(list(
      dat = tibble::tibble(), 
      plots = list(ddm = ggplot(data.frame(x = 0, y = 0), aes(x = x, y = y)) +
                     geom_text(label = "No results available", size = 15) +
                     theme(text = element_blank(), line = element_blank()))
    ))
  }
  
  ddm_data <- features_with_design %>%
    dplyr::filter(Age != "8 months") %>%
    dplyr::right_join(
      ddm_signif,
      by = c("groupId", "data_type", "data_modality", "feature_name"),
      relationship = "many-to-many"
      ) %>%
    dplyr::mutate(groupName = factor(groupName, levels = unique(feature_names)))
  
  if (length(unique(ddm_data$feature_name)) > length(feature_names)) {
    grouping_var_formula <- as.formula(~ feature_name)
  } else {
    grouping_var_formula <- as.formula(~ groupName)
  }
  
  ddm_plot <- ggplot(
    ddm_data,
    aes(x = is_ddm_sample, y = corrected_log2_fold_change)
  ) +
    geom_boxplot(outlier.shape = 1) +
    geom_hline(yintercept = 0, color = "RED") +
    ggbeeswarm::geom_beeswarm(size = 0.33) +
    ggtitle("Near death changes") +
    scale_x_discrete(
      "DDM Status (mouse dies <= 21 days after blood draw)",
      breaks = c("FALSE", "TRUE"),
      labels = c("non-DDM", "DDM")
      ) +
    scale_y_continuous(expression(log[2] ~ "fold change")) +
    facet_wrap(grouping_var_formula, scale = "free_y", ncol = facet_ncol) + 
    theme_minimal() + theme(
      text = element_text(size = 15),
      strip.background = element_rect(fill = "gray90"),
      panel.spacing = grid::unit(2, "lines")
      )
  
  return(list(
    dat = clean_plot_signif(ddm_signif), 
    plots = list(ddm = ddm_plot)
  ))
  
}

plot_aging <- function (
    features_with_design,
    model_signif,
    feature_names,
    facet_ncol = 3
  ) {
  
  age_signif <- model_signif %>%
    dplyr::filter(
      term %in% c("early age", "late age"),
      groupName %in% feature_names
      ) %>%
    dplyr::group_by(term, groupName, feature_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  stopifnot(nrow(age_signif %>% distinct(groupName)) == length(feature_names))
  
  age_data <- features_with_design %>%
    dplyr::filter(Age != "8 months") %>%
    dplyr::right_join(
      age_signif,
      by = c("groupId", "data_type", "data_modality", "feature_name"),
      relationship = "many-to-many"
      ) %>%
    dplyr::mutate(groupName = factor(groupName, levels = unique(feature_names))) %>%
    # drop undefined fold-changes
    dplyr::filter(!is.na(corrected_log2_fold_change))
  
  if (length(unique(age_data$feature_name)) > length(feature_names)) {
    grouping_var_formula <- as.formula(~ feature_name)
  } else {
    grouping_var_formula <- as.formula(~ groupName)
  }
  
  aging_plot <- ggplot(age_data, aes(x = Age, y = corrected_log2_fold_change, color = Age)) + 
    geom_boxplot(outlier.shape = 1) +
    geom_hline(yintercept = 0, color = "RED")  +
    ggbeeswarm::geom_beeswarm(size = 0.33) +
    ggtitle("Age-dependent molecules") +
    scale_y_continuous(expression(log[2] ~ "fold change")) +
    scale_color_manual("Chronological Age", values = c("8 months" = "goldenrod1", "14 months" = "darkorange2", "20 months" = "darkred")) +
    facet_wrap(grouping_var_formula, scale = "free_y", ncol = facet_ncol) + 
    theme_minimal() + theme(
      text = element_text(size = 15),
      strip.background = element_rect(fill = "gray90"),
      panel.spacing = grid::unit(2, "lines")
      )
  
  # stratify by generation, known batch effects and experimental sets  
  
  stratified_grobs <- purrr::map2(
    stratify_by$stratify_var,
    stratify_by$stratify_label,
    stratify_grob,
    grob = aging_plot,
    grouping_var_formula = grouping_var_formula
    )
  
  aging_strata <- (stratified_grobs[[1]] / stratified_grobs[[2]] / stratified_grobs[[3]]) +
    plot_layout(heights = c(3, 5, 3), guides = "collect")
  
  return(list(
    dat = clean_plot_signif(age_signif),
    plots = list(
      aging = aging_plot,
      aging_strata = aging_strata
    )
  ))
}

plot_sex <- function (
    features_with_design,
    model_signif,
    feature_names,
    facet_ncol = 3
    ) {
  
  female = intToUtf8(9792)
  male = intToUtf8(9794)
  
  sex_signif <- model_signif %>%
    dplyr::filter(
      term == "sex",
      groupName %in% feature_names
      )

  stopifnot(nrow(sex_signif %>% distinct(groupName)) == length(feature_names))
  
  sex_data <- features_with_design %>%
    dplyr::right_join(
      sex_signif,
      by = c("groupId", "data_type", "data_modality", "feature_name"),
      relationship = "many-to-many"
      ) %>%
    dplyr::mutate(groupName = factor(groupName, levels = unique(feature_names))) %>%
    dplyr::filter(!is_ddm_sample)
  
  if (length(unique(sex_data$feature_name)) > length(feature_names)) {
    grouping_var_formula <- as.formula(~ feature_name)
  } else {
    grouping_var_formula <- as.formula(~ groupName)
  }
  
  sex_plot <- ggplot(
    sex_data,
    aes(x = Sex, y = corrected_log2_abundance)
    ) +
    geom_boxplot(outlier.shape = 1) +
    ggbeeswarm::geom_beeswarm(size = 0.33) +
    ggtitle("Sex-dependent molecules") +
    scale_x_discrete("Sex", breaks = c("F", "M"), labels = c(female, male)) +
    scale_y_continuous(expression(log[2] ~ "abundance")) +
    facet_wrap(grouping_var_formula, scale = "free_y", ncol = facet_ncol) + 
    theme_minimal() + theme(
      text = element_text(size = 15),
      axis.text.x = element_text(size = 30),
      strip.background = element_rect(fill = "gray90"),
      panel.spacing = grid::unit(2, "lines")
      )
  
  return(list(
    dat = clean_plot_signif(sex_signif),
    plots = list(sex = sex_plot)
    ))
}

plot_lifespan_remaining <- function (
    features_with_design,
    model_signif,
    feature_names,
    facet_ncol = 3
    ) {
  
  lifespan_remaining_signif <- model_signif %>%
    dplyr::filter(
      model_name == "cross-sectional lm (lifespan-remaining), no DDM",
      term == 'lifespan remaining',
      groupName %in% feature_names
      )
  
  stopifnot(nrow(lifespan_remaining_signif) > 0)
  
  lifespan_remaining_data <- features_with_design %>%
    dplyr::right_join(
      lifespan_remaining_signif,
      by = c("groupId", "data_type", "data_modality", "feature_name"),
      relationship = "many-to-many"
      ) %>%
    dplyr::mutate(groupName = factor(groupName, levels = unique(feature_names))) %>%
    dplyr::filter(!is_ddm_sample)
  
  if (length(unique(lifespan_remaining_data$feature_name)) > length(feature_names)) {
    grouping_var_formula <- as.formula(~ feature_name)
  } else {
    grouping_var_formula <- as.formula(~ groupName)
  }
  
  lifespan_remaining_plot <- ggplot(
    lifespan_remaining_data, 
    aes(x = lifespan_remaining, y = corrected_log2_abundance)) + 
    geom_point(aes(color = Age)) +
    geom_smooth(aes(color = Age, fill = Age), method = "lm", formula = y ~ x, alpha = 0.25) +
    geom_smooth(method = "lm", formula = y ~ x) +
    ggtitle("Lifespan remaining-dependent molecules") +
    scale_x_continuous("Lifespan remaining (days)") +
    scale_y_continuous(expression(log[2] ~ "abundance")) +
    scale_color_manual("Chronological Age", values = c("8 months" = "goldenrod1", "14 months" = "darkorange2", "20 months" = "darkred")) +
    scale_fill_manual("Chronological Age", values = c("8 months" = "goldenrod1", "14 months" = "darkorange2", "20 months" = "darkred")) +
    facet_wrap(grouping_var_formula, scale = "free_y", ncol = facet_ncol) + 
    theme_minimal() + theme(
      text = element_text(size = 15),
      strip.background = element_rect(fill = "gray90"),
      panel.spacing = grid::unit(2, "lines")
    )
  
  stratified_grobs <- purrr::map2(
    stratify_by$stratify_var,
    stratify_by$stratify_label,
    stratify_grob,
    grob = lifespan_remaining_plot,
    grouping_var_formula = grouping_var_formula
  )
  
  lifespan_remaining_strata <- (stratified_grobs[[1]] / stratified_grobs[[2]] / stratified_grobs[[3]]) +
    plot_layout(heights = c(3, 5, 3), guides = "collect")
  
  return(list(
    dat = clean_plot_signif(lifespan_remaining_signif),
    plots = list(
      lifespan_remaining = lifespan_remaining_plot,
      lifespan_remaining_strata = lifespan_remaining_strata
    )))
}

plot_fraction_of_life_lived <- function (
    features_with_design,
    model_signif,
    feature_names,
    facet_ncol = 3
    ) {
  
  ffl_signif <- model_signif %>%
    dplyr::filter(
      term == "fraction of life lived",
      groupName %in% feature_names
    )
  
  stopifnot(nrow(ffl_signif %>% distinct(groupName)) == length(feature_names))
  
  ffl_data <- features_with_design %>%
    dplyr::right_join(
      ffl_signif,
      by = c("groupId", "data_type", "data_modality", "feature_name"),
      relationship = "many-to-many"
      ) %>%
    dplyr::mutate(groupName = factor(groupName, levels = unique(feature_names))) %>%
    dplyr::filter(!is_ddm_sample)
  
  if (length(unique(ffl_data$feature_name)) > length(feature_names)) {
    grouping_var_formula <- as.formula(~ feature_name)
  } else {
    grouping_var_formula <- as.formula(~ groupName)
  }
  
  ffl_plot <- ggplot(ffl_data,
    aes(x = fraction_of_life_lived, y = corrected_log2_abundance)) +
    geom_point(aes(color = Age)) +
    geom_smooth(aes(color = Age, fill = Age), method = "lm", formula = y ~ x, alpha = 0.25) +
    geom_smooth(method = "lm", formula = y ~ x) +
    ggtitle("Fraction of lifespan lived-dependent molecules") +
    scale_x_continuous("Fraction of life lived") +
    scale_y_continuous(expression(log[2] ~ "abundance")) +
    scale_color_manual("Chronological Age", values = c("8 months" = "goldenrod1", "14 months" = "darkorange2", "20 months" = "darkred")) +
    scale_fill_manual("Chronological Age", values = c("8 months" = "goldenrod1", "14 months" = "darkorange2", "20 months" = "darkred")) +
    facet_wrap(grouping_var_formula, scale = "free_y", ncol = facet_ncol) + 
    theme_minimal() + theme(
      text = element_text(size = 15),
      strip.background = element_rect(fill = "gray90"),
      panel.spacing = grid::unit(2, "lines")
    )
  
  stratified_grobs <- purrr::map2(
    stratify_by$stratify_var,
    stratify_by$stratify_label,
    stratify_grob,
    grob = ffl_plot,
    grouping_var_formula = grouping_var_formula
  )
  
  fraction_of_life_lived_strata <- (stratified_grobs[[1]] / stratified_grobs[[2]] / stratified_grobs[[3]]) +
    plot_layout(heights = c(3, 5, 3), guides = "collect")
  
  return(list(
    dat = clean_plot_signif(ffl_signif),
    plots = list(
      fraction_of_life_lived = ffl_plot,
      fraction_of_life_lived_strata = fraction_of_life_lived_strata
    )))
}

plot_lifespan <- function (
    features_with_design,
    model_signif,
    feature_names,
    facet_ncol = 3
    ) {
  
  lifespan_effects <- model_signif %>%
    dplyr::filter(
      term == 'lifespan',
      groupName %in% feature_names
      )
  
  stopifnot(nrow(lifespan_effects %>% distinct(groupName)) == length(feature_names))
  
  top_lifespan_data <- features_with_design %>%
    dplyr::right_join(
      lifespan_effects,
      by = c("groupId", "data_type", "data_modality", "feature_name"),
      relationship = "many-to-many"
      ) %>%
    dplyr::mutate(groupName = factor(groupName, levels = unique(feature_names))) %>%
    dplyr::filter(!is_ddm_sample)
  
  if (length(unique(top_lifespan_data$feature_name)) > length(feature_names)) {
    grouping_var_formula <- as.formula(~ feature_name)
  } else {
    grouping_var_formula <- as.formula(~ groupName)
  }
  
  lifespan_plot <- ggplot(top_lifespan_data, aes(x = lifespan, y = corrected_log2_abundance)) + 
    geom_point(aes(color = Age)) +
    geom_smooth(aes(group = groupId), method = "lm", se = TRUE, formula = "y ~ x", alpha = 0.25) +
    ggtitle("Lifespan-dependent molecules, controlling for ages") +
    scale_x_continuous("Lifespan (days)") +
    scale_y_continuous(expression(log[2] ~ "abundance")) +
    scale_color_manual("Chronological Age", values = c("8 months" = "goldenrod1", "14 months" = "darkorange2", "20 months" = "darkred")) +
    facet_wrap(grouping_var_formula, scale = "free_y", ncol = facet_ncol) + 
    theme_minimal() + theme(text = element_text(size = 15), strip.background = element_rect(fill = "gray90"), panel.spacing = grid::unit(2, "lines"))
  
  stratified_grobs <- purrr::map2(
    stratify_by$stratify_var,
    stratify_by$stratify_label,
    stratify_grob,
    grob = lifespan_plot,
    grouping_var_formula = grouping_var_formula
  )
  
  lifespan_strata <- (stratified_grobs[[1]] / stratified_grobs[[2]] / stratified_grobs[[3]]) +
    plot_layout(heights = c(3, 5, 3), guides = "collect")
  
  return(list(
    dat = clean_plot_signif(lifespan_effects),
    plots = list(
      lifespan = lifespan_plot,
      lifespan_strata = lifespan_strata
    )))
}


plot_age_by_lifespan <- function (
    features_with_design,
    model_signif,
    feature_names,
    facet_ncol = 3
    ) {
  
  age_by_lifespan_effects <- model_signif %>%
    dplyr::filter(
      term %in% c('early age x lifespan', 'late age x lifespan'),
      groupName %in% feature_names
    )
  
  stopifnot(nrow(age_by_lifespan_effects) > 0)
  
  age_by_lifespan_data <- features_with_design %>%
    dplyr::filter(Age != "8 months", !is_ddm_sample) %>%
    dplyr::right_join(
      age_by_lifespan_effects,
      by = c("groupId", "data_type", "data_modality", "feature_name"),
      relationship = "many-to-many"
      ) %>%
    dplyr::mutate(groupName = factor(groupName, levels = unique(feature_names))) %>%
    dplyr::filter(!is_ddm_sample) %>%
    dplyr::filter(!is.na(log2_fold_change))
  
  if (length(unique(age_by_lifespan_data$feature_name)) > length(feature_names)) {
    grouping_var_formula <- as.formula(~ feature_name)
  } else {
    grouping_var_formula <- as.formula(~ groupName)
  }
  
  age_by_lifespan_plot <- age_by_lifespan_data %>%
    ggplot(aes(x = lifespan, y = corrected_log2_fold_change)) + 
    geom_point(aes(color = Age)) +
    geom_smooth(aes(group = Age, color = Age, fill = Age), method = "lm", se = TRUE, formula = "y ~ x", alpha = 0.25) +
    ggtitle("Age x lifespan dependent molecules") +
    scale_x_continuous("Lifespan (days)") +
    scale_y_continuous(expression(log[2] ~ "fold change")) +
    scale_color_manual("Chronological Age", values = c("8 months" = "goldenrod1", "14 months" = "darkorange2", "20 months" = "darkred")) +
    scale_fill_manual("Chronological Age", values = c("8 months" = "goldenrod1", "14 months" = "darkorange2", "20 months" = "darkred")) +
    facet_wrap(grouping_var_formula, scale = "free_y", ncol = facet_ncol) + 
    theme_minimal() + theme(
      text = element_text(size = 15),
      strip.background = element_rect(fill = "gray90"),
      panel.spacing = grid::unit(2, "lines")
      )
  
  stratified_grobs <- purrr::map2(
    stratify_by$stratify_var,
    stratify_by$stratify_label,
    stratify_grob,
    grob = age_by_lifespan_plot,
    grouping_var_formula = grouping_var_formula
  )
  
  age_by_lifespan_strata <- (stratified_grobs[[1]] / stratified_grobs[[2]] / stratified_grobs[[3]]) +
    plot_layout(heights = c(3, 5, 3), guides = "collect")
  
  return(list(
    dat = clean_plot_signif(age_by_lifespan_effects),
    plots = list(
      age_by_lifespan = age_by_lifespan_plot,
      age_by_lifespan_strata = age_by_lifespan_strata
    )))
}

clean_plot_signif <- function (plot_signif) {
  
  # identify groupNames which match to multiple features
  degenerate_groups_present <- plot_signif %>%
    distinct(groupName, feature_name) %>%
    count(groupName) %>%
    {any(.$n > 1)}
  
  if (degenerate_groups_present) {
    features_defined_by = c("groupName", "feature_name")
  } else {
    features_defined_by <- "groupName"
  }

  plot_signif %>%
    dplyr::select(term, !!!rlang::syms(features_defined_by), data_modality, estimate:q.value) %>%
    dplyr::mutate(
      estimate = round(estimate, 5),
      std.error = round(std.error, 5),
      statistic = round(statistic, 3)
    )
}

#' Shiny ggplot Test
#'
#' Test the shiny ggplot module as a stand-alone application.
#'
#' @inheritParams tomic_to
#'
#' @returns A \code{shiny} app
#'
#' @examples
#'
#' plot_list <- plot_fraction_of_lifespan_lived(features_with_design, model_signif, feature_names)
#' #plot_list <- plot_aging(features_with_design, model_signif, feature_names)
#' #plot_list <- plot_lifespan_remaining(features_with_design, model_signif, feature_names)
#' test_shiny_signif(plot_list)
#'
#' @export
test_shiny_signif <- function(plot_list) {
  checkmate::assertList(plot_list)
  stopifnot(names(plot_list) == c("dat", "plots"))
  
  shiny::shinyApp(
    ui = shiny::fluidPage(
      
      # Sidebar with a slider input for the number of bins
      shiny::verticalLayout(
        signifOutput("signif_out")
      )
    ),
    server = function(input, output, session) {
      shiny::observe({
        signifServer("signif_out", plot_list)
      })
    }
  )
}

#' Significance Output
#'
#' UI components for the significance module.
#'
#' @inheritParams shiny::moduleServer
#'
#' @returns A \code{shiny} UI
#'
#' @export
signifOutput <- function(
  id
) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::verticalLayout(
      shiny::h3("Significance"),
      shiny::dataTableOutput(
        ns("data")
      ),
      shiny::downloadButton(ns("download_data"), "Download .tsv"),
      shiny::hr(),
      shiny::h3("Visualize Effects"),
      shiny::plotOutput(
        ns("ggplot1"),
        height = 800
      ),
      romic::plotsaverInput(ns("ggsave1"), ui_format = "wide"),
      shiny::hr(),
      shiny::uiOutput(ns("ggplot2_ui"))
    )
  )
}

#' Significance Server
#'
#' Server components for the significance module.
#'
#' @inheritParams shiny::moduleServer
#' @param plot_list List containing dat and plots
#'
#' @returns Nothing; used for side-effects.
#'
#' @export
signifServer <- function(id, plot_list) {
  
  checkmate::assertList(plot_list)
  stopifnot(names(plot_list) == c("dat", "plots"))
  
  n_plots <- length(plot_list$plots)
  plot1 <- plot_list$plots[[1]]
  checkmate::assertClass(plot1, "ggplot")
  
  shiny::moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      
      shiny::observe({
        ns <- session$ns
        
        # table
        output$data <- shiny::renderDataTable(plot_list$dat)
        plot_basename <- names(plot_list$plots)[1]
        
        output$download_data <- shiny::downloadHandler(
          filename = glue::glue("{plot_basename}_significance.tsv"),
          content = function(file) {
            vroom::vroom_write(plot_list$dat, file)
          }
        )
        
        # plots
        output$ggplot1 <- shiny::renderPlot(plot1)
        romic::plotsaverServer("ggsave1", plot1, filename = paste0(plot_basename, ".png"))
        
        # display a second plot if it exists 
        if (n_plots >= 2) {
          plot2 <- plot_list$plots[[2]]
          checkmate::assertClass(plot2, "ggplot")
          output$ggplot2 <- shiny::renderPlot(plot2)
          
          output$ggplot2_ui <- shiny::renderUI({
            shiny::verticalLayout(
              shiny::h3("Effect Consistency Across Biological and Technical Batches"),
              shiny::plotOutput(
                ns("ggplot2"),
                height = 2000
                )
              )
            })
          }
        })
    })
  }

