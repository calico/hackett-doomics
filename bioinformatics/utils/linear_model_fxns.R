# define models used in study
do_linear_models <- tibble::tribble(
  ~ model_name, ~ formula, ~ filter, ~ contrasts, ~ no_ddm, ~ n.svs, ~ term,

  "fold-change lm (age), no DDM", "log2_fold_change ~ Age + fc_batch + 0",
  quo(Age != "8 months"), list(Age = "contr.treatment"), TRUE, 0, c("early age", "late age", "late blood draw", "early G8"),
  
  "fold-change lm (age x sex), no DDM", "log2_fold_change ~ Age + Age:Sex + fc_batch + 0",
  quo(Age != "8 months"), list(Age = "contr.treatment"), TRUE, 0, c("early age x male", "late age x male"),
  
  "fold-change lm (age x lifespan), no DDM", "log2_fold_change ~ Age + Age:lifespan + fc_batch + 0",
  quo(Age != "8 months"), list(Age = "contr.treatment"), TRUE, 0, c("early age x lifespan", "late age x lifespan"),
  
  "fold-change lm (DDM)", "log2_fold_change ~ Age + is_ddm_sample + fc_batch + 0",
  quo(Age != "8 months"), list(Age = "contr.treatment"), FALSE, 0, "ddm",
  
  # age -> age x lifespan
  "cross-sectional lm (age), no DDM", "normalized_log2_abundance ~ Age + Sex + xs_batch", 
  NA, list(Age = "contr.helmert"), TRUE, 0, NA,
  
  "cross-sectional lm (int-age), no DDM", "normalized_log2_abundance ~ blood_draw_age_actual + Sex + xs_batch", 
  NA, NULL, TRUE, 0, NA,
  
  "cross-sectional lm (age + lifespan), no DDM", "normalized_log2_abundance ~ Age + lifespan + Sex + xs_batch", 
  NA, list(Age = "contr.helmert"), TRUE, 0, c("lifespan", "sex"),
  
  "cross-sectional lm (int-age + lifespan), no DDM", "normalized_log2_abundance ~ blood_draw_age_actual + lifespan + Sex + xs_batch", 
  NA, NULL, TRUE, 0, NA,
  
  "cross-sectional lm (age x lifespan), no DDM", "normalized_log2_abundance ~ Age*lifespan + Sex + xs_batch", 
  NA, list(Age = "contr.helmert"), TRUE, 0, NA,
  
  "cross-sectional lm (int-age x lifespan), no DDM", "normalized_log2_abundance ~ blood_draw_age_actual * lifespan + Sex + xs_batch", 
  NA, NULL, TRUE, 0, NA,
  
  # ffl -> ffl x lifespan
  "cross-sectional lm (fll), no DDM", "normalized_log2_abundance ~ fraction_of_life_lived + Sex + xs_batch",
  NA, NULL, TRUE, 0, "fraction of life lived",
  
  "cross-sectional lm (fll + lifespan), no DDM", "normalized_log2_abundance ~ fraction_of_life_lived + lifespan + Sex + xs_batch",
  NA, NULL, TRUE, 0, NA, 
  
  "cross-sectional lm (fll x lifespan), no DDM", "normalized_log2_abundance ~ fraction_of_life_lived*lifespan + Sex + xs_batch",
  NA, NULL, TRUE, 0, NA, 
  
  # lifespan
  "cross-sectional lm (lifespan), no DDM", "normalized_log2_abundance ~ lifespan + Sex + xs_batch", 
  NA, NULL, TRUE, 0, NA,
  
  # lifespan remaining
  "cross-sectional lm (lifespan-remaining), no DDM", "normalized_log2_abundance ~ lifespan_remaining + Sex + xs_batch",
  NA, NULL, TRUE, 0, "lifespan remaining",
  
  "cross-sectional lm (lifespan-remaining + fll), no DDM", "normalized_log2_abundance ~ lifespan_remaining + fraction_of_life_lived + Sex + xs_batch",
  NA, NULL, TRUE, 0, NA,
  
  # f(sex)
  "cross-sectional lm (fll x sex), no DDM", "normalized_log2_abundance ~ fraction_of_life_lived*Sex + xs_batch",
  NA, NULL, TRUE, 0, "fraction of life lived x sex",
  
  # f(age, ffl)
  "cross-sectional lm (fll + age), no DDM", "normalized_log2_abundance ~ fraction_of_life_lived + Age + Sex + xs_batch",
  NA, list(Age = "contr.helmert"), TRUE, 0, NA
  ) %>%
  tidyr::nest(model_params = -model_name)

# filter model data to only include molecular entities (groupIds)
# with sufficient observations (e.g. more than one level, can compare to Age8months / SexMale) to fit a regression
filter_data <- function(model_data, model_params){
  
  n_levels <- function(var, data_group){
    n_levels <- data_group %>% 
      dplyr::filter(!is.na(!!rlang::sym(var))) %>%
      dplyr::select(!!var) %>%
      dplyr::distinct() %>%
      dplyr::summarise(n = dplyr::n())
    unlist(n_levels)
  }
  
  sufficient_levels <- function(data_group, vars){
    all(purrr::map(vars, n_levels, data_group) >= 2)
  }
  
  formula_vars = all.vars(as.formula(as.character(model_params$formula)))
  dependent_var = formula_vars[1]
  
  model_data.list <- model_data %>%
    dplyr::filter(!is.na(!!rlang::sym(dependent_var))) %>%
    dplyr::group_split(groupId)
  
  can_model <- model_data.list %>%
    purrr::map_lgl(sufficient_levels, formula_vars)
  
  dplyr::bind_rows(model_data.list[can_model])
}

estimate_surrogate_vars <- function(model_data, model_params) {
  
  sample_metadata <- model_data %>%
    dplyr::distinct(sampleId, set, unique_batch, sample, Mouse.ID, Age, Sex,
                    lifespan, lifespan_remaining, is_ddm_sample)
  
  
  dependent_var <- all.vars(as.formula(model_params$formula))[1]
  
  Y <- model_data %>%
    reshape2::acast(groupId ~ sampleId, value.var = dependent_var)
  
  ordered_sample_metadata <- tibble::tibble(sampleId = colnames(Y)) %>%
    dplyr::left_join(sample_metadata, by = "sampleId")
  
  one_sided_formula <- stringr::str_split(model_params$formula, "~")[[1]][-1] %>%
    {paste0("~", .)} %>%
    as.formula()
  
  model_matrix = model.matrix(one_sided_formula, contrasts = model_params$contrasts[[1]], data = ordered_sample_metadata)
  
  sva_model = sva::sva(Y, model_matrix, n.sv = model_params$n.svs)
  
  surrogate_vars <- sva_model$sv
  rownames(surrogate_vars) <- ordered_sample_metadata$sampleId
  colnames(surrogate_vars) <- paste0("S", 1:ncol(surrogate_vars))
  
  surrogate_vars <- as.data.frame(surrogate_vars) %>%
    dplyr::mutate(sampleId = rownames(.)) %>%
    dplyr::tbl_df()
  
  return(surrogate_vars)
}

bootstrap_hypothesis_test <- function(lm_fit, model_params, data_group,
                                      weights, n_bootstraps) {
  
  # remove unused levels (e.g. 8 month measurements were filtered)
  # so they will not be included in the design matrix below
  data_group <- data_group %>% 
    purrr::modify_if(is.factor, factor)
  
  formula <- as.formula(model_params$formula)
  design_matrix <- model.matrix(formula, data = data_group, contrasts.arg = model_params$contrasts[[1]])
  response_var <- all.vars(formula)[1]
  n_samples <- nrow(design_matrix)
  
  # residuals of model fit adjusted for fitted degrees of freedom
  
  inflated_residuals = residuals(lm_fit) * sqrt(n_samples / (n_samples - length(coef(lm_fit))))
  prediction = predict(lm_fit)
  
  # generate bootstrapped response variables
  
  bootstrap_residuals <- matrix(sample(inflated_residuals, size = n_bootstraps * n_samples, replace = TRUE),
                                nrow = n_samples, ncol = n_bootstraps)
  bootstrap_predictions <- bootstrap_residuals + prediction
  
  # linear algebra implementation of OLS linear regression
  
  lsq_fit <- if(is.null(weights)) { 
    solve(t(design_matrix) %*% design_matrix) %*% t(design_matrix) %*% bootstrap_predictions
  } else {
    w <- weights[!is.na(data_group[ ,response_var])]
    w <- diag(w)
    solve(t(design_matrix) %*% w %*% design_matrix) %*% t(design_matrix) %*% w %*% bootstrap_predictions
  }
  
  # apply a two-sided test of coefficient significance:  beta != 0
  
  bootstrap_pvalues <- apply(lsq_fit, 1, function(estimate) {
    p_value <- 1 - 2*abs(0.5 - sum(estimate > 0)/n_bootstraps) + 1/n_bootstraps
    p_value <- pmin(p_value, 1)
  })
  
  bootstrap_pvalues
  
}  

lm_fitter_group <- function(data_group, model_params, surrogate_vars, n_bootstraps) {
  
  if (model_params$n.svs > 0) {
    # update formula to include surrogate variables
    
    surrogate_var_formula <- paste(setdiff(colnames(surrogate_vars), "sampleId"), collapse = " + ")
    model_params$formula <- paste0(model_params$formula, " + ", surrogate_var_formula)
    
    data_group <- data_group %>%
      dplyr::left_join(surrogate_vars, by = "sampleId")
  }
  
  # CompMS proteomics data provides variance values for each protein over it's peptides.
  # This values are useful for weight each protein observations by it's certainty in
  # a linear regression model.
  
  weights <- if(all(is.na(data_group$weight))){
    # no weights, non-compms dataset
    NULL
  } else if(all(!is.na(data_group$weight))){
    # available weights, compms dataset
    data_group$weight
  } else {
    # weights available for only some measurements, compms data is improperly formatted
    stop("variance values available for only some measurements")
  }
  
  lm_fit <- lm(
    data = data_group,
    formula = as.formula(model_params$formula),
    contrasts = model_params$contrasts[[1]],
    weights = weights
    )
  
  map_term <- function(term){
    dplyr::case_when(
      term == "SexM" ~ "sex",
      term == "is_ddm_sampleTRUE" ~ "ddm",
      term == "Age1" ~ "early age",
      term == "Age14 months" ~ "early age",
      term == "Age2" ~ "late age",
      term == "Age20 months" ~ "late age",
      term == "Age14 months:lifespan" ~ "early age x lifespan",
      term == "Age20 months:lifespan" ~ "late age x lifespan",
      term == "lifespan:Age1" ~ "early age x lifespan",
      term == "lifespan:Age2" ~ "late age x lifespan",
      term == "Age1:lifespan" ~ "early age x lifespan",
      term == "Age2:lifespan" ~ "late age x lifespan",
      term == "Age14 months:SexM" ~ "early age x male",
      term == "Age20 months:SexM" ~ "late age x male",
      term == "lifespan_remaining" ~ "lifespan remaining",
      term == "xs_batchLBD" ~ "late blood draw",
      term == "fc_batchLBD" ~ "late blood draw",
      term == "xs_batchEG8" ~ "early G8",
      term == "fc_batchEG8" ~ "early G8",
      term == "fraction_of_life_lived" ~ "fraction of life lived",
      term == "fraction_of_life_lived:lifespan" ~ "fraction of life lived x lifespan",
      term == "fraction_of_life_lived:SexM" ~ "fraction of life lived x sex",
      term == "blood_draw_age_actual" ~ "age",
      term == "blood_draw_age_actual:lifespan" ~ "age x lifespan",
      TRUE ~ term
      )
  }
  
  params <- broom::tidy(lm_fit) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(term = map_term(term))
  
  # Fit linear regression to a regression model predictions and bootstrapped 
  # residuals to form an empirical null distribution of p-values
  if(n_bootstraps > 0){
    bootstrap_pvalues <- bootstrap_hypothesis_test(lm_fit, model_params, data_group, weights, n_bootstraps)
    names(bootstrap_pvalues) <- map_term(names(bootstrap_pvalues))
    params <- params %>%
      dplyr::rename(pvalue_ols = p.value) %>%
      dplyr::mutate(pvalue_bs = bootstrap_pvalues[term])
  }
  
  augmented_df <- broom::augment(lm_fit)
  desired_augment_params <- c(".fitted", ".resid", ".std.resid")
  # raw residuals may be missing for some models on some R versions
  augmented_df <- augmented_df %>%
    dplyr::select(!!!rlang::syms(intersect(desired_augment_params, colnames(augmented_df))))
  
  fits <- data_group %>%
    dplyr::select(
      sampleId, set, unique_batch, sample, Mouse.ID, Age, Sex,
      blood_draw_age_actual, lifespan, lifespan_remaining,
      fraction_of_life_lived, is_ddm_sample
      ) %>%
    dplyr::bind_cols(augmented_df)

  model_summary <- broom::glance(lm_fit)
  
  list(params = params, fits = fits, model_summary = model_summary)
}

lm_fitter <- function(model_data, model_params, n_bootstraps = 100000) {
  checkmate::assertLogical(model_params$no_ddm, len = 1)
  
  if (model_params$no_ddm) {
    model_data <- model_data %>%
      dplyr::filter(!is_ddm_sample)
  }
  
  if (!is.na(model_params$filter)) {
    # filter with a quosure if there is one
    model_data <- model_data %>%
      dplyr::filter(!!!model_params$filter)
  }
  
  if (model_params$n.svs != 0) {
    surrogate_vars <- estimate_surrogate_vars(model_data, model_params)
  } else {
    surrogate_vars <- tibble(sampleId = NA_character_) %>% dplyr::slice(-1)
  }
  
  filtered_data <- filter_data(model_data, model_params)
  
  # some datasets / data_types will not have any values
  # for the relevant model e.g. olink data was only collected on set3,
  # which did not include any ddm samples
  if(nrow(filtered_data) == 0){
    return(NULL)
  }
  
  filtered_data %>%
    tidyr::nest(data_group = -groupId) %>%
    dplyr::mutate(fit_group = purrr::map(
      data_group,
      lm_fitter_group,
      model_params,
      surrogate_vars,
      n_bootstraps
      )) %>%
    dplyr::select(-data_group) %>%
    tidyr::unnest_wider(fit_group)
}
