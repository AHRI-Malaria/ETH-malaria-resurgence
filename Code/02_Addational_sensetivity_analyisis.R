################################################################################
#  ADDITIONAL SENSITIVITY 

# ------------------------------------------------------------------------------
# A. High species-completeness subset of the all-HF data:
#      >=80% of months observed for BOTH species in EACH calendar year.
#    This directly addresses concern that primary 934-woreda missingness is
#    geographically/non-randomly distributed.
#
# B. Elevation omitted from the 934-woreda primary model.
#    This assesses dependence on completed/imputed population-weighted elevation.
#
# C. Exclude woredas whose elevation was newly extension-imputed (75 expected).
#    This retains the wider historical elevation lookup but removes the new
#    75-woreda extension step.
#
# D. IDP functional-form/context sensitivity:
#      - log1p(IDP) instead of linear standardised IDP
#      - returning IDP added
#      - conflict events added
#      - conflict-related fatalities added
#
# E. An. stephensi coding:
#      woreda-only detection instead of woreda-or-neighbour detection.
#
# F. Intervention targeting/resistance:
#      add between-woreda mean LLIN access, mean IRS coverage, and zonal
#      deltamethrin resistance to the within-woreda intervention terms.
#
# G. Exploratory surveillance indicators:
#      HRP2 deletion and ACT-resistance indicators (Pf only).
#
################################################################################


# ==============================================================================
# 0. PACKAGES
# ==============================================================================


library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(writexl)

set.seed(20260830)


# ==============================================================================
# 1. PATHS AND SWITCHES
# ==============================================================================

BASE_DIR <- "result_species_only_reanalysis_934_primary"

MODEL_DIR <- file.path(
  BASE_DIR,
  "models"
)

PREP_DIR <- file.path(
  BASE_DIR,
  "prepared_data"
)

TABLE_DIR <- file.path(
  BASE_DIR,
  "tables"
)

SENS_DIR <- file.path(
  BASE_DIR,
  "additional_sensitivity"
)

dir.create(
  SENS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

PRIMARY_DATA_PATH <- file.path(
  PREP_DIR,
  "PRIMARY_allHF_noTigray_934_model_data.rds"
)

PRIMARY_GRAPH_PATH <- file.path(
  PREP_DIR,
  "woreda_all_HF_no_Tigray_primary.adj"
)

SHAPE_CANDIDATES <- c(
  "data/shp/eth_adm3_clean_id.shp",
  "eth_adm3_clean_id.shp"
)


RUN_HIGH_COMPLETENESS_80 <- TRUE
RUN_NO_ELEVATION <- TRUE
RUN_EXCLUDE_NEW_ELEVATION_EXTENSION <- TRUE

RUN_LOG_IDP <- TRUE
RUN_RETURNING_IDP <- TRUE
RUN_CONFLICT <- TRUE
RUN_FATALITIES <- TRUE

RUN_STEPHENSI_WOREDA_ONLY <- TRUE
RUN_PROGRAMME_TARGETING_RESISTANCE <- TRUE

# These two are explicitly exploratory because surveillance indicators may
# represent surveillance geography as much as biological prevalence/resistance.
RUN_EXPLORATORY_HRP2_ACT <- TRUE


# ==============================================================================
# 2. HELPERS
# ==============================================================================

first_existing <- function(
    paths,
    label
) {

  hit <- paths[
    file.exists(paths)
  ]

  if (length(hit) == 0) {
    stop(
      "Could not find ",
      label,
      ". Checked:\n",
      paste(paths, collapse = "\n")
    )
  }

  hit[1]
}


safe_scalar <- function(
    x,
    element,
    nested
) {

  z <- tryCatch(
    x[[element]],
    error = function(e) NULL
  )

  if (is.null(z)) {
    return(NA_real_)
  }

  # Standard INLA structure: e.g. model[["waic"]][["waic"]].
  if (is.list(z)) {

    zz <- tryCatch(
      z[[nested]],
      error = function(e) NULL
    )

    if (!is.null(zz) &&
        length(zz) > 0) {
      return(
        as.numeric(
          zz[1]
        )
      )
    }
  }

  # Some INLA versions/components can be named atomic vectors.
  if (is.atomic(z)) {

    if (!is.null(names(z)) &&
        nested %in% names(z)) {
      return(
        as.numeric(
          z[[nested]]
        )
      )
    }

    if (length(z) == 1) {
      return(
        as.numeric(
          z[1]
        )
      )
    }
  }

  NA_real_
}


extract_irrs <- function(
    model,
    label
) {

  if (!inherits(model, "inla")) {
    stop(
      label,
      ": not a valid INLA model."
    )
  }

  model$summary.fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      "term"
    ) %>%
    transmute(
      model = label,
      term,
      beta = mean,
      beta_sd = sd,
      beta_lcl = `0.025quant`,
      beta_ucl = `0.975quant`,
      IRR = exp(mean),
      IRR_lcl = exp(`0.025quant`),
      IRR_ucl = exp(`0.975quant`),
      `IRR (95% CrI)` =
        sprintf(
          "%.3f (%.3f–%.3f)",
          IRR,
          IRR_lcl,
          IRR_ucl
        )
    )
}


extract_fit <- function(
    model,
    label
) {

  if (!inherits(model, "inla")) {
    stop(
      label,
      ": not a valid INLA model."
    )
  }

  tibble(
    model = label,
    WAIC = safe_scalar(
      model,
      "waic",
      "waic"
    ),
    effective_parameters = safe_scalar(
      model,
      "waic",
      "p.eff"
    ),
    DIC = safe_scalar(
      model,
      "dic",
      "dic"
    )
  )
}


z_primary <- function(x) {

  s <- sd(
    x,
    na.rm = TRUE
  )

  if (!is.finite(s) ||
      s <= 0) {
    stop(
      "Cannot standardise sensitivity covariate."
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) /
    s
}


PRIMARY_FIXED <- c(
  "rainfall_lag2_std",
  "max_temp_lag2_std",
  "min_temp_night_lag2_std",
  "evi_lag2_std",
  "elev_pop_weighted_std",
  "idp_ind_std",
  "stephensi_in_neighboring_woredas",
  "llin_within_std",
  "irs_within_std"
)


make_formula <- function(
    fixed_terms,
    graph_obj
) {

  rhs_fixed <- paste(
    fixed_terms,
    collapse = " + "
  )

  rhs_random <- paste(
    c(
      "region:time_linear",

      paste0(
        "f(woreda_idx, model='bym2', graph=graph_obj, ",
        "scale.model=TRUE, constr=TRUE, adjust.for.con.comp=TRUE)"
      ),

      paste0(
        "f(time_month, model='rw1', ",
        "hyper=list(prec=list(prior='pc.prec', param=c(0.3, 0.01))))"
      ),

      paste0(
        "f(season_month, model='rw1', cyclic=TRUE, ",
        "hyper=list(prec=list(prior='pc.prec', param=c(0.2, 0.01))))"
      ),

      paste0(
        "f(st_id, model='iid', ",
        "hyper=list(prec=list(prior='pc.prec', param=c(1, 0.01))))"
      )
    ),
    collapse = " + "
  )

  fml <- as.formula(
    paste0(
      ".outcome ~ ",
      rhs_fixed,
      " + ",
      rhs_random,
      " + offset(log(woreda_population))"
    )
  )

  environment(fml) <- environment()

  fml
}


fit_sensitivity <- function(
    d,
    outcome,
    graph_obj,
    model_name,
    label,
    fixed_terms = PRIMARY_FIXED
) {

  summary_path <- file.path(
    SENS_DIR,
    paste0(
      model_name,
      "_summary.rds"
    )
  )

  if (file.exists(summary_path)) {
    message(
      "Loading completed sensitivity summary: ",
      model_name
    )
    return(
      readRDS(
        summary_path
      )
    )
  }

  message("\n============================================================")
  message("Sensitivity model: ", label)
  message("============================================================")

  model_dat <- d %>%
    mutate(
      .outcome =
        as.numeric(
          .data[[outcome]]
        )
    )

  non_integer <- is.finite(
    model_dat$.outcome
  ) &
    abs(
      model_dat$.outcome -
        round(
          model_dat$.outcome
        )
    ) >
      1e-8

  if (any(non_integer)) {
    model_dat$.outcome[
      non_integer
    ] <- round(
      model_dat$.outcome[
        non_integer
      ]
    )
  }

  if (any(
    model_dat$.outcome < 0,
    na.rm = TRUE
  )) {
    stop(
      "Negative outcome in ",
      model_name
    )
  }

  fml <- make_formula(
    fixed_terms,
    graph_obj
  )

  fit <- INLA::inla(
    formula = fml,
    family = "nbinomial",
    data = model_dat,
    control.predictor = list(
      compute = FALSE
    ),
    control.compute = list(
      waic = TRUE,
      dic = TRUE,
      cpo = FALSE,
      config = FALSE
    ),
    verbose = FALSE
  )

  if (!inherits(fit, "inla")) {
    stop(
      model_name,
      ": invalid INLA fit."
    )
  }

  ans <- list(
    irr = extract_irrs(
      fit,
      label
    ),
    fit = extract_fit(
      fit,
      label
    )
  )

  saveRDS(
    ans,
    summary_path
  )

  rm(
    fit,
    model_dat
  )

  invisible(
    gc()
  )

  ans
}


# ==============================================================================
# 3. LOAD PRIMARY DATA AND GRAPH
# ==============================================================================

if (!file.exists(PRIMARY_DATA_PATH)) {
  stop(
    "Run 01_REANALYSIS_934_PRIMARY_620_SENSITIVITY_FINAL.R first."
  )
}

if (!file.exists(PRIMARY_GRAPH_PATH)) {
  stop(
    "Primary 934-woreda graph not found: ",
    PRIMARY_GRAPH_PATH
  )
}

dat <- readRDS(
  PRIMARY_DATA_PATH
)

graph_primary <- INLA::inla.read.graph(
  PRIMARY_GRAPH_PATH
)

stopifnot(
  n_distinct(
    dat$id_1082
  ) ==
    934
)

SHAPE_FILE <- first_existing(
  SHAPE_CANDIDATES,
  "woreda shapefile"
)


# ==============================================================================
# 4. SUBSET-GRAPH HELPER
# ==============================================================================

prepare_subset <- function(
    d,
    ids_keep,
    label
) {

  dd <- d %>%
    filter(
      id_1082 %in%
        ids_keep
    )

  shp <- sf::st_read(
    SHAPE_FILE,
    quiet = TRUE
  ) %>%
    sf::st_make_valid() %>%
    mutate(
      id_1082 =
        as.integer(
          id_1082
        )
    ) %>%
    filter(
      id_1082 %in%
        ids_keep
    ) %>%
    arrange(
      id_1082
    )

  if (nrow(shp) !=
      length(
        unique(
          ids_keep
        )
      )) {
    stop(
      label,
      ": shapefile/subset mismatch."
    )
  }

  lookup <- shp %>%
    st_drop_geometry() %>%
    transmute(
      id_1082,
      woreda_idx_new =
        row_number()
    )

  shp <- shp %>%
    left_join(
      lookup,
      by = "id_1082"
    ) %>%
    arrange(
      woreda_idx_new
    )

  nb <- spdep::poly2nb(
    shp,
    queen = TRUE,
    row.names =
      shp$woreda_idx_new
  )

  safe_label <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    label
  )

  graph_path <- file.path(
    SENS_DIR,
    paste0(
      "graph_",
      safe_label,
      ".adj"
    )
  )

  spdep::nb2INLA(
    graph_path,
    nb
  )

  graph <- INLA::inla.read.graph(
    graph_path
  )

  dd <- dd %>%
    select(
      -woreda_idx
    ) %>%
    left_join(
      lookup,
      by = "id_1082"
    ) %>%
    rename(
      woreda_idx =
        woreda_idx_new
    ) %>%
    arrange(
      woreda_idx,
      date
    ) %>%
    mutate(
      st_id =
        as.integer(
          interaction(
            woreda_idx,
            time_month,
            drop = TRUE
          )
        )
    )

  list(
    data = dd,
    graph = graph,
    nb = nb,
    n_woredas =
      n_distinct(
        dd$id_1082
      )
  )
}


# ==============================================================================
# 5. ADDITIONAL COVARIATES ON THE 934-PRIMARY SCALE
# ==============================================================================

dat <- dat %>%
  group_by(
    id_1082
  ) %>%
  mutate(
    mean_llin =
      mean(
        prop_with_llin_access2,
        na.rm = TRUE
      ),
    mean_irs =
      mean(
        prop_with_irs_cov,
        na.rm = TRUE
      )
  ) %>%
  ungroup()


dat$idp_log_std <-
  z_primary(
    log1p(
      dat$idp_ind
    )
  )

dat$returning_idp_log_std <-
  z_primary(
    log1p(
      dat$returning_idp_ind
    )
  )

dat$conflict_log_std <-
  z_primary(
    log1p(
      dat$conflict
    )
  )

dat$fatalities_log_std <-
  z_primary(
    log1p(
      dat$fatalities
    )
  )

dat$mean_llin_std <-
  z_primary(
    dat$mean_llin
  )

dat$mean_irs_std <-
  z_primary(
    dat$mean_irs
  )

dat$ins_res_zone_std <-
  z_primary(
    dat$ins_res_in_zone
  )

dat$hrp2_prop_zone_std <-
  z_primary(
    dat$hrp2_prop_in_zone
  )

dat$act_res_zone_std <-
  z_primary(
    dat$act_res_in_zone
  )

dat$stephensi_woreda_binary <-
  factor(
    as.integer(
      dat$stephensi_in_woreda >
        0
    ),
    levels = c(
      0,
      1
    )
  )


# ==============================================================================
# 6. HIGH-COMPLETENESS >=80% EACH YEAR
# ==============================================================================

results <- list()

completeness <- dat %>%
  group_by(
    id_1082,
    year
  ) %>%
  summarise(
    both_species_observed_pct =
      100 *
      mean(
        !is.na(
          pf_cases_model
        ) &
          !is.na(
            pv_cases_model
          )
      ),
    .groups = "drop"
  )


high_complete_ids <- completeness %>%
  group_by(
    id_1082
  ) %>%
  summarise(
    qualifies =
      all(
        both_species_observed_pct >=
          80
      ),
    .groups = "drop"
  ) %>%
  filter(
    qualifies
  ) %>%
  pull(
    id_1082
  )


high_complete_audit <- completeness %>%
  filter(
    id_1082 %in%
      high_complete_ids
  )


message(
  "High-completeness sensitivity woredas: ",
  length(
    high_complete_ids
  )
)


if (RUN_HIGH_COMPLETENESS_80) {

  high80 <- prepare_subset(
    dat,
    high_complete_ids,
    "high_completeness_80_each_year"
  )

  results[["high80_pf"]] <-
    fit_sensitivity(
      high80$data,
      "pf_cases_model",
      high80$graph,
      "SENS_high80_each_year_Pf",
      paste0(
        "P. falciparum — >=80% species reporting each year (",
        high80$n_woredas,
        " woredas)"
      )
    )

  results[["high80_pv"]] <-
    fit_sensitivity(
      high80$data,
      "pv_cases_model",
      high80$graph,
      "SENS_high80_each_year_Pv",
      paste0(
        "P. vivax — >=80% species reporting each year (",
        high80$n_woredas,
        " woredas)"
      )
    )

  rm(
    high80
  )

  invisible(
    gc()
  )
}


# ==============================================================================
# 7. ELEVATION SENSITIVITY
# ==============================================================================

if (RUN_NO_ELEVATION) {

  no_elev_terms <- setdiff(
    PRIMARY_FIXED,
    "elev_pop_weighted_std"
  )

  results[["noelev_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_PRIMARY934_no_elevation_Pf",
      "P. falciparum — elevation omitted",
      fixed_terms =
        no_elev_terms
    )

  results[["noelev_pv"]] <-
    fit_sensitivity(
      dat,
      "pv_cases_model",
      graph_primary,
      "SENS_PRIMARY934_no_elevation_Pv",
      "P. vivax — elevation omitted",
      fixed_terms =
        no_elev_terms
    )
}


if (RUN_EXCLUDE_NEW_ELEVATION_EXTENSION) {

  ids_no_extension <- dat %>%
    group_by(
      id_1082
    ) %>%
    summarise(
      newly_extended =
        any(
          new_extension_imputed ==
            1,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) %>%
    filter(
      !newly_extended
    ) %>%
    pull(
      id_1082
    )

  noext <- prepare_subset(
    dat,
    ids_no_extension,
    "exclude_new_75_elevation_extension"
  )

  message(
    "Elevation-extension exclusion sensitivity woredas: ",
    noext$n_woredas
  )

  results[["noext_pf"]] <-
    fit_sensitivity(
      noext$data,
      "pf_cases_model",
      noext$graph,
      "SENS_exclude_new_elevation_extension_Pf",
      paste0(
        "P. falciparum — exclude newly extension-imputed elevation woredas (",
        noext$n_woredas,
        ")"
      )
    )

  results[["noext_pv"]] <-
    fit_sensitivity(
      noext$data,
      "pv_cases_model",
      noext$graph,
      "SENS_exclude_new_elevation_extension_Pv",
      paste0(
        "P. vivax — exclude newly extension-imputed elevation woredas (",
        noext$n_woredas,
        ")"
      )
    )

  rm(
    noext
  )

  invisible(
    gc()
  )
}


# ==============================================================================
# 8. IDP / CONFLICT SENSITIVITIES
# ==============================================================================

if (RUN_LOG_IDP) {

  fixed_log_idp <- PRIMARY_FIXED
  fixed_log_idp[
    fixed_log_idp ==
      "idp_ind_std"
  ] <- "idp_log_std"

  results[["logidp_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_log_IDP_Pf",
      "P. falciparum — log-transformed IDP",
      fixed_terms =
        fixed_log_idp
    )

  results[["logidp_pv"]] <-
    fit_sensitivity(
      dat,
      "pv_cases_model",
      graph_primary,
      "SENS_log_IDP_Pv",
      "P. vivax — log-transformed IDP",
      fixed_terms =
        fixed_log_idp
    )
}


if (RUN_RETURNING_IDP) {

  results[["returning_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_returning_IDP_added_Pf",
      "P. falciparum — returning IDP added",
      fixed_terms =
        c(
          PRIMARY_FIXED,
          "returning_idp_log_std"
        )
    )

  results[["returning_pv"]] <-
    fit_sensitivity(
      dat,
      "pv_cases_model",
      graph_primary,
      "SENS_returning_IDP_added_Pv",
      "P. vivax — returning IDP added",
      fixed_terms =
        c(
          PRIMARY_FIXED,
          "returning_idp_log_std"
        )
    )
}


if (RUN_CONFLICT) {

  results[["conflict_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_conflict_added_Pf",
      "P. falciparum — conflict events added",
      fixed_terms =
        c(
          PRIMARY_FIXED,
          "conflict_log_std"
        )
    )

  results[["conflict_pv"]] <-
    fit_sensitivity(
      dat,
      "pv_cases_model",
      graph_primary,
      "SENS_conflict_added_Pv",
      "P. vivax — conflict events added",
      fixed_terms =
        c(
          PRIMARY_FIXED,
          "conflict_log_std"
        )
    )
}


if (RUN_FATALITIES) {

  results[["fatalities_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_fatalities_added_Pf",
      "P. falciparum — conflict-related fatalities added",
      fixed_terms =
        c(
          PRIMARY_FIXED,
          "fatalities_log_std"
        )
    )

  results[["fatalities_pv"]] <-
    fit_sensitivity(
      dat,
      "pv_cases_model",
      graph_primary,
      "SENS_fatalities_added_Pv",
      "P. vivax — conflict-related fatalities added",
      fixed_terms =
        c(
          PRIMARY_FIXED,
          "fatalities_log_std"
        )
    )
}


# ==============================================================================
# 9. AN. STEPHENSI CODING
# ==============================================================================

if (RUN_STEPHENSI_WOREDA_ONLY) {

  fixed_steph <- PRIMARY_FIXED
  fixed_steph[
    fixed_steph ==
      "stephensi_in_neighboring_woredas"
  ] <- "stephensi_woreda_binary"

  results[["steph_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_stephensi_woreda_only_Pf",
      "P. falciparum — An. stephensi in woreda only",
      fixed_terms =
        fixed_steph
    )

  results[["steph_pv"]] <-
    fit_sensitivity(
      dat,
      "pv_cases_model",
      graph_primary,
      "SENS_stephensi_woreda_only_Pv",
      "P. vivax — An. stephensi in woreda only",
      fixed_terms =
        fixed_steph
    )
}


# ==============================================================================
# 10. INTERVENTION TARGETING + DELTAMETHRIN RESISTANCE
# ==============================================================================

if (RUN_PROGRAMME_TARGETING_RESISTANCE) {

  fixed_programme <- unique(
    c(
      PRIMARY_FIXED,
      "mean_llin_std",
      "mean_irs_std",
      "ins_res_zone_std"
    )
  )

  results[["programme_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_programme_targeting_resistance_Pf",
      "P. falciparum — average LLIN/IRS + deltamethrin resistance added",
      fixed_terms =
        fixed_programme
    )

  results[["programme_pv"]] <-
    fit_sensitivity(
      dat,
      "pv_cases_model",
      graph_primary,
      "SENS_programme_targeting_resistance_Pv",
      "P. vivax — average LLIN/IRS + deltamethrin resistance added",
      fixed_terms =
        fixed_programme
    )
}


# ==============================================================================
# 11. EXPLORATORY HRP2 / ACT INDICATORS
# ==============================================================================

if (RUN_EXPLORATORY_HRP2_ACT) {

  results[["hrp2_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_HRP2_indicator_exploratory_Pf",
      "P. falciparum — HRP2 indicator exploratory",
      fixed_terms =
        c(
          PRIMARY_FIXED,
          "hrp2_prop_zone_std"
        )
    )

  results[["act_pf"]] <-
    fit_sensitivity(
      dat,
      "pf_cases_model",
      graph_primary,
      "SENS_ACT_resistance_indicator_exploratory_Pf",
      "P. falciparum — ACT-resistance indicator exploratory",
      fixed_terms =
        c(
          PRIMARY_FIXED,
          "act_res_zone_std"
        )
    )
}


# ==============================================================================
# 12. EXPORT
# ==============================================================================

all_irrs <- bind_rows(
  purrr::map(
    results,
    "irr"
  )
)

all_fit <- bind_rows(
  purrr::map(
    results,
    "fit"
  )
)


# Retain selected sensitivity coefficients that are most directly interpretable
# in the supplement/reviewer response.
selected_terms <- c(
  "rainfall_lag2_std",
  "max_temp_lag2_std",
  "min_temp_night_lag2_std",
  "evi_lag2_std",
  "elev_pop_weighted_std",
  "idp_ind_std",
  "idp_log_std",
  "returning_idp_log_std",
  "conflict_log_std",
  "fatalities_log_std",
  "stephensi_in_neighboring_woredas1",
  "stephensi_woreda_binary1",
  "llin_within_std",
  "irs_within_std",
  "mean_llin_std",
  "mean_irs_std",
  "ins_res_zone_std",
  "hrp2_prop_zone_std",
  "act_res_zone_std"
)

selected_irrs <- all_irrs %>%
  filter(
    term %in%
      selected_terms
  )


coverage_audit <- tibble(
  sensitivity = c(
    "Primary all-HF",
    "High completeness >=80% each year",
    "Exclude newly extension-imputed elevation"
  ),
  n_woredas = c(
    n_distinct(dat$id_1082),
    length(high_complete_ids),
    n_distinct(
      dat$id_1082[
        dat$new_extension_imputed !=
          1 |
          is.na(
            dat$new_extension_imputed
          )
      ]
    )
  )
)


writexl::write_xlsx(
  list(
    selected_adjusted_IRRs =
      selected_irrs,
    all_adjusted_IRRs =
      all_irrs,
    model_fit =
      all_fit,
    high_completeness_audit =
      high_complete_audit,
    coverage_audit =
      coverage_audit
  ),
  file.path(
    TABLE_DIR,
    "06_ADDITIONAL_sensitivity_models_PRIMARY934.xlsx"
  )
)


capture.output(
  sessionInfo(),
  file =
    file.path(
      SENS_DIR,
      "sessionInfo_additional_sensitivity_934_primary.txt"
    )
)


