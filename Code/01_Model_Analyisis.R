################################################################################
# PRIMARY ANALYSIS

# ------------------------------------------------------------------------------
# All 934 informative non-Tigray woredas from all reporting health facilities.
# Species-specific Pf/Pv/mixed observations are used directly.
# 
# REPORTING-QUALITY SENSITIVITY
# ------------------------------------------------------------------------------
# 620 woredas containing >=1 retained health facility meeting the prespecified
# >=50% annual reporting criterion. The same full model specification and the
# same 934-primary covariate scaling are used.
#
# WHY THIS HIERARCHY
# ------------------------------------------------------------------------------
# The 934-woreda panel maximises geographic/population coverage after excluding
# Tigray, while annual species missingness is modest and declines over time.
# The 620-woreda panel provides a stringent reporting-quality sensitivity.
#
# Final MODEL
# ------------------------------------------------------------------------------
# Negative-binomial INLA model with:
#   - environmental/programmatic fixed effects
#   - region-specific linear monthly trend
#   - BYM2 woreda spatial effect
#   - RW1 long-term monthly trend
#   - cyclic RW1 seasonality
#   - IID woreda-month interaction
#   - log(population) offset
#

################################################################################

# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(tidyverse)
library(lubridate)
library(INLA)
library(sf)
library(spdep)
library(Matrix)
library(writexl)

set.seed(20260828)


# ==============================================================================
### raw data all health facilty exculding tigray and 620 woredas
PATH_RETAINED_HF <-
  "processed_output/eth_resurgence_retainedHF_620.csv"

PATH_ALL_HF <-
  "processed_output/eth_resurgence_allHF_noTigray_934.csv"


PATH_SHAPE_CANDIDATES <- c(
  "data/shp/eth_adm3_clean_id.shp",
  "eth_adm3_clean_id.shp"
)


OUT_DIR <-
  "result_species_only_reanalysis_934_primary"

TABLE_DIR <-
  file.path(
    OUT_DIR,
    "tables"
  )

MODEL_DIR <-
  file.path(
    OUT_DIR,
    "models"
  )

DATA_DIR <-
  file.path(
    OUT_DIR,
    "prepared_data"
  )

VALID_DIR <-
  file.path(
    OUT_DIR,
    "validation"
  )

purrr::walk(
  c(
    OUT_DIR,
    TABLE_DIR,
    MODEL_DIR,
    DATA_DIR,
    VALID_DIR
  ),
  ~dir.create(
    .x,
    recursive = TRUE,
    showWarnings = FALSE
  )
)


# ==============================================================================
# 2.  primary model fit 
# ==============================================================================

# Pf and Pv primary models 
RUN_PRIMARY_TOTAL_MODEL <- TRUE

# Retained-HF sensitivity uses the same full specification as the 934-woreda primary model.
RUN_RETAINED_HF_SENSITIVITY <- TRUE


# ==============================================================================
# 3.  HELPERS function
# ==============================================================================

resolve_first_existing <- function(
    candidates,
    label
) {

  hit <-
    candidates[
      file.exists(
        candidates
      )
    ]

  if (length(hit) == 0) {
    stop(
      "Could not find ",
      label,
      ". Checked:\n",
      paste(
        candidates,
        collapse = "\n"
      )
    )
  }

  normalizePath(
    hit[1],
    mustWork = TRUE
  )
}


safe_sum <- function(x) {

  if (all(is.na(x))) {
    return(NA_real_)
  }

  sum(
    x,
    na.rm = TRUE
  )
}


safe_mean <- function(x) {

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(
    x,
    na.rm = TRUE
  )
}


fill_small_gap <- function(x) {

  x <-
    as.numeric(
      x
    )

  idx <-
    seq_along(
      x
    )

  ok <-
    is.finite(
      x
    )

  if (sum(ok) == 0) {
    return(x)
  }

  if (sum(ok) == 1) {
    return(
      rep(
        x[ok][1],
        length(x)
      )
    )
  }

  approx(
    x =
      idx[ok],
    y =
      x[ok],
    xout =
      idx,
    method =
      "linear",
    rule =
      2
  )$y
}


extract_irrs <- function(
    model,
    model_name
) {

  model$summary.fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      "term"
    ) %>%
    transmute(
      model =
        model_name,
      term,
      beta =
        mean,
      beta_sd =
        sd,
      beta_lcl =
        `0.025quant`,
      beta_ucl =
        `0.975quant`,
      IRR =
        exp(
          mean
        ),
      IRR_lcl =
        exp(
          `0.025quant`
        ),
      IRR_ucl =
        exp(
          `0.975quant`
        ),
      `IRR (95% CrI)` =
        sprintf(
          "%.3f (%.3f–%.3f)",
          IRR,
          IRR_lcl,
          IRR_ucl
        )
    )
}


extract_hyperparameters <- function(
    model,
    model_name
) {

  model$summary.hyperpar %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      "hyperparameter"
    ) %>%
    mutate(
      model =
        model_name
    ) %>%
    relocate(
      model,
      hyperparameter
    )
}


extract_fit_statistics <- function(
    model,
    model_name
) {

  if (
    !inherits(
      model,
      "inla"
    )
  ) {

    stop(
      model_name,
      ": expected an INLA model object, but received class: ",
      paste(
        class(
          model
        ),
        collapse =
          ", "
      ),
      ". Refit the model using fit_inla_model()."
    )
  }


  waic_value <-
    if (
      is.list(
        model$waic
      ) &&
      !is.null(
        model$waic$waic
      )
    ) {
      as.numeric(
        model$waic$waic
      )
    } else {
      NA_real_
    }


  p_eff_value <-
    if (
      is.list(
        model$waic
      ) &&
      !is.null(
        model$waic$p.eff
      )
    ) {
      as.numeric(
        model$waic$p.eff
      )
    } else {
      NA_real_
    }


  dic_value <-
    if (
      is.list(
        model$dic
      ) &&
      !is.null(
        model$dic$dic
      )
    ) {
      as.numeric(
        model$dic$dic
      )
    } else {
      NA_real_
    }


  tibble(
    model =
      model_name,
    WAIC =
      waic_value,
    effective_parameters =
      p_eff_value,
    DIC =
      dic_value
  )
}


# ==============================================================================
# 4. READ THE TWO NEW SPECIES-ONLY DATASETS
# ==============================================================================

if (!file.exists(PATH_RETAINED_HF)) {
  stop(
    "Retained-HF dataset not found: ",
    PATH_RETAINED_HF
  )
}

if (!file.exists(PATH_ALL_HF)) {
  stop(
    "All-HF dataset not found: ",
    PATH_ALL_HF
  )
}


message(
  "Reading retained-HF primary dataset: ",
  normalizePath(
    PATH_RETAINED_HF,
    mustWork = TRUE
  )
)

retained_raw <- readr::read_csv(
  PATH_RETAINED_HF,
  col_types =
    readr::cols(
      .default =
        readr::col_guess(),
      date =
        readr::col_character()
    ),
  show_col_types =
    FALSE
)


message(
  "Reading all-HF sensitivity dataset: ",
  normalizePath(
    PATH_ALL_HF,
    mustWork = TRUE
  )
)

all_hf_raw <- readr::read_csv(
  PATH_ALL_HF,
  col_types =
    readr::cols(
      .default =
        readr::col_guess(),
      date =
        readr::col_character()
    ),
  show_col_types =
    FALSE
)


required_new_cols <- c(
  "id_1082",
  "date",
  "mixed_confirmed",
  "pf_confirmed",
  "pv_confirmed",
  "region",
  "zone",
  "woreda",
  "woreda_population",
  "prop_with_irs_cov",
  "prop_with_llin_access2",
  "conflict",
  "fatalities",
  "idp_ind",
  "returning_idp_ind",
  "rainfall",
  "rainfall_lag1",
  "rainfall_lag2",
  "rainfall_lag3",
  "max_temp_lag1",
  "max_temp_lag2",
  "max_temp_lag3",
  "min_temp_night_lag1",
  "min_temp_night_lag2",
  "min_temp_night_lag3",
  "evi",
  "evi_lag1",
  "evi_lag2",
  "evi_lag3",
  "hrp2_prop_in_zone",
  "ins_res_in_zone",
  "stephensi_in_woreda",
  "act_res_in_zone",
  "elev_pop_weghted"
  
)


for (
  object_name in
  c(
    "retained_raw",
    "all_hf_raw"
  )
) {

  d <-
    get(
      object_name
    )

  missing_cols <-
    setdiff(
      required_new_cols,
      names(
        d
      )
    )

  if (length(missing_cols) > 0) {
    stop(
      object_name,
      " is missing required variables: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
}


# ==============================================================================
# 6. CLEAN DATA WITHOUT IMPUTING SPECIES
# ==============================================================================

clean_species_panel <- function(
    d,
    dataset_label,
    exclude_tigray = FALSE
) {
  
  d <- d %>%
    mutate(
      id_1082 =
        as.integer(
          id_1082
        ),
      date =
        lubridate::ymd(
          date
        ),
      year =
        lubridate::year(
          date
        ),
      month =
        lubridate::month(
          date
        )
    )
  
  if (exclude_tigray) {
    
    d <- d %>%
      filter(
        region !=
          "Tigray"
      )
  }
  
  
  duplicated_panel <- d %>%
    count(
      id_1082,
      date,
      name =
        "n"
    ) %>%
    filter(
      n >
        1
    )
  
  
  if (nrow(duplicated_panel) > 0) {
    stop(
      dataset_label,
      " contains duplicated woreda-month rows."
    )
  }
  
  
  if (
    any(
      d$woreda_population <=
      0,
      na.rm =
      TRUE
    )
  ) {
    
    stop(
      dataset_label,
      " contains non-positive population values."
    )
  }
  

  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # Use species-specific case fields exactly as supplied.
  # Do not derive Pf/Pv from total positives.
  # Do not replace missing species outcomes with zero.
  # ---------------------------------------------------------------------------

  d <- d %>%
    mutate(
      pf_cases =
        if_else(
          is.na(
            pf_confirmed
          ) |
            is.na(
              mixed_confirmed
            ),
          NA_real_,
          as.numeric(
            pf_confirmed
          ) +
            as.numeric(
              mixed_confirmed
            )
        ),

      pv_cases =
        if_else(
          is.na(
            pv_confirmed
          ) |
            is.na(
              mixed_confirmed
            ),
          NA_real_,
          as.numeric(
            pv_confirmed
          ) +
            as.numeric(
              mixed_confirmed
            )
        ),

      total_cases =
        if_else(
          is.na(
            pf_confirmed
          ) |
            is.na(
              pv_confirmed
            ) |
            is.na(
              mixed_confirmed
            ),
          NA_real_,
          as.numeric(
            pf_confirmed
          ) +
            as.numeric(
              pv_confirmed
            ) +
            as.numeric(
              mixed_confirmed
            )
        ),

      # The negative-binomial likelihood is discrete.
      # Descriptive totals above remain unchanged.
      pf_cases_model =
        if_else(
          is.na(
            pf_cases
          ),
          NA_integer_,
          as.integer(
            round(
              pf_cases
            )
          )
        ),

      pv_cases_model =
        if_else(
          is.na(
            pv_cases
          ),
          NA_integer_,
          as.integer(
            round(
              pv_cases
            )
          )
        ),

      total_cases_model =
        if_else(
          is.na(
            total_cases
          ),
          NA_integer_,
          as.integer(
            round(
              total_cases
            )
          )
        ),

      pf_api =
        1000 *
        pf_cases /
        woreda_population,

      pv_api =
        1000 *
        pv_cases /
        woreda_population,

      total_api =
        1000 *
        total_cases /
        woreda_population,

      dataset =
        dataset_label
    )


  d <- d %>%
    mutate(
      elev_pop_weghted =
        as.numeric(
          elev_pop_weghted
        ),

      elevation_imputed =
        as.integer(
          elevation_imputed
        ),

      new_extension_imputed =
        as.integer(
          new_extension_imputed
        )
    )


  if (
    anyNA(
      d$elev_pop_weghted
    )
  ) {

    stop(
      dataset_label,
      ": elevation is missing in the merged analytical input."
    )
  }


  d <- d %>%
    mutate(
      max_temp_lag2_original_missing =
        as.integer(
          is.na(
            max_temp_lag2
          )
        ),

      min_temp_night_lag2_original_missing =
        as.integer(
          is.na(
            min_temp_night_lag2
          )
        )
    ) %>%
    group_by(
      id_1082
    ) %>%
    arrange(
      date,
      .by_group =
        TRUE
    ) %>%
    mutate(
      max_temp_lag1_model =
        fill_small_gap(
          max_temp_lag1
        ),

      max_temp_lag2_model =
        fill_small_gap(
          max_temp_lag2
        ),

      max_temp_lag3_model =
        fill_small_gap(
          max_temp_lag3
        ),

      min_temp_night_lag1_model =
        fill_small_gap(
          min_temp_night_lag1
        ),

      min_temp_night_lag2_model =
        fill_small_gap(
          min_temp_night_lag2
        ),

      min_temp_night_lag3_model =
        fill_small_gap(
          min_temp_night_lag3
        )
    ) %>%
    ungroup()


  if (
    anyNA(
      d$max_temp_lag2_model
    ) ||
    anyNA(
      d$min_temp_night_lag2_model
    )
  ) {

    stop(
      dataset_label,
      ": climate covariate gaps remain after interpolation. "
    )
  }


  # Within-woreda intervention deviations are the primary intervention measures.
  d <- d %>%
    group_by(
      id_1082
    ) %>%
    mutate(
      mean_llin =
        mean(
          prop_with_llin_access2,
          na.rm =
            TRUE
        ),

      mean_irs =
        mean(
          prop_with_irs_cov,
          na.rm =
            TRUE
        ),

      llin_within =
        prop_with_llin_access2 -
        mean_llin,

      irs_within =
        prop_with_irs_cov -
        mean_irs
    ) %>%
    ungroup()


  d
}



# Primary : all 934 informative non-Tigray woredas.
all_hf_base <- clean_species_panel(
  all_hf_raw,
  dataset_label = "All HFs, Tigray excluded (PRIMARY)",
  exclude_tigray = FALSE
)

# Reporting-quality sensitivity: retained HFs after the >=50% criterion.
retained_base <- clean_species_panel(
  retained_raw,
  dataset_label = "Retained HFs >=50% (SENSITIVITY)",
  exclude_tigray = FALSE
)

if (any(all_hf_base$region == "Tigray", na.rm = TRUE)) {
  stop("Tigray is unexpectedly present in the 934-woreda primary input.")
}

if (any(retained_base$region == "Tigray", na.rm = TRUE)) {
  stop("Tigray is unexpectedly present in the retained-HF sensitivity input.")
}

# Stable analytical aliases used below.
primary_base <- all_hf_base
sensitivity_base <- retained_base


# ==============================================================================
# 7. DATA  BEFORE MODELLING
# ==============================================================================

panel_summary <- bind_rows(
  primary_base %>%
    summarise(
      dataset = "PRIMARY: all HFs, Tigray excluded",
      rows = n(),
      woredas = n_distinct(id_1082),
      months = n_distinct(date),
      regions = n_distinct(region),
      pf_missing_rows = sum(is.na(pf_cases)),
      pv_missing_rows = sum(is.na(pv_cases)),
      elevation_missing_rows = sum(is.na(elev_pop_weghted)),
      newly_extended_elevation_woredas =
        n_distinct(id_1082[new_extension_imputed == 1])
    ),

  sensitivity_base %>%
    summarise(
      dataset = "SENSITIVITY: retained HFs >=50%",
      rows = n(),
      woredas = n_distinct(id_1082),
      months = n_distinct(date),
      regions = n_distinct(region),
      pf_missing_rows = sum(is.na(pf_cases)),
      pv_missing_rows = sum(is.na(pv_cases)),
      elevation_missing_rows = sum(is.na(elev_pop_weghted)),
      newly_extended_elevation_woredas =
        n_distinct(id_1082[new_extension_imputed == 1])
    )
)

print(panel_summary)

stopifnot(
  nrow(primary_base) == 934 * 60,
  n_distinct(primary_base$id_1082) == 934,
  nrow(sensitivity_base) == 620 * 60,
  n_distinct(sensitivity_base$id_1082) == 620,
  !anyNA(primary_base$elev_pop_weghted),
  !anyNA(sensitivity_base$elev_pop_weghted)
)

primary_species_coverage_by_woreda <- primary_base %>%
  group_by(id_1082, region, zone, woreda) %>%
  summarise(
    n_pf_observed_months = sum(!is.na(pf_cases)),
    n_pv_observed_months = sum(!is.na(pv_cases)),
    n_both_species_observed_months =
      sum(!is.na(pf_cases) & !is.na(pv_cases)),
    .groups = "drop"
  )

no_species_information <- primary_species_coverage_by_woreda %>%
  filter(
    n_pf_observed_months == 0 &
      n_pv_observed_months == 0
  )

if (nrow(no_species_information) > 0) {
  stop("The 934-woreda primary file contains a woreda with no species information.")
}

primary_elevation_audit <- primary_base %>%
  distinct(
    id_1082, region, zone, woreda,
    elev_pop_weghted, elevation_source,
    elevation_imputed, new_extension_imputed
  ) %>%
  arrange(id_1082)

message("PRIMARY all-HF woredas: ", n_distinct(primary_base$id_1082))
message("SENSITIVITY retained-HF woredas: ", n_distinct(sensitivity_base$id_1082))
message(
  "Newly extension-imputed elevation woredas in primary: ",
  sum(primary_elevation_audit$new_extension_imputed == 1, na.rm = TRUE)
)


# ==============================================================================
# 8. DESCRIPTIVE TOTALS AND SPECIES MISSINGNESS
# ==============================================================================

annual_summary <- function(d) {
  d %>%
    group_by(year) %>%
    summarise(
      n_woredas = n_distinct(id_1082),
      n_woreda_months = n(),
      pf_observed_woreda_months = sum(!is.na(pf_cases)),
      pv_observed_woreda_months = sum(!is.na(pv_cases)),
      pf_missing_woreda_months = sum(is.na(pf_cases)),
      pv_missing_woreda_months = sum(is.na(pv_cases)),
      pf_missing_pct = 100 * mean(is.na(pf_cases)),
      pv_missing_pct = 100 * mean(is.na(pv_cases)),
      pf_cases_reported = safe_sum(pf_cases),
      pv_cases_reported = safe_sum(pv_cases),
      total_cases_reported = safe_sum(total_cases),
      represented_population =
        sum(woreda_population[!duplicated(id_1082)], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year) %>%
    mutate(
      pf_year_on_year_pct = 100 * (pf_cases_reported / lag(pf_cases_reported) - 1),
      pv_year_on_year_pct = 100 * (pv_cases_reported / lag(pv_cases_reported) - 1),
      total_year_on_year_pct = 100 * (total_cases_reported / lag(total_cases_reported) - 1),
      pf_change_from_2020_pct = 100 * (pf_cases_reported / first(pf_cases_reported) - 1),
      pv_change_from_2020_pct = 100 * (pv_cases_reported / first(pv_cases_reported) - 1),
      total_change_from_2020_pct = 100 * (total_cases_reported / first(total_cases_reported) - 1)
    )
}

primary_annual <- annual_summary(primary_base) %>%
  mutate(dataset = "PRIMARY: all HFs, Tigray excluded") %>%
  relocate(dataset)

sensitivity_annual <- annual_summary(sensitivity_base) %>%
  mutate(dataset = "SENSITIVITY: retained HFs >=50%") %>%
  relocate(dataset)

annual_cases_all <- bind_rows(
  primary_annual,
  sensitivity_annual
)

missingness_year <- bind_rows(
  primary_base %>%
    group_by(year) %>%
    summarise(
      dataset = "PRIMARY: all HFs, Tigray excluded",
      rows = n(),
      pf_missing = sum(is.na(pf_cases)),
      pv_missing = sum(is.na(pv_cases)),
      mixed_missing = sum(is.na(mixed_confirmed)),
      pf_missing_pct = 100 * pf_missing / rows,
      pv_missing_pct = 100 * pv_missing / rows,
      .groups = "drop"
    ),

  sensitivity_base %>%
    group_by(year) %>%
    summarise(
      dataset = "SENSITIVITY: retained HFs >=50%",
      rows = n(),
      pf_missing = sum(is.na(pf_cases)),
      pv_missing = sum(is.na(pv_cases)),
      mixed_missing = sum(is.na(mixed_confirmed)),
      pf_missing_pct = 100 * pf_missing / rows,
      pv_missing_pct = 100 * pv_missing / rows,
      .groups = "drop"
    )
)

missingness_region_year <- bind_rows(
  primary_base %>%
    group_by(region, year) %>%
    summarise(
      dataset = "PRIMARY: all HFs, Tigray excluded",
      rows = n(),
      pf_missing = sum(is.na(pf_cases)),
      pv_missing = sum(is.na(pv_cases)),
      pf_missing_pct = 100 * pf_missing / rows,
      pv_missing_pct = 100 * pv_missing / rows,
      .groups = "drop"
    ),

  sensitivity_base %>%
    group_by(region, year) %>%
    summarise(
      dataset = "SENSITIVITY: retained HFs >=50%",
      rows = n(),
      pf_missing = sum(is.na(pf_cases)),
      pv_missing = sum(is.na(pv_cases)),
      pf_missing_pct = 100 * pf_missing / rows,
      pv_missing_pct = 100 * pv_missing / rows,
      .groups = "drop"
    )
)

writexl::write_xlsx(
  list(
    panel_summary = panel_summary,
    annual_cases = annual_cases_all,
    missingness_by_year = missingness_year,
    missingness_region_year = missingness_region_year,
    primary_woreda_species_coverage = primary_species_coverage_by_woreda,
    primary_elevation_audit = primary_elevation_audit
  ),
  file.path(
    TABLE_DIR,
    "01_PRIMARY934_species_only_data_audit.xlsx"
  )
)


# ==============================================================================
# 9. STANDARDISATION — USE THE 934-WOREDA PRIMARY SCALE
# ==============================================================================
# Continuous covariates in BOTH the primary and 620-woreda sensitivity models
# are standardised with the 934-primary mean and SD. Therefore a one-SD IRR is
# directly comparable across the two analytical panels.

STANDARDISE_MAP <- c(
  rainfall_lag2_std = "rainfall_lag2",
  max_temp_lag2_std = "max_temp_lag2_model",
  min_temp_night_lag2_std = "min_temp_night_lag2_model",
  evi_lag2_std = "evi_lag2",
  elev_pop_weighted_std = "elev_pop_weghted",
  idp_ind_std = "idp_ind",
  llin_within_std = "llin_within",
  irs_within_std = "irs_within"
)

primary_scaler <- purrr::map_dfr(
  names(STANDARDISE_MAP),
  function(std_name) {
    raw_name <- STANDARDISE_MAP[[std_name]]
    mu <- mean(primary_base[[raw_name]], na.rm = TRUE)
    sigma <- sd(primary_base[[raw_name]], na.rm = TRUE)

    if (!is.finite(sigma) || sigma <= 0) {
      stop("Cannot standardise ", raw_name, ": undefined/non-positive SD.")
    }

    tibble(
      standardized_variable = std_name,
      source_variable = raw_name,
      primary_mean = mu,
      primary_sd = sigma
    )
  }
)

apply_primary_scaling <- function(d) {
  for (j in seq_len(nrow(primary_scaler))) {
    std_name <- primary_scaler$standardized_variable[j]
    raw_name <- primary_scaler$source_variable[j]
    mu <- primary_scaler$primary_mean[j]
    sigma <- primary_scaler$primary_sd[j]

    d[[std_name]] <- (d[[raw_name]] - mu) / sigma
  }
  d
}

primary_scaled <- apply_primary_scaling(primary_base)
sensitivity_scaled <- apply_primary_scaling(sensitivity_base)

writexl::write_xlsx(
  list(primary_scaling = primary_scaler),
  file.path(
    TABLE_DIR,
    "02_PRIMARY934_covariate_scaling_parameters.xlsx"
  )
)

# ==============================================================================
# 10. READ SHAPEFILE AND BUILD THE FULL GEOGRAPHIC NEIGHBOUR STRUCTURE
# ==============================================================================

PATH_SHAPE <-
  resolve_first_existing(
    PATH_SHAPE_CANDIDATES,
    "woreda shapefile"
  )


woredas_all <- sf::st_read(
  PATH_SHAPE,
  quiet =
    TRUE
) %>%
  sf::st_make_valid()


if (
  !"id_1082" %in%
  names(
    woredas_all
  )
) {

  stop(
    "The shapefile must contain the stable woreda identifier `id_1082`."
  )
}


woredas_all <- woredas_all %>%
  mutate(
    id_1082 =
      as.integer(
        id_1082
      )
  ) %>%
  arrange(
    id_1082
  )


if (
  anyDuplicated(
    woredas_all$id_1082
  ) >
    0
) {

  stop(
    "The shapefile contains duplicated `id_1082` values."
  )
}


# ------------------------------------------------------------------------------
# 10.1 Full geographic adjacency
# ------------------------------------------------------------------------------

woredas_full_geo <- woredas_all %>%
  mutate(
    full_geo_idx =
      row_number()
  )


nb_full_geo <- spdep::poly2nb(
  woredas_full_geo,
  queen =
    TRUE,
  row.names =
    woredas_full_geo$full_geo_idx
)


full_geo_components <-
  spdep::n.comp.nb(
    nb_full_geo
  )


full_geo_audit <- woredas_full_geo %>%
  sf::st_drop_geometry() %>%
  transmute(
    id_1082,
    full_geo_idx,
    spatial_component =
      full_geo_components$comp.id,
    n_touching_neighbours =
      spdep::card(
        nb_full_geo
      )
  )


readr::write_csv(
  full_geo_audit,
  file.path(
    DATA_DIR,
    "FULL_shapefile_geographic_adjacency_audit.csv"
  )
)


message(
  "Full shapefile adjacency: ",
  full_geo_components$nc,
  " connected component(s); ",
  sum(
    full_geo_audit$n_touching_neighbours ==
      0
  ),
  " singleton polygon(s)."
)


# Binary adjacency + self for direct geographic exposure.
A_full_geo <- spdep::nb2mat(
  nb_full_geo,
  style =
    "B",
  zero.policy =
    TRUE
)


A_full_geo_self <- Matrix::Matrix(
  A_full_geo +
    diag(
      nrow(
        A_full_geo
      )
    ),
  sparse =
    TRUE
)


full_geo_lookup <- woredas_full_geo %>%
  sf::st_drop_geometry() %>%
  select(
    id_1082,
    full_geo_idx
  )


# ------------------------------------------------------------------------------
# 10.2 Full-geography An. stephensi exposure

stephensi_source <- all_hf_raw %>%
  transmute(
    id_1082 =
      as.integer(
        id_1082
      ),

    date =
      lubridate::ymd(
        date
      ),

    stephensi_in_woreda =
      as.numeric(
        stephensi_in_woreda
      )
  ) %>%
  group_by(
    id_1082,
    date
  ) %>%
  summarise(
    stephensi_in_woreda =
      as.integer(
        any(
          stephensi_in_woreda >
            0,
          na.rm =
            TRUE
        )
      ),
    .groups =
      "drop"
  )


analysis_dates <-
  sort(
    unique(
      lubridate::ymd(
        all_hf_raw$date
      )
    )
  )


stephensi_full_grid <- tidyr::expand_grid(
  date =
    analysis_dates,
  full_geo_idx =
    seq_len(
      nrow(
        woredas_full_geo
      )
    )
) %>%
  left_join(
    full_geo_lookup,
    by =
      "full_geo_idx"
  ) %>%
  left_join(
    stephensi_source,
    by =
      c(
        "id_1082",
        "date"
      )
  ) %>%
  mutate(
    stephensi_in_woreda =
      tidyr::replace_na(
        stephensi_in_woreda,
        0L
      )
  ) %>%
  arrange(
    date,
    full_geo_idx
  )


stephensi_full_neighbor_lookup <- stephensi_full_grid %>%
  group_by(
    date
  ) %>%
  group_modify(
    ~{

      .x <-
        .x %>%
        arrange(
          full_geo_idx
        )

      x <-
        as.numeric(
          .x$stephensi_in_woreda
        )


      if (
        length(
          x
        ) !=
          nrow(
            A_full_geo_self
          )
      ) {

        stop(
          "Full-geography An. stephensi vector does not match the shapefile ",
          "adjacency matrix."
        )
      }


      .x$stephensi_in_neighboring_woredas <-
        as.integer(
          as.vector(
            A_full_geo_self %*%
              x
          ) >
            0
        )


      .x
    }
  ) %>%
  ungroup() %>%
  select(
    id_1082,
    date,
    stephensi_in_woreda_full_source =
      stephensi_in_woreda,
    stephensi_in_neighboring_woredas
  )


readr::write_csv(
  stephensi_full_neighbor_lookup,
  file.path(
    DATA_DIR,
    "FULL_geographic_An_stephensi_neighbor_exposure.csv"
  )
)


# ==============================================================================
# 11. SPATIAL/TEMPORAL PREPARATION FUNCTION
# ==============================================================================

prepare_model_panel <- function(
    d,
    analysis_name
) {

  missing_shape_ids <-
    setdiff(
      unique(
        d$id_1082
      ),
      unique(
        woredas_all$id_1082
      )
    )


  if (
    length(
      missing_shape_ids
    ) >
      0
  ) {

    stop(
      analysis_name,
      ": data IDs absent from shapefile: ",
      paste(
        missing_shape_ids,
        collapse =
          ", "
      )
    )
  }


  woredas_model <- woredas_all %>%
    filter(
      id_1082 %in%
        unique(
          d$id_1082
        )
    ) %>%
    arrange(
      id_1082
    ) %>%
    mutate(
      woreda_idx =
        row_number()
    )


  woreda_lookup <- woredas_model %>%
    sf::st_drop_geometry() %>%
    select(
      id_1082,
      woreda_idx
    )


  d <- d %>%
    select(
      -any_of(
        c(
          "woreda_idx",
          "stephensi_in_neighboring_woredas",
          "stephensi_in_woreda_full_source",
          "time_month",
          "time_linear",
          "season_month",
          "st_id",
          "zone_idx"
        )
      )
    ) %>%
    left_join(
      woreda_lookup,
      by =
        "id_1082"
    ) %>%
    left_join(
      stephensi_full_neighbor_lookup,
      by =
        c(
          "id_1082",
          "date"
        )
    )


  if (
    anyNA(
      d$woreda_idx
    )
  ) {

    stop(
      analysis_name,
      ": woreda-index assignment failed."
    )
  }


  if (
    anyNA(
      d$stephensi_in_neighboring_woredas
    )
  ) {

    stop(
      analysis_name,
      ": full-geography An. stephensi exposure could not be assigned to all ",
      "woreda-month rows."
    )
  }


  # ---------------------------------------------------------------------------
  # Analytical BYM2 adjacency.
  # ---------------------------------------------------------------------------

  nb_woreda <- spdep::poly2nb(
    woredas_model,
    queen =
      TRUE,
    row.names =
      woredas_model$woreda_idx
  )


  component_info <-
    spdep::n.comp.nb(
      nb_woreda
    )


  neighbour_count <-
    spdep::card(
      nb_woreda
    )


  woreda_names <- d %>%
    distinct(
      id_1082,
      region,
      zone,
      woreda
    )


  spatial_graph_audit <- tibble(
    id_1082 =
      woredas_model$id_1082,

    woreda_idx =
      woredas_model$woreda_idx,

    spatial_component =
      component_info$comp.id,

    n_touching_neighbours =
      neighbour_count
  ) %>%
    left_join(
      woreda_names,
      by =
        "id_1082"
    ) %>%
    relocate(
      id_1082,
      region,
      zone,
      woreda,
      woreda_idx,
      spatial_component,
      n_touching_neighbours
    )


  component_summary <- spatial_graph_audit %>%
    count(
      spatial_component,
      name =
        "n_woredas"
    ) %>%
    arrange(
      desc(
        n_woredas
      )
    )


  singleton_woredas <- spatial_graph_audit %>%
    filter(
      n_touching_neighbours ==
        0
    )


  safe_analysis_name <-
    gsub(
      "[^A-Za-z0-9]+",
      "_",
      analysis_name
    )


  readr::write_csv(
    spatial_graph_audit,
    file.path(
      DATA_DIR,
      paste0(
        "spatial_graph_audit_",
        safe_analysis_name,
        ".csv"
      )
    )
  )


  readr::write_csv(
    component_summary,
    file.path(
      DATA_DIR,
      paste0(
        "spatial_component_summary_",
        safe_analysis_name,
        ".csv"
      )
    )
  )


  readr::write_csv(
    singleton_woredas,
    file.path(
      DATA_DIR,
      paste0(
        "spatial_singletons_",
        safe_analysis_name,
        ".csv"
      )
    )
  )


  message(
    analysis_name,
    ": ",
    component_info$nc,
    " connected spatial component(s); ",
    nrow(
      singleton_woredas
    ),
    " singleton woreda(s)."
  )


  if (
    nrow(
      singleton_woredas
    ) >
      0
  ) {

    message(
      "Singleton IDs saved to: ",
      file.path(
        DATA_DIR,
        paste0(
          "spatial_singletons_",
          safe_analysis_name,
          ".csv"
        )
      )
    )
  }


  
  # INLA's BYM2 term below uses adjust.for.con.comp = TRUE.

  adjacency_filename <-
    paste0(
      "woreda_",
      safe_analysis_name,
      ".adj"
    )


  adjacency_path <-
    file.path(
      DATA_DIR,
      adjacency_filename
    )


  spdep::nb2INLA(
    adjacency_path,
    nb_woreda
  )


  graph_obj <-
    INLA::inla.read.graph(
      adjacency_path
    )


  # ---------------------------------------------------------------------------
  # Temporal indices and factors.
  # ---------------------------------------------------------------------------

  d <- d %>%
    mutate(
      time_month =
        as.integer(
          factor(
            date,
            levels =
              sort(
                unique(
                  date
                )
              )
          )
        ),

      time_linear =
        time_month -
        mean(
          time_month
        ),

      season_month =
        month,

      region =
        factor(
          region
        ),

      stephensi_in_neighboring_woredas =
        factor(
          stephensi_in_neighboring_woredas,
          levels =
            c(
              0,
              1
            )
        )
    ) %>%
    mutate(
      st_id =
        as.integer(
          interaction(
            woreda_idx,
            time_month,
            drop =
              TRUE
          )
        )
    )


  zone_lookup <- d %>%
    distinct(
      zone
    ) %>%
    arrange(
      zone
    ) %>%
    mutate(
      zone_idx =
        row_number()
    )


  d <- d %>%
    left_join(
      zone_lookup,
      by =
        "zone"
    )


  list(
    data =
      d,
    woredas =
      woredas_model,
    nb =
      nb_woreda,
    graph =
      graph_obj,
    adjacency_path =
      adjacency_path,
    graph_audit =
      spatial_graph_audit,
    component_summary =
      component_summary,
    singleton_woredas =
      singleton_woredas
  )
}



# ------------------------------------------------------------------------------
# Prepare the two analytical spatial structures.
# ------------------------------------------------------------------------------

primary_spatial <- prepare_model_panel(
  primary_scaled,
  "all_HF_no_Tigray_primary"
)

sensitivity_spatial <- prepare_model_panel(
  sensitivity_scaled,
  "retained_HF_sensitivity"
)

# Save prepared datasets for every downstream script.
saveRDS(
  primary_spatial$data,
  file.path(
    DATA_DIR,
    "PRIMARY_allHF_noTigray_934_model_data.rds"
  )
)

saveRDS(
  sensitivity_spatial$data,
  file.path(
    DATA_DIR,
    "SENS_retained_HF_620_model_data.rds"
  )
)

graph_summary <- bind_rows(
  tibble(
    analysis = "PRIMARY: all HFs, Tigray excluded (934)",
    n_woredas = nrow(primary_spatial$woredas),
    n_components = nrow(primary_spatial$component_summary),
    n_singletons = nrow(primary_spatial$singleton_woredas)
  ),
  tibble(
    analysis = "SENSITIVITY: retained HFs >=50% (620)",
    n_woredas = nrow(sensitivity_spatial$woredas),
    n_components = nrow(sensitivity_spatial$component_summary),
    n_singletons = nrow(sensitivity_spatial$singleton_woredas)
  )
)

writexl::write_xlsx(
  list(
    graph_summary = graph_summary,
    primary_graph = primary_spatial$graph_audit,
    sensitivity_graph = sensitivity_spatial$graph_audit
  ),
  file.path(
    TABLE_DIR,
    "02b_PRIMARY934_SPATIAL_graph_audit.xlsx"
  )
)

# ==============================================================================
# 12. PRIMARY MODEL SPECIFICATION
# ==============================================================================

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
    fixed_terms
) {

  rhs_fixed <-
    paste(
      fixed_terms,
      collapse =
        " + "
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
    collapse =
      " + "
  )


  as.formula(
    paste0(
      ".outcome ~ ",
      rhs_fixed,
      " + ",
      rhs_random,
      " + offset(log(woreda_population))"
    )
  )
}


# ==============================================================================
# 13. MODEL-FITTING FUNCTION
# ==============================================================================

fit_inla_model <- function(
    d,
    outcome,
    graph_obj,
    model_name,
    fixed_terms =
      PRIMARY_FIXED,
    compute_cpo =
      FALSE,
    compute_config =
      FALSE
) {

  message("\n")
  message("============================================================")
  message("FITTING MODEL: ", model_name)
  message("============================================================")


  if (
    !outcome %in%
      names(
        d
      )
  ) {

    stop(
      model_name,
      ": outcome variable `",
      outcome,
      "` is absent from the model dataset."
    )
  }


  model_dat <- d %>%
    mutate(
      .outcome =
        as.numeric(
          .data[[
            outcome
          ]]
        )
    )


  # Negative-binomial likelihood requires non-negative integer counts.
  non_integer <-
    is.finite(
      model_dat$.outcome
    ) &
    abs(
      model_dat$.outcome -
        round(
          model_dat$.outcome
        )
    ) >
      1e-8


  if (
    any(
      non_integer
    )
  ) {

    message(
      "Rounding ",
      sum(
        non_integer
      ),
      " non-integer outcomes for the negative-binomial likelihood."
    )

    model_dat$.outcome[
      non_integer
    ] <-
      round(
        model_dat$.outcome[
          non_integer
        ]
      )
  }


  if (
    any(
      model_dat$.outcome <
        0,
      na.rm =
        TRUE
    )
  ) {

    stop(
      model_name,
      ": negative outcome values detected."
    )
  }


  n_likelihood <-
    sum(
      !is.na(
        model_dat$.outcome
      )
    )


  if (
    n_likelihood ==
      0
  ) {

    stop(
      model_name,
      ": no observed outcome rows are available for the likelihood."
    )
  }


  message(
    "Rows in model frame: ",
    nrow(
      model_dat
    )
  )

  message(
    "Observed likelihood rows: ",
    n_likelihood
  )

  message(
    "Woredas: ",
    n_distinct(
      model_dat$woreda_idx
    )
  )

  message(
    "Connected components handled with adjust.for.con.comp = TRUE."
  )

  message(
    "CPO: ",
    compute_cpo,
    "; posterior config: ",
    compute_config
  )


  fml <-
    make_formula(
      fixed_terms =
        fixed_terms
    )


  environment(
    fml
  ) <-
    environment()


  fit <- INLA::inla(
    formula =
      fml,

    family =
      "nbinomial",

    data =
      model_dat,

    control.predictor =
      list(
        compute =
          TRUE,
        link =
          1
      ),

   
    control.compute =
      list(
        waic =
          TRUE,
        dic =
          TRUE,
        cpo =
          compute_cpo,
        config =
          compute_config
      ),

    verbose =
      FALSE
  )


  if (
    !inherits(
      fit,
      "inla"
    )
  ) {

    stop(
      model_name,
      ": INLA did not return a valid `inla` model object."
    )
  }


  saveRDS(
    fit,
    file.path(
      MODEL_DIR,
      paste0(
        model_name,
        ".rds"
      )
    )
  )


  message(
    "Model completed and saved: ",
    file.path(
      MODEL_DIR,
      paste0(
        model_name,
        ".rds"
      )
    )
  )


  fit
}


# ==============================================================================
# 14.  PRIMARY 934-WOREDA ALL-HF MODELS
# ==============================================================================

primary_dat <- primary_spatial$data

mod_pf_primary <- fit_inla_model(
  d = primary_dat,
  outcome = "pf_cases_model",
  graph_obj = primary_spatial$graph,
  model_name = "PRIMARY_allHF_noTigray_Pf_full_934"
)

mod_pv_primary <- fit_inla_model(
  d = primary_dat,
  outcome = "pv_cases_model",
  graph_obj = primary_spatial$graph,
  model_name = "PRIMARY_allHF_noTigray_Pv_full_934"
)

stopifnot(
  inherits(mod_pf_primary, "inla"),
  inherits(mod_pv_primary, "inla")
)

if (RUN_PRIMARY_TOTAL_MODEL) {
  mod_total_primary <- fit_inla_model(
    d = primary_dat,
    outcome = "total_cases_model",
    graph_obj = primary_spatial$graph,
    model_name = "PRIMARY_allHF_noTigray_Total_full_934"
  )
}

# Stable aliases used by all downstream scripts.
saveRDS(
  mod_pf_primary,
  file.path(
    MODEL_DIR,
    "mod_pf_PRIMARY_934_species_only.rds"
  )
)

saveRDS(
  mod_pv_primary,
  file.path(
    MODEL_DIR,
    "mod_pv_PRIMARY_934_species_only.rds"
  )
)


# ==============================================================================
# 15. FIT 620-WOREDA REPORTING-QUALITY SENSITIVITY
# ==============================================================================

if (RUN_RETAINED_HF_SENSITIVITY) {

  sensitivity_dat <- sensitivity_spatial$data

  mod_pf_sens620 <- fit_inla_model(
    d = sensitivity_dat,
    outcome = "pf_cases_model",
    graph_obj = sensitivity_spatial$graph,
    model_name = "SENS_retainedHF_Pf_full_620",
    fixed_terms = PRIMARY_FIXED
  )

  mod_pv_sens620 <- fit_inla_model(
    d = sensitivity_dat,
    outcome = "pv_cases_model",
    graph_obj = sensitivity_spatial$graph,
    model_name = "SENS_retainedHF_Pv_full_620",
    fixed_terms = PRIMARY_FIXED
  )

  stopifnot(
    inherits(mod_pf_sens620, "inla"),
    inherits(mod_pv_sens620, "inla")
  )

  saveRDS(
    mod_pf_sens620,
    file.path(
      MODEL_DIR,
      "mod_pf_SENS_retainedHF_620.rds"
    )
  )

  saveRDS(
    mod_pv_sens620,
    file.path(
      MODEL_DIR,
      "mod_pv_SENS_retainedHF_620.rds"
    )
  )
}


# ==============================================================================
# 16. PRIMARY MODEL RESULTS
# ==============================================================================

primary_irrs <- bind_rows(
  extract_irrs(
    mod_pf_primary,
    "P. falciparum — PRIMARY all HFs, Tigray excluded (934)"
  ),
  extract_irrs(
    mod_pv_primary,
    "P. vivax — PRIMARY all HFs, Tigray excluded (934)"
  )
)

primary_hyperparameters <- bind_rows(
  extract_hyperparameters(
    mod_pf_primary,
    "P. falciparum — PRIMARY all HFs, Tigray excluded (934)"
  ),
  extract_hyperparameters(
    mod_pv_primary,
    "P. vivax — PRIMARY all HFs, Tigray excluded (934)"
  )
)

primary_fit <- bind_rows(
  extract_fit_statistics(
    mod_pf_primary,
    "P. falciparum — PRIMARY all HFs, Tigray excluded (934)"
  ),
  extract_fit_statistics(
    mod_pv_primary,
    "P. vivax — PRIMARY all HFs, Tigray excluded (934)"
  )
)

if (RUN_PRIMARY_TOTAL_MODEL) {
  primary_irrs <- bind_rows(
    primary_irrs,
    extract_irrs(
      mod_total_primary,
      "Total malaria — PRIMARY all HFs, Tigray excluded (934)"
    )
  )

  primary_hyperparameters <- bind_rows(
    primary_hyperparameters,
    extract_hyperparameters(
      mod_total_primary,
      "Total malaria — PRIMARY all HFs, Tigray excluded (934)"
    )
  )

  primary_fit <- bind_rows(
    primary_fit,
    extract_fit_statistics(
      mod_total_primary,
      "Total malaria — PRIMARY all HFs, Tigray excluded (934)"
    )
  )
}

writexl::write_xlsx(
  list(
    adjusted_IRRs = primary_irrs,
    hyperparameters = primary_hyperparameters,
    model_fit = primary_fit,
    annual_cases = primary_annual,
    scaling_parameters = primary_scaler,
    missingness_by_year =
      missingness_year %>%
      filter(dataset == "PRIMARY: all HFs, Tigray excluded"),
    elevation_audit = primary_elevation_audit
  ),
  file.path(
    TABLE_DIR,
    "03_PRIMARY_allHF_noTigray_934_INLA_results.xlsx"
  )
)


# ==============================================================================
# 17. RETAINED-HF 620-WOREDA SENSITIVITY RESULTS
# ==============================================================================

if (RUN_RETAINED_HF_SENSITIVITY) {

  sensitivity_irrs <- bind_rows(
    extract_irrs(
      mod_pf_sens620,
      "P. falciparum — SENSITIVITY retained HFs >=50% (620)"
    ),
    extract_irrs(
      mod_pv_sens620,
      "P. vivax — SENSITIVITY retained HFs >=50% (620)"
    )
  )

  sensitivity_fit <- bind_rows(
    extract_fit_statistics(
      mod_pf_sens620,
      "P. falciparum — SENSITIVITY retained HFs >=50% (620)"
    ),
    extract_fit_statistics(
      mod_pv_sens620,
      "P. vivax — SENSITIVITY retained HFs >=50% (620)"
    )
  )

  sensitivity_hyper <- bind_rows(
    extract_hyperparameters(
      mod_pf_sens620,
      "P. falciparum — SENSITIVITY retained HFs >=50% (620)"
    ),
    extract_hyperparameters(
      mod_pv_sens620,
      "P. vivax — SENSITIVITY retained HFs >=50% (620)"
    )
  )

  writexl::write_xlsx(
    list(
      adjusted_IRRs = sensitivity_irrs,
      model_fit = sensitivity_fit,
      hyperparameters = sensitivity_hyper,
      annual_cases = sensitivity_annual,
      missingness_by_year =
        missingness_year %>%
        filter(dataset == "SENSITIVITY: retained HFs >=50%")
    ),
    file.path(
      TABLE_DIR,
      "04_SENS_retainedHF_620_full_spec_results.xlsx"
    )
  )
}


# ==============================================================================
# 18. DIRECT PRIMARY 934 VS SENSITIVITY 620 COEFFICIENT COMPARISON
# ==============================================================================

main_terms <- c(
  "rainfall_lag2_std",
  "max_temp_lag2_std",
  "min_temp_night_lag2_std",
  "evi_lag2_std",
  "elev_pop_weighted_std",
  "idp_ind_std",
  "stephensi_in_neighboring_woredas1",
  "llin_within_std",
  "irs_within_std"
)

comparison_irrs <- bind_rows(
  primary_irrs %>%
    mutate(analysis = "Primary: all HFs, Tigray excluded (934)"),
  sensitivity_irrs %>%
    mutate(analysis = "Sensitivity: retained HFs >=50% (620)")
) %>%
  filter(term %in% main_terms) %>%
  mutate(
    species = case_when(
      str_detect(model, "falciparum") ~ "P. falciparum",
      str_detect(model, "vivax") ~ "P. vivax",
      TRUE ~ "Other"
    )
  ) %>%
  select(
    species, analysis, term,
    beta, beta_lcl, beta_ucl,
    IRR, IRR_lcl, IRR_ucl,
    `IRR (95% CrI)`
  ) %>%
  arrange(species, term, analysis)

writexl::write_xlsx(
  list(primary934_vs_sens620 = comparison_irrs),
  file.path(
    TABLE_DIR,
    "05_PRIMARY934_vs_SENS620_coefficient_comparison.xlsx"
  )
)


# ==============================================================================
# 19. PRIMARY 934 MODEL VALIDATION
# ==============================================================================

validation_one_model <- function(
    d,
    model,
    outcome,
    model_name
) {

  observed <- as.numeric(d[[outcome]])
  predicted <- as.numeric(model$summary.fitted.values[["mean"]])

  keep <- is.finite(observed) & is.finite(predicted)

  tibble(
    model = model_name,
    n_observed = sum(keep),
    observed_cases = sum(observed[keep], na.rm = TRUE),
    predicted_cases = sum(predicted[keep], na.rm = TRUE),
    observed_to_predicted = observed_cases / predicted_cases,
    MAE = mean(abs(observed[keep] - predicted[keep])),
    RMSE = sqrt(mean((observed[keep] - predicted[keep])^2)),
    log1p_RMSE =
      sqrt(mean((log1p(observed[keep]) - log1p(predicted[keep]))^2)),
    Spearman_r =
      suppressWarnings(
        cor(observed[keep], predicted[keep], method = "spearman")
      )
  )
}


primary_validation <- bind_rows(
  validation_one_model(
    primary_dat,
    mod_pf_primary,
    "pf_cases_model",
    "P. falciparum"
  ),
  validation_one_model(
    primary_dat,
    mod_pv_primary,
    "pv_cases_model",
    "P. vivax"
  )
)


residual_moran <- function(
    d,
    model,
    outcome,
    nb,
    model_name
) {

  temp <- d %>%
    transmute(
      woreda_idx,
      observed = as.numeric(.data[[outcome]]),
      predicted = as.numeric(model$summary.fitted.values[["mean"]])
    ) %>%
    filter(
      is.finite(observed),
      is.finite(predicted)
    ) %>%
    mutate(
      log_residual = log1p(observed) - log1p(predicted)
    ) %>%
    group_by(woreda_idx) %>%
    summarise(
      residual = mean(log_residual, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(woreda_idx)

  all_indices <- seq_along(nb)

  if (!all(all_indices %in% temp$woreda_idx)) {
    keep_idx <- all_indices %in% temp$woreda_idx

    nb_use <- spdep::subset.nb(
      nb,
      subset = keep_idx
    )

    temp <- temp %>%
      filter(
        woreda_idx %in% all_indices[keep_idx]
      ) %>%
      arrange(woreda_idx)

  } else {
    nb_use <- nb
  }

  listw <- spdep::nb2listw(
    nb_use,
    style = "W",
    zero.policy = TRUE
  )

  mt <- spdep::moran.test(
    temp$residual,
    listw,
    zero.policy = TRUE
  )

  tibble(
    model = model_name,
    Moran_I = unname(mt$estimate[["Moran I statistic"]]),
    expected_I = unname(mt$estimate[["Expectation"]]),
    p_value = mt$p.value
  )
}


primary_moran <- bind_rows(
  residual_moran(
    primary_dat,
    mod_pf_primary,
    "pf_cases_model",
    primary_spatial$nb,
    "P. falciparum"
  ),
  residual_moran(
    primary_dat,
    mod_pv_primary,
    "pv_cases_model",
    primary_spatial$nb,
    "P. vivax"
  )
)


writexl::write_xlsx(
  list(
    validation_metrics = primary_validation,
    residual_Morans_I = primary_moran
  ),
  file.path(
    VALID_DIR,
    "06_PRIMARY934_model_validation.xlsx"
  )
)


# ==============================================================================
# 20. SAVE PRIMARY 934 FITTED VALUES FOR MAIN/SUPPLEMENTARY FIGURES
# ==============================================================================

pf_primary_fitted <- primary_dat %>%
  mutate(
    fitted_mean =
      as.numeric(mod_pf_primary$summary.fitted.values[["mean"]]),
    fitted_median =
      as.numeric(mod_pf_primary$summary.fitted.values[["0.5quant"]])
  )

pv_primary_fitted <- primary_dat %>%
  mutate(
    fitted_mean =
      as.numeric(mod_pv_primary$summary.fitted.values[["mean"]]),
    fitted_median =
      as.numeric(mod_pv_primary$summary.fitted.values[["0.5quant"]])
  )

saveRDS(
  pf_primary_fitted,
  file.path(
    MODEL_DIR,
    "PRIMARY934_Pf_fitted_values.rds"
  )
)

saveRDS(
  pv_primary_fitted,
  file.path(
    MODEL_DIR,
    "PRIMARY934_Pv_fitted_values.rds"
  )
)


predicted_2024 <- bind_rows(
  pf_primary_fitted %>%
    filter(year == 2024) %>%
    group_by(id_1082, region, zone, woreda) %>%
    summarise(
      species = "P. falciparum",
      observed_cases = safe_sum(pf_cases_model),
      predicted_cases = sum(fitted_mean, na.rm = TRUE),
      population = mean(woreda_population, na.rm = TRUE),
      observed_api = 1000 * observed_cases / population,
      predicted_api = 1000 * predicted_cases / population,
      .groups = "drop"
    ),

  pv_primary_fitted %>%
    filter(year == 2024) %>%
    group_by(id_1082, region, zone, woreda) %>%
    summarise(
      species = "P. vivax",
      observed_cases = safe_sum(pv_cases_model),
      predicted_cases = sum(fitted_mean, na.rm = TRUE),
      population = mean(woreda_population, na.rm = TRUE),
      observed_api = 1000 * observed_cases / population,
      predicted_api = 1000 * predicted_cases / population,
      .groups = "drop"
    )
)

writexl::write_xlsx(
  list(observed_predicted_2024 = predicted_2024),
  file.path(
    VALID_DIR,
    "07_PRIMARY934_observed_predicted_2024.xlsx"
  )
)


# ==============================================================================
# 21. PRIMARY 934 ANNUAL/REGIONAL DESCRIPTIVE RESULTS
# ==============================================================================

primary_regional_annual <- primary_base %>%
  group_by(region, year) %>%
  summarise(
    n_woredas = n_distinct(id_1082),
    pf_cases = safe_sum(pf_cases),
    pv_cases = safe_sum(pv_cases),
    total_cases = safe_sum(total_cases),
    pf_missing_months = sum(is.na(pf_cases)),
    pv_missing_months = sum(is.na(pv_cases)),
    .groups = "drop"
  ) %>%
  arrange(region, year) %>%
  group_by(region) %>%
  mutate(
    pf_yoy_pct = 100 * (pf_cases / lag(pf_cases) - 1),
    pv_yoy_pct = 100 * (pv_cases / lag(pv_cases) - 1),
    total_yoy_pct = 100 * (total_cases / lag(total_cases) - 1)
  ) %>%
  ungroup()

writexl::write_xlsx(
  list(
    annual_primary934 = primary_annual,
    annual_regional_primary934 = primary_regional_annual,
    annual_sensitivity620 = sensitivity_annual
  ),
  file.path(
    TABLE_DIR,
    "08_PRIMARY934_descriptive_results_for_manuscript.xlsx"
  )
)

