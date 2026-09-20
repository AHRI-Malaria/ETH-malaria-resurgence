################################################################################
# COUNTERFACTUAL analyise
 
# ==============================================================================
# 0. PACKAGES
# ==============================================================================


library(tidyverse)
library(lubridate)
library(INLA)
library(sf)
library(writexl)
library(readxl)

set.seed(20260828)


# ==============================================================================
# 1. PATHS
# ==============================================================================

BASE_DIR <-
  "result_species_only_reanalysis_934_primary"

MODEL_DIR <-
  file.path(
    BASE_DIR,
    "models"
  )

PREP_DIR <-
  file.path(
    BASE_DIR,
    "prepared_data"
  )

OUT_DIR <-
  file.path(
    BASE_DIR,
    "final_after_model_lock"
  )

TABLE_DIR <-
  file.path(
    OUT_DIR,
    "tables"
  )

VALID_DIR <-
  file.path(
    OUT_DIR,
    "validation"
  )

POST_DIR <-
  file.path(
    OUT_DIR,
    "posterior"
  )

purrr::walk(
  c(
    OUT_DIR,
    TABLE_DIR,
    VALID_DIR,
    POST_DIR
  ),
  ~dir.create(
    .x,
    recursive = TRUE,
    showWarnings = FALSE
  )
)


PRIMARY_DATA_PATH <-
  file.path(
    PREP_DIR,
    "PRIMARY_allHF_noTigray_934_model_data.rds"
  )

PRIMARY_GRAPH_PATH <-
  file.path(
    PREP_DIR,
    "woreda_all_HF_no_Tigray_primary.adj"
  )

PF_MODEL_CANDIDATES <- c(
  file.path(
    MODEL_DIR,
    "PRIMARY_allHF_noTigray_Pf_full_934.rds"
  ),
  file.path(
    MODEL_DIR,
    "mod_pf_PRIMARY_934_species_only.rds"
  )
)

PV_MODEL_CANDIDATES <- c(
  file.path(
    MODEL_DIR,
    "PRIMARY_allHF_noTigray_Pv_full_934.rds"
  ),
  file.path(
    MODEL_DIR,
    "mod_pv_PRIMARY_934_species_only.rds"
  )
)

SHAPE_CANDIDATES <- c(
  "data/shp/eth_adm3_clean_id.shp",
  "eth_adm3_clean_id.shp"
)


# ==============================================================================
# 2. RUN SWITCHES
# ==============================================================================

RUN_COUNTERFACTUAL <- TRUE

RUN_GEOGRAPHIC_CV <- TRUE

RUN_FINAL_TABLES <- TRUE

N_POSTERIOR_SAMPLES <- 1000L

POSTERIOR_BATCH_SIZE <- 10L

N_CV_FOLDS <- 5L

CV_SEED <- 20260828L

RESUME_COMPLETED_FILES <- TRUE


# ==============================================================================
# 3. FILE HELPERS
# ==============================================================================

first_existing <- function(
    paths,
    label
) {

  hit <-
    paths[
      file.exists(
        paths
      )
    ]

  if (length(hit) == 0) {
    stop(
      "Could not find ",
      label,
      ". Checked:\n",
      paste(
        paths,
        collapse = "\n"
      )
    )
  }

  hit[1]
}


if (!file.exists(PRIMARY_DATA_PATH)) {
  stop(
    "Primary prepared data not found: ",
    PRIMARY_DATA_PATH
  )
}

if (!file.exists(PRIMARY_GRAPH_PATH)) {
  stop(
    "Primary adjacency graph not found: ",
    PRIMARY_GRAPH_PATH
  )
}


PF_MODEL_PATH <-
  first_existing(
    PF_MODEL_CANDIDATES,
    "saved primary Pf INLA model"
  )

PV_MODEL_PATH <-
  first_existing(
    PV_MODEL_CANDIDATES,
    "saved primary Pv INLA model"
  )

SHAPE_PATH <-
  first_existing(
    SHAPE_CANDIDATES,
    "woreda shapefile"
  )


# ==============================================================================
# 4. LOAD LOCKED PRIMARY DATA / MODELS / GRAPH
# ==============================================================================

primary_dat <-
  readRDS(
    PRIMARY_DATA_PATH
  )

graph_primary <-
  INLA::inla.read.graph(
    PRIMARY_GRAPH_PATH
  )

mod_pf_primary <-
  readRDS(
    PF_MODEL_PATH
  )

mod_pv_primary <-
  readRDS(
    PV_MODEL_PATH
  )


if (!inherits(mod_pf_primary, "inla")) {
  stop(
    "Saved Pf object is not a valid INLA model: ",
    PF_MODEL_PATH
  )
}

if (!inherits(mod_pv_primary, "inla")) {
  stop(
    "Saved Pv object is not a valid INLA model: ",
    PV_MODEL_PATH
  )
}


required_data_cols <- c(
  "id_1082",
  "date",
  "year",
  "month",
  "region",
  "zone",
  "woreda",
  "woreda_population",
  "woreda_idx",
  "time_month",
  "time_linear",
  "season_month",
  "st_id",
  "pf_cases_model",
  "pv_cases_model",
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

missing_data_cols <-
  setdiff(
    required_data_cols,
    names(
      primary_dat
    )
  )

if (length(missing_data_cols) > 0) {
  stop(
    "Prepared primary model data are missing: ",
    paste(
      missing_data_cols,
      collapse = ", "
    )
  )
}


stopifnot(
  n_distinct(
    primary_dat$id_1082
  ) ==
    934
)


message(
  "Primary 934-woreda data loaded: ",
  nrow(primary_dat),
  " rows; ",
  n_distinct(primary_dat$id_1082),
  " woredas."
)


# ==============================================================================
# 5. LOCKED PRIMARY MODEL FORMULA
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
    fixed_terms =
      PRIMARY_FIXED
) {

  rhs_fixed <-
    paste(
      fixed_terms,
      collapse =
        " + "
    )


  rhs_random <-
    paste(
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
# 6. TABLE 1 HELPER — PRIMARY ADJUSTED IRRs
# ==============================================================================

get_fixed_row <- function(
    model,
    exact = NULL,
    prefix = NULL
) {

  sfixed <-
    model$summary.fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      "term"
    )


  if (!is.null(exact)) {

    hit <-
      sfixed %>%
      filter(
        term ==
          exact
      )

    if (nrow(hit) == 1) {
      return(hit)
    }
  }


  if (!is.null(prefix)) {

    hit <-
      sfixed %>%
      filter(
        startsWith(
          term,
          prefix
        )
      )

    if (nrow(hit) == 1) {
      return(hit)
    }
  }


  stop(
    "Could not identify fixed-effect row: ",
    ifelse(
      is.null(exact),
      prefix,
      exact
    )
  )
}


term_dictionary <- tribble(
  ~display_order, ~label, ~exact, ~prefix,
  1L, "Rainfall (2-month lag)", "rainfall_lag2_std", NA_character_,
  2L, "Maximum temperature (2-month lag)", "max_temp_lag2_std", NA_character_,
  3L, "Minimum night-time temperature (2-month lag)", "min_temp_night_lag2_std", NA_character_,
  4L, "Enhanced Vegetation Index (2-month lag)", "evi_lag2_std", NA_character_,
  5L, "Population-weighted elevation", "elev_pop_weighted_std", NA_character_,
  6L, "Internally displaced persons", "idp_ind_std", NA_character_,
  7L, "An. stephensi in woreda or neighbouring woredas", NA_character_, "stephensi_in_neighboring_woredas",
  8L, "LLIN access, within-woreda change", "llin_within_std", NA_character_,
  9L, "IRS coverage, within-woreda change", "irs_within_std", NA_character_
)


extract_primary_table1 <- function(
    model,
    species
) {

  purrr::pmap_dfr(
    term_dictionary,
    function(
        display_order,
        label,
        exact,
        prefix
    ) {

      row <-
        get_fixed_row(
          model,
          exact =
            if (
              is.na(
                exact
              )
            ) {
              NULL
            } else {
              exact
            },
          prefix =
            if (
              is.na(
                prefix
              )
            ) {
              NULL
            } else {
              prefix
            }
        )


      tibble(
        display_order =
          display_order,
        covariate =
          label,
        species =
          species,
        beta =
          row$mean,
        beta_sd =
          row$sd,
        beta_lcl =
          row$`0.025quant`,
        beta_ucl =
          row$`0.975quant`,
        IRR =
          exp(
            row$mean
          ),
        IRR_lcl =
          exp(
            row$`0.025quant`
          ),
        IRR_ucl =
          exp(
            row$`0.975quant`
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
  )
}


table1_long <- bind_rows(
  extract_primary_table1(
    mod_pf_primary,
    "P. falciparum"
  ),
  extract_primary_table1(
    mod_pv_primary,
    "P. vivax"
  )
)


table1_wide <- table1_long %>%
  select(
    display_order,
    covariate,
    species,
    `IRR (95% CrI)`
  ) %>%
  pivot_wider(
    names_from =
      species,
    values_from =
      `IRR (95% CrI)`
  ) %>%
  arrange(
    display_order
  ) %>%
  select(
    -display_order
  )


# ==============================================================================
# 7. POSTERIOR COUNTERFACTUAL — CONFIG=TRUE MODEL
# ==============================================================================

fit_scenario_config_model <- function(
    outcome,
    model_name
) {

  model_path <-
    file.path(
      POST_DIR,
      paste0(
        model_name,
        "_configTRUE.rds"
      )
    )


  if (
    RESUME_COMPLETED_FILES &&
    file.exists(
      model_path
    )
  ) {

    message(
      "Loading completed config=TRUE model: ",
      model_path
    )

    fit <-
      readRDS(
        model_path
      )

    if (!inherits(fit, "inla")) {
      stop(
        "Saved scenario model is not a valid INLA object: ",
        model_path
      )
    }

    return(
      fit
    )
  }


  d <-
    primary_dat %>%
    mutate(
      .outcome =
        as.numeric(
          .data[[
            outcome
          ]]
        )
    )


  graph_obj <-
    graph_primary


  fml <-
    make_formula(
      PRIMARY_FIXED
    )


  environment(
    fml
  ) <-
    environment()


  message(
    "\n============================================================\n",
    "Fitting posterior-sampling model: ",
    model_name,
    "\n============================================================"
  )


  fit <-
    INLA::inla(
      formula =
        fml,

      family =
        "nbinomial",

      data =
        d,

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
            FALSE,
          dic =
            FALSE,
          cpo =
            FALSE,
          config =
            TRUE
        ),

      verbose =
        FALSE
    )


  if (!inherits(fit, "inla")) {
    stop(
      model_name,
      ": config=TRUE fit did not return a valid INLA object."
    )
  }


  saveRDS(
    fit,
    model_path
  )


  fit
}


latent_vector <- function(
    sample_object
) {

  x <-
    as.numeric(
      sample_object$latent
    )

  names(x) <-
    rownames(
      sample_object$latent
    )

  x
}


get_latent_scalar <- function(
    lat,
    exact_term = NULL,
    prefix = NULL
) {

  nm <-
    names(
      lat
    )


  if (!is.null(exact_term)) {

    candidates <- c(
      exact_term,
      paste0(
        exact_term,
        ":1"
      )
    )

    hit <-
      candidates[
        candidates %in%
          nm
      ]

    if (length(hit) == 1) {
      return(
        as.numeric(
          lat[
            hit
          ]
        )
      )
    }
  }


  if (!is.null(prefix)) {

    hit <-
      nm[
        startsWith(
          nm,
          prefix
        )
      ]


    hit <-
      hit[
        !grepl(
          "Theta|Log precision|Precision",
          hit,
          ignore.case =
            TRUE
        )
      ]


    if (length(hit) >= 1) {

      hit1 <-
        hit[
          grepl(
            ":1$",
            hit
          )
        ]


      if (length(hit1) == 1) {
        return(
          as.numeric(
            lat[
              hit1
            ]
          )
        )
      }


      if (length(hit) == 1) {
        return(
          as.numeric(
            lat[
              hit
            ]
          )
        )
      }
    }
  }


  stop(
    "Could not identify posterior coefficient: ",
    ifelse(
      is.null(
        exact_term
      ),
      prefix,
      exact_term
    )
  )
}


get_predictor_sample <- function(
    lat,
    n_expected
) {

  nm <-
    names(
      lat
    )[
      grepl(
        "^Predictor:",
        names(
          lat
        )
      )
    ]


  if (length(nm) != n_expected) {
    stop(
      "Posterior draw contains ",
      length(nm),
      " Predictor entries; expected ",
      n_expected,
      "."
    )
  }


  index <-
    suppressWarnings(
      as.integer(
        sub(
          "^Predictor:",
          "",
          nm
        )
      )
    )


  if (anyNA(index)) {
    stop(
      "Could not parse posterior Predictor indices."
    )
  }


  nm <-
    nm[
      order(
        index
      )
    ]


  as.numeric(
    lat[
      nm
    ]
  )
}


scenario_reference <- function(
    data,
    ref_years
) {

  ref <-
    data %>%
    filter(
      year %in%
        ref_years
    )


  steph <-
    as.integer(
      as.character(
        ref$stephensi_in_neighboring_woredas
      )
    )


  list(
    rainfall =
      mean(
        ref$rainfall_lag2_std,
        na.rm =
          TRUE
      ),

    max_temp =
      mean(
        ref$max_temp_lag2_std,
        na.rm =
          TRUE
      ),

    min_temp =
      mean(
        ref$min_temp_night_lag2_std,
        na.rm =
          TRUE
      ),

    evi =
      mean(
        ref$evi_lag2_std,
        na.rm =
          TRUE
      ),

    idp =
      mean(
        ref$idp_ind_std,
        na.rm =
          TRUE
      ),

    llin =
      mean(
        ref$llin_within_std,
        na.rm =
          TRUE
      ),

    irs =
      mean(
        ref$irs_within_std,
        na.rm =
          TRUE
      ),

    stephensi =
      as.integer(
        mean(
          steph,
          na.rm =
            TRUE
        ) >=
          0.5
      )
  )
}


process_scenario_draw <- function(
    sample_object,
    model_data,
    outcome_label,
    ref_years,
    ref_label,
    draw_id
) {

  lat <-
    latent_vector(
      sample_object
    )


  eta_obs <-
    get_predictor_sample(
      lat,
      n_expected =
        nrow(
          model_data
        )
    )


  b_rain <-
    get_latent_scalar(
      lat,
      exact_term =
        "rainfall_lag2_std"
    )

  b_tmax <-
    get_latent_scalar(
      lat,
      exact_term =
        "max_temp_lag2_std"
    )

  b_tmin <-
    get_latent_scalar(
      lat,
      exact_term =
        "min_temp_night_lag2_std"
    )

  b_evi <-
    get_latent_scalar(
      lat,
      exact_term =
        "evi_lag2_std"
    )

  b_idp <-
    get_latent_scalar(
      lat,
      exact_term =
        "idp_ind_std"
    )

  b_llin <-
    get_latent_scalar(
      lat,
      exact_term =
        "llin_within_std"
    )

  b_irs <-
    get_latent_scalar(
      lat,
      exact_term =
        "irs_within_std"
    )

  b_steph <-
    get_latent_scalar(
      lat,
      prefix =
        "stephensi_in_neighboring_woredas"
    )


  ref <-
    scenario_reference(
      model_data,
      ref_years
    )


  idx_future <-
    which(
      model_data$year %in%
        2022:2024
    )


  x_steph <-
    as.integer(
      as.character(
        model_data$stephensi_in_neighboring_woredas
      )
    )


  eta_rain <-
    eta_obs

  eta_rain[
    idx_future
  ] <-
    eta_obs[
      idx_future
    ] +
    b_rain *
    (
      ref$rainfall -
      model_data$rainfall_lag2_std[
        idx_future
      ]
    )


  eta_temp <-
    eta_obs

  eta_temp[
    idx_future
  ] <-
    eta_obs[
      idx_future
    ] +
    b_tmax *
    (
      ref$max_temp -
      model_data$max_temp_lag2_std[
        idx_future
      ]
    ) +
    b_tmin *
    (
      ref$min_temp -
      model_data$min_temp_night_lag2_std[
        idx_future
      ]
    )


  eta_evi <-
    eta_obs

  eta_evi[
    idx_future
  ] <-
    eta_obs[
      idx_future
    ] +
    b_evi *
    (
      ref$evi -
      model_data$evi_lag2_std[
        idx_future
      ]
    )


  eta_idp <-
    eta_obs

  eta_idp[
    idx_future
  ] <-
    eta_obs[
      idx_future
    ] +
    b_idp *
    (
      ref$idp -
      model_data$idp_ind_std[
        idx_future
      ]
    )


  eta_steph <-
    eta_obs

  eta_steph[
    idx_future
  ] <-
    eta_obs[
      idx_future
    ] +
    b_steph *
    (
      ref$stephensi -
      x_steph[
        idx_future
      ]
    )


  eta_programme <-
    eta_obs

  eta_programme[
    idx_future
  ] <-
    eta_obs[
      idx_future
    ] +
    b_llin *
    (
      ref$llin -
      model_data$llin_within_std[
        idx_future
      ]
    ) +
    b_irs *
    (
      ref$irs -
      model_data$irs_within_std[
        idx_future
      ]
    )


  eta_all <-
    eta_obs

  eta_all[
    idx_future
  ] <-
    eta_obs[
      idx_future
    ] +
    b_rain *
    (
      ref$rainfall -
      model_data$rainfall_lag2_std[
        idx_future
      ]
    ) +
    b_tmax *
    (
      ref$max_temp -
      model_data$max_temp_lag2_std[
        idx_future
      ]
    ) +
    b_tmin *
    (
      ref$min_temp -
      model_data$min_temp_night_lag2_std[
        idx_future
      ]
    ) +
    b_evi *
    (
      ref$evi -
      model_data$evi_lag2_std[
        idx_future
      ]
    ) +
    b_idp *
    (
      ref$idp -
      model_data$idp_ind_std[
        idx_future
      ]
    ) +
    b_steph *
    (
      ref$stephensi -
      x_steph[
        idx_future
      ]
    ) +
    b_llin *
    (
      ref$llin -
      model_data$llin_within_std[
        idx_future
      ]
    ) +
    b_irs *
    (
      ref$irs -
      model_data$irs_within_std[
        idx_future
      ]
    )


  predicted_2224 <-
    sum(
      exp(
        eta_obs[
          idx_future
        ]
      )
    )


  ref_modelled_by_year <-
    purrr::map_dbl(
      ref_years,
      function(
          yr
      ) {

        sum(
          exp(
            eta_obs[
              model_data$year ==
                yr
            ]
          )
        )
      }
    )


  reference_3yr <-
    mean(
      ref_modelled_by_year
    ) *
    3


  modelled_increase <-
    predicted_2224 -
    reference_3yr


  counterfactual <- c(
    Rainfall =
      sum(
        exp(
          eta_rain[
            idx_future
          ]
        )
      ),

    Temperature =
      sum(
        exp(
          eta_temp[
            idx_future
          ]
        )
      ),

    `Enhanced Vegetation Index` =
      sum(
        exp(
          eta_evi[
            idx_future
          ]
        )
      ),

    `Internally displaced persons` =
      sum(
        exp(
          eta_idp[
            idx_future
          ]
        )
      ),

    `An. stephensi` =
      sum(
        exp(
          eta_steph[
            idx_future
          ]
        )
      ),

    `LLIN and IRS` =
      sum(
        exp(
          eta_programme[
            idx_future
          ]
        )
      ),

    `All examined time-varying drivers` =
      sum(
        exp(
          eta_all[
            idx_future
          ]
        )
      )
  )


  tibble(
    outcome =
      outcome_label,

    reference =
      ref_label,

    draw =
      draw_id,

    factor =
      names(
        counterfactual
      ),

    predicted_cases_2022_2024 =
      predicted_2224,

    modelled_reference_cases_3yr =
      reference_3yr,

    modelled_increase =
      modelled_increase,

    counterfactual_cases_2022_2024 =
      as.numeric(
        counterfactual
      ),

    model_based_case_difference =
      predicted_cases_2022_2024 -
      counterfactual_cases_2022_2024,

    contribution_to_modelled_increase_pct =
      ifelse(
        abs(
          modelled_increase
        ) <
          1e-12,
        NA_real_,
        100 *
          model_based_case_difference /
          modelled_increase
      )
  )
}


run_posterior_scenarios <- function(
    model,
    model_data,
    outcome_label,
    n_samples =
      N_POSTERIOR_SAMPLES,
    batch_size =
      POSTERIOR_BATCH_SIZE
) {

  refs <- list(
    "2020" =
      2020,

    "2021" =
      2021,

    "2020–2021 mean" =
      c(
        2020,
        2021
      )
  )


  output_chunks <- list()

  draw_counter <- 0L

  chunk_counter <- 0L


  while (
    draw_counter <
      n_samples
  ) {

    n_this <-
      min(
        batch_size,
        n_samples -
          draw_counter
      )


    chunk_counter <-
      chunk_counter +
      1L


    message(
      outcome_label,
      ": posterior draws ",
      draw_counter + 1L,
      "–",
      draw_counter + n_this,
      " of ",
      n_samples
    )


    samples <-
      INLA::inla.posterior.sample(
        n =
          n_this,

        result =
          model,

        seed =
          20260828L +
          chunk_counter
      )


    chunk_result <-
      purrr::map_dfr(
        seq_along(
          samples
        ),
        function(
            j
        ) {

          current_draw <-
            draw_counter +
            j


          purrr::imap_dfr(
            refs,
            function(
                years_ref,
                label_ref
            ) {

              process_scenario_draw(
                samples[[
                  j
                ]],
                model_data,
                outcome_label,
                years_ref,
                label_ref,
                current_draw
              )
            }
          )
        }
      )


    output_chunks[[chunk_counter]] <-
      chunk_result


    draw_counter <-
      draw_counter +
      n_this


    rm(
      samples
    )

    invisible(
      gc()
    )
  }


  bind_rows(
    output_chunks
  )
}


summarise_counterfactual <- function(
    draws
) {

  draws %>%
    group_by(
      outcome,
      reference,
      factor
    ) %>%
    summarise(
      predicted_cases =
        mean(
          predicted_cases_2022_2024,
          na.rm =
            TRUE
        ),

      predicted_cases_lcl =
        quantile(
          predicted_cases_2022_2024,
          0.025,
          na.rm =
            TRUE
        ),

      predicted_cases_ucl =
        quantile(
          predicted_cases_2022_2024,
          0.975,
          na.rm =
            TRUE
        ),

      counterfactual_cases =
        mean(
          counterfactual_cases_2022_2024,
          na.rm =
            TRUE
        ),

      counterfactual_lcl =
        quantile(
          counterfactual_cases_2022_2024,
          0.025,
          na.rm =
            TRUE
        ),

      counterfactual_ucl =
        quantile(
          counterfactual_cases_2022_2024,
          0.975,
          na.rm =
            TRUE
        ),

      case_difference =
        mean(
          model_based_case_difference,
          na.rm =
            TRUE
        ),

      case_difference_lcl =
        quantile(
          model_based_case_difference,
          0.025,
          na.rm =
            TRUE
        ),

      case_difference_ucl =
        quantile(
          model_based_case_difference,
          0.975,
          na.rm =
            TRUE
        ),

      contribution_pct =
        mean(
          contribution_to_modelled_increase_pct,
          na.rm =
            TRUE
        ),

      contribution_lcl =
        quantile(
          contribution_to_modelled_increase_pct,
          0.025,
          na.rm =
            TRUE
        ),

      contribution_ucl =
        quantile(
          contribution_to_modelled_increase_pct,
          0.975,
          na.rm =
            TRUE
        ),

      `Contribution % (95% CrI)` =
        sprintf(
          "%.1f (%.1f–%.1f)",
          contribution_pct,
          contribution_lcl,
          contribution_ucl
        ),

      .groups =
        "drop"
    )
}


COUNTERFACTUAL_XLSX <-
  file.path(
    TABLE_DIR,
    "11_FINAL_posterior_counterfactual_PRIMARY934.xlsx"
  )


if (RUN_COUNTERFACTUAL) {

  pf_draw_path <-
    file.path(
      POST_DIR,
      "pf_counterfactual_draws.csv.gz"
    )

  pv_draw_path <-
    file.path(
      POST_DIR,
      "pv_counterfactual_draws.csv.gz"
    )


  if (
    RESUME_COMPLETED_FILES &&
    file.exists(
      pf_draw_path
    )
  ) {

    message(
      "Loading completed Pf posterior draws."
    )

    pf_counterfactual_draws <-
      readr::read_csv(
        pf_draw_path,
        show_col_types =
          FALSE
      )

  } else {

    mod_pf_scenario <-
      fit_scenario_config_model(
        "pf_cases_model",
        "PRIMARY934_scenario_Pf"
      )


    pf_counterfactual_draws <-
      run_posterior_scenarios(
        mod_pf_scenario,
        primary_dat,
        "P. falciparum"
      )


    readr::write_csv(
      pf_counterfactual_draws,
      pf_draw_path
    )


    rm(
      mod_pf_scenario
    )

    invisible(
      gc()
    )
  }


  if (
    RESUME_COMPLETED_FILES &&
    file.exists(
      pv_draw_path
    )
  ) {

    message(
      "Loading completed Pv posterior draws."
    )

    pv_counterfactual_draws <-
      readr::read_csv(
        pv_draw_path,
        show_col_types =
          FALSE
      )

  } else {

    mod_pv_scenario <-
      fit_scenario_config_model(
        "pv_cases_model",
        "PRIMARY934_scenario_Pv"
      )


    pv_counterfactual_draws <-
      run_posterior_scenarios(
        mod_pv_scenario,
        primary_dat,
        "P. vivax"
      )


    readr::write_csv(
      pv_counterfactual_draws,
      pv_draw_path
    )


    rm(
      mod_pv_scenario
    )

    invisible(
      gc()
    )
  }


  counterfactual_draws <-
    bind_rows(
      pf_counterfactual_draws,
      pv_counterfactual_draws
    )


  counterfactual_summary <-
    summarise_counterfactual(
      counterfactual_draws
    )


  writexl::write_xlsx(
    list(
      posterior_summary =
        counterfactual_summary,
      main_reference_2020_2021 =
        counterfactual_summary %>%
        filter(
          reference ==
            "2020–2021 mean"
        ),
      reference_sensitivity =
        counterfactual_summary %>%
        filter(
          reference %in%
            c(
              "2020",
              "2021"
            )
        )
    ),
    COUNTERFACTUAL_XLSX
  )

} else {

  if (!file.exists(COUNTERFACTUAL_XLSX)) {
    stop(
      "RUN_COUNTERFACTUAL=FALSE but counterfactual workbook does not exist."
    )
  }

  counterfactual_summary <-
    readxl::read_xlsx(
      COUNTERFACTUAL_XLSX,
      sheet =
        "posterior_summary"
    )
}


# ==============================================================================
# 8. CORRECTED GEOGRAPHIC 5-FOLD CROSS-VALIDATION
# ==============================================================================

make_primary934_spatial_sf <- function() {

  shp <-
    sf::st_read(
      SHAPE_PATH,
      quiet =
        TRUE
    ) %>%
    sf::st_make_valid() %>%
    mutate(
      id_1082 =
        as.integer(
          id_1082
        )
    )


  if (!"id_1082" %in% names(shp)) {
    stop(
      "Shapefile is missing id_1082."
    )
  }


  lookup <-
    primary_dat %>%
    distinct(
      id_1082,
      woreda_idx
    )


  sf_primary <-
    shp %>%
    filter(
      id_1082 %in%
        lookup$id_1082
    ) %>%
    select(
      id_1082,
      geometry
    ) %>%
    left_join(
      lookup,
      by =
        "id_1082"
    ) %>%
    arrange(
      woreda_idx
    )


  if (
    nrow(
      sf_primary
    ) !=
      934 ||
    anyNA(
      sf_primary$woreda_idx
    )
  ) {

    stop(
      "Could not reconstruct all 934 primary spatial polygons."
    )
  }


  sf_primary
}


woredas_primary934_sf <-
  make_primary934_spatial_sf()


make_geographic_folds <- function(
    sf_woredas,
    k =
      5L,
    seed =
      20260828L
) {

  set.seed(
    seed
  )


  cent <-
    sf::st_point_on_surface(
      sf::st_transform(
        sf_woredas,
        3857
      )
    )


  xy <-
    sf::st_coordinates(
      cent
    )


  km <-
    stats::kmeans(
      scale(
        xy
      ),
      centers =
        k,
      nstart =
        100,
      iter.max =
        100
    )


  tibble(
    woreda_idx =
      sf_woredas$woreda_idx,
    id_1082 =
      sf_woredas$id_1082,
    cv_fold =
      as.integer(
        km$cluster
      )
  )
}


cv_folds <-
  make_geographic_folds(
    woredas_primary934_sf,
    k =
      N_CV_FOLDS,
    seed =
      CV_SEED
  )


readr::write_csv(
  cv_folds,
  file.path(
    VALID_DIR,
    "10_PRIMARY934_geographic_CV_fold_allocation.csv"
  )
)


dat_cv <-
  primary_dat %>%
  select(
    -any_of(
      "cv_fold"
    )
  ) %>%
  left_join(
    cv_folds %>%
      select(
        woreda_idx,
        cv_fold
      ),
    by =
      "woreda_idx"
  )


if (anyNA(dat_cv$cv_fold)) {
  stop(
    "Some primary records were not allocated to a geographic CV fold."
  )
}


fit_one_cv_fold <- function(
    outcome,
    fold,
    species_prefix
) {

  pred_file <-
    file.path(
      VALID_DIR,
      paste0(
        "CV_",
        species_prefix,
        "_fold",
        fold,
        "_predictions.csv.gz"
      )
    )


  if (
    RESUME_COMPLETED_FILES &&
    file.exists(
      pred_file
    )
  ) {

    message(
      "Loading completed CV ",
      species_prefix,
      " fold ",
      fold
    )


    return(
      readr::read_csv(
        pred_file,
        show_col_types =
          FALSE
      )
    )
  }


  message(
    "\n============================================================\n",
    "Geographic CV ",
    species_prefix,
    " fold ",
    fold,
    "/",
    N_CV_FOLDS,
    "\n============================================================"
  )


  d <-
    dat_cv %>%
    mutate(
      .outcome =
        as.numeric(
          .data[[
            outcome
          ]]
        ),

      observed_cv =
        as.numeric(
          .data[[
            outcome
          ]]
        )
    )


  heldout <-
    d$cv_fold ==
      fold


  # Remove outcomes from the held-out geographic block. Primary missing outcomes
  # remain NA as well; no missing species report is converted to zero.
  d$.outcome[
    heldout
  ] <-
    NA_real_


  # IMPORTANT:
  # The IID woreda-month interaction cannot be estimated for an entirely held-out
  # observation. Setting st_id to NA removes this unpredictable residual term
  # from the held-out point prediction.
  d$st_id[
    heldout
  ] <-
    NA_integer_


  graph_obj <-
    graph_primary


  fml <-
    make_formula(
      PRIMARY_FIXED
    )


  environment(
    fml
  ) <-
    environment()


  fit <-
    INLA::inla(
      formula =
        fml,

      family =
        "nbinomial",

      data =
        d,

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
            FALSE,
          dic =
            FALSE,
          cpo =
            FALSE,
          config =
            FALSE
        ),

      verbose =
        FALSE
    )


  fv <-
    fit$summary.fitted.values %>%
    as.data.frame()


  if (!"0.5quant" %in% names(fv)) {
    stop(
      "INLA fitted values do not contain posterior median `0.5quant`."
    )
  }


  prediction_median <-
    as.numeric(
      fv[["0.5quant"]]
    )


  prediction_mean <-
    if (
      "mean" %in%
        names(
          fv
        )
    ) {
      as.numeric(
        fv[["mean"]]
      )
    } else {
      rep(
        NA_real_,
        nrow(
          d
        )
      )
    }


  ans <-
    d %>%
    transmute(
      id_1082,
      date,
      region,
      zone,
      woreda,
      woreda_idx,
      cv_fold,
      observed =
        observed_cv,
      predicted_median =
        prediction_median,
      predicted_mean =
        prediction_mean
    ) %>%
    filter(
      cv_fold ==
        fold
    )


  readr::write_csv(
    ans,
    pred_file
  )


  rm(
    fit,
    d,
    fv,
    prediction_median,
    prediction_mean
  )


  invisible(
    gc()
  )


  ans
}


cv_metrics <- function(
    data,
    species_label
) {

  d <-
    data %>%
    filter(
      is.finite(
        observed
      ),
      is.finite(
        predicted_median
      ),
      predicted_median >=
        0
    )


  fold_metrics <-
    d %>%
    group_by(
      cv_fold
    ) %>%
    summarise(
      species =
        species_label,
      n =
        n(),
      n_woredas =
        n_distinct(
          woreda_idx
        ),
      observed_cases =
        sum(
          observed
        ),
      predicted_cases =
        sum(
          predicted_median
        ),
      O_to_P =
        observed_cases /
        predicted_cases,
      MAE =
        mean(
          abs(
            observed -
              predicted_median
          )
        ),
      RMSE =
        sqrt(
          mean(
            (
              observed -
                predicted_median
            )^2
          )
        ),
      log1p_RMSE =
        sqrt(
          mean(
            (
              log1p(
                observed
              ) -
                log1p(
                  predicted_median
                )
            )^2
          )
        ),
      Spearman_r =
        suppressWarnings(
          cor(
            observed,
            predicted_median,
            method =
              "spearman"
          )
        ),
      .groups =
        "drop"
    )


  overall <-
    d %>%
    summarise(
      species =
        species_label,
      n =
        n(),
      n_woredas =
        n_distinct(
          woreda_idx
        ),
      observed_cases =
        sum(
          observed
        ),
      predicted_cases =
        sum(
          predicted_median
        ),
      O_to_P =
        observed_cases /
        predicted_cases,
      MAE =
        mean(
          abs(
            observed -
              predicted_median
          )
        ),
      RMSE =
        sqrt(
          mean(
            (
              observed -
                predicted_median
            )^2
          )
        ),
      log1p_RMSE =
        sqrt(
          mean(
            (
              log1p(
                observed
              ) -
                log1p(
                  predicted_median
                )
            )^2
          )
        ),
      Spearman_r =
        suppressWarnings(
          cor(
            observed,
            predicted_median,
            method =
              "spearman"
          )
        )
    )


  list(
    fold =
      fold_metrics,
    overall =
      overall
  )
}


calibration_deciles <- function(
    data,
    species_label
) {

  data %>%
    filter(
      is.finite(
        observed
      ),
      is.finite(
        predicted_median
      )
    ) %>%
    mutate(
      prediction_decile =
        ntile(
          predicted_median,
          10
        )
    ) %>%
    group_by(
      prediction_decile
    ) %>%
    summarise(
      species =
        species_label,
      n =
        n(),
      observed_mean =
        mean(
          observed
        ),
      predicted_mean =
        mean(
          predicted_median
        ),
      observed_sum =
        sum(
          observed
        ),
      predicted_sum =
        sum(
          predicted_median
        ),
      O_to_P =
        observed_sum /
        predicted_sum,
      .groups =
        "drop"
    )
}


CV_XLSX <-
  file.path(
    VALID_DIR,
    "10_FINAL_geographic_5fold_cross_validation_PRIMARY934.xlsx"
  )


if (RUN_GEOGRAPHIC_CV) {

  cv_pf <-
    purrr::map_dfr(
      seq_len(
        N_CV_FOLDS
      ),
      ~fit_one_cv_fold(
        outcome =
          "pf_cases_model",
        fold =
          .x,
        species_prefix =
          "Pf"
      )
    )


  cv_pv <-
    purrr::map_dfr(
      seq_len(
        N_CV_FOLDS
      ),
      ~fit_one_cv_fold(
        outcome =
          "pv_cases_model",
        fold =
          .x,
        species_prefix =
          "Pv"
      )
    )


  m_pf <-
    cv_metrics(
      cv_pf,
      "P. falciparum"
    )


  m_pv <-
    cv_metrics(
      cv_pv,
      "P. vivax"
    )


  cv_calibration <-
    bind_rows(
      calibration_deciles(
        cv_pf,
        "P. falciparum"
      ),
      calibration_deciles(
        cv_pv,
        "P. vivax"
      )
    )


  writexl::write_xlsx(
    list(
      overall_metrics =
        bind_rows(
          m_pf$overall,
          m_pv$overall
        ),

      fold_metrics =
        bind_rows(
          m_pf$fold,
          m_pv$fold
        ),

      calibration_deciles =
        cv_calibration,

      fold_allocation =
        cv_folds
    ),
    CV_XLSX
  )
}


# ==============================================================================
# 9. FINAL TABLE 2 FROM MEAN 2020–2021 REFERENCE
# ==============================================================================

factor_order <- c(
  "Rainfall",
  "Temperature",
  "Enhanced Vegetation Index",
  "Internally displaced persons",
  "An. stephensi",
  "LLIN and IRS",
  "All examined time-varying drivers"
)


table2_long <-
  counterfactual_summary %>%
  filter(
    reference ==
      "2020–2021 mean"
  ) %>%
  mutate(
    factor =
      factor(
        factor,
        levels =
          factor_order
      ),

    `Model-based case difference, thousands (95% CrI)` =
      sprintf(
        "%.1f (%.1f–%.1f)",
        case_difference /
          1000,
        case_difference_lcl /
          1000,
        case_difference_ucl /
          1000
      ),

    `Contribution to modelled increase, % (95% CrI)` =
      sprintf(
        "%.1f (%.1f–%.1f)",
        contribution_pct,
        contribution_lcl,
        contribution_ucl
      )
  ) %>%
  arrange(
    factor,
    outcome
  )


table2_wide <-
  table2_long %>%
  transmute(
    factor =
      as.character(
        factor
      ),

    outcome,

    contribution =
      `Contribution to modelled increase, % (95% CrI)`,

    case_difference =
      `Model-based case difference, thousands (95% CrI)`
  ) %>%
  pivot_wider(
    names_from =
      outcome,
    values_from =
      c(
        contribution,
        case_difference
      )
  ) %>%
  mutate(
    factor =
      factor(
        factor,
        levels =
          factor_order
      )
  ) %>%
  arrange(
    factor
  ) %>%
  mutate(
    factor =
      as.character(
        factor
      )
  )


# ==============================================================================
# 10. PRIMARY DESCRIPTIVE ANNUAL TOTALS FOR FINAL MANUSCRIPT CHECK
# ==============================================================================

annual_primary <-
  primary_dat %>%
  group_by(
    year
  ) %>%
  summarise(
    n_woredas =
      n_distinct(
        id_1082
      ),

    Pf_cases =
      sum(
        pf_cases,
        na.rm =
          TRUE
      ),

    Pv_cases =
      sum(
        pv_cases,
        na.rm =
          TRUE
      ),

    total_cases =
      sum(
        total_cases,
        na.rm =
          TRUE
      ),

    .groups =
      "drop"
  ) %>%
  arrange(
    year
  ) %>%
  mutate(
    Pf_change_from_2020_pct =
      100 *
      (
        Pf_cases /
          first(
            Pf_cases
          ) -
          1
      ),

    Pv_change_from_2020_pct =
      100 *
      (
        Pv_cases /
          first(
            Pv_cases
          ) -
          1
      ),

    total_change_from_2020_pct =
      100 *
      (
        total_cases /
          first(
            total_cases
          ) -
          1
      )
  )


# ==============================================================================
# 11. WRITE FINAL MANUSCRIPT SOURCE WORKBOOK
# ==============================================================================

FINAL_TABLES_XLSX <-
  file.path(
    TABLE_DIR,
    "12_UPDATED_MAIN_TABLES_PRIMARY934.xlsx"
  )


if (RUN_FINAL_TABLES) {

  writexl::write_xlsx(
    list(
      Table1_display =
        table1_wide,

      Table1_numeric =
        table1_long,

      Table2_display =
        table2_wide,

      Table2_numeric =
        table2_long,

      annual_primary =
        annual_primary,

      counterfactual_reference_sensitivity =
        counterfactual_summary %>%
        arrange(
          outcome,
          factor,
          reference
        )
    ),
    FINAL_TABLES_XLSX
  )
}


# ==============================================================================
# 12. SAVE SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file =
    file.path(
      OUT_DIR,
      "sessionInfo_final_after_model_lock.txt"
    )
)


