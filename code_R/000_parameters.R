START_YEAR <- 2000
END_YEAR   <- 2040

config <- list()
config$run_days_lost       <- FALSE
config$run_days_lost_lever <- FALSE
config$run_projection      <- FALSE
config$lever_change_start_at      <- 1860

if(FALSE)
{
  config$run_days_lost       <- TRUE
  config$run_days_lost_lever <- TRUE

}

config$run_projection      <- TRUE

scenarios <- list(    pediatric_incidence_plus_10_perc = TRUE
                    , pediatric_incidence_minus_10_perc = TRUE
                    , pediatric_incidence_plus_25_perc = TRUE
                    , pediatric_incidence_minus_25_perc = TRUE

                    , adult_incidence_all_studies = TRUE
                    # , adult_incidence_17_studies = TRUE

                    , iot_no_growth = TRUE
                    , iot_global_curve_for_all = TRUE

                    , smr_age_cruve_flat_average = TRUE

                    , smr_plus_10_perc  = TRUE
                    , smr_minus_10_perc = TRUE

                    , diagnosis_rate_plus_25_pp = TRUE
                    , diagnosis_rate_minus_25_pp = TRUE
                    , adult_onset_zero = TRUE
                    , diagnosis_rate_left = TRUE
                    , diagnosis_rate_right = TRUE
                    , lever_change_at_2022 = TRUE

)
scenarios[] <- FALSE
