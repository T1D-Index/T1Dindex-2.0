# Mortality calculation functions
#
# Confidential
# Copyright (c) JDRF 2020, All rights reserved
#
# This file calculates mortality rates for individuals in the model. Mortality
# is higher for people with T1D. We express this using standardized mortality
# ratios (SMRs) applied to the mortality rates given by UNWPP life tables.
#
# Currently, we use an SMR relationship derived from EDC study data, which gives
# an upward-sloping, exponential relationship between HbA1c and SMR. Assumed
# HbA1c for each country/year is applied to obtain an SMR, which is then applied
# to the background mortality rate from the UNWPP.

#' Get function that returns mortality by age for years
#'
#' This is background mortality, formally the probability of death at age x,
#' $$q_x$$ in life tables, for the given country. Data are as provided by
#' UN Population Prospects.
#'
#' @param country_wb_name (scalar character value) World Bank name for country
#' @return function of year that returns mortality rates (qx) for ages 0-99
#' @export
#' @references DESA, UN. "World Population Prospects 2019: Highlights."
#'             New York (US): United Nations Department for Economic and
#'             Social Affairs (2019).
#' @importFrom rlang .data
#' @importFrom DBI dbSendStatement dbFetch
#' @importFrom dplyr filter select mutate
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map_dfr simplify set_names
background_mortality_function <- function(country_wb_name,yr) {

  if(FALSE)
  {


    loc_id <- get_loc_id(country_wb_name)
    con <- get_database_connection(read_only=TRUE)
    all_data <- dbGetQuery(
      con,
      'SELECT start_year, age_grp_start, qx
      FROM life_table
      WHERE loc_id = :loc_id
      ORDER BY start_year;',
      list(loc_id=loc_id)
    )
    dbDisconnect(con)
    # Returns a vector of age-specific background mortality rates for the given year
    draw <- 1

    clamp_yr <- pmin(pmax(yr, 1950), 2040)
    yr5_low <- 5*floor(min(clamp_yr)/5)
    yr5_high <- 5*ceiling(max(clamp_yr)/5)
    res <- filter(all_data, start_year >= yr5_low & start_year <= yr5_high)

    divided <- split(res$qx, res$start_year)
    qx_for <- function(y) {
      base_yr <- 5*floor(y/5)

        weight1=base_yr + 5-y
        qx5_1=divided[[as.character(base_yr)]]
        qx5_2=divided[[as.character(base_yr + 5)]]

        #' Interpolate 5x5 year life table qx to single year
        #'
        #' Recall that in 5x5 year data, qx is given for 0 year olds,
        #' 1-4 year olds, then 5-9, 10-14, ..., 95-99. We trim out 100+
        #' at the data loading stage.
        #'
        #' @param weight1 weight of pop5_1 time five
        #' @param qx5_1 vector of length 21, 5-year-binned qx with weight weight1/5
        #' @param qx5_2 vector of length 21, 5-year-binned qx with weight (5-weight1)/5. If weight1 is 5 this parameter ignored.
        #' @return 100-vector of single-year qx
        stopifnot(weight1 <= 5 & weight1 > 0 & length(qx5_1) == 21)
        # 1 year, 4 year, then 5 year bucketed qxs, but interpolated to match target year
        if (weight1 == 5) {
          qx5 <- qx5_1
        } else {
          stopifnot(length(qx5_2) == 21)
          qx5 <- weight1/5*qx5_1 + (5-weight1)/5*qx5_2
        }
        # linearly interpolated to single-year ages
        qx <- qx5[seq(2,21.8,0.2)]
        qx[1] <- qx5[1]

        # convert 4- and 5-year bucket qx probabilities to 1-year probabilities
        qx[2:5] <- (1+qx[2:5])^(1/4)-1
        qx[-(1:5)] <- (1+qx[-(1:5)])^(1/5)-1

        set_names(qx, AGES)

    }

    res <- t(qx_for( y=set_names(clamp_yr, yr) ))%>% as.matrix

    # res <- map_dfr(set_names(clamp_yr, yr), qx_for) %>% as.matrix
    rownames(res) <- yr
    res


  }

  if(TRUE)
  {


    loc_id <- get_loc_id(country_wb_name)
    con <- get_database_connection(read_only=TRUE)
    all_data <- dbGetQuery(
      con,
      'SELECT start_year, age_grp_start, qB
      FROM life_table_single_year
      WHERE loc_id = :loc_id
      ORDER BY start_year ,age_grp_start asc ;',
      list(loc_id=loc_id)


    )
    dbDisconnect(con)

    res <- filter(all_data, start_year >= min(yr) & start_year <= max(yr))

    data_wide_age           <- (spread(res, age_grp_start, qB))
    rownames(data_wide_age) <- data_wide_age$start_year
    res <- select(data_wide_age,-start_year)%>% as.matrix
  }

res


}

#' Probabilistic standardised T1D mortality ratio
#'
#' @param country_wb_name (scalar character value) World Bank name for country.
#' @param log_adjust (function) Adjustment made to link function for mortality
#'   probability. Creat function using `t1d_log_adjust`.
#' @return function of year that yields $$q_x$$ vectors for ages 0..99 and for
#'   each draw.
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate as_tibble rename
#' @importFrom tidyr pivot_longer spread
#' @importFrom purrr map_dfr set_names
#' @importFrom DBI dbGetQuery
#' @importFrom stats model.matrix
#' @export
t1d_mortality_function <- function(country_wb_name, years, log_adjust = NULL) {

  if(FALSE)
  {
    # old approach , smr based on hba1c -----------------------------------------------------------------
    loc_id <- get_loc_id(country_wb_name)
    con <- get_database_connection()
    sql <- "SELECT year, smr FROM hba1c WHERE loc_id = :loc_id ORDER BY year"
    smrs <- dbGetQuery(con, sql, list(loc_id=loc_id))
    dbDisconnect(con)


    draw <- 1

    stopifnot(is.numeric(years), all(years > 1800))
    clamp_yr_covar <- pmin(pmax(years, min(smrs$year)), max(smrs$year))

    smr_vec <- with(smrs, approx(x=year, y=smr, xout=years))$y



    smr_matrix <- matrix(rep(smr_vec, each=MAX_AGE), ncol=MAX_AGE, byrow=TRUE,
                         dimnames=list(years, AGES))


  }

  if(FALSE)
  {
    # old approach , smr based on lnear regression -----------------------------------------------------------------

    loc_id <- get_loc_id(country_wb_name)

    data_long <- run_query_df(paste0("SELECT year,age, value as smr FROM classical_model_smr WHERE loc_id = ",loc_id," ORDER BY year") )

    data_long <- data_long[ data_long$year <= max(years) & data_long$year >= min(years),]
    data_wide <- spread(data_long, age, smr)
    rownames(data_wide) <- data_wide[,1]
    data_wide <- data_wide[,-1]
    smr_matrix  <- data_wide %>% as.matrix
    smr_matrix

  }
  if(TRUE)
  {
    # machine learning, with standard of care  -----------------------------------------------------------------
    loc_id <- get_loc_id(country_wb_name)

    data_long <- run_query_df(paste0("SELECT *  FROM machine_learning_model_smr WHERE loc_id = ",loc_id," ORDER BY year") )
    data_long <- data_long[ data_long$year <= max(years) & data_long$year >= min(years),]

    data_wide <- spread(select(data_long,year,age,smr= value_smr_non_minimal_care), age, smr)
    rownames(data_wide) <- data_wide[,1]
    data_wide <- data_wide[,-1]
    matrix_smr_non_minimal_care  <- data_wide %>% as.matrix


    data_wide <- spread(select(data_long,year,age,smr= value_smr_minimal_care), age, smr)
    rownames(data_wide) <- data_wide[,1]
    data_wide <- data_wide[,-1]
    matrix_smr_minimal_care  <- data_wide %>% as.matrix


    data_wide <- spread(select(data_long,year,age,smr= value_percent_non_minimal_care), age, smr)
    rownames(data_wide) <- data_wide[,1]
    data_wide <- data_wide[,-1]
    matrix_smr_percent_non_minimal_care  <- data_wide %>% as.matrix

  return(list(matrix_smr_non_minimal_care=matrix_smr_non_minimal_care
              , matrix_smr_minimal_care=matrix_smr_minimal_care
              , matrix_smr_percent_non_minimal_care=matrix_smr_percent_non_minimal_care ))

  }
}


#' Interpolate 5x5 year life table qx to single year
#'
#' Recall that in 5x5 year data, qx is given for 0 year olds,
#' 1-4 year olds, then 5-9, 10-14, ..., 95-99. We trim out 100+
#' at the data loading stage.
#'
#' @param weight1 weight of pop5_1 time five
#' @param qx5_1 vector of length 21, 5-year-binned qx with weight weight1/5
#' @param qx5_2 vector of length 21, 5-year-binned qx with weight (5-weight1)/5. If weight1 is 5 this parameter ignored.
#' @return 100-vector of single-year qx
interp_qx_5x5 <- function(weight1, qx5_1, qx5_2) {
  stopifnot(weight1 <= 5 & weight1 > 0 & length(qx5_1) == 21)
  # 1 year, 4 year, then 5 year bucketed qxs, but interpolated to match target year
  if (weight1 == 5) {
    qx5 <- qx5_1
  } else {
    stopifnot(length(qx5_2) == 21)
    qx5 <- weight1/5*qx5_1 + (5-weight1)/5*qx5_2
  }
  # linearly interpolated to single-year ages
  qx <- qx5[seq(2,21.8,0.2)]
  qx[1] <- qx5[1]

  # convert 4- and 5-year bucket qx probabilities to 1-year probabilities
  qx[2:5] <- (1+qx[2:5])^(1/4)-1
  qx[-(1:5)] <- (1+qx[-(1:5)])^(1/5)-1

  set_names(qx, AGES)
}


#' Create T1D mortality rate log adjustment
#'
#' @details Creates a function which gives the log adjustment for T1D mortality
#' rate based on hba1c changes. The formula is a simplification of
#'
#' \deqn{\log{\frac{SMR_2(y)}{SMR_1(y)}} = 0.3545 (hba1c_2(yr) - hba1c_1(yr)),}
#'
#' where each SMR is calculated using the formula \eqn{\exp\left\{-1.5274 + 0.3545 hba1c(yr)\right\}}.
#'
#' @param hba1c1 (function) hba1c function with original hba1c levels.
#' @param hba1c2 (function) hba1c function with adjusted hba1c levels.
#'
#' @return (function) function that returns t1d mortality probability log
#'   adjustment for each year
#' @export
t1d_log_adjust <- function(hba1c1, hba1c2) {
  function(yr) {
    0.3545 * (hba1c2(yr) - hba1c1(yr))
  }
}
