# Calculate prevalence of T1D-related complications
#
# Confidential
# Copyright (c) JDRF 2020, All rights reserved
#
# Complications have either constant risk (ie they can re-occur) or
# modeled as accelerated failure time time-to-event models, with a Weibull
# hazard function and HbA1c as a linear predictor. See below for the
# Weibull parameterization.
#
# There are several different risk models given below; we are still figuring
# out which one we will use. The data are only fitted up to 30 years with the
# disease. Projecting the Weibull curve out-of-sample is inappropriate because
# the nonlinearity of the Weibull curve usually appears beyond the 30 year data
# limit. Therefore, beyond 30 years we want to be conservative, and the
# various hazard models below provide options for doing that.
#
# In each case, we discretize hazard/risk annually.

#' Weibull survival function with variable HbA1c
#'
#' Calculates survival proportions for time 0, 1, ..., T
#' given HbA1c values for T=0, T=1, etc. The Weibull hazard
#' function is given by
#' \deqn{\lambda(t; b; m; s) = \frac{t^{\left(\frac{1-s}{s}\right)}}{s}\exp\left\{-\frac{b+s h}{s}\right\}.}
#'
#' We use numerical integration to compute
#' \deqn{S(t)=\exp\left\{-\Lambda(t)\right\}=\exp\left\{-\int_0^t\lambda(\tau) d\tau\right\}}.
#'
#' @param hba1c Vector of HbA1c values as a percentage for time periods starting with 0
#' @param intercept Weibull intercept parameter b
#' @param slope Weibull slope parameter (regression coefficient for HbA1c) m
#' @param scale Weibull scale parameter s
#' @export
weib_survival <- function(hba1c, intercept, slope, scale) {

  T <- length(hba1c)-1 # total time length (years) incl. T=0
  hs_odd <- 0.5 * (hba1c[-1] + hba1c[-(T+1)]) # interpolate at half-year points
  ts <- pmin(0:T, 30)  # constant hazard after 30 years

  # Simpson's rule
  l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c)/scale) / scale
  l_odd <- (ts[-1] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale
  Lambda <- 1/6 * (l_even + cumsum(2 * c(0, 0, l_even[-c(1, T+1)]) + 4 * c(0, l_odd[-(T+1)])))

  exp(-Lambda) # survival function

if(FALSE)
{ # testing out vector operation
  system.time({

    # hba1c <- hba1c_full
    # hba1c <- hba1c_full[1,,drop=FALSE]
    # hba1c <- hba1c_full[1,]


    hba1c <- rbind(hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full
                   ,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full
                   ,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full
                   ,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full
                   ,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full,hba1c_full)
    # T <- length(hba1c)-1 # total time length (years) incl. T=0
    T <- ncol(hba1c)-1 # total time length (years) incl. T=0
    hs_odd <- 0.5 * (hba1c[,-1] + hba1c[,-(T+1)]) # interpolate at half-year points
    ts <- pmin(0:T, 30)  # constant hazard after 30 years
    ts <- matrix(ts, nrow = nrow(hba1c), ncol = ncol(hba1c), byrow = TRUE)
    # Simpson's rule
    # l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c)/scale) / scale
    # l_odd <- (ts[-1] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale
    l_even <-  exp(-(intercept + slope * hba1c)/scale) * (ts^((1 - scale)/scale))  / scale
    l_odd  <-  exp(-(intercept + slope * hs_odd)/scale) * ((ts[,-1] - 0.5)^((1 - scale)/scale))  / scale
    # Lambda <- 1/6 * (l_even + t(apply(2 * cbind(0, 0, l_even[,-c(1, T+1)]) + 4 * cbind(0, l_odd[,-(T+1)]), 1, cumsum)))
    Lambda <- 1/6 * (l_even + rowCumsums(2 * cbind(0, 0, l_even[,-c(1, T+1)]) + 4 * cbind(0, l_odd[,-(T+1)])))
    # rowCumsums(my_matrix)
    exp(-Lambda) # survival function
  })

}


}

#' Prevalence of T1D complications
#'
#' This function uses calculated T1D prevalence to compute prevalence of T1D complications.
#' Complications are modeled as a function of T1D prevalence and average ha1bc.
#'
#' The model of complications varies by age, time since T1D diagnosis (in whole years),
#' and average hba1c. Here, we work with multidimensional arrays with dimensions
#' time x complication x age x cohort, where cohorts are indexed by age at diagnosis.
#'
#' There are three complication models:
#' \itemize{
#'   \item Time-to-event models, structured as an accelerated failure time model
#'   in hba1c. These complications are presumed to be one-way transitions: once
#'   they occur they are permanent.
#'   \item Constant risk complications have a constant risk of occurrence
#'   (conditional on hba1c level), for all individuals with the disease.
#'   \item End-stage kidney disease (aka renal failure aka RF) is modelled
#'   separately. Incidence follows a time-to-event model, as above, but
#'   conditional on having RF we model kidney transplant and subsequent
#'   transplant failure. For more details of this sub-model, see documentation
#'   for \link{kidney_disease}().
#' }
#'
#' Two types of results are returned:
#' \itemize{
#' \item **levels** are expressed as numbers of individuals in the reference
#' population, and
#' \item **risks** (`year_comp_age_cohort_prev_risk`) are expressed as the share
#' of individuals with T1DM, ie a probability of having that complication,
#' conditional on having T1DM.
#' }
#'
#' Returns a list containing `year_comp_age_cohort` (4D array) complication
#' prevalence with dimensions  year x complication x age x cohort
#'
#' @param prev prevalence
#' @param hba1c_f hba1c function of year, which yields 100-vectors or matrices
#' @return a list with elements as described above
#' @export
complication_prevalence <- function(years,P_cohorts_level,smr_matrix_n, smr_matrix_m,qT1D_percent_n, hba1c_f=NULL) {

 #  years=seq(start_year, end_year); P_cohorts_level=P_cohorts_level;smr_matrix_n=  matrices_list$smr_matrix_n; smr_matrix_m= matrices_list$smr_matrix_m; qT1D_percent_n =matrices_list$qT1D_percent_n;hba1c_f=NULL

  start_year      <- min(years)
  end_year        <- max(years)

  # nyears          <- length(prev$years)
  data_start_year <- start_year-MAX_AGE+1
  all_years       <- seq(data_start_year, end_year)

  # if (is.null(hba1c_f)) {
  #   hba1c_f <- hba1c_function(prev$country)  # dimensions: year x age
  # }
  # hba1c_matrix <- hba1c_f(all_years)

  smr <- smr_matrix_n * qT1D_percent_n + smr_matrix_m * (1- qT1D_percent_n)
  hba1c_matrix <- (log(smr) +  1.5274 )/ 0.3545

  hba1c_matrix <- hba1c_matrix[-1,]

  complication_parameters <- get_complication_parameters()
  comp_names              <- with(complication_parameters, c(weibull$abbrev, constant$complication))
  nweib                   <- nrow(complication_parameters$weibull)
  nconst                  <- nrow(complication_parameters$constant)

  # 4D comp_prev array tabulates complication prevalence by time, complication,
  # age, and age at diagnosis. Dimensions are `year` x `complication` x `age` x
  # `age at diagnosis`. Values are numbers of individuals in the reference
  # population. That is, a value of 1,000 would indicate 1,000 individuals in
  # the reference country, at that age and cohort.
  comp_prev <- array(NA,
                     dim=list(length(comp_names), length(years), MAX_AGE, MAX_AGE),
                     dimnames=list(comp_names, years, AGES, AGES))

  # Weibull time-to-event complications
  risk <- array(0,
               dim=c(length(comp_names), length(years), MAX_AGE, MAX_AGE),
               dimnames = list(comp=comp_names, year=years, age=AGES, cohort=AGES))
  can_reuse_S <- any(('uniform_ages' %in% attributes(hba1c_f)) & c(attr(hba1c_f, 'uniform_ages'), FALSE))

  if (can_reuse_S) {
    # the case where we can assume hba1cs are identical for all ages and
    # therefore only integrate once per incidence year
    for (comp_i in seq_len(nweib)) {
      intercept <- complication_parameters$weibull$intercept[[comp_i]]
      slope <- complication_parameters$weibull$slope[[comp_i]]
      scale <- complication_parameters$weibull$scale[[comp_i]]
      for (incidence_year in seq(end_year, data_start_year)) {
        # maximum years an individual born (age 0) in incidence_year can be in the data
        # hba1c array index usually starts before start_year because integration starts
        # from age of incidence
        h_length <- min(1 + end_year - incidence_year, 100)
        h_age_idx <- seq(1, h_length)
        h_year_idx <- seq(incidence_year - data_start_year + 1, incidence_year - data_start_year + h_length)
        # survivor curve starting in incidence_year at age 0
        hba1cs <- hba1c_matrix[cbind(h_year_idx, h_age_idx)]
        S <- weib_survival(hba1cs, intercept, slope, scale)[1:h_length]
        for (incidence_age in seq(0, 99 - max(0, start_year - incidence_year))) {
          # risk indexes are potentially shorter, only starting at start_year
          first_risk_age <- incidence_age + max(start_year - incidence_year, 0)
          first_risk_year <- max(start_year, incidence_year)
          risk_length <- min(h_length - (first_risk_age - incidence_age), 100 - first_risk_age)
          keep_risk_idx <- seq(1 + first_risk_year - incidence_year, first_risk_year - incidence_year + risk_length)
          cohort_idx <- incidence_age + 1
          risk_year_idx <- seq(max(0, incidence_year - start_year) + 1, max(0, incidence_year - start_year) + risk_length)
          risk_age_idx <- seq(first_risk_age + 1, first_risk_age + risk_length)
          # risk is the complement of survival
          risk[cbind(comp=comp_i, year=risk_year_idx, age=risk_age_idx, cohort=cohort_idx)] <- (1 - S)[keep_risk_idx]
        }
      }
    }
  } else {
    # the case where we can't assume hba1cs are uniform across ages, e.g. under children's programs
    # of course this case is up to 100x slower, which justifies keeping code for both this and simple case
      for (incidence_year in seq(end_year, data_start_year)) {
        # incidence_year <- 1980
        # print(incidence_year)
        for (incidence_age in seq(0, 99 - max(0, start_year - incidence_year))) {
          # incidence_age <- 20
          max_age <- min(incidence_age + end_year - incidence_year, 99)
          # hba1c array indexes - potentially starts before start_year because
          # integration starts from age of incidence
          h_length <- max_age - incidence_age + 1
          h_age_idx <- seq(incidence_age + 1, max_age + 1)
          h_year_idx <- seq(incidence_year - data_start_year + 1, incidence_year - data_start_year + h_length)
          # survivor curve starting in incidence_year at age incidence_age
          hba1cs <- hba1c_matrix[cbind(h_year_idx, h_age_idx)]
          # risk indexes are potentially shorter, only starting at start_year
          first_risk_age <- incidence_age + max(start_year - incidence_year, 0)
          first_risk_year <- max(start_year, incidence_year)
          risk_length <- h_length - (first_risk_age - incidence_age)
          keep_risk_idx <- seq(1 + first_risk_year - incidence_year, first_risk_year - incidence_year + risk_length)
          cohort_idx <- incidence_age + 1
          risk_year_idx <- seq(max(0, incidence_year - start_year) + 1, max(0, incidence_year - start_year) + risk_length)
          risk_age_idx <- seq(first_risk_age + 1, first_risk_age + risk_length)
          # risk is the complement of survival
          inputs <- list()
          intercept <- complication_parameters$weibull$intercept
          slope <- complication_parameters$weibull$slope
          scale <- complication_parameters$weibull$scale

          for (comp_i in seq_len(nweib)) {
            # comp_i <- 1; print(comp_i)
            T <- length(hba1cs)-1 # total time length (years) incl. T=0
            hs_odd <- 0.5 * (hba1cs[-1] + hba1cs[-(T+1)]) # interpolate at half-year points
            ts <- pmin(0:T, 30)  # constant hazard after 30 years
            inputs[[comp_i]] <- list(hba1c=hba1cs,intercept=intercept[comp_i],slope=slope[comp_i],scale=scale[comp_i],T=T,hs_odd=hs_odd,ts=ts)
          }

          weib_survival <- function(inputs) {
            # # inputs <- inputs[[2]]
              hba1c <- inputs$hba1c ; intercept<- inputs$intercept; slope<- inputs$slope; scale<- inputs$scale
              # hba1c <- hba1cs
              T      <- inputs$T # total time length (years) incl. T=0
              hs_odd <- inputs$hs_odd # interpolate at half-year points
              ts     <- inputs$ts  # constant hazard after 30 years
              # T <- length(hba1c)-1 # total time length (years) incl. T=0
              # hs_odd <- 0.5 * (hba1c[-1] + hba1c[-(T+1)]) # interpolate at half-year points
              # ts <- pmin(0:T, 30)  # constant hazard after 30 years
              l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c)/scale) / scale
              l_odd <- (ts[-1] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale
              Lambda <- 1/6 * (l_even + cumsum(2 * c(0, 0, l_even[-c(1, T+1)]) + 4 * c(0, l_odd[-(T+1)])))
              exp(-Lambda) # survival function
            # exp(-inputs$hba1c )
          }
           S <-  lapply( inputs, weib_survival )
           for (comp_i in seq_len(nweib)) {
             risk[cbind(comp=comp_i, year=risk_year_idx, age=risk_age_idx, cohort=cohort_idx)] <- (1 - S[[comp_i]])[keep_risk_idx]
           }

        }
      }

  }
  # This is Graham's and Gabriel's 'journeys' adjustments.
  # Reduces ON and PR by survivor shares of RF and BL, respectively,
  # effectively making this a multistate model (rather than a pure
  # time-to-event model)
  risk['ON',,,] <- risk['ON',,,] * (1 - risk['RF',,,])
  risk['PR',,,] <- risk['PR',,,] * (1 - risk['BL',,,])
  for (comp_i in seq_len(nweib)) {
    comp_prev[comp_i,,,] <- P_cohorts_level * risk[comp_i,,,]
  }

  # constant risk complications
  for (comp_i in seq_len(nconst)) { # comp_i <- 1
    const_risk <- complication_parameters$constant[comp_i,]
    lt8 <- 1e-2 * const_risk$hba1c_lt_8
    gt8lt9 <- 1e-2 * const_risk$hba1c_8_to_9
    gt9 <- 1e-2 * const_risk$hba1c_gteq_9
    cr <- function(h) ifelse(h < lt8, lt8, ifelse(h < 9, gt8lt9, gt9))
    rsk <- apply(hba1c_matrix[-(1:99),], MARGIN=1:2, FUN=cr)
    for (cohort_i in seq_len(MAX_AGE)) {
      risk[nweib + comp_i,,,cohort_i] <- rsk
    }
    comp_prev[nweib + comp_i,,,] <- risk[nweib + comp_i,,,] * P_cohorts_level
  }
  structure(list(
    year_comp_age_cohort=comp_prev,
    year_comp_age_cohort_prev_risk=risk
  ), class='complications')
}

#' Convert complications into tidy format
#'
#' @param comp complication prevalence, output of
#'        \link{complication_prevalence}().
#' @return tibble with columns year, age, and a column per complication
#' @importFrom purrr set_names map_dfc as_vector
#' @importFrom tidyr crossing
#' @importFrom dplyr bind_cols mutate
tidy_complications <- function(comp) {
  if (class(comp) != 'complications') {
    stop('Complications object has incorrect type')
  }
  # comp$year_comp_age has dimensions `year` x `complication` x `age`
  year <- dimnames(comp$year_comp_age_cohort)[[2]]
  comps <- dimnames(comp$year_comp_age_cohort)[[1]]
  stopifnot(all(MAX_AGE == dim(comp$year_comp_age_cohort)[3]))
  stack_cols <- function(mtrx) as.numeric(t(mtrx))
  year_comp_age <- apply(comp$year_comp_age_cohort, c(1,2,3), sum)
  bind_cols(
    crossing(year=as.numeric(year), age=AGES),
    map_dfc(set_names(comps, comps),
            function(cmp) stack_cols(year_comp_age[cmp,,]))
  )
}
