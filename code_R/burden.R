# Disease burden calculations
#
# Confidential
# Copyright (c) JDRF and affiliates 2019-2020. All rights reserved.
#
# Disease burden is a measure of disability, calculated as a discount to
# life-years lived. We arrive at the overall discount (called
# disability-adjusted life years, or DALYs) by multiplicatively combining
# individual disease weights (DWs) according to complications given by the
# model. We follow the method used in Gabriel and Graham's paper. The
# description provided by Emma is as follows.
#
# These were obtained from the WHO Global Burden of Disease Estimates (GBD)
# (2000-2015 and 2000-2011) to calculate DALYs for the complication listed.
#
# Gregory, G. A., Guo, J., Klatman, E. L., Ahmadov, G. A., Besançon, S.,
#   Gomez, E. D., ... & Ogle, G. D. (2020). Costs and outcomes of “intermediate”
#   vs “minimal” care for youth‐onset type 1 diabetes in six countries.
#   Pediatric Diabetes, 21(4), 628-636. https://doi.org/10.1111/pedi.12988
#


#' Disease burden as disability-adjusted-life years (DALYs)
#'
#' `calculate_dalys` is a generic function that calculates disease burden. The function invokes methods which depend on the class of `prev`.
#'
#' Note that disease burden compounds, so that we calculate DALYs as
#' the product of (1 - discount) for each disability.
#'
#' Disabilities are identified by abbreviations, which are as follows:
#' \itemize{
#'   \item DSP - distal symmetric polyneuropathy
#'   \item PR - proliferative retinopathy
#'   \item BL - blindness
#'   \item ON - overt nephropathy
#'   \item RF_transplant - renal failure (on transplant)
#'   \item RF_dialysis - renal failure (on dialysis)
#'   \item UoA - foot ulcer or amputation
#'   \item HoM - Hypertension or Microalbuminuria
#'   \item nfMI - non fatal MI
#'   \item nfCBVD - non fatal stroke
#' }
#'
#' Disease weights are stored in the global datafile, `R/sysdata.rda`.
#' Source is `data-raw/disease_weights.csv`.
#'
#' Returned list summarizes DALYs in different dimensions, all scaled to
#' population so that units are life-years:
#' \itemize{
#'   \item `year` - total DALYs by year
#'   \item `year_age` - matrix of year by age
#'   \item `year_cohort` - matrix of year and cohort of T1D onset
#'   \item `year_age_cohort` - 3D array with margins year, age, cohort of T1D onset
#' }
#'
#' @param prev prevalence, output of \link{prevalence_and_ghost_pop}()
#' @param comp complications, output of \link{complication_prevalence}()
#'
#' @return Either a list of arrays of DALYs class or, if draws are present, a list of DALYS objects corresponding to each parameter draw.
#' @export
calculate_dalys <- function(P_cohorts_level, comp) {

  # Complication probabilities by time, age, and cohort
  props <- comp$year_comp_age_cohort_prev_risk
  # Calculate disability weights to apply to prevalence. Since 100% of prevalent
  # cases have T1D, we apply the T1D disability weight directly. For all other
  # complications, we scale by the modeled complication probability.
  # Terminology note: disease weight == disability weight
  disease_weights <- get_disease_weights()
  t1d <- disease_weights$T1D
  dsp <- props["DSP",,,] * disease_weights$DSP
  pr <- props["PR",,,] * disease_weights$PR + props["BL",,,] * disease_weights$BL
  on <- props["ON",,,] * disease_weights$ON
  #rf <- (props["RF_transplant",,,] * disease_weights$Transplant
  #  + props["RF_dialysis",,,] * disease_weights$Dialysis)
  rf <- props["RF",,,] * disease_weights$Dialysis
  uoa <- props["UoA",,,] * disease_weights$UoA
  hom <- props["HoM",,,] * disease_weights$HoM
  nfMI <- props["nfMI",,,] * disease_weights$nfMI
  nfCBVD <- props["nfCBVD",,,] * disease_weights$nfCBVD
  # disability weights are combined multiplicatively and applied to prevalence level
  disability_wts <- (
    1 - (1 - t1d) * (1 - dsp) * (1 - pr) * (1 - on - rf) * (1 - uoa) *
      (1 - hom) * (1 - nfMI) * (1 - nfCBVD))
  disability <- disability_wts * P_cohorts_level
  structure(
    list(
      year=apply(disability, 1, sum),
      year_age=apply(disability, c(1, 2), sum),
      year_cohort=apply(disability, c(1, 3), sum),
      year_age_cohort = disability
    ),
    class='DALYs')
}

#' Calculate complications and DALYS
#'
#' Helper function to calculate complications and burden at the same time. Useful when running for large number of draws as it avoids memory issues caused by the large year_comp_age_cohort and year_comp_age_cohort_prev_risk elements in the complications objects. Complications are discarded after being used to calculate burden for each draw. Only burden results are returned.
#'
#' @param prev (prevalence_draws) Prevalence draws object.
#' @param hazard_model hazard model to use, either 'weibull', 'piecewise', or
#'   'truncated'
#' @param msm_adjustment if TRUE, apply 'journeys' adjustment to ON and PR for
#'   blindness and RF
#' @param dalys_names Vector of names of dalys elements to return.
#'
#' @return  list of DALYS objects corresponding to each parameter draw.
#' @export
calculate_comp_and_dalys <- function(prev,
                                     hazard_model='compromise',
                                     msm_adjustment=FALSE,
                                     dalys_names="year_age") {
  dalys_list <- list()
  for (i_draw in 1:prev$draws) {
    comp <- complication_prevalence(prev$sims[[i_draw]])
    dalys <- calculate_dalys(prev$sims[[i_draw]], comp)
    dalys_list[[i_draw]] <- dalys[dalys_names]
  }
  structure(
    list(
      sims = dalys_list,
      draws = prev$draws
    ),
    class = "DALYs_draws"
  )
}

#' Pivot DALYs into tidy tibble
#'
#' We only convert the year-age DALYs. If you want cohorts or other summaries
#' then calculate them yourself!
#'
#' @param prev prevalence - output of \link{prevalence_and_ghost_pop}()
#' @param dalys dalys - output of \link{calculate_dalys}()
#' @return tibble in tidy format
#' @export
#' @importFrom dplyr bind_cols inner_join select arrange mutate
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
tidy_dalys <- function(prev, dalys) {
  if (class(prev) != 'prevalence') {
    stop('Prevalence object is of incorrect type')
  }
  if (class(dalys) != 'DALYs') {
    stop('DALYs object is of incorrect type')
  }
  tidy_up <- function(data, colname) {
    bind_cols(age=AGES, as_tibble(t(data))) %>%
      pivot_longer(-.data$age, names_to='year', values_to=colname) %>%
      mutate(year = as.integer(.data$year))
  }
  dis <- tidy_up(dalys$year_age, 'disability')
  ddx <- tidy_up(prev$ghost_ddx_level, 'ddx')
  dexcess <- tidy_up(prev$ghost_hba1c_level, 'excess')
  dis %>%
    inner_join(ddx, by=c('year','age')) %>%
    inner_join(dexcess, by=c('year', 'age')) %>%
    select(.data$year, .data$age, .data$disability, .data$ddx, .data$excess) %>%
    arrange(.data$year, .data$age)
}

#' Summarise DALYs (and death flows) for writing to a workbook
#'
#' @param prev prevalence - output of \link{prevalence_and_ghost_pop}()
#' @param dalys dalys - output of \link{calculate_dalys}()
#' @return list of data frames
#' @export
dalys_summary <- function(prev, dalys) {
  if (class(prev) != 'prevalence') {
    stop('Prevalence object is of incorrect type')
  }
  if (class(dalys) != 'DALYs') {
    stop('DALYs object is of incorrect type')
  }
  summary <- data.frame(
    year=prev$years,
    DALYs=apply(dalys$year_age, 1, sum),
    DDx=apply(prev$DDx_flow, 1, sum),
    DExcess=apply(prev$DT1D_flow, 1, sum)
  )
  df <- function(mtx) {
    bind_cols(
      data.frame(year=rownames(mtx)),
      as.data.frame(mtx))
  }
  list(
    summary=summary,
    DALYs_age=df(dalys$year_age),
    DALYs_cohort=df(dalys$year_cohort),
    DDx_age=df(prev$DDx_flow),
    DExcess_age=df(prev$DT1D_flow)
  )
}
