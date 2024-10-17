MAX_AGE <- 100

#' Vector of ages to use in the model
AGES <- seq(0, MAX_AGE-1)

# Prevalence calculation functions
#
# Confidential
# Copyright (c) JDRF 2020, All rights reserved
#
# Prevalence is calculated using an illness-death Markov model as described in
# the paper.
#
# Most quantities in this file (and indeed elsewhere in the model) are expressed
# as arrays with dimension MAX_AGE (100), corresponding to ages 0,1,...,99, ie
# up until an individual's 100th birthday. This approach allows us to benefit
# from vectorization and is what lets us calculate model results interactively.
#
# Most matrices in the file are indexed by {year}x{age}, ie if the model runs
# from 2005-2009 then each matrix will have dimension 5x100, and the value for
# age 0 in 2006 will be entry (1, 2). Pythonists, remember that R is barbaric
# and indexes from 1.


#' Calculate T1D prevalence
#'
#' @description Inputs `i`, `qB`, `qT1D`, and `dDx` are data matrixes of shape
#'   `nyears` x `MAX_AGE`, where each value relates to a single year and age.
#'
#' Outputs a list whose values are compartment values (`S`, `P`, `Pcohorts`,
#' `D`) expressed as a proportion of 1, as well as flows (`I`, `DDx`, `DT1D`) in
#' the same units.
#'
#' Equation numbers in code comments refer to the working paper.
#'
#' @param i Observed incidence (excludes deaths on diagnosis)
#' @param qB Background mortality from life tables
#' @param qT1D Type-1 diabetes mortality function
#' @param dDx Death on diagnosis rate
#' @param years Vector of consecutive years to model, including 100 years of warmup
#' @return structure as described above
#' @keywords internal
calculate_prevalence <- function(i, qB,  qT1D_n, qT1D_m, qT1D_percent_n, dDx, years) {

  # Data matrixes: years on rows, ages on columns. Memory footprint of storing
  # all this information really is tiny, a few MB. Similarly, a 3D array: {year}x{age}x{onset age}
  # Model compartments, all annual cohorts as a proportion of 1.
  S <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # S - susceptible (non-T1D)
  P <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # P - prevalence (T1D)
  D <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # D - deaths (absorbing state)

  # incidence flows are reflected the following year, for individuals one year older
  i_shift <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))
  i_shift[, 1] <- 0
  i_shift[,-1] <- i[,-MAX_AGE]
  i_all <- i / (1-dDx)              # includes deaths at onset
  i_all_shift <- i_shift / (1-dDx)  # shifted & includes deaths at onset

  # Track P cohorts and in/outflows in 3D array: {year}x{age}x{onset age}
  Pcohorts  <- array(NA, dim=list(length(years), MAX_AGE, MAX_AGE), dimnames=list(years, AGES, AGES))
  Icohorts  <- array(0, dim=list(length(years), MAX_AGE, MAX_AGE), dimnames=list(years, AGES, AGES))
  PDcohorts <- array(0, dim=list(length(years), MAX_AGE, MAX_AGE), dimnames=list(years, AGES, AGES))

  # The model proceeds yearwise, populating successive matrix rows.
  # Equation references below refer to the model summary.
  S[1,1] <- 1
  P[1,1] <- Pcohorts[1,1,1] <- D[1,1] <- 0
  for (t in seq_along(years[-1])) {
    # for (t in 1:90) {
    # bump previous period's data along one cell (ie one year older)
    Sshift <- c(1, S[t, -MAX_AGE])
    Pshift <- c(0, P[t, -MAX_AGE])
    Dshift <- c(0, D[t, -MAX_AGE])

    # equation (3) - susceptible compartment, S
    S[t+1,] <- Sshift * (1 - i_all_shift[t,]) * (1 - qB[t,])

    # equation (4) - prevalence compartment P
    # P[t+1,] <- Pshift * (1 - qT1D[t,]) + i_shift[t,] * Sshift

    P[t+1,] <-    qT1D_percent_n[t,]  *  Pshift  * (1 - qT1D_n[t,]) + i_shift[t,] * Sshift +
      (1-qT1D_percent_n[t,]) *  Pshift  * (1 - qT1D_m[t,])


    if(config$run_days_lost)  # do not run levers for paper stats
    {
      PCshift <- array(0, dim=c(MAX_AGE, MAX_AGE))
      PCshift[1,] <- 0
      PCshift[-1,] <- Pcohorts[t,-MAX_AGE,]

      Ishift <- matrix(NA, nrow = MAX_AGE, ncol = MAX_AGE)
      Ishift[,-MAX_AGE] <- diag(i_shift[t,] * Sshift)[,-1] # Only want to shift age, not onset age
      Ishift[,MAX_AGE] <- 0  # NB: half-cycle adjustment populates any incidence for MAX_AGE

      # equation (4') - prevalence compartment P by cohort
      Icohorts[t,,] <- diag(i[t,] * S[t,])
      # PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D[t,], c(MAX_AGE, MAX_AGE))
      PDcohorts[t,,] <- qT1D_percent_n[t,] * Pcohorts[t,,] * array(qT1D_n[t,], c(MAX_AGE, MAX_AGE)) +
        (1 - qT1D_percent_n[t,]) * Pcohorts[t,,] * array(qT1D_m[t,], c(MAX_AGE, MAX_AGE))

      # Pcohorts[t+1,,] <- PCshift * array(1 - qT1D[t,], c(MAX_AGE, MAX_AGE)) + Ishift
      Pcohorts[t+1,,] <- qT1D_percent_n[t,] * PCshift * array(1 - qT1D_n[t,], c(MAX_AGE, MAX_AGE)) + Ishift +
        (1 - qT1D_percent_n[t,]) * PCshift * array(1 - qT1D_m[t,], c(MAX_AGE, MAX_AGE))
    }
    # equation (5) - death compartment D
    D[t+1,] <- (Dshift
                + i_all_shift[t,] * dDx[t,] * Sshift
                # + qT1D[t,] * Pshift
                + qT1D_n[t,] * Pshift *       qT1D_percent_n[t,]
                + qT1D_m[t,] * Pshift *  (1 - qT1D_percent_n[t,])
                + Sshift * qB[t,] * (1 - i_all_shift[t,]))

  }
  # final year cohort flows - same as in the loop but for final period
  t <- length(years)
  Icohorts[t,,] <- diag(i[t,] * S[t,])

  # PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D[t,], c(MAX_AGE, MAX_AGE))
  PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D_n[t,], c(MAX_AGE, MAX_AGE)) *      qT1D_percent_n[t,]  +
    Pcohorts[t,,] * array(qT1D_m[t,], c(MAX_AGE, MAX_AGE)) * (1 - qT1D_percent_n[t,])

  # flows: based on unshifted versions of incidence
  I <- i * S                    # T1D incidence
  DDx <- i_all * dDx * S        # deaths at T1D onset
  # DDx <- i_all * S - I        # deaths at T1D onset
  # DT1D <- (qT1D - qB) * P       # T1D-cause mortality
  if(config$lever_change_start_at==1860)
  {
    DT1D <- (qT1D_n - qB) * P  *          qT1D_percent_n[t,]  +
      (qT1D_m - qB) * P  *     (1 - qT1D_percent_n[t,])   # T1D-cause mortality
  }else
  {
    DT1D <- (qT1D_n - qB) * P  *          qT1D_percent_n  +
      (qT1D_m - qB) * P  *     (1 - qT1D_percent_n)   # T1D-cause mortality
  }

  DBGP <- P * qB                # background mortality for people w/ T1D
  DBGS <- S * qB * (1 - i_all)  # background mortality for susceptible pop'n
  # DBGS <- S * qB  # background mortality for susceptible pop'n

  # Half-cycle adjustments - apply half of each flow in the reference
  # year for that flow. We expect that on average, 1/2 of the flow
  # occurs part-way through the year.
  S <- S - 0.5 * I - 0.5 * DDx - 0.5 * DBGS
  P <- P + 0.5 * I - 0.5 * (DT1D + DBGP)
  D <- D + 0.5 * (DDx + DBGS + DT1D + DBGP)
  Pcohorts <- Pcohorts + 0.5 * Icohorts - 0.5 * PDcohorts
  list(S=S, P=P, Pcohorts=Pcohorts, D=D, I=I, DDx=DDx, DT1D=DT1D,DBGP=DBGP, DBGS=DBGS)
}



data_long_2_matrices <- function(data_long,data_long_default_run,smr_scale_factor=1,incidence_scale_factor=1)
{

  # sensitivity analysis
  if(scenarios$pediatric_incidence_plus_10_perc) { data_long$incidence_rate[data_long$age <= 18 ] <- data_long$incidence_rate[data_long$age <= 18 ]  * 1.1  }
  if(scenarios$pediatric_incidence_minus_10_perc){ data_long$incidence_rate[data_long$age <= 18 ] <- data_long$incidence_rate[data_long$age <= 18 ]  * 0.9  }
  if(scenarios$pediatric_incidence_plus_25_perc) { data_long$incidence_rate[data_long$age <= 18 ] <- data_long$incidence_rate[data_long$age <= 18 ]  * 1.25  }
  if(scenarios$pediatric_incidence_minus_25_perc){ data_long$incidence_rate[data_long$age <= 18 ] <- data_long$incidence_rate[data_long$age <= 18 ]  * 0.75  }

  if(scenarios$smr_plus_10_perc) { data_long$value_smr_minimal_care     <- data_long$value_smr_minimal_care  * 1.1  }
  if(scenarios$smr_plus_10_perc) { data_long$value_smr_non_minimal_care <- data_long$value_smr_non_minimal_care  * 1.1  }
  if(scenarios$smr_minus_10_perc){ data_long$value_smr_minimal_care     <- data_long$value_smr_minimal_care  * 0.75  }
  if(scenarios$smr_minus_10_perc){ data_long$value_smr_non_minimal_care <- data_long$value_smr_non_minimal_care  * 0.75  }

  if(scenarios$diagnosis_rate_plus_25_pp){data_long$mortality_undiagnosed_rate  <- sapply( data_long$mortality_undiagnosed_rate  - 0.25,max,0)}
  index_diagnosis <- data_long$income_class != "HIC" & data_long$age <= 24
  if(scenarios$diagnosis_rate_minus_25_pp & sum(index_diagnosis)){  data_long$mortality_undiagnosed_rate[index_diagnosis] <- sapply( data_long$mortality_undiagnosed_rate[index_diagnosis]  + 0.25,min,0.99999999)  }

  if(scenarios$diagnosis_rate_left){ data_long$mortality_undiagnosed_rate <- data_long$mortality_undiagnosed_rate_right }
  if(scenarios$diagnosis_rate_right){ data_long$mortality_undiagnosed_rate <- data_long$mortality_undiagnosed_rate_left }



  # write.csv(data_long,"temp/input_all.csv")
  # data_long <- data_long[ data_long$year <= max(years) & data_long$year >= min(years),]

  data_wide <- spread(dplyr::select(data_long,year,age,value=background_mortality_rate ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  qB                  <- data_wide

  data_wide <- spread(dplyr::select(data_long,year,age,value=background_population ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  pop                  <- data_wide

  data_wide <- spread(dplyr::select(data_long,year,age,value=mortality_undiagnosed_rate), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  dDx                 <- data_wide
  dDx[dDx >= 1] <- 0.95  # Ddx can't be 100% or i_all is undefined

  data_wide <- spread(dplyr::select(data_long,year,age,value=incidence_rate ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  i                   <- incidence_scale_factor * data_wide / 1e5

  data_wide <- spread(dplyr::select(data_long_default_run,year,age,value=incidence_rate ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  i_default                   <- incidence_scale_factor * data_wide / 1e5


  # t1d_mortality_function  --------------------------------------------------------------------


  data_wide <- spread(dplyr::select(data_long,year,age,value=value_smr_non_minimal_care ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  smr_matrix_n                  <- smr_scale_factor * data_wide

  data_wide <- spread(dplyr::select(data_long_default_run,year,age,value=value_smr_non_minimal_care ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  smr_matrix_n_default                  <- smr_scale_factor * data_wide

  data_wide <- spread(dplyr::select(data_long,year,age,value=value_smr_minimal_care ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  smr_matrix_m                  <- smr_scale_factor * data_wide

  data_wide <- spread(dplyr::select(data_long_default_run,year,age,value=value_smr_minimal_care ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  smr_matrix_m_default                  <- smr_scale_factor * data_wide



  data_wide <- spread(dplyr::select(data_long,year,age,value=value_percent_non_minimal_care ), age, value)
  rownames(data_wide) <- data_wide$year
  data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  qT1D_percent_n                  <- data_wide

  # apply growth rate after 2020 ------------------
  apply_Growth_Rate <- function(matrix_rate,matrix_rate_default=NULL,year_start=2021, rate=NA)
  {
    # matrix_rate <- smr_matrix_m
    # data.frame(x=rownames(matrix_rate),age_1= matrix_rate[,1],age_25=matrix_rate[,25],age_50=matrix_rate[,50],age_75=matrix_rate[,75],age_99=matrix_rate[,99]) %>%
    #   e_charts(x) %>%e_line(age_1)%>%e_line(age_25)%>%e_line(age_50)%>%e_line(age_75)%>%e_line(age_99)%>%e_tooltip(trigger = "item")%>% e_x_axis(min='1980')
    if(is.null(matrix_rate_default)){ matrix_rate_default <-matrix_rate }
    rate_avg_5 <- (    matrix_rate_default[as.character(year_start),]/matrix_rate_default[as.character(year_start-1),] +
                         matrix_rate_default[as.character(year_start-1),]/matrix_rate_default[as.character(year_start-2),] +
                         matrix_rate_default[as.character(year_start-2),]/matrix_rate_default[as.character(year_start-3),] +
                         matrix_rate_default[as.character(year_start-3),]/matrix_rate_default[as.character(year_start-4),] +
                         matrix_rate_default[as.character(year_start-4),]/matrix_rate_default[as.character(year_start-5),]
                       +matrix_rate_default[as.character(year_start-5),]/matrix_rate_default[as.character(year_start-6),]
                       +matrix_rate_default[as.character(year_start-6),]/matrix_rate_default[as.character(year_start-7),]
                       +matrix_rate_default[as.character(year_start-7),]/matrix_rate_default[as.character(year_start-8),]
                       +matrix_rate_default[as.character(year_start-8),]/matrix_rate_default[as.character(year_start-9),]
                       +matrix_rate_default[as.character(year_start-9),]/matrix_rate_default[as.character(year_start-10),]
    )/10
    if(!is.na(rate))
    {rate_avg_5[!is.nan(rate_avg_5)] <- 1+rate/100 }

    rate_avg_5[is.nan(rate_avg_5)] <- 0

    for(year in (year_start+1): max(rownames(matrix_rate)))
    {#  year <- 2022
      matrix_rate[as.character(year),] <- matrix_rate[as.character(year-1),] *  rate_avg_5
    };return(matrix_rate)
  }
  # if(TRUE)
  print("  print(config$run_projection)       ------------------")
  print(config$run_projection)
  if(config$run_projection)
  {
    # i   <- apply_Growth_Rate(i,year_start=2021,rate = 1.5)
    i   <- apply_Growth_Rate(matrix_rate=i,matrix_rate_default=i_default,year_start=2021)
    dDx <- apply_Growth_Rate(matrix_rate=dDx,year_start=2021)
    smr_matrix_n <- apply_Growth_Rate(matrix_rate=smr_matrix_n,matrix_rate_default=smr_matrix_n_default,year_start=2021)
    smr_matrix_m <- apply_Growth_Rate(matrix_rate=smr_matrix_m,matrix_rate_default=smr_matrix_m_default,year_start=2021)
  }

  i               <- pre_fill_constant_1860_1900(i)  ;    # mean(i[as.character(2010:2019),19]) *100000
  qB              <- pre_fill_constant_1860_1900(qB)
  # qT1D_n          <- pre_fill_constant_1860_1900(qT1D_n, years)
  # qT1D_m          <- pre_fill_constant_1860_1900(qT1D_m, years)
  qT1D_percent_n  <- pre_fill_constant_1860_1900(qT1D_percent_n)
  smr_matrix_m    <- pre_fill_constant_1860_1900(smr_matrix_m)
  smr_matrix_n    <- pre_fill_constant_1860_1900(smr_matrix_n)
  dDx             <- pre_fill_constant_1860_1900(dDx)
  pop             <- pre_fill_constant_1860_1900(pop)
  return(list(i=i,qB=qB,pop=pop, dDx=dDx,smr_matrix_n=smr_matrix_n, smr_matrix_m=smr_matrix_m,qT1D_percent_n=qT1D_percent_n))
}

pre_fill_constant_1860_1900 <- function (rate)
{ # rate <- i ;
  years <- 1860
  rate_pre_fill <- rate[1,,drop=FALSE]
  rate_pre_fill <- rate_pre_fill[rep(1, times = 1900-years), ]
  if(nrow(rate_pre_fill))
  { rownames(rate_pre_fill) <- as.character(years:1899)
  rate <- rbind(rate_pre_fill,rate)
  }; return(rate)
}

#' Simplified T1D prevalence calculator based on incidence levels
#'
#' Inputs `I`, `qB`, `qT1D`, and `dDx` are data matrixes of shape
#'   `nyears` x `MAX_AGE`, where each value relates to a single year and age.
#'
#' Outputs a list whose values are compartment values (`P`) expressed as a number
#' of individuals, as well as flows (`I`, `DDx`, `DT1D`) in the same units.
#'
#' There is no half-cycle adjustment. Incidence is contemporaneous (2020 incidence
#' will appear as 2020 prevalence). This is for simplicity, and to allow the outputs
#' to line up with the inputs.
#'
#' @param I Incidence, expressed as a number of individuals
#' @param qB Background mortality from life tables, as an annualized probability
#' @param qT1D Type-1 diabetes mortality function, as an annualized probability
#' @param years Vector of consecutive years to model. Warmup periods are not necessary.
#' @return structure as described above
#' @keywords internal
prevalence_from_incidence_level <- function(I, qB, qT1D, years) {
  stopifnot(
    all(dim(I) == dim(qB)),
    all(dim(I) == dim(qT1D)),
    all(nrow(I) == length(years))
  )
  # T1D prevalence
  P <- matrix(0, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))
  Pcohorts <- array(0, dim=list(length(years), MAX_AGE, MAX_AGE),
                    dimnames=list(years, AGES, AGES))
  # The model proceeds year-wise, populating successive matrix rows.
  # Unlike the other model, incidence is contemporaneous, so mind the
  # different indexes
  P[1,] <- I[1,]
  Pcohorts[1,,] <- diag(I[1,])
  for (t in seq_along(years[-1])) {
    Pshift <- c(0, P[t, -MAX_AGE])
    PCshift <- array(0, dim=c(MAX_AGE, MAX_AGE))
    PCshift[1,] <- 0
    PCshift[-1,] <- Pcohorts[t,-MAX_AGE,]
    # equation (4) - prevalence compartment P
    P[t+1,] <- Pshift * (1 - qT1D[t,]) + I[t+1,]
    Pcohorts[t+1,,] <- PCshift * array(1 - qT1D[t,], c(MAX_AGE, MAX_AGE)) + diag(I[t+1,])
  }
  # T1D-cause mortality (flow)
  DT1D <- (qT1D - qB) * P
  # sanity check
  stopifnot(all(apply(Pcohorts, 1, sum) - apply(P, 1, sum) < 1e-10))

  list(P=P, Pcohorts=Pcohorts, DT1D=DT1D)
}


#' Calculate prevalence and ghost population with direct incidence
#'
#' @param country_name world bank name of country
#' @param years vector of years to model; additional history not required
#' @param incidence_fn callback function returning vector of incidence, measured in persons, of same length as years
#' @param hba1c_fn hba1c callback function
#' @param bg_mort_fn background mortality callback function
#' @param t1d_mort_fn t1d mortality function
prevalence_and_ghost_pop_from_inc_level <- function(
    country_name,
    years,
    incidence_fn,
    hba1c_fn,
    bg_mort_fn,
    t1d_mort_fn
) {
  start_time <- Sys.time()

  # Track model parameters in matrixes X[year-1, age]
  mtrx_from_function <- function(f) {
    matrix(sapply(years, f), nrow=length(years), byrow=TRUE, dimnames=list(years, AGES))
  }

  qB <- t(mapply(background_mortality_function,country_wb_name,years))
  dimnames(qB) <- list(years, AGES)

  qT1D <- mtrx_from_function(t1d_mortality_function(country_name))
  I <- mtrx_from_function(incidence_fn)

  # T1D prevalence
  prev <- prevalence_from_incidence_level(I, qB, qT1D, years)
  # counterfactual, same but using background as T1D mortality
  cfac <- prevalence_from_incidence_level(I, qB, qB, years)

  # by assumption, ddx is zero
  ghost_ddx_level <- prev$P * 0
  # reduced hba1c mort rate shows up in prevalence
  ghost_hba1c_level <- cfac$P - prev$P
  # ghost population is the sum of these two
  ghost_level <- ghost_hba1c_level + ghost_ddx_level

  # result structure mirrors the usual prevalence_and_ghost_pop() function
  time_taken <- Sys.time() - start_time
  with(prev,
       structure(
         list(
           P_level=P,
           S_level=NA,
           P_cohorts_level=Pcohorts,
           I_flow=I,
           DDx_flow=I*0,
           DT1D_flow=DT1D,
           ghost_level=ghost_level,
           ghost_ddx_level=ghost_ddx_level,
           ghost_hba1c_level=ghost_hba1c_level,
           pop=NA,
           country=country_name,
           time=time_taken,
           hba1c=hba1c_fn,
           years=years,
           pop_scale_factor=1
         ),
         class='prevalence'))
}


#' Convert several matrixes to a tibble
#'
#' The matrixes must all have the same dimension.
#'
#' @param years optional vector of years spanned by data
#' @param ages optional vector of ages spanned by data
#' @param ... matrices of the same dimension, with age on columns and year on rows
#' @return tibble
#' @importFrom tidyr crossing
#' @importFrom dplyr bind_cols
#' @importFrom purrr map_dfc
#' @export
matrixes_to_long_format <- function(years=NULL, ages=NULL, ...) {
  ms <- list(...)
  if(length(ms) == 0) {
    stop('No matrices specified.')
  }
  to_long <- function(X) as.vector(t(X))
  if (is.null(years)) {
    years <- as.numeric(rownames(ms[[1]]))
  }
  if (is.null(ages)) {
    ages <- as.numeric(colnames(ms[[1]]))
  }
  index <- crossing(year=years, age=ages)
  bind_cols(index, map_dfc(ms, to_long))
}

#' Summarize prevalence for writing to a workbook
#'
#' @param prev prevalence, output of \link{prevalence_and_ghost_pop}().
#' @return list of matrixes and data frames
#' @importFrom tibble tibble
prevalence_summary <- function(prev) {
  stopifnot('prevalence' %in% class(prev))
  ghost_summary <- with(prev,
                        data.frame(year=years,
                                   ghost=apply(ghost_level, 1, sum),
                                   ghost_ddx=apply(ghost_ddx_level, 1, sum),
                                   ghost_hba1c=apply(ghost_hba1c_level, 1, sum)))
  df <- function(mtx) bind_cols(data.frame(year=rownames(mtx)), as.data.frame(mtx))
  with(prev, list(
    prevalence_age=df(apply(P_cohorts_level, c(1, 2), sum)),
    prevalence_cohort=df(apply(P_cohorts_level, c(1, 3), sum)),
    ghost=ghost_summary,
    ghost_age=df(ghost_level),
    ghost_ddx_age=df(ghost_ddx_level),
    ghost_hba1c_age=df(ghost_hba1c_level)
  ))
}

#' Tabulate prevalence by year of onset
#'
#' We currently track prevalence cohorts by age at onset. This function tabulates
#' prevalence by *year* of onset. The total overall prevalence is unchanges; this
#' function just provides an alternative breakdown.
#'
#' @param prev prevalence, output of \link{prevalence_and_ghost_pop}().
#' @return prevalence matrix, year and year of onset
#' @export
prevalence_by_cohort <- function(prev) {
  if (class(prev) != 'prevalence') {
    stop('Prevalence object is of incorrect type')
  }
  P_cohorts_level <- prev$P_cohorts_level
  years <- prev$years
  cohort_years <- seq(min(years)-MAX_AGE, max(years))
  # prevalence by reference year, and year of onset of T1D
  year_prev <- matrix(0, nrow=length(cohort_years), ncol=length(years),
                      dimnames=list(sort(cohort_years, decreasing=TRUE), years))
  # The model currently tracks *age* of onset of T1D, not year of onset, so
  # to construct this matrix we need to sum the diagonals of each year-slice
  # of the P_cohorts_level 3D array. This is a lower-triangular matrix, where
  # the main diagonal is those with onset in the reference (ie current) year
  # because (age of onset) == (age).
  # i indexes reference year
  COHORTS <- AGES  # 'years ago' of onset, starting with 0 for current year
  nyears <- length(years)
  for (i in seq_along(years)) {
    # the year is presented in reverse order
    # j indexes cohorts, ie age of onset
    for (j in COHORTS) {
      # start with main diagonal, work our way down
      year_prev[j+nyears-i+1, i] <- ifelse(
        j == 99,
        sum(P_cohorts_level[i, (1+j):MAX_AGE, 1:(MAX_AGE-j)]),
        sum(diag(P_cohorts_level[i, (1+j):MAX_AGE, 1:(MAX_AGE-j)]))
      )
    }
  }
  year_prev
}




Apply_levers_to_input_matrices <- function(matrices_list, lever=1, lever_year_range)
{
  i                 <- matrices_list$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list$dDx
  smr_matrix_n      <- matrices_list$smr_matrix_n
  smr_matrix_m      <- matrices_list$smr_matrix_m
  qT1D_percent_n    <- matrices_list$qT1D_percent_n

  # lever=1; lever_year_range = as.character(config$lever_change_start_at:max(years))
  if(lever>=1)
  { # full diagnosis
    i_new <- i
    i_new[lever_year_range,]   <- (i_new/(1-dDx))[lever_year_range,]

    dDx_new <- dDx
    dDx_new[lever_year_range,] <- 0

    smr_matrix_n_new = smr_matrix_n
    smr_matrix_m_new = smr_matrix_m
    qT1D_percent_n_new = qT1D_percent_n

  }

  if(lever>=2)
  { # basic care, insulin, strip, education.  SMR  median( c(3.7, 4.4)) , 4.05,  across age 20-50 -----
    ratio        <- 4.05/(rowSums(smr_matrix_n[,as.character(20:50)])/length(20:50))
    ratio_matrix <- t(matrix(ratio, nrow=100, ncol=length(ratio), byrow=TRUE))
    ratio_matrix[ratio_matrix>=1] <- 1
    smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- smr_matrix_n[smr_matrix_adjusted<1]

    smr_matrix_n_new <- smr_matrix_n
    smr_matrix_n_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]
    smr_matrix_m_new <- smr_matrix_m

    qT1D_percent_n_new <- qT1D_percent_n
    qT1D_percent_n_new[as.character(lever_year_range),] <- 1

  }

  if(lever>=3)
  { # best care Pumps and FGMs/CGMs SMR median( c(2.2, 2.6)) , median 2.4--------------------------------------------------------------------
    ratio        <- 2.4/(rowSums(smr_matrix_n[,as.character(20:50)])/length(20:50))
    ratio_matrix <- t(matrix(ratio, nrow=100, ncol=length(ratio), byrow=TRUE))
    ratio_matrix[ratio_matrix>=1] <- 1
    smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- smr_matrix_n[smr_matrix_adjusted<1]

    smr_matrix_n_new <- smr_matrix_n
    smr_matrix_n_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]
    smr_matrix_m_new <- smr_matrix_m

    qT1D_percent_n_new <- qT1D_percent_n
    qT1D_percent_n_new[as.character(lever_year_range),] <- 1

  }

  if(lever>=4)
  { # cure
    smr_matrix_n_new <- smr_matrix_n
    smr_matrix_n_new[as.character(lever_year_range),] <- 1
    smr_matrix_m_new <- smr_matrix_m
    smr_matrix_m_new[as.character(lever_year_range),] <- 1

    qT1D_percent_n_new <- qT1D_percent_n
    qT1D_percent_n_new[as.character(lever_year_range),] <- 1

  }
  return(list(i=i_new,dDx=dDx_new,smr_matrix_n=smr_matrix_n_new,smr_matrix_m=smr_matrix_m_new,qT1D_percent_n=qT1D_percent_n_new,qB=matrices_list$qB ,pop=matrices_list$pop  ))
}




Apply_smr_to_input_matrices <- function(matrices_list, lever=1, smr_low,smr_high,lever_year_range)
{   #  smr_low <- 5.3 ;  smr_high <- 6.3
  i                 <- matrices_list$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list$dDx
  smr_matrix_n      <- matrices_list$smr_matrix_n
  smr_matrix_m      <- matrices_list$smr_matrix_m
  qT1D_percent_n    <- matrices_list$qT1D_percent_n

  # lever=1; lever_year_range = as.character(config$lever_change_start_at:max(years))
  if(TRUE)
  { # full diagnosis
    i_new <- i
    i_new[lever_year_range,]   <- (i_new/(1-dDx))[lever_year_range,]

    dDx_new <- dDx
    dDx_new[lever_year_range,] <- 0

    smr_matrix_n_new = smr_matrix_n
    smr_matrix_m_new = smr_matrix_m
    qT1D_percent_n_new = qT1D_percent_n

  }

  if(TRUE)
  { # basic care, insulin, strip, education.  SMR  median( c(3.7, 4.4)) , 4.05,  across age 20-50 -----
    ratio        <- mean(smr_low,smr_high)/(rowSums(smr_matrix_n[,as.character(20:50)])/length(20:50))
    ratio_matrix <- t(matrix(ratio, nrow=100, ncol=length(ratio), byrow=TRUE))
    # ratio_matrix[ratio_matrix>=1] <- 1
    smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- smr_matrix_n[smr_matrix_adjusted<1]

    smr_matrix_n_new <- smr_matrix_n
    smr_matrix_n_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]
    smr_matrix_m_new <- smr_matrix_m

    qT1D_percent_n_new <- qT1D_percent_n
    qT1D_percent_n_new[as.character(lever_year_range),] <- 1

  }


  return(list(i=i_new,dDx=dDx_new,smr_matrix_n=smr_matrix_n_new,smr_matrix_m=smr_matrix_m_new,qT1D_percent_n=qT1D_percent_n_new,qB=matrices_list$qB ,pop=matrices_list$pop  ))
}







Agg_country <- function(country_data_merge ,key_list)
{

  agg_life_expectancy_median <- function(year,frequency)
  { #
    #     years <- 1960 ; ages= 98
    #
    # # print(years)
    # year <- country_data_merge$`Life expectency (6 t1d cure)` [country_data_merge$Year==years & country_data_merge$Age==ages ]
    # frequency <- (country_data_merge$Prevalence + country_data_merge$Ghosts)[country_data_merge$Year==years& country_data_merge$Age==ages]


    # df <- data.frame(year=year,frequency= frequency)
    # df <- setDT(df)[,list(frequency=sum(frequency)),by="year"]
    # df <- df[order(df$year)]
    # df$frequency <- df$frequency + 0.01  #  for extrem cases
    # df$frequency_cumsum <- cumsum(df$frequency)
    # index_same_value_index <- df$frequency_cumsum== df$frequency_cumsum[nrow(df)]/2
    # if(sum(index_same_value_index)>0)
    # {
    #   median <- (df$year[which(index_same_value_index) ] + df$year[which(index_same_value_index)+1]) /2
    #   median
    # }else{
    #   index <- min(which(df$frequency_cumsum >  df$frequency_cumsum[nrow(df)]/2) )
    #   median <- df$year[index]
    #   median
    # }

    median <- sum(year * frequency )/ sum(frequency)
    median
  }

  country_data <- setDT(country_data_merge)[,list(
    `Ann. background mortality` = sum( `Ann. background mortality` )
    ,`Ann. background population` = sum( `Ann. background population` )
    ,`Ann. days lost (1 base)` = sum( `Ann. days lost (1 base)` )
    ,`Ann. days lost (2 diagnosis)`=sum(`Ann. days lost (2 diagnosis)`)
    ,`Ann. days lost (3 basic care)`  = sum(`Ann. days lost (3 basic care)`)
    ,`Ann. days lost (4 best care)`  = sum(`Ann. days lost (4 best care)`)
    ,`Ann. days lost (5 cure)`  = sum(`Ann. days lost (5 cure)`)
    ,`Ann. early deaths`  = sum(`Ann. early deaths`)
    ,`Ann. early deaths (2 diagnosis)`  = sum(`Ann. early deaths (2 diagnosis)`)
    ,`Ann. early deaths (3 basic care)`  = sum(`Ann. early deaths (3 basic care)`)
    ,`Ann. early deaths (4 best care)`  = sum(`Ann. early deaths (4 best care)`)
    ,`Ann. early deaths (5 cure)`  = sum(`Ann. early deaths (5 cure)`)
    ,`Ann. onset deaths`  = sum(`Ann. onset deaths`)
    ,`Ghosts`  = sum(`Ghosts`)
    ,`Ghosts (delta basic care)`  = sum(`Ghosts (delta basic care)`)
    ,`Ghosts (delta best care)`   = sum(`Ghosts (delta best care)`)
    ,`Ghosts (delta cure)`   = sum(`Ghosts (delta cure)`)
    ,`Ghosts (early death)`  = sum(`Ghosts (early death)`)
    ,`Ghosts (onset death)`  = sum(`Ghosts (onset death)`)
    ,`Incidence (1 base)`  = sum(`Incidence (1 base)`)
    ,`Incidence (2 diagnosis)`  = sum(`Incidence (2 diagnosis)`)
    ,`Prevalence`  = sum(`Prevalence`)
    ,`1 in x families`  = sum(`1 in x families` *`Ann. background population` )/sum(`Ann. background population`)

  ),by=key_list]


  # country_data_age10 <- setDT(country_data_merge)[Age==10,list(
  country_data_age10 <- setDT(country_data_merge)[,list(

    `Life expectency (1 background)`      =  agg_life_expectancy_median(`Life expectency (1 background)`      , `Ann. background population`  )
    ,`Life expectency (2 t1d base)`       =  agg_life_expectancy_median(`Life expectency (2 t1d base)`        , Prevalence+Ghosts  )
    ,`Life expectency (3 t1d diagnosis)`  =  agg_life_expectancy_median(`Life expectency (3 t1d diagnosis)`   , Prevalence+Ghosts  )
    ,`Life expectency (4 t1d basic care)` =  agg_life_expectancy_median(`Life expectency (4 t1d basic care)`  , Prevalence+Ghosts  )
    ,`Life expectency (5 t1d best care)`  =  agg_life_expectancy_median(`Life expectency (5 t1d best care)`   , Prevalence+Ghosts  )
    ,`Life expectency (6 t1d cure)`       =  agg_life_expectancy_median(`Life expectency (6 t1d cure)`        , Prevalence+Ghosts  )

    ,`Lifetime years lost (2 t1d base) (complication)`        = agg_life_expectancy_median(`Lifetime years lost (2 t1d base) (complication)`              , Prevalence+Ghosts  )
    ,`Lifetime years lost (3 t1d diagnosis) (complication)`   = agg_life_expectancy_median(`Lifetime years lost (3 t1d diagnosis) (complication)`         , Prevalence+Ghosts  )
    ,`Lifetime years lost (4 t1d basic care) (complication)`  = agg_life_expectancy_median(`Lifetime years lost (4 t1d basic care) (complication)`        , Prevalence+Ghosts  )
    ,`Lifetime years lost (5 t1d best care) (complication)`   = agg_life_expectancy_median(`Lifetime years lost (5 t1d best care) (complication)`         , Prevalence+Ghosts  )
    ,`Lifetime years lost (6 t1d cure) (complication)`        = agg_life_expectancy_median(`Lifetime years lost (6 t1d cure) (complication)`              , Prevalence+Ghosts  )

    ,`Lifetime years lost (2 t1d base) (treatment)`        = agg_life_expectancy_median(`Lifetime years lost (2 t1d base) (treatment)`              , Prevalence+Ghosts  )
    ,`Lifetime years lost (3 t1d diagnosis) (treatment)`   = agg_life_expectancy_median(`Lifetime years lost (3 t1d diagnosis) (treatment)`         , Prevalence+Ghosts  )
    ,`Lifetime years lost (4 t1d basic care) (treatment)`  = agg_life_expectancy_median(`Lifetime years lost (4 t1d basic care) (treatment)`        , Prevalence+Ghosts  )
    ,`Lifetime years lost (5 t1d best care) (treatment)`   = agg_life_expectancy_median(`Lifetime years lost (5 t1d best care) (treatment)`         , Prevalence+Ghosts  )
    ,`Lifetime years lost (6 t1d cure) (treatment)`        = agg_life_expectancy_median(`Lifetime years lost (6 t1d cure) (treatment)`              , Prevalence+Ghosts  )
    ,`Life expectency (strip low)`   =   agg_life_expectancy_median(`Life expectency (strip low)`              , Prevalence+Ghosts  )
    ,`Life expectency (strip hig)`   =   agg_life_expectancy_median(`Life expectency (strip hig)`              , Prevalence+Ghosts  )
    ,`Lifetime years lost (strip low)`   =   agg_life_expectancy_median(`Lifetime years lost (strip low)`              , Prevalence+Ghosts  )
    ,`Lifetime years lost (strip hig)`   =   agg_life_expectancy_median(`Lifetime years lost (strip hig)`              , Prevalence+Ghosts  )
    ,`Life expectency (sensor low)`   =   agg_life_expectancy_median(`Life expectency (sensor low)`              , Prevalence+Ghosts  )
    ,`Life expectency (sensor hig)`   =   agg_life_expectancy_median(`Life expectency (sensor hig)`              , Prevalence+Ghosts  )
    ,`Lifetime years lost (sensor low)`   =   agg_life_expectancy_median(`Lifetime years lost (sensor low)`              , Prevalence+Ghosts  )
    ,`Lifetime years lost (sensor hig)`   =   agg_life_expectancy_median(`Lifetime years lost (sensor hig)`              , Prevalence+Ghosts  )
    ,`% Odds living to`   =   agg_life_expectancy_median(`% Odds living to`              , Prevalence+Ghosts  )

  ),by=key_list]

  country_data_age10$`Lifetime years lost (2 t1d base)`       <- country_data_age10$`Lifetime years lost (2 t1d base) (complication)` + country_data_age10$`Lifetime years lost (2 t1d base) (treatment)`
  country_data_age10$`Lifetime years lost (3 t1d diagnosis)`  <- country_data_age10$`Lifetime years lost (3 t1d diagnosis) (complication)` + country_data_age10$`Lifetime years lost (3 t1d diagnosis) (treatment)`
  country_data_age10$`Lifetime years lost (4 t1d basic care)` <- country_data_age10$`Lifetime years lost (4 t1d basic care) (complication)` + country_data_age10$`Lifetime years lost (4 t1d basic care) (treatment)`
  country_data_age10$`Lifetime years lost (5 t1d best care)`  <- country_data_age10$`Lifetime years lost (5 t1d best care) (complication)` + country_data_age10$`Lifetime years lost (5 t1d best care) (treatment)`
  country_data_age10$`Lifetime years lost (6 t1d cure)`       <- country_data_age10$`Lifetime years lost (6 t1d cure) (complication)` + country_data_age10$`Lifetime years lost (6 t1d cure) (treatment)`

  country_data_merge_global <- country_data %>% left_join(country_data_age10, by=key_list)
  country_data_merge_global <- country_data_merge_global %>% mutate_if(is.numeric, round, digits=2)

  country_data_merge_global
}

#' Construct prevalence and ghost population
#'
#' This function calculates prevalence for the given country, as well as the
#' 'ghost population'. For full details, see the working paper.
#'
#' Prevalence is calculated using \link{calculate_prevalence}(), where the
#' parameter matrices are derived from the callback functions given as
#' arguments. Where callbacks are not provided, defaults for
#' \code{country_wb_name} are used. The ghost population is calculated by taking
#' the difference between a counterfactual population with zero deaths at onset
#' and the T1D mortality rate the same as background mortality.
#'
#' The return value is a list with the following elements, where `N` is the
#' number of years of output:
#' \itemize{
#'   \item `P_level` - `N` x `MAX_AGE` matrix of prevalence (T1D) compartments,
#'         by year and age. Values are scaled by country population, so units
#'         are persons.
#'   \item `S_level` - `N` x `MAX_AGE` matrix of susceptible (healthy) compartments,
#'         by year and age. Values are scaled by country population, so units
#'         are persons.
#'   \item `P_cohorts_level` - `N` x `MAX_AGE` x `MAX_AGE` array of prevalence(T1D)
#'         compartments by year, age, and age at diagnosis. This is a
#'         more detailed breakdown of `P_level`, which we often refer to as a
#'         'cohort' breakdown. Values are scaled by country population, so units
#'         are persons.
#'   \item `I_flow` - `N` x `MAX_AGE` matrix of incidence flows, by year and age
#'         i.e. transitions from `S` (healthy/susceptible) to `P` (T1D/prevalent).
#'         Values are scaled by country population, so units are persons.
#'   \item `DDx_flow` - `N` x `MAX_AGE` matrix of death at onset flows by year and
#'         age. These are transitions directly from `S` (healthy/susceptible) to
#'         `D` (death). Terminology: this flow was historically called
#'         'death on diagnosis', hence the abbrevation DDx, although this was
#'         abandoned for reasons of accuracy. Values are scaled by country
#'         population, so units are persons.
#'   \item `DT1D_flow` - `N` x `MAX_AGE` matrix of excess deaths due to T1D by
#'         year and age. This is a subset of the flows from `P` (prevalent/T1D)
#'         to `D` over and above background deaths. Values are scaled by country
#'         population, so units are persons.
#'   \item `ghost_level` - `N` x `MAX_AGE` matrix of total ghost population,
#'         in units of persons.
#'   \item `ghost_ddx_level` - `N` x `MAX_AGE` matrix of ghost population due to
#'         deaths at onset, in units of persons.
#'   \item `ghost_hba1c_level` - `N` x `MAX_AGE` matrix of ghost population due to
#'         inadequate care (excess deaths due to T1D), in units of persons.
#'   \item `pop` - `N` x `MAX_AGE` matrix of ghost population
#'   \item `country` - World Bank name of country
#'   \item `time` - Elapsed time for the function to run
#'   \item `hba1c` - mean hba1c level
#'   \item `years` - vector of >=2 consecutive years to run the model for
#'   \item `pop_scale_factor` - population scaling factors used
#'}
#'
#' @param country_wb_name world bank name of country
#' @param hba1c average hba1c to model
#' @param incidence_fn Observed incidence function (excludes deaths on diagnosis)
#' @param bg_mort_fn Background mortality function
#' @param t1d_mort_fn Type-1 diabetes mortality function
#' @param population_fn Population function
#' @param death_on_diagnosis_fn Death on diagnosis rate function
#' @param start_year First year to model
#' @param end_year Last year to model
#' @param draw (integer) number of population trajectories to return
#' @param pop_scale_factor optional population scaling factor matrix
#' @param prev_reference optional prevalence reference level, a point estimate
#' @param prev_reference_year optional prevalence reference year
#' @return structure of class `prevalence` as described above
#' @export
prevalence_and_ghost_pop_1_0_0 <- function(
    country_wb_name,
    # hba1c=make_age_function(7.5),
    start_year,
    end_year,
    data_dir,
    smr_scale_factor=1,
    incidence_scale_factor=1,
    data_long = NULL,
    matrices_list = NULL,
    log_timing=FALSE
) {

  # smr_scale_factor=1 ;   incidence_scale_factor=1;data_long = NULL;matrices_list = NULL;log_timing=FALSE; start_year=1960; end_year=2040


  if(log_timing){ sink("log.txt",append=TRUE);cat(paste0(Sys.time()," prevalence_and_ghost_pop() \n") );sink()}

  ptm <- proc.time()
  # function to drop warm up periods
  dwu <- function(X) X[-(1:MAX_AGE),]
  # warm up period to seed prevalence: start model MAX_AGE years earlier
  years <- (start_year - MAX_AGE):end_year

  # Take draws from prevalence posterior distribution by simulating
  # Load rates from input_rates_combined -----------------------------------------------------------------
  # country_wb_name <- "Morocco"
  country_region <- Data_Run_query_return_df('SELECT loc_id,world_bank_name  FROM index_parameters.country_region')
  loc_id <- country_region$loc_id[country_region$world_bank_name==country_wb_name]

  # loc_id      <- get_loc_id(country_wb_name)

  # data_long   <- run_query_df (paste0("SELECT * FROM input_rates_combined  WHERE loc_id = '",loc_id,"'" ) )
  if(is.null(data_long))             {data_long   <- readRDS(paste0(data_dir,"/",loc_id,".Rds"))}
  data_long_default_run   <- readRDS(paste0(strsplit(data_dir,"_mcmc")[[1]][1],"/",loc_id,".Rds"))

  if(is.null(matrices_list))
  {
    matrices_list <- data_long_2_matrices (data_long=data_long,data_long_default_run=data_long_default_run)
    i                 <- matrices_list$i   # mean(i[as.character(2010:2019),19]) *100000
    qB                <- matrices_list$qB
    pop               <- matrices_list$pop
    dDx               <- matrices_list$dDx
    smr_matrix_n      <- matrices_list$smr_matrix_n
    smr_matrix_m      <- matrices_list$smr_matrix_m
    qT1D_percent_n    <- matrices_list$qT1D_percent_n
  }else
  { # Shiny input
    i                 <- matrices_list$i   # mean(i[as.character(2010:2019),19]) *100000
    qB                <- matrices_list$qB
    pop               <- matrices_list$pop
    dDx               <- matrices_list$dDx
    smr_matrix_n      <- matrices_list$smr_matrix_n
    smr_matrix_m      <- matrices_list$smr_matrix_m
    qT1D_percent_n    <- matrices_list$qT1D_percent_n
  }
  # set up parameters
  zero <- matrix_from_function(make_age_function(0), years)
  lever_year_range        <- as.character(config$lever_change_start_at:max(years))  #  only apply lever to years after lever_change_at


  # life expectency , background ---------------------------------------------------------------------
  life_table_background    <- dplyr::select(data_long, year,age,qx =background_mortality_rate)
  life_table_background    <- calculate_ex(life_table_background)

  # make sure t1d mortality ( smr * backgroun_mortality  ) is always <=1 -------
  Get_qT1D_from_smr_matrix <- function(qB, smr_matrix)
  {
    diab_odds  <- qB / (1 - qB) * smr_matrix
    qT1D       <- diab_odds / (1 + diab_odds)
    qT1D
  }
  #----------- life expectancy t1ds -----------------------------------------------------------------------
  Get_life_expectancy_t1d <- function(qB, smr_matrix_m,smr_matrix_n, qT1D_percent_n)
  {
    smr_matrix <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n
    smr_matrix_long    <- as.data.frame.table(smr_matrix,stringsAsFactors = F)%>% mutate_all(as.numeric) # when shiny run simulation, apply percent back to data_long
    colnames(smr_matrix_long) <- c("year",  "age",   "smr")

    t1d_mortality_matrix <- Get_qT1D_from_smr_matrix (qB, smr_matrix)

    t1d_mortality_long <- as.data.frame.table(t1d_mortality_matrix,stringsAsFactors = F)%>% mutate_all(as.numeric) # when shiny run simulation, apply percent back to data_long
    colnames(t1d_mortality_long) <- c("year",  "age",   "Value")
    t1d_mortality_long$smr <- smr_matrix_long$smr
    t1d_mortality_long <- t1d_mortality_long[t1d_mortality_long$year >= 1900,]
    t1d_mortality_long <- t1d_mortality_long[with(t1d_mortality_long, order( year,age)), ]

    t1d_mortality_long    <- dplyr::select(t1d_mortality_long, year,age,smr,qx =Value)
    life_table_t1d        <- calculate_ex_lifetime_years_lost(t1d_mortality_long)
    life_table_t1d
  }

  life_table_base_scenario <- Get_life_expectancy_t1d (qB, smr_matrix_m,smr_matrix_n, qT1D_percent_n)

  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)

  # prev base scenario  --------------------------------------------------------------------
  prev                    <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,  years)
  # prev                    <- calculate_prevalence_gpu(i        , qB, qT1D_n, qT1D_m, qT1D_percent_n, dDx ,  years)


  matrices_list_lever1 <- Apply_levers_to_input_matrices (matrices_list, lever=1, lever_year_range=lever_year_range)

  i                 <- matrices_list_lever1$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list_lever1$dDx
  smr_matrix_n      <- matrices_list_lever1$smr_matrix_n
  smr_matrix_m      <- matrices_list_lever1$smr_matrix_m
  qT1D_percent_n    <- matrices_list_lever1$qT1D_percent_n

  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)

  prev_100d               <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,  years) # deaths on diagnosis are converted to higer incidence
  # prev lever 2: basic care Insulin, strips and education SMR  median( c(3.7, 4.4)) , 4.05--------------------------------------------------------------------

  matrices_list_lever2 <- Apply_levers_to_input_matrices (matrices_list, lever=2, lever_year_range=lever_year_range)

  i                 <- matrices_list_lever2$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list_lever2$dDx
  smr_matrix_n      <- matrices_list_lever2$smr_matrix_n
  smr_matrix_m      <- matrices_list_lever2$smr_matrix_m
  qT1D_percent_n    <- matrices_list_lever2$qT1D_percent_n

  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
  prev_100d_basic_care    <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,  years) # smr all non-minimal

  life_table_lever_2     <- Get_life_expectancy_t1d (qB, smr_matrix_m,smr_matrix_n, qT1D_percent_n)


  # # prev best care Pumps and FGMs/CGMs SMR median( c(2.2, 2.6)) , median 2.4--------------------------------------------------------------------

  matrices_list_lever3 <- Apply_levers_to_input_matrices (matrices_list, lever=3, lever_year_range=lever_year_range)

  i                 <- matrices_list_lever3$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list_lever3$dDx
  smr_matrix_n      <- matrices_list_lever3$smr_matrix_n
  smr_matrix_m      <- matrices_list_lever3$smr_matrix_m
  qT1D_percent_n    <- matrices_list_lever3$qT1D_percent_n

  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)

  prev_100d_best_care     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,  years) # deaths on diagnosis are converted to higer incidence

  life_table_lever_3     <- Get_life_expectancy_t1d (qB, smr_matrix_m,smr_matrix_n, qT1D_percent_n)


  # counterfactual: ghost pop is diff between prevalence and counterfactual where deaths on diagnosis are converted to incidence, and hba1c is fully controlled (4.3)
  # prev cure care --------------------------------------------------------------------

  matrices_list_lever4 <- Apply_levers_to_input_matrices (matrices_list, lever=4, lever_year_range=lever_year_range)

  i                 <- matrices_list_lever4$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list_lever4$dDx
  smr_matrix_n      <- matrices_list_lever4$smr_matrix_n
  smr_matrix_m      <- matrices_list_lever4$smr_matrix_m
  qT1D_percent_n    <- matrices_list_lever4$qT1D_percent_n

  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)

  prev_100d_cure          <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,  years)


  life_table_lever_4     <- Get_life_expectancy_t1d (qB, smr_matrix_m,smr_matrix_n, qT1D_percent_n)

  # calculating years gained --------
  matrices_list_strips_low  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 5.3,6.3,lever_year_range)
  matrices_list_strips_hig <- Apply_smr_to_input_matrices(matrices_list, lever=1, 3.1,3.7,lever_year_range)
  matrices_list_sensor_low  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 4.0,4.7,lever_year_range)
  matrices_list_sensor_hig <- Apply_smr_to_input_matrices(matrices_list, lever=1, 2.2,2.8,lever_year_range)

  life_table_strips_low     <- Get_life_expectancy_t1d (qB, matrices_list_strips_low$smr_matrix_m,matrices_list_strips_low$smr_matrix_n, matrices_list_strips_low$qT1D_percent_n)
  life_table_strips_hig     <- Get_life_expectancy_t1d (qB, matrices_list_strips_hig$smr_matrix_m,matrices_list_strips_hig$smr_matrix_n, matrices_list_strips_hig$qT1D_percent_n)
  life_table_sensor_low     <- Get_life_expectancy_t1d (qB, matrices_list_sensor_low$smr_matrix_m,matrices_list_sensor_low$smr_matrix_n, matrices_list_sensor_low$qT1D_percent_n)
  life_table_sensor_hig     <- Get_life_expectancy_t1d (qB, matrices_list_sensor_hig$smr_matrix_m,matrices_list_sensor_hig$smr_matrix_n, matrices_list_sensor_hig$qT1D_percent_n)

  # write.csv(life_table_sensor_low[life_table_sensor_low$year==2022],paste0("temp/life_table_sensor_low_",country_wb_name,".csv") )
  # write.csv(life_table_sensor_hig[life_table_sensor_hig$year==2022],paste0("temp/life_table_sensor_hig_",country_wb_name,".csv") )

  # ----------------------------------------------------------------------------------------------------------------------------------------
  # translate proportions to population levels
  # P + S + D = pop_scale_factor
  # start simulation from 1800 ----------------------------------------------------

  pop_scale_factor       <- pop / (prev$P + prev$S)   # exclude the dead from population
  pop_scale_factor_100d  <- pop / (prev_100d$P + prev_100d$S)   # exclude the dead from population
  pop_scale_factor_100d_bacare  <- pop / (prev_100d_basic_care$P + prev_100d_basic_care$S)   # exclude the dead from population
  pop_scale_factor_100d_becare  <- pop / (prev_100d_best_care$P + prev_100d_best_care$S)   # exclude the dead from population
  pop_scale_factor_100d_cure    <- pop / (prev_100d_cure$P + prev_100d_cure$S)   # exclude the dead from population

  P_level         <- pop_scale_factor * prev$P  #  sum(P_level["2021",]) ; sum( (pop_scale_factor_100d * prev_100d$P) ["2021",])
  S_level         <- pop_scale_factor * prev$S
  I_flow          <- pop_scale_factor * prev$I

  I_flow_diagnosis           <- pop_scale_factor * prev_100d$I
  I_flow_basic_care          <- pop_scale_factor * prev_100d_basic_care$I
  I_flow_best_care           <- pop_scale_factor * prev_100d_best_care$I
  I_flow_cure                <- pop_scale_factor * prev_100d_cure$I

  DDx_flow                   <- pop_scale_factor * prev$DDx   # sum(DDx_flow["2021",])

  DT1D_flow                  <- pop_scale_factor * prev$DT1D
  DT1D_flow_100d             <- pop_scale_factor * prev_100d$DT1D
  DT1D_flow_basic_care       <- pop_scale_factor * prev_100d_basic_care$DT1D
  DT1D_flow_best_care        <- pop_scale_factor * prev_100d_best_care$DT1D
  DT1D_flow_cure             <- pop_scale_factor * prev_100d_cure$DT1D

  BD_flow         <- pop * qB

  # Calculate missing prevalence  and delta of levers. --------------------------------------------------------------------------------
  ghost_ddx_level   <- (prev_100d$P              - prev$P) * pop_scale_factor #
  ghost_hba1c_level <- (prev_100d_cure$P         - prev_100d$P) * pop_scale_factor #

  ghost_basic_care  <- (prev_100d_basic_care$P   - prev_100d$P) * pop_scale_factor #  delta
  ghost_best_care   <- (prev_100d_best_care$P    - prev_100d_basic_care$P) * pop_scale_factor #
  ghost_cure        <- (prev_100d_cure$P         - prev_100d_best_care$P) * pop_scale_factor #  sum(ghost_hba1c_level["2021",])

  # # Calculate missing prevalence  and delta of levers. --------------------------------------------------------------------------------
  # ghost_ddx_level   <- prev_100d$P *pop_scale_factor_100d             - prev$P* pop_scale_factor   #
  # ghost_hba1c_level <- prev_100d_cure$P   *pop_scale_factor_100d_cure       - prev_100d$P *pop_scale_factor_100d
  #
  # ghost_basic_care  <- prev_100d_basic_care$P *pop_scale_factor_100d_bacare   - prev_100d$P *pop_scale_factor_100d
  # ghost_best_care   <- prev_100d_best_care$P *pop_scale_factor_100d_becare    - prev_100d_basic_care$P *pop_scale_factor_100d_bacare
  # ghost_cure        <- prev_100d_cure$P *pop_scale_factor_100d_cure         - prev_100d_best_care$P *pop_scale_factor_100d_becare

  # ghost population is the sum of these two
  ghost_level <- ghost_hba1c_level + ghost_ddx_level #   sum(ghost_level["2021",])
  print(paste0("ghost onset death: ",round(sum(ghost_ddx_level["2021",]),2),"; ghost early death: ",round(sum(ghost_hba1c_level["2021",]),2) ) )

  # df <- data.frame(x=rownames(ghost_level)
  #                  ,diagnosed           = (rowSums(P_level))
  #                  ,ghost_undiagnosed   = (rowSums(ghost_ddx_level))
  #                  ,ghost_diagnosed     = (rowSums(ghost_hba1c_level))
  #                  ,ghost_basic_care    = (rowSums(ghost_basic_care))
  #                  ,ghost_best_care     = (rowSums(ghost_best_care))
  #                  ,ghost_cure          = (rowSums(ghost_cure))
  #                  )
  # df[df$x>=1960 & df$x<=2020,] %>%
  #   e_charts(x) %>%
  #   e_bar(diagnosed, stack = "grp") %>%
  #   e_bar(ghost_undiagnosed, stack = "grp") %>%
  #   # e_bar(ghost_diagnosed, stack = "grp") %>%
  #   e_bar(ghost_basic_care, stack = "grp") %>%
  #   e_bar(ghost_best_care, stack = "grp") %>%
  #   e_bar(ghost_cure, stack = "grp") %>%
  #   e_title("T1D population",country_wb_name)%>%
  #   e_legend(bottom=0)
  # complications , burden

  dalys <- list()
  dalys$year_age <-  dwu(P_level)

  if(config$run_days_lost) # do not run burden,  days lost same as prevalance----
  {
    P_cohorts_level <- prev      $Pcohorts * array(pop_scale_factor, c(length(years),MAX_AGE,MAX_AGE))
    P_cohorts_level <- P_cohorts_level[-(1:MAX_AGE),,]
    comp  <- complication_prevalence(years=seq(start_year, end_year),P_cohorts_level,matrices_list$smr_matrix_n, matrices_list$smr_matrix_m,matrices_list$qT1D_percent_n)
    dalys <- calculate_dalys(P_cohorts_level, comp)

  }

  # test  <- apply(P_cohorts_level, c(1, 2), sum) # check if cohort is adds up to prevalence.




  dalys_100d              <- dalys
  dalys_100d_basic_care   <- dalys
  dalys_100d_best_care    <- dalys
  dalys_100d_cure         <- dalys

  if(config$run_days_lost_lever)  # do not run levers for paper stats
  {
    P_cohorts_level <- prev_100d$Pcohorts * array(pop_scale_factor, c(length(years),MAX_AGE,MAX_AGE))
    P_cohorts_level <- P_cohorts_level[-(1:MAX_AGE),,]
    # comp  <- complication_prevalence(start_year, end_year,P_cohorts_level,smr_matrix_n, smr_matrix_m,qT1D_percent_n)
    dalys_100d <- calculate_dalys(P_cohorts_level, comp) # same complication , as smr did not change

    P_cohorts_level <- prev_100d_basic_care$Pcohorts * array(pop_scale_factor, c(length(years),MAX_AGE,MAX_AGE))
    P_cohorts_level <- P_cohorts_level[-(1:MAX_AGE),,]
    comp  <- complication_prevalence(years=seq(start_year, end_year),P_cohorts_level,matrices_list_lever2$smr_matrix_n, matrices_list_lever2$smr_matrix_m,matrices_list_lever2$qT1D_percent_n)
    dalys_100d_basic_care <- calculate_dalys(P_cohorts_level, comp)


    P_cohorts_level <- prev_100d_best_care$Pcohorts * array(pop_scale_factor, c(length(years),MAX_AGE,MAX_AGE))
    P_cohorts_level <- P_cohorts_level[-(1:MAX_AGE),,]
    comp  <- complication_prevalence(years=seq(start_year, end_year),P_cohorts_level,matrices_list_lever3$smr_matrix_n, matrices_list_lever3$smr_matrix_m,matrices_list_lever3$qT1D_percent_n)
    dalys_100d_best_care <- calculate_dalys(P_cohorts_level, comp)

    P_cohorts_level <- prev_100d_cure$Pcohorts * array(pop_scale_factor, c(length(years),MAX_AGE,MAX_AGE))
    P_cohorts_level <- P_cohorts_level[-(1:MAX_AGE),,]
    comp  <- complication_prevalence(years=seq(start_year, end_year),P_cohorts_level,matrices_list_lever4$smr_matrix_n, matrices_list_lever4$smr_matrix_m,matrices_list_lever4$qT1D_percent_n)
    dalys_100d_cure <- calculate_dalys(P_cohorts_level, comp)

  }

  # merge all outputs ------------------------------------------------------------------------------
  prev_merge <-                     cbind(Value_type="Ann. background mortality" , Year = rownames(dwu(pop))      , as.data.frame(dwu(BD_flow) ))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. background population", Year = rownames(dwu(pop))      , as.data.frame(dwu(pop) )))

  prev_merge <- rbind( prev_merge,  cbind(Value_type="Prevalence"                , Year = rownames(dwu(pop))      , as.data.frame(dwu(P_level))))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Incidence (1 base)"        , Year = rownames(dwu(pop))      , as.data.frame(dwu(I_flow))))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Incidence (2 diagnosis)"   , Year = rownames(dwu(pop))      , as.data.frame(dwu(I_flow_diagnosis))))

  # not needed, same value,
  # prev_merge <- rbind( prev_merge,  cbind(Value_type="Incidence (3 basic care)"   , Year = rownames(dwu(pop))      , as.data.frame(dwu(I_flow_basic_care))))
  # prev_merge <- rbind( prev_merge,  cbind(Value_type="Incidence (4 best care)"   , Year = rownames(dwu(pop))      , as.data.frame(dwu(I_flow_best_care))))
  # prev_merge <- rbind( prev_merge,  cbind(Value_type="Incidence (5 cure)"   , Year = rownames(dwu(pop))      , as.data.frame(dwu(I_flow_cure))))


  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ghosts"                    , Year = rownames(dwu(pop))      , as.data.frame(dwu(ghost_level)) ))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ghosts (onset death)"      , Year = rownames(dwu(pop))      , as.data.frame(dwu(ghost_ddx_level) ) ))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ghosts (early death)"      , Year = rownames(dwu(pop))      , as.data.frame(dwu(ghost_hba1c_level) ) ))

  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ghosts (delta basic care)" , Year = rownames(dwu(pop))      , as.data.frame(dwu(ghost_basic_care) )))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ghosts (delta best care)"  , Year = rownames(dwu(pop))      , as.data.frame(dwu(ghost_best_care) )))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ghosts (delta cure)"       , Year = rownames(dwu(pop))      , as.data.frame(dwu(ghost_cure) )))

  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. onset deaths"         , Year = rownames(dwu(pop))      , as.data.frame(dwu(DDx_flow) )))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. early deaths"         , Year = rownames(dwu(pop))      , as.data.frame(dwu(DT1D_flow) )))

  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. early deaths (2 diagnosis)"         , Year = rownames(dwu(pop))      , as.data.frame(dwu(DT1D_flow_100d) )))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. early deaths (3 basic care)"         , Year = rownames(dwu(pop))      , as.data.frame(dwu(DT1D_flow_basic_care) )))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. early deaths (4 best care)"         , Year = rownames(dwu(pop))      , as.data.frame(dwu(DT1D_flow_best_care) )))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. early deaths (5 cure)"         , Year = rownames(dwu(pop))      , as.data.frame(dwu(DT1D_flow_cure) )))

  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. days lost (1 base)"       , Year = rownames(dwu(pop))      , as.data.frame(dalys$year_age) ))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. days lost (2 diagnosis)"  , Year = rownames(dwu(pop))      , as.data.frame(dalys_100d$year_age) ))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. days lost (3 basic care)" , Year = rownames(dwu(pop))      , as.data.frame(dalys_100d_basic_care$year_age) ))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. days lost (4 best care)"  , Year = rownames(dwu(pop))      , as.data.frame(dalys_100d_best_care$year_age) ))
  prev_merge <- rbind( prev_merge,  cbind(Value_type="Ann. days lost (5 cure)"       , Year = rownames(dwu(pop))      , as.data.frame(dalys_100d_cure$year_age) ))

  prev_merge_long      <- gather(prev_merge ,key = "Age", value = "Value",-Value_type,-Year)
  prev_merge_long$Year <- as.numeric(prev_merge_long$Year)
  prev_merge_long$Age  <- as.numeric(prev_merge_long$Age)
  prev_merge_long      <- cbind(Country= country_wb_name,Type="Country",prev_merge_long)

  # saveRDS(prev_merge_long,"temp/test.Rds")
  prev_merge_wide   <- spread(prev_merge_long, Value_type, Value )

  prev_merge_wide$`Life expectency (1 background)`      <- life_table_background$ex[life_table_background$year>=1960]

  # for t1d base life expectancy, include diagnosis rate effect in the begining
  prev_merge_wide$`Life expectency (2 t1d base)`        <- life_table_base_scenario$ex[life_table_base_scenario$year>=1960]
  diagnosis_rate <- (prev_merge_wide$`Incidence (1 base)`+0.001) / (prev_merge_wide$`Ann. onset deaths`  + prev_merge_wide$`Incidence (1 base)`+0.001)
  prev_merge_wide$`Life expectency (2 t1d base)`        <- prev_merge_wide$`Life expectency (2 t1d base)` * diagnosis_rate +  (1-diagnosis_rate)*0.5

  prev_merge_wide$`Life expectency (3 t1d diagnosis)`   <- life_table_base_scenario$ex[life_table_base_scenario$year>=1960]
  prev_merge_wide$`Life expectency (4 t1d basic care)`  <- life_table_lever_2$ex[life_table_lever_2$year>=1960]
  prev_merge_wide$`Life expectency (5 t1d best care)`   <- life_table_lever_3$ex[life_table_lever_3$year>=1960]
  prev_merge_wide$`Life expectency (6 t1d cure)`        <- life_table_lever_4$ex[life_table_background$year>=1960]

  prev_merge_wide$`Lifetime years lost (2 t1d base) (complication)`        <- life_table_base_scenario$ex_complication[life_table_base_scenario$year>=1960] *  prev_merge_wide$`Life expectency (2 t1d base)`/ life_table_base_scenario$ex[life_table_base_scenario$year>=1960]
  prev_merge_wide$`Lifetime years lost (3 t1d diagnosis) (complication)`   <- life_table_base_scenario$ex_complication[life_table_base_scenario$year>=1960]
  prev_merge_wide$`Lifetime years lost (4 t1d basic care) (complication)`  <- life_table_lever_2$ex_complication[life_table_lever_2$year>=1960]
  prev_merge_wide$`Lifetime years lost (5 t1d best care) (complication)`   <- life_table_lever_3$ex_complication[life_table_lever_3$year>=1960]
  prev_merge_wide$`Lifetime years lost (6 t1d cure) (complication)`        <- life_table_lever_4$ex_complication[life_table_background$year>=1960]

  prev_merge_wide$`Lifetime years lost (2 t1d base) (treatment)`        <- life_table_base_scenario$ex_treatment[life_table_base_scenario$year>=1960] *  prev_merge_wide$`Life expectency (2 t1d base)`/ life_table_base_scenario$ex[life_table_base_scenario$year>=1960]
  prev_merge_wide$`Lifetime years lost (3 t1d diagnosis) (treatment)`   <- life_table_base_scenario$ex_treatment[life_table_base_scenario$year>=1960]
  prev_merge_wide$`Lifetime years lost (4 t1d basic care) (treatment)`  <- life_table_lever_2$ex_treatment[life_table_lever_2$year>=1960]
  prev_merge_wide$`Lifetime years lost (5 t1d best care) (treatment)`   <- life_table_lever_3$ex_treatment[life_table_lever_3$year>=1960]
  prev_merge_wide$`Lifetime years lost (6 t1d cure) (treatment)`        <- life_table_lever_4$ex_treatment[life_table_background$year>=1960]


  prev_merge_wide$`Life expectency (strip low)`        <- life_table_strips_low$ex[life_table_strips_low$year>=1960]
  prev_merge_wide$`Life expectency (strip hig)`        <- life_table_strips_hig$ex[life_table_strips_hig$year>=1960]
  prev_merge_wide$`Lifetime years lost (strip low)`        <- (life_table_strips_low$ex_complication + life_table_strips_low$ex_treatment)  [life_table_strips_low$year>=1960]
  prev_merge_wide$`Lifetime years lost (strip hig)`        <- (life_table_strips_hig$ex_complication + life_table_strips_hig$ex_treatment)  [life_table_strips_hig$year>=1960]


  prev_merge_wide$`Life expectency (sensor low)`        <- life_table_sensor_low$ex[life_table_sensor_low$year>=1960]
  prev_merge_wide$`Life expectency (sensor hig)`        <- life_table_sensor_hig$ex[life_table_sensor_hig$year>=1960]
  prev_merge_wide$`Lifetime years lost (sensor low)`        <- (life_table_sensor_low$ex_complication + life_table_sensor_low$ex_treatment)  [life_table_sensor_low$year>=1960]
  prev_merge_wide$`Lifetime years lost (sensor hig)`        <- (life_table_sensor_hig$ex_complication + life_table_sensor_hig$ex_treatment)  [life_table_sensor_hig$year>=1960]

  prev_merge_wide$`% Odds living to`  <- prev_merge_wide$Prevalence / (prev_merge_wide$Prevalence + prev_merge_wide$Ghosts) * 100
  prev_merge_wide$`% Odds living to`[is.na(prev_merge_wide$`% Odds living to`)] <- 0

  JDRF_calculations <- extract_for_purpose(prev_merge_wide,0,99,2022)
  days_ <- JDRF_calculations$days_1
  prev_merge_wide   <- prev_merge_wide %>% dplyr::inner_join(dplyr::select(days_,Country,Year,`1 in x families` = `1 in x families`
                                                                           # ,`Lifetime years lost (2 t1d base)`=lifetime_years_lost_complication_and_treatment
                                                                           # ,`Lifetime years lost (2 t1d base) (complication)`=lifetime_years_lost_complication
                                                                           # ,`Lifetime years lost (2 t1d base) (treatment)`   =lifetime_years_lost_treatment
  ), by=c("Country","Year"))


  # days_ <- JDRF_calculations$days_2
  # prev_merge_wide   <- prev_merge_wide %>% dplyr::inner_join(dplyr::select(days_,Country,Year
  #                                                                          ,`Lifetime years lost (3 t1d diagnosis)`=lifetime_years_lost_complication_and_treatment
  #                                                                          ,`Lifetime years lost (3 t1d diagnosis) (complication)`=lifetime_years_lost_complication
  #                                                                          ,`Lifetime years lost (3 t1d diagnosis) (treatment)`   =lifetime_years_lost_treatment
  #                                                                          ), by=c("Country","Year"))
  #
  #
  # days_ <- JDRF_calculations$days_3
  # prev_merge_wide   <- prev_merge_wide %>% dplyr::inner_join(dplyr::select(days_,Country,Year
  #                                                                          ,`Lifetime years lost (4 t1d basic care)`=lifetime_years_lost_complication_and_treatment
  #                                                                          ,`Lifetime years lost (4 t1d basic care) (complication)`=lifetime_years_lost_complication
  #                                                                          ,`Lifetime years lost (4 t1d basic care) (treatment)`   =lifetime_years_lost_treatment
  #                                                                          ), by=c("Country","Year"))
  #
  # days_ <- JDRF_calculations$days_4
  # prev_merge_wide   <- prev_merge_wide %>% dplyr::inner_join(dplyr::select(days_,Country,Year
  #                                                                          ,`Lifetime years lost (5 t1d best care)`=lifetime_years_lost_complication_and_treatment
  #                                                                          ,`Lifetime years lost (5 t1d best care) (complication)`=lifetime_years_lost_complication
  #                                                                          ,`Lifetime years lost (5 t1d best care) (treatment)`   =lifetime_years_lost_treatment
  #                                                                          ), by=c("Country","Year"))
  #
  # days_ <- JDRF_calculations$days_5
  # prev_merge_wide   <- prev_merge_wide %>% dplyr::inner_join(dplyr::select(days_,Country,Year
  #                                                                          ,`Lifetime years lost (6 t1d cure)`=lifetime_years_lost_complication_and_treatment
  #                                                                          ,`Lifetime years lost (6 t1d cure) (complication)`=lifetime_years_lost_complication
  #                                                                          ,`Lifetime years lost (6 t1d cure) (treatment)`   =lifetime_years_lost_treatment
  #                                                                          ), by=c("Country","Year"))

  prev_merge_wide <- prev_merge_wide %>% mutate_if(is.numeric, round, digits=2)





  # if(FALSE)
  # {
  #   prev_merge_wide_old <- read_parquet( "C:/DropboxT1D/reruns_scenario_2/rerun_0.4.15_0/AFG.binary")
  #   prev_merge_wide_old2 <- read_parquet( "C:/DropboxT1D/reruns_scenario_2/rerun_0.4.15.lever2023_0/AFG.binary")
  # }

  # for(i in 1:48)
  # {
  #   print(i)
  #  print( sum((prev_merge_wide_old!=prev_merge_wide_old2)[,i]) )
  # }


  time_taken <- proc.time() - ptm

  structure(with(prev,
                 list(
                   prev_merge_wide=prev_merge_wide
                 )),
            class='prevalence')
}
