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
calculate_prevalence <- function(i, qB,  qT1D_n, qT1D_m, qT1D_percent_n, dDx, smr_matrix, years, run_complications=FALSE) {
  # run_complications <- TRUE  ; run_complications <- FALSE
  weibull <- data.frame()
  if(run_complications)  # not run complications
  {

    hba1c_full   <- (log(smr_matrix) +  1.5274 )/ 0.3545
    weibull      <- get_complication_parameters()$weibull
    # weibull <- weibull[1:2,]
    weib_survival_props <- function(hba1c_year,age_pair, weibull) {
      tail <- as.numeric(colnames(age_pair))
      head <- round(age_pair,0)
      props_full     <- array(0,dim=c(nrow(weibull)+1, length(age_pair)), dimnames = list(comp=c(weibull$abbrev,"AtLeast1C"), age=1:length(age_pair)))

      for(a in 1:length( age_pair) )
      { # a <- 1
        if(tail[a] >  head[a])
        {
          hba1c <- hba1c_year
          hba1c <- hba1c[head[a]:tail[a] ]
          # weibull=complication_parameters$weibull  ;  hba1c <- c(10,10,10,10,10,10)
          props      <- array(0,dim=c(nrow(weibull)+1, length(hba1c)), dimnames = list(comp=c(weibull$abbrev,"AtLeast1C"), age=1:length(hba1c)))
          for(i_ in 1:nrow(weibull))
          {
            intercept <- weibull$intercept[i_]
            slope     <- weibull$slope[i_]
            scale     <- weibull$scale[i_]
            T <- length(hba1c)-1 # total time length (years) incl. T=0
            hs_odd <- 0.5 * (hba1c[-1] + hba1c[-(T+1)]) # interpolate at half-year points
            ts <- pmin(0:T, 30)  # constant hazard after 30 years
            # Simpson's rule
            l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c)/scale) / scale
            l_odd <- (ts[-1] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale
            Lambda <- 1/6 * (l_even + cumsum(2 * c(0, 0, l_even[-c(1, T+1)]) + 4 * c(0, l_odd[-(T+1)])))
            exp(-Lambda) # survival function
            # print( exp(-Lambda) )
            props[weibull$abbrev[i_],] <- 1- exp(-Lambda)
          }
          # This is Graham's and Gabriel's 'journeys' adjustments. Reduces ON and PR by survivor shares of RF and BL, respectively,
          # effectively making this a multistate model (rather than a pure time-to-event model)
          props['ON',] <- props['ON',] * (1 - props['RF',])
          props['PR',] <- props['PR',] * (1 - props['BL',])
          props['AtLeast1C',] <- 1-(apply((1-props), 2, prod))
          props[,length(hba1c),drop=FALSE]
          props_full[,a] <-  props[,length(hba1c),drop=FALSE]
        }
      }
      props_full
    }
  }
  # Data matrixes: years on rows, ages on columns. Memory footprint of storing
  # all this information really is tiny, a few MB. Similarly, a 3D array: {year}x{age}x{onset age}
  # Model compartments, all annual cohorts as a proportion of 1.

  S <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # S - susceptible (non-T1D) ,equals SS0+SS1+SS2
  P <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # P - prevalence (T1D)
  D <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # D - deaths (absorbing state)


  O    <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # O , mean age of Onset
  C_I  <-  array(0, dim=list(nrow(weibull)+1, length(years), MAX_AGE), dimnames=list(c(weibull$abbrev,"AtLeast1C"),years, AGES))  #  complication incidence  probability
  C_P  <-  array(NA, dim=list(nrow(weibull)+1, length(years), MAX_AGE), dimnames=list(c(weibull$abbrev,"AtLeast1C"),years, AGES)) #  complication prevalence cases


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

  O[,1] <- 0
  C_P[,1,1] <- 0

  t_loop <- seq_along(years[-1])

  for (t in t_loop  ) {  # t <- 131   ;   t <- 1
    # for (t in 1:129) {  # t <- 135   ;   t <- 2
      # for (t in 1:90) {
      # bump previous period's data along one cell (ie one year older)
    Sshift     <- c(1, S[t, -MAX_AGE]) ; names(Sshift) <- as.character(0:99)
    Pshift     <- c(0, P[t, -MAX_AGE]) ; names(Pshift) <- as.character(0:99)
    Dshift     <- c(0, D[t, -MAX_AGE]) ; names(Dshift) <- as.character(0:99)

    # equation (3) - susceptible compartment, S
    S[t+1,] <- Sshift * (1 - i_all_shift[t,]) * (1 - qB[t,])

    # print(S[1:15,1:15])
    # equation (4) - prevalence compartment P
    # P[t+1,] <- Pshift * (1 - qT1D[t,]) + i_shift[t,] * Sshift
    P[t+1,] <-    qT1D_percent_n[t,]  *  Pshift  * (1 - qT1D_n[t,])  +
               (1-qT1D_percent_n[t,]) *  Pshift  * (1 - qT1D_m[t,])  + i_shift[t,] * Sshift


    if(run_complications)  # not run complications
    {
      Oshift <- c(0, O[t, -MAX_AGE]) ; names(Oshift) <- as.character(0:99) # mean age of onset plus 1 with shifting

      C_P_shift <- cbind(rep(0, nrow(weibull)+1), C_P[,t, -MAX_AGE]) ; colnames(C_P_shift) <- as.character(0:99) #

      index_not_na <- !is.na(Pshift) ; index_not_na[1] <- FALSE

      if(sum(index_not_na))
      {
        O[t+1,index_not_na] <- ( Oshift[index_not_na] * Pshift[index_not_na] +  c(0:99) [index_not_na]* i_shift[t,index_not_na]* Sshift[index_not_na] + 1e-20 ) / ( Pshift[index_not_na] + i_shift[t,index_not_na]* Sshift[index_not_na] + 1e-20 )

        C_I[,t+1,index_not_na]  <- weib_survival_props (hba1c_year= hba1c_full[t,], age_pair = O[t+1,index_not_na,drop=FALSE], weibull)

        C_P[,t+1,]  <-  qT1D_percent_n[t,]  *  C_P_shift  * (1 - qT1D_n[t,])  +
          (1-qT1D_percent_n[t,]) *  C_P_shift  * (1 - qT1D_m[t,])  +   C_I[,t+1,] * (matrix(rep(P[t+1,], 11), nrow=11, byrow=TRUE) - C_P_shift  )
      }
    }
    # sweep(C_I[,t+1,], MARGIN = 2, STATS = P[t+1,], FUN = "*")
    # equation (5) - death compartment D
    D[t+1,] <- (Dshift
                + i_all_shift[t,] * dDx[t,] * Sshift
                # + qT1D[t,] * Pshift
                + qT1D_n[t,] * Pshift *       qT1D_percent_n[t,]
                + qT1D_m[t,] * Pshift *  (1 - qT1D_percent_n[t,])
                + Sshift * qB[t,] * (1 - i_all_shift[t,]))
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

  }
  C_P[,"1960","99"] <- 0 # Replace NA becasuse of 1 year lagging due to simulation start years

  # final year cohort flows - same as in the loop but for final period
  t <- length(years)
  Icohorts[t,,] <- diag(i[t,] * S[t,])

  # PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D[t,], c(MAX_AGE, MAX_AGE))
  PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D_n[t,], c(MAX_AGE, MAX_AGE)) *      qT1D_percent_n[t,]  +
                    Pcohorts[t,,] * array(qT1D_m[t,], c(MAX_AGE, MAX_AGE)) * (1 - qT1D_percent_n[t,])

  # flows: based on unshifted versions of incidence
  I   <- i * S                    # T1D incidence
  DDx <- i_all * dDx * S        # deaths at T1D onset
  # DDx <- i_all * S - I        # deaths at T1D onset
  # DT1D <- (qT1D - qB) * P       # T1D-cause mortalitys
  DT1D <- (qT1D_n - qB) * P  *          qT1D_percent_n  +
          (qT1D_m - qB) * P  *     (1 - qT1D_percent_n)   # T1D-cause mortality


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

  list(S=S, P=P, Pcohorts=Pcohorts, D=D, I=I, DDx=DDx, DT1D=DT1D,DBGP=DBGP, DBGS=DBGS,C_P=C_P)
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

  # data_wide <- spread(dplyr::select(data_long,year,age,value=gdp_base_2010 ), age, value)
  # rownames(data_wide) <- data_wide$year
  # data_wide           <- dplyr::select(data_wide,-year)%>% as.matrix
  # gdp                  <- data_wide

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
    rate_avg_5[is.infinite(rate_avg_5)] <- 0

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
    # if( !grepl("ward", config$version_no_input) ) # do not project for ward simulation---
    # {
      i   <- apply_Growth_Rate(matrix_rate=i,matrix_rate_default=i_default,year_start=2021)
      dDx <- apply_Growth_Rate(matrix_rate=dDx,year_start=2021)

    # }

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

pre_fill_constant_1860_1900 <- function (rate,year_start=1860)
{ # rate <- i ;
  years <- year_start
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
  # lever <- 2 ; lever_year_range = as.character(2023:max(years))
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
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1

    smr_matrix_n_new <- smr_matrix_n
    smr_matrix_n_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]
    # smr_matrix_m_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]

    qT1D_percent_n_new <- qT1D_percent_n
    qT1D_percent_n_new[as.character(lever_year_range),] <- 1

  }

  if(lever>=2.5)
  { # best care Pumps and FGMs/CGMs SMR median( c(2.2, 2.6)) , median 2.4----every one get basic care . 20 get devices. 100 diagnosis ,----------------------------------------------------------------

    smr_matrix_m_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]

    ratio        <- 2.4/(rowSums(smr_matrix_n[,as.character(20:50)])/length(20:50))
    ratio_matrix <- t(matrix(ratio, nrow=100, ncol=length(ratio), byrow=TRUE))
    ratio_matrix[ratio_matrix>=1] <- 1
    smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1

    smr_matrix_n_new <- smr_matrix_n
    smr_matrix_n_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]
    # smr_matrix_m_new <- smr_matrix_m

    qT1D_percent_n_new <- qT1D_percent_n
    qT1D_percent_n_new[as.character(lever_year_range),] <- 0.2

  }

  if(lever>=3)
  { # best care Pumps and FGMs/CGMs SMR median( c(2.2, 2.6)) , median 2.4--------------------------------------------------------------------
    ratio        <- 2.4/(rowSums(smr_matrix_n[,as.character(20:50)])/length(20:50))
    ratio_matrix <- t(matrix(ratio, nrow=100, ncol=length(ratio), byrow=TRUE))
    ratio_matrix[ratio_matrix>=1] <- 1
    smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1

    smr_matrix_n_new <- smr_matrix_n
    smr_matrix_n_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]
    # smr_matrix_m_new <- smr_matrix_m

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
    # smr_matrix_adjusted[smr_matrix_adjusted<1] <- smr_matrix_n[smr_matrix_adjusted<1]
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1

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

    ,`Ghosts lever2023`  = sum(`Ghosts lever2023`)
    ,`Ghosts (delta basic care) lever2023`  = sum(`Ghosts (delta basic care) lever2023`)
    ,`Ghosts (delta best care) lever2023`   = sum(`Ghosts (delta best care) lever2023`)
    ,`Ghosts (delta cure) lever2023`   = sum(`Ghosts (delta cure) lever2023`)
    ,`Ghosts (early death) lever2023`  = sum(`Ghosts (early death) lever2023`)
    ,`Ghosts (onset death) lever2023`  = sum(`Ghosts (onset death) lever2023`)

    ,`Incidence (1 base)`       = sum(`Incidence (1 base)`)
    ,`true_incidence_age_10`       = sum(`true_incidence_age_10`)

    ,`Incidence (2 diagnosis)`  = sum(`Incidence (2 diagnosis)`)
    ,`Prevalence`  = sum(`Prevalence`)
    ,`1 in x families`  = sum(`1 in x families` *`Ann. background population` )/sum(`Ann. background population`)

  ),by=key_list]


  # country_data_age10 <- setDT(country_data_merge)[Age==10,list(
    country_data_age10 <- setDT(country_data_merge)[,list(

    `Life expectency (1 background)`      =  agg_life_expectancy_median(`Life expectency (1 background)`      , `prevalence plus missing (age 10)`  )
    ,`Life expectency (2 t1d base)`       =  agg_life_expectancy_median(`Life expectency (2 t1d base)`        , `prevalence plus missing (age 10)`  )
    ,`Life expectency (3 t1d diagnosis)`  =  agg_life_expectancy_median(`Life expectency (3 t1d diagnosis)`   , `prevalence plus missing (age 10)`  )
    ,`Life expectency (4 t1d basic care)` =  agg_life_expectancy_median(`Life expectency (4 t1d basic care)`  , `prevalence plus missing (age 10)`  )
    ,`Life expectency (5 t1d best care)`  =  agg_life_expectancy_median(`Life expectency (5 t1d best care)`   , `prevalence plus missing (age 10)`  )
    ,`Life expectency (6 t1d cure)`       =  agg_life_expectancy_median(`Life expectency (6 t1d cure)`        , `prevalence plus missing (age 10)`  )

    ,`Lifetime years lost (2 t1d base) (complication)`        = agg_life_expectancy_median(`Lifetime years lost (2 t1d base) (complication)`              , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (3 t1d diagnosis) (complication)`   = agg_life_expectancy_median(`Lifetime years lost (3 t1d diagnosis) (complication)`         , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (4 t1d basic care) (complication)`  = agg_life_expectancy_median(`Lifetime years lost (4 t1d basic care) (complication)`        , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (5 t1d best care) (complication)`   = agg_life_expectancy_median(`Lifetime years lost (5 t1d best care) (complication)`         , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (6 t1d cure) (complication)`        = agg_life_expectancy_median(`Lifetime years lost (6 t1d cure) (complication)`              , `prevalence plus missing (age 10)`  )

    ,`Lifetime years lost (2 t1d base) (treatment)`        = agg_life_expectancy_median(`Lifetime years lost (2 t1d base) (treatment)`              , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (3 t1d diagnosis) (treatment)`   = agg_life_expectancy_median(`Lifetime years lost (3 t1d diagnosis) (treatment)`         , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (4 t1d basic care) (treatment)`  = agg_life_expectancy_median(`Lifetime years lost (4 t1d basic care) (treatment)`        , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (5 t1d best care) (treatment)`   = agg_life_expectancy_median(`Lifetime years lost (5 t1d best care) (treatment)`         , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (6 t1d cure) (treatment)`        = agg_life_expectancy_median(`Lifetime years lost (6 t1d cure) (treatment)`              , `prevalence plus missing (age 10)`  )
    ,`Life expectency (strip low)`   =   agg_life_expectancy_median(`Life expectency (strip low)`              , `prevalence plus missing (age 10)`  )
    ,`Life expectency (strip hig)`   =   agg_life_expectancy_median(`Life expectency (strip hig)`              , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (strip low)`   =   agg_life_expectancy_median(`Lifetime years lost (strip low)`              , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (strip hig)`   =   agg_life_expectancy_median(`Lifetime years lost (strip hig)`              , `prevalence plus missing (age 10)`  )
    ,`Life expectency (sensor low)`   =   agg_life_expectancy_median(`Life expectency (sensor low)`              , `prevalence plus missing (age 10)`  )
    ,`Life expectency (sensor hig)`   =   agg_life_expectancy_median(`Life expectency (sensor hig)`              , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (sensor low)`   =   agg_life_expectancy_median(`Lifetime years lost (sensor low)`              , `prevalence plus missing (age 10)`  )
    ,`Lifetime years lost (sensor hig)`   =   agg_life_expectancy_median(`Lifetime years lost (sensor hig)`              , `prevalence plus missing (age 10)`  )
    ,`% Odds living to`   =   agg_life_expectancy_median(`% Odds living to`              , `prevalence plus missing (age 10)`  )

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






HLY_breakdown_standard_of_care <- function()
{
  lever_year_range <- as.character(1860:2040)
  source("code_R_data_prep/DATA_CONSTANTS.R")
  avg_life_expe_at_10 <- Data_Run_query_return_df('SELECT * FROM main_0_4_15_paper_median WHERE "Age"=10  AND "Country" in (select world_bank_name from index_parameters.country) ')
  country_list        <- unique(avg_life_expe_at_10$Country)
  country_region <- Data_Run_query_return_df('SELECT loc_id,world_bank_name  FROM country_region')
  country        <- Data_Run_query_return_df('SELECT *  FROM index_parameters.country')
  result_df_all <- data.frame()



  for(i in 1:length(country_list))
  {  # i <- 84
    print(i)
    year_at <- 2022
    country_wb_name <- country_list[i]
    loc_id <- country_region$loc_id[country_region$world_bank_name==country_wb_name]

    data_long  <- Data_Run_query_return_df (paste0("SELECT * FROM index_input_0_4_15  WHERE loc_id = '",loc_id,"'" ) )
    # data_long   <- readRDS(paste0(data_dir,"/",loc_id,".Rds"))
    matrices_list <- data_long_2_matrices (data_long=data_long,data_long_default_run=data_long)

    # get different standard of care for india of 1a, 2b1, 3c.
    matrices_list_temp  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 31.1,31.1,lever_year_range)
    life_table_0A     <- Get_life_expectancy_t1d (qB, matrices_list_temp$smr_matrix_m,matrices_list_temp$smr_matrix_n, matrices_list_temp$qT1D_percent_n)
    life_table_0A$level_of_care <- "0A"
    matrices_list_temp  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 15.3,31.1,lever_year_range)
    life_table_1A     <- Get_life_expectancy_t1d (qB, matrices_list_temp$smr_matrix_m,matrices_list_temp$smr_matrix_n, matrices_list_temp$qT1D_percent_n)
    life_table_1A$level_of_care <- "1A"

    matrices_list_temp  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 3.7,4.4,lever_year_range)
    life_table_2B1     <- Get_life_expectancy_t1d (qB, matrices_list_temp$smr_matrix_m,matrices_list_temp$smr_matrix_n, matrices_list_temp$qT1D_percent_n)
    life_table_2B1$level_of_care <- "2B1"

    matrices_list_temp  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 2.2,2.6,lever_year_range)
    life_table_3C       <- Get_life_expectancy_t1d (qB, matrices_list_temp$smr_matrix_m,matrices_list_temp$smr_matrix_n, matrices_list_temp$qT1D_percent_n)
    life_table_3C$level_of_care <- "3C"

    # matrices_list_temp  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 1,1,lever_year_range)
    # matrices_list_temp$smr_matrix_n <- 1  # this is the real cure s

    life_table_Cure     <- Get_life_expectancy_t1d (qB, matrices_list_temp$smr_matrix_m,matrices_list_temp$smr_matrix_n, matrices_list_temp$qT1D_percent_n)

    result_df <- select(life_table_0A[year==year_at&age==10],level_of_care,ex,ex_complication ,ex_treatment)

    result_df <- rbind(result_df, select(life_table_1A[year==year_at&age==10],level_of_care,ex ,ex_complication ,ex_treatment) )
    result_df <- rbind(result_df, select(life_table_2B1[year==year_at&age==10],level_of_care,ex ,ex_complication ,ex_treatment) )
    result_df <- rbind(result_df, select(life_table_3C[year==year_at&age==10],level_of_care,ex,ex_complication ,ex_treatment ) )
    # result_df <- rbind(result_df, select(life_table_Cure[year==year_at&age==10],ex,ex_complication ,ex_treatment ) )
    # result_df$avg_life_expe_at_10 <- 73.46
    result_df$avg_life_expe_at_10    <- avg_life_expe_at_10$`Life expectency (1 background)`[avg_life_expe_at_10$Year==year_at&avg_life_expe_at_10$Country==country_wb_name&avg_life_expe_at_10$age_bracket=="00_99"] + 10
    result_df$prevalence_and_missing <- (avg_life_expe_at_10$Prevalence + avg_life_expe_at_10$Ghosts)[avg_life_expe_at_10$Year==year_at&avg_life_expe_at_10$Country==country_wb_name&avg_life_expe_at_10$age_bracket=="00_99"]

    result_df$value_1 <- result_df$avg_life_expe_at_10 - (result_df$ex +10 )
    result_df$value_2 <- result_df$ex_treatment
    result_df$value_3 <- result_df$ex_complication
    result_df$value_4 <- result_df$ex - result_df$ex_treatment - result_df$ex_complication + 10

    result_df <- cbind(Country=country_wb_name,Year=year_at, result_df)

    result_df_all <- rbind(result_df_all,result_df)
  }

  result_df_all_full <- result_df_all %>% left_join(select(country,Country= world_bank_name,wd_region ))




  result_df_all_full_region <- setDT(result_df_all_full)[,list(
     value_1= sum(value_1 * prevalence_and_missing)/sum(prevalence_and_missing)
    ,value_2= sum(value_2 * prevalence_and_missing)/sum(prevalence_and_missing)
    ,value_3= sum(value_3 * prevalence_and_missing)/sum(prevalence_and_missing)
    ,value_4= sum(value_4 * prevalence_and_missing)/sum(prevalence_and_missing)
  ),by=c("wd_region","Year","level_of_care")]



  # write.csv(result_df_all_full_region[result_df_all_full_region$wd_region=="Middle East & North Africa",],"temp/standard_of_care_MENA.csv")



  result_df_all_full_global <- setDT(result_df_all_full)[,list(
     value_1= sum(value_1 * prevalence_and_missing)/sum(prevalence_and_missing)
    ,value_2= sum(value_2 * prevalence_and_missing)/sum(prevalence_and_missing)
    ,value_3= sum(value_3 * prevalence_and_missing)/sum(prevalence_and_missing)
    ,value_4= sum(value_4 * prevalence_and_missing)/sum(prevalence_and_missing)
  ),by=c("Year","level_of_care")]

  # write.csv(result_df_all_full_global,"temp/standard_of_care_GLOBAL.csv")



}

Get_qT1D_from_smr_matrix <- function(qB, smr_matrix)
{
  diab_odds  <- qB / (1 - qB) * smr_matrix
  qT1D       <- diab_odds / (1 + diab_odds)
  qT1D
}
