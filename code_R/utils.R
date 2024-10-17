# Utility functions
#
# Confidential
# Copyright (c) JDRF 2020, All rights reserved
#

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Convert a constant to univariate function, if needed
#'
#' This variant broadcasts the numeric to width \code{MAX_AGE}.
#'
#' @param f function or numeric constant
#' @return a function, possibly converting f to a function
#' @export
make_age_function <- function(f) {
  if(is.function(f)) {
    f
  } else {
    function(year, draw=1) {
      res <- rep(f, MAX_AGE)
      if (length(year) > 1) {
        matrix(rep(res, length(year)), ncol=100, byrow=FALSE,
               dimnames = list(year, AGES))
      } else {
        res
      }
    }
  }
}


#' Returns an array based on the function provided
#'
#' Internal function that returns an array based on the closure provided for
#' each year.
#'
#' @param f function. Can be incidence, mortality, death or population.
#' @param years numeric vector containing years to calculate statistics for.
#' @param ... additional parameters to pass through to f
#'
#' @return array
#' @export
matrix_from_function <- function(f, years, ...) {
  a <- t(sapply(years, f, simplify='matrix', ...))
  dimnames(a) <- list(years, AGES)
  a
}

inpute_to_1900_2040 <- function(states_population_proportion, year_start=1899, year_end=2040)
{
  year_end_input <- year_end
  year_end <- 2040
  # states_population_proportion <- states_profile
  # impute to 1900-2040   # states_population_proportion <- states_lifetable_ratio
  year_columns <- colnames(states_population_proportion)[colnames(states_population_proportion) %in% as.character(year_start:year_end)]
  non_year_columns <- colnames(states_population_proportion)[!colnames(states_population_proportion) %in% as.character(year_start:year_end)]


  year_range <- range(year_columns)
  previous_years <- data.frame(rep(states_population_proportion[,year_range[1],drop=FALSE], as.numeric(year_range[1]) - year_start ))
  colnames(previous_years) <- as.character(year_start:(as.numeric(year_range[1])-1))

  future_years <- data.frame(rep(states_population_proportion[,year_range[2],drop=FALSE], year_end- as.numeric(year_range[2]) ))

  if(ncol(future_years))
  {
    colnames(future_years) <- as.character((as.numeric(year_range[2])+1):year_end )
    population_mapping  <- cbind(states_population_proportion[,non_year_columns,drop=FALSE],previous_years,states_population_proportion[,year_columns,drop=FALSE], future_years)
  }else
  {
    population_mapping  <- cbind(states_population_proportion[,non_year_columns,drop=FALSE],previous_years,states_population_proportion[,year_columns,drop=FALSE])
  }


  population_mapping <- population_mapping[,c(non_year_columns,as.character((year_start+1):year_end_input) )]
  population_mapping

}

#' @param data_long  with column (loc_id, year,age,value).  imputate to year :1900 - 2040, keep constant outside,  age : 0-99
#' @export
imputation_year_age <- function(data_long)
{
  # data_long <- model_smr_minimal_care
  # data_long <- model_smr_minimal_care
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))


  # data_long <- od_rates [od_rates$loc_id==462,]

  data_wide <-  spread(data_long, year, "value")
  data_wide_prefix    <-data_wide[,c(1,2),drop=F]
  data_wide           <-data_wide[,c(-1,-2)]


  # Impute within the years
  new_cname_expand <- range(colnames(data_wide))[1] : range(colnames(data_wide))[2]
  data_wide <- data.frame(t(apply(data_wide, 1, function(x_row) spline(x=colnames(data_wide), y=log(x_row), method="natural",xout=new_cname_expand  )$y)))
  data_wide <- exp(data_wide)
  colnames(data_wide) <- new_cname_expand
  data_wide[is.nan(data_wide) ] <- 0
  # data.frame(x=colnames(data_wide), y = c(t(data_wide[1,]) )) %>% e_chart(x) %>% e_scatter(y)


  # Keep values before and after fit cuve constant --------------
  data_wide <- inpute_to_1900_2040(data_wide, year_start=1899, year_end=2040)


  # data.frame(x=colnames(data_wide), y = c(t(data_wide[1,]) )) %>% e_chart(x) %>% e_scatter(y)

  data_wide <- cbind(data_wide_prefix,data_wide)

  data_long <-  gather(data_wide, key = year, value ="value" , -loc_id,-age)

  data_wide <-  spread(data_long, age, value)

  data_wide_prefix    <-data_wide[,c(1,2),drop=F]
  data_wide           <-data_wide[,c(-1,-2),drop=FALSE]

  new_cname_expand <- seq(0, 99, 1)

  data_wide <- data.frame(t(apply(data_wide, 1, function(x) spline(1:ncol(data_wide), x, length(new_cname_expand))$y)))

  colnames(data_wide) <- new_cname_expand

  data_wide <- cbind(data_wide_prefix,data_wide)


  data_long  <-  gather(data_wide, key = age, value ="value" , -loc_id,-year)

  data_long$year <- as.numeric(data_long$year)
  data_long$age <- as.numeric(data_long$age)
  data_long
}


# calculate ex  life expectency at age x ---------------
calculate_ex <- function(life_table,add_year=30)
{
  # life_table <- dplyr::select(life_table_background, loc_id, year,age,qx =background_mortality_rate)
  # life_table <- states_lifetable
  # life_table <- life_table_t_all
  # life_table <- data.frame(life_table_background_)
  # add another 30 years to be close to real situation.
  # add_year <- 30

  if(! is.na(add_year) )
  {
    life_table_add_year_x <- life_table[life_table$age>=100-add_year,]
    life_table_add_year_x$age<- life_table_add_year_x$age + add_year
    for(i in range(life_table_add_year_x$age )[1] :range(life_table_add_year_x$age )[2]   )
    {  # i <- 100
      life_table_add_year_x[life_table_add_year_x$age==i,]$qx  <-  sapply(life_table[life_table$age==99,]$qx +  (life_table[life_table$age==99,]$qx - life_table[life_table$age==98,]$qx)*(i-99)  ,min   ,1)
    }
    life_table <- rbind(life_table,life_table_add_year_x )
    # print(colnames(life_table))
    if("loc_id" %in%  colnames(life_table))
    {
      life_table <- life_table %>% dplyr::arrange(loc_id,year, age)
    }else
    {
      life_table <- life_table %>% dplyr::arrange(year, age)
    }
  }

  life_table$lx <- NA
  life_table$dx <- NA
  life_table$lx[life_table$age==0] <- 100000
  life_table$dx[life_table$age==0] <- 100000 * life_table$qx[life_table$age==0]
  for(i in 1:max(life_table$age))
  {# i <- 1
    life_table$lx[life_table$age==i] <- life_table$lx[life_table$age==(i-1)] -  life_table$dx[life_table$age==(i-1)]
    life_table$dx[life_table$age==i] <- life_table$lx[life_table$age==i] * life_table$qx[life_table$age==i]
  }

  life_table$Lx <- c(life_table$lx[-1],0)
  life_table$Lx[life_table$age==max(life_table$age)] <- 0
  life_table$Lx <- life_table$Lx + 0.5 * life_table$dx

  life_table$Tx <- 0
  life_table$Tx[life_table$age==max(life_table$age)] <- life_table$Lx[life_table$age==max(life_table$age)]

  for(i in (max(life_table$age)-1):0)
  { # i <- 1
    life_table$Tx[life_table$age==i] <- life_table$Tx[life_table$age==(i+1) ]+ life_table$Lx[life_table$age==i ]
  }

  life_table$ex <- ifelse(life_table$Tx==0 , 0 ,life_table$Tx / life_table$lx )
  life_table <- life_table[life_table$age <=99,]
  return(life_table)
}

# calculate ex  life expectency at age x ---------------
calculate_ex_matrix <- function(qB_full,add_year=30)
{
  # life_table <- dplyr::select(life_table_background, loc_id, year,age,qx =background_mortality_rate)
  # life_table <- data.matrix(life_table_background_)
  life_table_matrix <- qB_full
  dim_size <- dim(life_table_matrix)     ;dim_size[3] <- 130
  dim_name <- dimnames(life_table_matrix);dim_name[[3]] <- as.character(0:129)
  # Model compartments, all annual cohorts as a proportion of 1.

  qx <-  array(NA, dim=list(dim_size[1], dim_size[2], dim_size[3]), dimnames=dim_name)# S - susceptible (non-T1D)
  lx <-  array(NA, dim=list(dim_size[1], dim_size[2], dim_size[3]), dimnames=dim_name)# S - susceptible (non-T1D)
  dx <-  array(NA, dim=list(dim_size[1], dim_size[2], dim_size[3]), dimnames=dim_name)# S - susceptible (non-T1D)
  Lx <-  array(NA, dim=list(dim_size[1], dim_size[2], dim_size[3]), dimnames=dim_name)# S - susceptible (non-T1D)
  Tx <-  array(NA, dim=list(dim_size[1], dim_size[2], dim_size[3]), dimnames=dim_name)# S - susceptible (non-T1D)
  ex <-  array(NA, dim=list(dim_size[1], dim_size[2], dim_size[3]), dimnames=dim_name)# S - susceptible (non-T1D)

  qx[,,1:100] <- life_table_matrix[,,]

  for(i in 100:129)
  { # i <- 100
    qx[,,as.character(i)]  <- qx[,,"99",drop=FALSE] + (qx[,,"99",drop=FALSE] - qx[,,"98",drop=FALSE])*(i-99)
  }
  qx[qx>=1] <- 1

  lx[,,"0"] <- 100000
  dx[,,"0"] <- 100000 * qx[,,"0",drop=FALSE]

  for(i in 1:129)
  {
    # i <- 1   ;
    i_minus <- as.character(i-1)
    i       <-  as.character(i)
    lx[,,i] <- lx[,,i_minus] -  dx[,,i_minus]
    dx[,,i] <- lx[,,i] * qx[,,i]
  }

  Lx[,,as.character(0:128) ] <- lx[,,-1]
  Lx[,,as.character(129) ] <- 0
  Lx <- Lx + 0.5 * dx

  Tx[] <- 0
  Tx[,,"129"] <- Lx[,,"129",drop=FALSE]

  for(i in 128:0)
  { # i <- 128
    Tx[,, as.character(i)] <- Tx[,, as.character(i+1),drop=FALSE]+ Lx[,, as.character(i),drop=FALSE]
  }
  ex <- Tx / lx
  ex[is.nan(ex)] <- 0

  ex <- ex[,,1:100,drop=FALSE]
  lx <- lx[,,1:100,drop=FALSE]
  Lx <- Lx[,,1:100,drop=FALSE]

  return( list(ex=ex,Lx=Lx,lx=lx) )

}


# calculate ex  life expectency at age x ---------------
calculate_ex_lifetime_years_lost_matrix <- function(qB,smr_matrix,Enable_DKA_Hypo=FALSE,years_range=as.character(1960:2040))
{  # in the case of smr is 1 :  life_table$smr <- 2
  # life_table <- t1d_mortality_long
  qB         <- qB[,years_range,,drop=FALSE]
  smr_matrix <- smr_matrix[,years_range,,drop=FALSE]
  # life_table1 <-  calculate_ex(life_table,add_year = 1)
  diab_odds    <- qB / (1 - qB) * smr_matrix
  qT1D         <- diab_odds / (1 + diab_odds)

  ex_list  <-  calculate_ex_matrix (qT1D,add_year=30)
  ex       <-  ex_list$ex
  lx       <-  ex_list$lx
  Lx       <-  ex_list$Lx

  #Add complications :  Life time years lost factor ---------------------------------------------------
  hba1c   <- (log(smr_matrix) +  1.5274 )/ 0.3545
  # View(life_table %>% mutate_if(is.numeric, round, digits=2))
  weib_survival_burden_matrix <- function(hba1c) {
    # load parameters
    complication_parameters <- get_complication_parameters()
    weibull                 <- complication_parameters$weibull
    disease_weights         <- get_disease_weights()

    intercept <- complication_parameters$weibull$intercept
    slope <- complication_parameters$weibull$slope
    scale <- complication_parameters$weibull$scale

    # weibull=complication_parameters$weibull  ;  hba1c <- c(10,10,10,10,10,10)
    props <- array(0,
                   dim=c(nrow(weibull), dim(hba1c)[[1]], dim(hba1c)[[2]], dim(hba1c)[[3]]-10 )
                   ,dimnames = list(comp=weibull$abbrev,dimnames(hba1c)[[1]],dimnames(hba1c)[[2]],dimnames(hba1c)[[3]][11:100]  )
                   )

    for(i in 1:nrow(weibull))
    { # i <- 1

      hba1c_age_cut_10 <-  hba1c[,,11:100,drop=FALSE]

      intercept <- weibull$intercept[i]
      slope     <- weibull$slope[i]
      scale     <- weibull$scale[i]

      T <- dim(hba1c_age_cut_10)[[3]]-1 # total time length (years) incl. T=0
      hs_odd <- 0.5 * (hba1c_age_cut_10[,,-1,drop=FALSE] + hba1c_age_cut_10[,,-(T+1),drop=FALSE]) # interpolate at half-year points
      ts <- pmin(0:T, 30)  # constant hazard after 30 years
      ts <- array(ts,dim=c(length(ts), dim(hba1c_age_cut_10)[[2]], dim(hba1c_age_cut_10)[[1]] ))
      ts <- aperm(ts,c(3,2,1))
      # Simpson's rule

      l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c_age_cut_10)/scale) / scale
      l_odd  <- (ts[,,-1,drop=FALSE] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale

      # l_even <- (expm1((1 - scale)/scale * log1p(ts - 1)) + 1) * exp(-(intercept + slope * hba1c_age_cut_10)/scale) / scale
      # l_odd  <- (expm1(((1 - scale)/scale) * log1p((ts[,,-1] - 0.5) - 1)) + 1)  * exp(-(intercept + slope * hs_odd)/scale) / scale

      l_even_cumsum <- abind( array(0, dim=c(dim(hba1c_age_cut_10)[[1]], dim(hba1c_age_cut_10)[[2]], 2 )), l_even[,,-c(1, T+1),drop=FALSE],along=3)
      l_odd_cumsum  <- abind( array(0, dim=c(dim(hba1c_age_cut_10)[[1]], dim(hba1c_age_cut_10)[[2]], 1 )),  l_odd[,,-(T+1),drop=FALSE]    ,along=3)

      cumsum_matrix <-  2 * l_even_cumsum + 4 * l_odd_cumsum
      for(j in dim(hba1c_age_cut_10)[[3]]:2)
      { # j <- 100
        cumsum_matrix[,,j] <-    rowSums(cumsum_matrix[,,1:j,drop=FALSE], dims = 2)
      }
      Lambda <- 1/6 * (l_even + cumsum_matrix)

      props[weibull$abbrev[i],,,] <- 1- exp(-Lambda)
    }
    # This is Graham's and Gabriel's 'journeys' adjustments. Reduces ON and PR by survivor shares of RF and BL, respectively,
    # effectively making this a multistate model (rather than a pure time-to-event model)
    props['ON',,,] <- props['ON',,,] * (1 - props['RF',,,])
    props['PR',,,] <- props['PR',,,] * (1 - props['BL',,,])

    # Calculate disability weights to apply to prevalence. Since 100% of prevalent
    # cases have T1D, we apply the T1D disability weight directly. For all other
    # complications, we scale by the modeled complication probability.
    # Terminology note: disease weight == disability weight
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
    # disability_wts <- (
    #   1 - (1 - t1d) * (1 - dsp) * (1 - pr) * (1 - on - rf) * (1 - uoa) *
    #     (1 - hom) * (1 - nfMI) * (1 - nfCBVD))
    disability_wts <- (
      1 -  (1 - dsp) * (1 - pr) * (1 - on - rf) * (1 - uoa) *
        (1 - hom) * (1 - nfMI) * (1 - nfCBVD))
    props <- (
      1 -  (1 - props["DSP",,,]) * (1 - props["PR",,,]) * (1 - props["ON",,,] - props["RF",,,]) * (1 - props["UoA",,,]) *
        (1 - props["HoM",,,]) * (1 - props["nfMI",,,]) * (1 - props["nfCBVD",,,]))

    list(disability_wts=disability_wts,props=props,disease_weights=disease_weights )
    # disability_wts
  }

  risk      <-   hba1c ; risk[] <- 0
  risk_list <-   weib_survival_burden_matrix (hba1c)
  risk[,,as.character(10:99)] <-   risk_list$disability_wts
  disease_weights             <-   risk_list$disease_weights

  Lx_original <- Lx

  # T1D complication lifetime years -----------------------------------------------------------------------------------------------
  Lx <- Lx_original * risk
  Tx <-   hba1c ; Tx[] <- 0
  Tx[,,"99"] <- Lx[,,"99"]

  for(i in 98:0)
  { # i <- 128
    Tx[,, as.character(i)] <- Tx[,, as.character(i+1)]+ Lx[,, as.character(i)]
  }
  ex_complication <- Tx / lx
  ex_complication[is.nan(ex_complication)] <- 0

  # T1D treatment lifetime years ---------------------------------------------------------------------------------------------------

  Lx <- Lx_original * disease_weights$T1D
  Tx <-   hba1c ; Tx[] <- 0
  Tx[,,"99"] <- Lx[,,"99"]

  for(i in 98:0)
  { # i <- 128
    Tx[,, as.character(i)] <- Tx[,, as.character(i+1)]+ Lx[,, as.character(i)]
  }
  ex_treatment <- Tx / lx
  ex_treatment[is.nan(ex_treatment)] <- 0

  return(list(ex=ex, ex_complication=ex_complication,ex_treatment=ex_treatment ))
}


# calculate ex  life expectency at age x ---------------
calculate_ex_lifetime_years_lost <- function(life_table,Enable_DKA_Hypo=FALSE)
{  # in the case of smr is 1 :  life_table$smr <- 2
  # life_table <- t1d_mortality_long

  # life_table1 <-  calculate_ex(life_table,add_year = 1)
  life_table <-  calculate_ex(life_table)

  #Add complications :  Life time years lost factor ---------------------------------------------------
  life_table$hba1c   <- (log(life_table$smr) +  1.5274 )/ 0.3545
  # life_table$hba1c  <- 4
  # View(life_table %>% mutate_if(is.numeric, round, digits=2))

  complication_parameters <- get_complication_parameters()
  weibull                 <- complication_parameters$weibull
  disease_weights         <- get_disease_weights()

  inputs <- list()
  intercept <- complication_parameters$weibull$intercept
  slope <- complication_parameters$weibull$slope
  scale <- complication_parameters$weibull$scale

  weib_survival_burden <- function(hba1c, weibull,disease_weights) {
    # hba1c <- hba1c[1,1,]

    # weibull=complication_parameters$weibull  ;  hba1c <- c(10,10,10,10,10,10)
    props <- array(0,
                   dim=c(nrow(weibull), length(hba1c)),
                   dimnames = list(comp=weibull$abbrev, age=1:length(hba1c)))

    for(i in 1:nrow(weibull))
    {
      intercept <- weibull$intercept[i]
      slope     <- weibull$slope[i]
      scale     <- weibull$scale[i]
      T <- length(hba1c)-1 # total time length (years) incl. T=0
      hs_odd <- 0.5 * (hba1c[-1] + hba1c[-(T+1)]) # interpolate at half-year points
      ts <- pmin(0:T, 30)  # constant hazard after 30 years

      # Simpson's rule
      l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c)/scale) / scale
      l_odd <- (ts[-1] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale
      Lambda <- 1/6 * (l_even + cumsum(2 * c(0, 0, l_even[-c(1, T+1)]) + 4 * c(0, l_odd[-(T+1)])))

      exp(-Lambda) # survival function
      # print( exp(-Lambda) )
      props[weibull$abbrev[i],] <- 1- exp(-Lambda)

    }
    # This is Graham's and Gabriel's 'journeys' adjustments. Reduces ON and PR by survivor shares of RF and BL, respectively,
    # effectively making this a multistate model (rather than a pure time-to-event model)
    props['ON',] <- props['ON',] * (1 - props['RF',])
    props['PR',] <- props['PR',] * (1 - props['BL',])

    # Calculate disability weights to apply to prevalence. Since 100% of prevalent
    # cases have T1D, we apply the T1D disability weight directly. For all other
    # complications, we scale by the modeled complication probability.
    # Terminology note: disease weight == disability weight
    t1d <- disease_weights$T1D
    dsp <- props["DSP",] * disease_weights$DSP
    pr <- props["PR",] * disease_weights$PR + props["BL",] * disease_weights$BL
    on <- props["ON",] * disease_weights$ON
    #rf <- (props["RF_transplant",,,] * disease_weights$Transplant
    #  + props["RF_dialysis",,,] * disease_weights$Dialysis)
    rf <- props["RF",] * disease_weights$Dialysis
    uoa <- props["UoA",] * disease_weights$UoA
    hom <- props["HoM",] * disease_weights$HoM
    nfMI <- props["nfMI",] * disease_weights$nfMI
    nfCBVD <- props["nfCBVD",] * disease_weights$nfCBVD

    # disability_wts <- (
    #   1 - (1 - t1d) * (1 - dsp) * (1 - pr) * (1 - on - rf) * (1 - uoa) *
    #     (1 - hom) * (1 - nfMI) * (1 - nfCBVD))
    disability_wts <- (
      1 -  (1 - dsp) * (1 - pr) * (1 - on - rf) * (1 - uoa) *
        (1 - hom) * (1 - nfMI) * (1 - nfCBVD))
    props <- (
      1 -  (1 - props["DSP",]) * (1 - props["PR",]) * (1 - props["ON",] - props["RF",]) * (1 - props["UoA",]) *
        (1 - props["HoM",]) * (1 - props["nfMI",]) * (1 - props["nfCBVD",]))

    list(disability_wts=disability_wts,props=props )
    # disability_wts

  }

  weib_survival_burden_matrix <- function(hba1c, weibull,disease_weights) {
    # hba1c <- hba1c[1,1,]
    hba1c <-  (log(smr_matrix_n) +  1.5274 )/ 0.3545

    hba1c <- hba1c[,101:181,11:100]

    # weibull=complication_parameters$weibull  ;  hba1c <- c(10,10,10,10,10,10)
    props <- array(0,
                   dim=c(nrow(weibull), length(hba1c)),
                   dimnames = list(comp=weibull$abbrev, age=1:length(hba1c)))

    for(i in 1:nrow(weibull))
    { # i <- 1
      intercept <- weibull$intercept[i]
      slope     <- weibull$slope[i]
      scale     <- weibull$scale[i]

      T <- dim(hba1c)[[3]]-1 # total time length (years) incl. T=0
      hs_odd <- 0.5 * (hba1c[,,-1] + hba1c[,,-(T+1)]) # interpolate at half-year points
      ts <- pmin(0:T, 30)  # constant hazard after 30 years
      ts <- array(ts,dim=c(length(ts), dim(hba1c)[[2]], dim(hba1c)[[1]] ))
      ts <- aperm(ts,c(3,2,1))
      # Simpson's rule

      # asdf <- exp(ts * log((1 - scale)/scale))
      # asdf <- (expm1((1 - scale)/scale * log1p(ts - 1)) + 1)  # Adjusting for b close to 1
      # asdf2 <- ts^((1 - scale)/scale)  # Adjusting for b close to 1

      l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c)/scale) / scale
      l_odd  <- (ts[,,-1] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale

      # l_even <- (expm1((1 - scale)/scale * log1p(ts - 1)) + 1) * exp(-(intercept + slope * hba1c)/scale) / scale
      # l_odd  <- (expm1(((1 - scale)/scale) * log1p((ts[,,-1] - 0.5) - 1)) + 1)  * exp(-(intercept + slope * hs_odd)/scale) / scale

      l_even_cumsum <- abind( array(0, dim=c(dim(hba1c)[[1]], dim(hba1c)[[2]], 2 )), l_even[,,-c(1, T+1)],along=3)
      l_odd_cumsum  <- abind( array(0, dim=c(dim(hba1c)[[1]], dim(hba1c)[[2]], 1 )),  l_odd[,,-(T+1)]    ,along=3)

      cumsum_matrix <-  2 * l_even_cumsum + 4 * l_odd_cumsum
      for(j in dim(hba1c)[[3]]:2)
      { # j <- 100
        cumsum_matrix[,,j] <-    rowSums(cumsum_matrix[,,1:j], dims = 2)
      }
      Lambda <- 1/6 * (l_even + cumsum_matrix)

      # exp(-Lambda) # survival function
      # print( exp(-Lambda) )
      # props[weibull$abbrev[i],] <- 1- exp(-Lambda)

    }

    # This is Graham's and Gabriel's 'journeys' adjustments. Reduces ON and PR by survivor shares of RF and BL, respectively,
    # effectively making this a multistate model (rather than a pure time-to-event model)
    props['ON',] <- props['ON',] * (1 - props['RF',])
    props['PR',] <- props['PR',] * (1 - props['BL',])

    # Calculate disability weights to apply to prevalence. Since 100% of prevalent
    # cases have T1D, we apply the T1D disability weight directly. For all other
    # complications, we scale by the modeled complication probability.
    # Terminology note: disease weight == disability weight
    t1d <- disease_weights$T1D
    dsp <- props["DSP",] * disease_weights$DSP
    pr <- props["PR",] * disease_weights$PR + props["BL",] * disease_weights$BL
    on <- props["ON",] * disease_weights$ON
    #rf <- (props["RF_transplant",,,] * disease_weights$Transplant
    #  + props["RF_dialysis",,,] * disease_weights$Dialysis)
    rf <- props["RF",] * disease_weights$Dialysis
    uoa <- props["UoA",] * disease_weights$UoA
    hom <- props["HoM",] * disease_weights$HoM
    nfMI <- props["nfMI",] * disease_weights$nfMI
    nfCBVD <- props["nfCBVD",] * disease_weights$nfCBVD

    # disability_wts <- (
    #   1 - (1 - t1d) * (1 - dsp) * (1 - pr) * (1 - on - rf) * (1 - uoa) *
    #     (1 - hom) * (1 - nfMI) * (1 - nfCBVD))
    disability_wts <- (
      1 -  (1 - dsp) * (1 - pr) * (1 - on - rf) * (1 - uoa) *
        (1 - hom) * (1 - nfMI) * (1 - nfCBVD))
    props <- (
      1 -  (1 - props["DSP",]) * (1 - props["PR",]) * (1 - props["ON",] - props["RF",]) * (1 - props["UoA",]) *
        (1 - props["HoM",]) * (1 - props["nfMI",]) * (1 - props["nfCBVD",]))

    list(disability_wts=disability_wts,props=props )
    # disability_wts

  }



  setDT(life_table)[age>=10,risk:=weib_survival_burden(hba1c=hba1c[age>=10], weibull=weibull,disease_weights)$disability_wts,by=c("year")]
  # setDT(life_table)[age>=10,prob:=weib_survival_burden(hba1c=hba1c[age>=10], weibull=weibull,disease_weights)$props,by=c("year")]
  life_table$Lx_original <- life_table$Lx

  life_table$risk[is.na(life_table$risk)] <- 0

  # T1D complication lifetime years -----------------------------------------------------------------------------------------------
  life_table$Lx <- life_table$Lx_original * life_table$risk
  life_table$Tx <- 0
  life_table$Tx[life_table$age==max(life_table$age)] <- life_table$Lx[life_table$age==max(life_table$age)]
  for(i in (max(life_table$age)-1):0)
  { # i <- 1
    life_table$Tx[life_table$age==i] <- life_table$Tx[life_table$age==(i+1) ]+ life_table$Lx[life_table$age==i ]
  }
  life_table$ex_complication <- ifelse(life_table$Tx==0 , 0 ,life_table$Tx / life_table$lx )

  # T1D treatment lifetime years ---------------------------------------------------------------------------------------------------
  life_table$Lx <-  life_table$Lx_original * disease_weights$T1D
  life_table$Tx <- 0
  life_table$Tx[life_table$age==max(life_table$age)] <- life_table$Lx[life_table$age==max(life_table$age)]
  for(i in (max(life_table$age)-1):0)
  { # i <- 1
    life_table$Tx[life_table$age==i] <- life_table$Tx[life_table$age==(i+1) ]+ life_table$Lx[life_table$age==i ]
  }
  life_table$ex_treatment <- ifelse(life_table$Tx==0 , 0 ,life_table$Tx / life_table$lx )


  if(Enable_DKA_Hypo)
  {


    # Estimated DKA counts
    constant                  <- complication_parameters$constant
    cr <- function(constant,hba1c_vector,dka_hypo)
    {
      const_risk <- constant[constant$complication== dka_hypo,]
      lt8 <- 1e-2 * const_risk$hba1c_lt_8
      gt8lt9 <- 1e-2 * const_risk$hba1c_8_to_9
      gt9 <- 1e-2 * const_risk$hba1c_gteq_9
      ifelse(hba1c_vector < 8, lt8, ifelse(hba1c_vector < 9, gt8lt9, gt9))
    }

    # Estimated DKA counts
    setDT(life_table)[age>=10,prob:=cr(constant,hba1c_vector=hba1c[age>=10],dka_hypo="DKA"),by=c("year")]
    life_table$prob[is.na(life_table$prob)] <- 0

    # T1D complication lifetime years -----------------------------------------------------------------------------------------------
    life_table$Lx <- life_table$Lx_original * life_table$prob
    life_table$Tx <- 0
    life_table$Tx[life_table$age==max(life_table$age)] <- life_table$Lx[life_table$age==max(life_table$age)]
    for(i in (max(life_table$age)-1):0)
    { # i <- 1
      life_table$Tx[life_table$age==i] <- life_table$Tx[life_table$age==(i+1) ]+ life_table$Lx[life_table$age==i ]
    }
    life_table$events_DKA <- ifelse(life_table$Tx==0 , 0 ,life_table$Tx / life_table$lx )

    # Estimated Hypo counts
    setDT(life_table)[age>=10,prob:=cr(constant,hba1c_vector=hba1c[age>=10],dka_hypo="Hypo"),by=c("year")]
    life_table$prob[is.na(life_table$prob)] <- 0

    # T1D complication lifetime years -----------------------------------------------------------------------------------------------
    life_table$Lx <- life_table$Lx_original * life_table$prob
    life_table$Tx <- 0
    life_table$Tx[life_table$age==max(life_table$age)] <- life_table$Lx[life_table$age==max(life_table$age)]
    for(i in (max(life_table$age)-1):0)
    { # i <- 1
      life_table$Tx[life_table$age==i] <- life_table$Tx[life_table$age==(i+1) ]+ life_table$Lx[life_table$age==i ]
    }
    life_table$events_Hypo <- ifelse(life_table$Tx==0 , 0 ,life_table$Tx / life_table$lx )

  }


  return(life_table)
}




random_gen <- function(mean,smr_ci_lower, smr_ci_upper)
{
  smr <- NULL
  smr_ci_upper <- sapply(smr_ci_upper, max,0)  +1
  smr_ci_lower <- sapply(smr_ci_lower, max,0)+1
  mean         <- sapply(mean, max,0)+1
  # # normal distribution
  #  sd <-                    (smr_ci_upper - smr_ci_lower)/(1.96*2)
  #  smr1                      <- rnorm(length(mean), mean=mean, sd= sd )/1000000
  #  smr1[smr1<0] <- 0
  #  smr <- smr1
  #  # log normal distribution
  sd <-                    (log(smr_ci_upper) - log(smr_ci_lower) )/(1.96*2)
  smr3                      <- round(exp(rnorm(length(mean), mean=log(mean), sd= sd ))-1,4)
  smr <- smr3
  # # gamma
  # mean[mean==0] <-0.0001
  # smr2   <- rgamma(length(mean),shape=mean^2/sd^2,scale=sd^2/mean)/1000000
  # smr2[is.nan(smr2)] <- 0
  # # hist(smr,breaks = 200)
  # smr <- smr2
  return(smr)
}

assert <- function(condition,message)
{
  try(if(!condition) stop(message))
}


get_income_class_gni_year <- function(df_year_gni)
  # dfget_income_class_gni_year <- function(df_year_gni)
{
  # https://datahelpdesk.worldbank.org/knowledgebase/articles/906519#High_income
  # df_year_gni<- country_indicator_input_remoteness
  df_year_gni$income_class <- NULL
  df_year_gni$income_class <- ""
  country_incomeclass <- readxl::read_xlsx("data_internal/OGHIST (FY2024).xlsx",sheet ="Country Analytical History",skip = 5)
  colnames(country_incomeclass)[1:2] <- c("iso3c","Country")
  country_incomeclass <- country_incomeclass[1:4,c(-1,-2)]


  # country_incomeclass$`2022` <- c("<= 1,085","1,086 - 4,255","4,256 - 13,205","> 13,205"    )
  # country_incomeclass_pre <- country_incomeclass[,1:(1987-1960)]
  # colnames(country_incomeclass_pre) <- as.character(1960:1986)
  # country_incomeclass_pre[,] <- country_incomeclass_pre[,1]
  # country_incomeclass <- cbind(country_incomeclass_pre,country_incomeclass)

  country_incomeclass[1,] <- t(trimws(gsub("<=","",country_incomeclass[1,] )))
  country_incomeclass[4,] <- t(trimws(gsub(">","",country_incomeclass[4,] )))
  # country_incomeclass <- country_incomeclass[country_incomeclass$iso3c%in% country_indicator$iso3c ,]

  for(i in 1: nrow(df_year_gni) )
  { # i <- 1
    df_year_gni_ <- df_year_gni[i,]
    if(df_year_gni_$year>=1987)
    {
      if(df_year_gni_$gni_pc <= as.numeric( gsub(",","", country_incomeclass[1,as.character(df_year_gni_$year)] ) ) )
      {
        df_year_gni_$income_class <- "LIC"
      } else if (df_year_gni_$gni_pc >  as.numeric( gsub(",","", country_incomeclass[4,as.character(df_year_gni_$year)] ) ) )
      {
        df_year_gni_$income_class <- "HIC"
      }else if( df_year_gni_$gni_pc < as.numeric( gsub(",","", strsplit(  unlist(country_incomeclass[3,as.character(df_year_gni_$year)]),"-")[[1]][1] ) )  )
      {
        df_year_gni_$income_class <- "LMIC"
      }else
      {
        df_year_gni_$income_class <- "UMIC"
      }
      df_year_gni[i,]$income_class <- df_year_gni_$income_class
    }
  }
  for(i in 1: nrow(df_year_gni) )
  { # i <- 1
    if(df_year_gni$year[i]<1987)
    {
      df_year_gni$income_class[i] <- df_year_gni$income_class[ df_year_gni$year== 1987 & df_year_gni$country== df_year_gni$country[i] ]
    }
  }
  return(df_year_gni$income_class)
}





world_bank_name_convert_site <- function(world_bank_name)
{  # world_bank_name <- country_parquet_list$country_name
  code_site_country_name <- read.csv("data_internal/code_site_country_name.csv")

  for(i in 1:length(world_bank_name))
  { # i <- 1
    if(world_bank_name[i] %in% code_site_country_name$code)
    {
      # print(country_parquet_list$country_name_display[i] )
      # print(code_site_country_name$site[code_site_country_name$code == country_parquet_list$country_name_display[i]] )
      world_bank_name[i] <- code_site_country_name$site[code_site_country_name$code ==world_bank_name[i]   ]
    }
  }
  world_bank_name
}


