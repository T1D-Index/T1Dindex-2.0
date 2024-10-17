
refresh_country_files_tensor <- function(data_long=NULL
                                         ,sim_enable= FALSE
                                         ,sim_parameters=list(sim_start_year=2024,sim_min_diag_rates=0,sim_min_non_minimal_care_perc=0,sim_min_non_minimal_care_level="default_care")
                                         ,switches      =list(run_site_data=TRUE ,run_delay_onset=TRUE,run_lever2023=TRUE,run_stages_1_2=TRUE)
                                         )
{ # sim_enable= FALSE
  # data_long                   <- setDF(arrow::read_feather(paste0(data_dir,"/",partition_id)))
  # sim_parameters=list(sim_start_year=2023,sim_min_diag_rates=0,sim_min_non_minimal_care_perc=70,sim_min_non_minimal_care_level="default_care") ; switches=list(run_site_data=TRUE ,run_delay_onset=TRUE,run_lever2023=TRUE,run_stages_1_2=TRUE)
    prev_merge <- data.frame()
  # tryCatch({
    MAX_AGE <- 100
    AGES    <- 0: 99
    years   <- (1960 - 100):2040
    sim_lever_year_range <-  as.character(sim_parameters$sim_start_year:max(years))

    # library(abind)
    loc_ids        <- unique(data_long$loc_id)

    #  pre fill from 1860 to 1899-------------------------------------
    data_long_1860_1900 <- data_long[data_long$year==1900,] %>% slice(rep(1:n(), 1900-1860))
    data_long_1860_1900$year <- sort(rep(1860:1899,length(loc_ids)*100))
    data_long <- rbind(data_long,data_long_1860_1900 )
    # data_long <- data_long[ with(data_long, order( loc_id,year,age)), ]   # order properly for transforming to matrix
    #--------------
    data_long <- data_long[ with(data_long, order(age,year, loc_id)), ]   # order properly for transforming to matrix
    #--------------

    data_long_to_array <- function(data_long,run_projection=FALSE)
    { # run_projection <- TRUE
      df <- dplyr::select(data_long,loc_id,year,age,value )
      matrix_rate <- array(data = df$value,
                  dim= c(length(unique(df$loc_id)), length(unique(df$year)), length(unique(df$age))),
                  dimnames=list(unique(df$loc_id)        , unique(df$year)        , unique(df$age)))
      matrix_rate
      if(run_projection)
      {
        year_start <- 2021
        matrix_rate_default <- matrix_rate
        rate_avg_5 <- (    matrix_rate_default[,as.character(year_start),]/matrix_rate_default[,as.character(year_start-1),] +
                             matrix_rate_default[,as.character(year_start-1),]/matrix_rate_default[,as.character(year_start-2),] +
                             matrix_rate_default[,as.character(year_start-2),]/matrix_rate_default[,as.character(year_start-3),] +
                             matrix_rate_default[,as.character(year_start-3),]/matrix_rate_default[,as.character(year_start-4),] +
                             matrix_rate_default[,as.character(year_start-4),]/matrix_rate_default[,as.character(year_start-5),]
                           +matrix_rate_default[,as.character(year_start-5),]/matrix_rate_default[,as.character(year_start-6),]
                           +matrix_rate_default[,as.character(year_start-6),]/matrix_rate_default[,as.character(year_start-7),]
                           +matrix_rate_default[,as.character(year_start-7),]/matrix_rate_default[,as.character(year_start-8),]
                           +matrix_rate_default[,as.character(year_start-8),]/matrix_rate_default[,as.character(year_start-9),]
                           +matrix_rate_default[,as.character(year_start-9),]/matrix_rate_default[,as.character(year_start-10),]
        )/10
        rate_avg_5[is.nan(rate_avg_5)] <- 0
        rate_avg_5[is.infinite(rate_avg_5)] <- 0

        for(year in (year_start+1): max(dimnames(matrix_rate)[[2]]))
        {#  year <- 2022
          matrix_rate[,as.character(year),] <- matrix_rate[,as.character(year-1),] *  rate_avg_5
        }
      }
      return(matrix_rate)
    }

    data_long$value <- data_long$background_mortality_rate
    qB_full         <- data_long_to_array(data_long)


    ex_background       <- calculate_ex_matrix(qB_full)$ex

    data_long$value <- data_long$background_population
    pop_full        <- data_long_to_array(data_long )

    data_long$value <- data_long$mortality_undiagnosed_rate
    dDx_full        <- data_long_to_array(data_long ,run_projection=TRUE)

    data_long$value <- data_long$incidence_rate
    i_full          <- data_long_to_array(data_long ,run_projection=TRUE)/ 1e5

    data_long$value <- data_long$value_smr_non_minimal_care
    smr_matrix_n    <- data_long_to_array(data_long, run_projection=TRUE )

    data_long$value <- data_long$value_smr_minimal_care
    smr_matrix_m    <- data_long_to_array(data_long, run_projection=TRUE )

    data_long$value       <- data_long$value_percent_non_minimal_care
    qT1D_percent_n_full   <- data_long_to_array(data_long )



    Smr_ratio_adjusting <- function(smr_matrix_n,smr_value,lever_year_range,as_value_or_sealevel="sealevel")
    {
      ratio        <- smr_value/( apply(smr_matrix_n[,,as.character(20:50),drop=FALSE], c(1, 2), sum)/length(20:50))
      ratio_matrix <- array(ratio, dim = c(dim(ratio), 100))
      if(as_value_or_sealevel=="sealevel")
      {
        ratio_matrix[ratio_matrix>=1] <- 1
      }

      smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
      smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1
      smr_matrix_n_new <- smr_matrix_n
      smr_matrix_n_new[,as.character(lever_year_range),] <- smr_matrix_adjusted[,as.character(lever_year_range),]
      # smr_matrix_m_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]
      return(smr_matrix_n_new)
    }


    if(sim_enable)
    {
      # sim  apply simulation/invervention =====================================================================================================================================================================
      # sim_min_diag_rates ----------------------
      i_full_true      <- i_full/(1-dDx_full)
      dDx_full_inverse <- (1-dDx_full)
      dDx_full_inverse_years <- dDx_full_inverse[,as.character(sim_lever_year_range),] # only apply to ages 0 to 24

      # dDx_full_inverse_years[dDx_full_inverse_years < sim_parameters$sim_min_diag_rates/100] <- sim_parameters$sim_min_diag_rates/100 # taka parameter as  sealevel
      # dDx_full_inverse_years[] <- dDx_full_inverse_years + sim_parameters$sim_min_diag_rates/100 # taka parameter as increment
      dDx_full_inverse_years[,,as.character(0:24)] <-  sim_parameters$sim_min_diag_rates/100 ; dDx_full_inverse_years[dDx_full_inverse_years ==0] <- 0.000000000000001 # taka parameter as value

      dDx_full_inverse_years[dDx_full_inverse_years >=1] <- 1
      #



      dDx_full_inverse[,as.character(sim_lever_year_range),] <- dDx_full_inverse_years
      dDx_full <- 1- dDx_full_inverse
      i_full   <- i_full_true * (1-dDx_full)

      # sim_min_non_minimal_care_perc ---------------------
      qT1D_percent_n_full_years <- qT1D_percent_n_full[,as.character(sim_lever_year_range),]

      # qT1D_percent_n_full_years[qT1D_percent_n_full_years < sim_parameters$sim_min_non_minimal_care_perc/100] <- sim_parameters$sim_min_non_minimal_care_perc/100 # taka parameter as  sealevel
      # qT1D_percent_n_full_years <- qT1D_percent_n_full_years + sim_parameters$sim_min_non_minimal_care_perc/100 # taka parameter as increment
      qT1D_percent_n_full_years[] <- sim_parameters$sim_min_non_minimal_care_perc/100 # taka parameter as value

      qT1D_percent_n_full_years[qT1D_percent_n_full_years >=1] <- 1
      # qT1D_percent_n_full_years[qT1D_percent_n_full_years ==0] <- 0.000000000000001


      qT1D_percent_n_full[,as.character(sim_lever_year_range),] <- qT1D_percent_n_full_years



      # sim_min_non_minimal_care_level ---------------------
      # c("default_care","basic_care","best_care", "cure") )
      as_value_or_sealevel <- "value"
      # as_value_or_sealevel <- "sealevel"
      if( sim_parameters$sim_min_non_minimal_care_level == "default_care" )
      {
        # 1A  Human premix  15.3 31.1
        # smr_matrix_n   <- Smr_ratio_adjusting (smr_matrix_n,smr_value=23.2,sim_lever_year_range)
        # smr_matrix_m   <- Smr_ratio_adjusting (smr_matrix_m,smr_value=23.2,sim_lever_year_range)
      }
      if( sim_parameters$sim_min_non_minimal_care_level == "minimal_care" )
      {
      #   # 1BHuman premix 1-2 tests Limited materials, diabetologist 6.3 15.3
        smr_matrix_n   <- Smr_ratio_adjusting (smr_matrix_n,smr_value=10.8,sim_lever_year_range,as_value_or_sealevel=as_value_or_sealevel)
      #   # smr_matrix_m   <- Smr_ratio_adjusting (smr_matrix_m,smr_value=23.2,sim_lever_year_range)
      }
      if( sim_parameters$sim_min_non_minimal_care_level == "basic_care" )
      {
        # 2B1 Human basal/bolus 4 tests   3.7 4.4
        smr_matrix_n   <- Smr_ratio_adjusting (smr_matrix_n,smr_value=4.05,sim_lever_year_range,as_value_or_sealevel=as_value_or_sealevel)
        # smr_matrix_m   <- Smr_ratio_adjusting (smr_matrix_m,smr_value=4.05,sim_lever_year_range)
        # qT1D_percent_n_full[,as.character(sim_lever_year_range),] <- 1

      }
      if( sim_parameters$sim_min_non_minimal_care_level == "best_care" )
      {
        # 3C Pump CGM   2.2 2.8
        smr_matrix_n   <- Smr_ratio_adjusting (smr_matrix_n,smr_value=2.5,sim_lever_year_range,as_value_or_sealevel=as_value_or_sealevel)
        # smr_matrix_m   <- Smr_ratio_adjusting (smr_matrix_m,smr_value=2.5,sim_lever_year_range)
        # qT1D_percent_n_full[,as.character(sim_lever_year_range),] <- 1

      }
      if( sim_parameters$sim_min_non_minimal_care_level == "cure" )
      {
        # smr = 1
        smr_matrix_n   <- Smr_ratio_adjusting (smr_matrix_n,smr_value=1,sim_lever_year_range,as_value_or_sealevel=as_value_or_sealevel)
        # smr_matrix_m   <- Smr_ratio_adjusting (smr_matrix_m,smr_value=1,sim_lever_year_range)
        # qT1D_percent_n_full[,as.character(sim_lever_year_range),] <- 1
      }

    }






    diab_odds    <- qB_full / (1 - qB_full) * smr_matrix_n
    qT1D_n_full  <- diab_odds / (1 + diab_odds)
    diab_odds    <- qB_full / (1 - qB_full) * smr_matrix_m
    qT1D_m_full  <- diab_odds / (1 + diab_odds)



    {
      # run prevalence =====================================================================================================================================================================

      lever_year_range        <- as.character(min(years):max(years))  #  only apply lever to years after lever_change_at

      # base scenario -----------------------------------------------------------------------------------------------------------------------------------------
      i <- i_full; qB <- qB_full; qT1D_n <- qT1D_n_full; qT1D_m <- qT1D_m_full; qT1D_percent_n <- qT1D_percent_n_full; dDx <- dDx_full

      smr_matrix_life            <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

      # prev                 <- calculate_prevalence_tensor(i,      qB, qT1D_n,       qT1D_m,           qT1D_percent_n,      dDx ,smr_matrix=   smr_matrix_life ,run_complications=TRUE)
      prev                 <- calculate_prevalence_tensor(i,      qB, qT1D_n,       qT1D_m,           qT1D_percent_n,      dDx ,smr_matrix=   smr_matrix_life ,run_complications=FALSE)

      life_table_base_scenario   <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=smr_matrix_life ,years_range=as.character(1960:2040))

      # full diag -------------------------------------------------------------------- ------------------------------------------------------------------------
      i_100d   <- i;       i_100d[,lever_year_range,]   <- (i/(1-dDx))[,lever_year_range,]
      dDx_100d <- dDx;     dDx_100d[,lever_year_range,] <- 0

      prev_100d               <- calculate_prevalence_tensor(i_100d, qB, qT1D_n,       qT1D_m,           qT1D_percent_n,      dDx_100d) # deaths on diagnosis are converted to higer incidence

      # basic care, insulin, strip, education.  SMR  median( c(3.7, 4.4)) , 4.05,  across age 20-50 ------------------------------------------------------------------
      smr_matrix_n_new   <- Smr_ratio_adjusting (smr_matrix_n,smr_value=4.05,lever_year_range)
      diab_odds  <- qB / (1 - qB) * smr_matrix_n_new;  qT1D_n_bacare <- diab_odds / (1 + diab_odds)
      qT1D_percent_n_100d <- qT1D_percent_n_full
      qT1D_percent_n_100d[,as.character(lever_year_range),] <- 1

      prev_100d_basic_care    <- calculate_prevalence_tensor(i_100d, qB, qT1D_n_bacare, qT1D_m,          qT1D_percent_n_100d, dDx_100d) # smr all non-minimal

      smr_matrix_life            <- smr_matrix_n_new
      life_table_lever_2         <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=smr_matrix_life ,years_range=as.character(1960:2040))

      # best care  --------------------------------------------------------------------  --------------------------------------------------------------------
      smr_matrix_n_new <- Smr_ratio_adjusting (smr_matrix_n,smr_value=2.5,lever_year_range)
      diab_odds      <- qB / (1 - qB) * smr_matrix_n_new;  qT1D_n_becare <- diab_odds / (1 + diab_odds)

      prev_100d_best_care     <- calculate_prevalence_tensor(i_100d, qB, qT1D_n_becare, qT1D_m,          qT1D_percent_n_100d, dDx_100d) # deaths on diagnosis are converted to higer incidence

      smr_matrix_life            <- smr_matrix_n_new
      life_table_lever_3         <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=smr_matrix_life ,years_range=as.character(1960:2040))

      # cure -------------------------------------------------------------------- --------------------------------------------------------------------
      qT1D_n_cure <- qT1D_n
      qT1D_n_cure[,lever_year_range,] <- qB[,lever_year_range,,drop=FALSE]
      qT1D_m_cure <- qT1D_m
      qT1D_m_cure[,lever_year_range,] <- qB[,lever_year_range,,drop=FALSE]

      prev_100d_cure          <- calculate_prevalence_tensor(i_100d, qB, qT1D_n_cure,    qT1D_m_cure,    qT1D_percent_n_100d, dDx_100d)

      smr_matrix_life            <- smr_matrix_n_new ;smr_matrix_life[] <- 1
      life_table_lever_4         <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=smr_matrix_life ,years_range=as.character(1960:2040))



      if(switches$run_lever2023)
      {
        # Lever start at sim_lever_year_range ===============================================================================================================================================================
        # Lever start at sim_lever_year_range ===============================================================================================================================================================

        # full diag -------------------------------------------------------------------- ------------------------------------------------------------------------
        i_100d   <- i;       i_100d[,sim_lever_year_range,]   <- (i/(1-dDx))[,sim_lever_year_range,]
        dDx_100d <- dDx;     dDx_100d[,sim_lever_year_range,] <- 0

        prev_100d_lever2023               <- calculate_prevalence_tensor(i_100d, qB, qT1D_n,        qT1D_m,         qT1D_percent_n,      dDx_100d) # deaths on diagnosis are converted to higer incidence

        # basic care, insulin, strip, education.  SMR  median( c(3.7, 4.4)) , 4.05,  across age 20-50 ------------------------------------------------------------------
        smr_matrix_n_new   <- Smr_ratio_adjusting (smr_matrix_n,smr_value=4.05,sim_lever_year_range)
        diab_odds  <- qB / (1 - qB) * smr_matrix_n_new;  qT1D_n_bacare <- diab_odds / (1 + diab_odds)
        qT1D_percent_n_100d <- qT1D_percent_n_full
        qT1D_percent_n_100d[,as.character(sim_lever_year_range),] <- 1

        prev_100d_basic_care_lever2023    <- calculate_prevalence_tensor(i_100d, qB, qT1D_n_bacare, qT1D_m,         qT1D_percent_n_100d, dDx_100d) # smr all non-minimal

        # best care  --------------------------------------------------------------------  --------------------------------------------------------------------
        smr_matrix_n_new <- Smr_ratio_adjusting (smr_matrix_n,smr_value=2.5,sim_lever_year_range)
        diab_odds      <- qB / (1 - qB) * smr_matrix_n_new;  qT1D_n_becare <- diab_odds / (1 + diab_odds)

        prev_100d_best_care_lever2023     <- calculate_prevalence_tensor(i_100d, qB, qT1D_n_becare, qT1D_m,         qT1D_percent_n_100d, dDx_100d) # deaths on diagnosis are converted to higer incidence

        # cure -------------------------------------------------------------------- ------------------------------------------------------------------------
        qT1D_n_cure <- qT1D_n
        qT1D_n_cure[,sim_lever_year_range,] <- qB[,sim_lever_year_range,,drop=FALSE]
        qT1D_m_cure <- qT1D_m
        qT1D_m_cure[,sim_lever_year_range,] <- qB[,sim_lever_year_range,,drop=FALSE]

        prev_100d_cure_lever2023          <- calculate_prevalence_tensor(i_100d, qB, qT1D_n_cure,   qT1D_m_cure,    qT1D_percent_n_100d, dDx_100d)
      }
    }

    # Delay onset by 1,3,5,8,13 years.-------------------------------------
    if(switches$run_delay_onset)
    {
      Delay_onset_by_years_incidence <- function(i_full,years,delay_years=1)
      {
        i                 <- i_full   # mean(i[as.character(2010:2019),19]) *100000
        for(j in 99:delay_years)
        {# j<- 99 #  shift incidence matrix to the right by delay_years
          i[,as.character(2023:max(years)),as.character(j)] <- i[,as.character(2023:max(years)),as.character(j-delay_years)]
        }
        for(j in (delay_years-1):0)
        {# j<- 99 #  fill 0 to  incidence matrix from age 0 to age delay_years
          i[,as.character(2023:max(years)),as.character(j)] <- 0
        }
        i[,,as.character(80:99)] <- 0
        return(i)
      }

      i  <- Delay_onset_by_years_incidence (i_full, years ,delay_years=1)
      prev_delay_year_1   <- calculate_prevalence_tensor(i,      qB, qT1D_n,       qT1D_m,           qT1D_percent_n,      dDx     )
      i  <- Delay_onset_by_years_incidence (i_full, years ,delay_years=3)
      prev_delay_year_3   <- calculate_prevalence_tensor(i,      qB, qT1D_n,       qT1D_m,           qT1D_percent_n,      dDx     )
      i  <- Delay_onset_by_years_incidence (i_full, years ,delay_years=5)
      prev_delay_year_5   <- calculate_prevalence_tensor(i,      qB, qT1D_n,       qT1D_m,           qT1D_percent_n,      dDx     )
      i  <- Delay_onset_by_years_incidence (i_full, years ,delay_years=8)
      prev_delay_year_8   <- calculate_prevalence_tensor(i,      qB, qT1D_n,       qT1D_m,           qT1D_percent_n,      dDx     )
      i  <- Delay_onset_by_years_incidence (i_full, years ,delay_years=13)
      prev_delay_year_13   <- calculate_prevalence_tensor(i,      qB, qT1D_n,       qT1D_m,           qT1D_percent_n,      dDx     )

      Delay_onset_by_years_smr <- function(smr_matrix,years,delay_years=1)
      {
        smr_matrix[,as.character(2023:max(years)),as.character(10:(10+delay_years-1))] <- 1
        return(smr_matrix)
      }

      smr_matrix            <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

      life_table_delay_year_1   <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=Delay_onset_by_years_smr(smr_matrix,years,delay_years=1) ,years_range=as.character(1960:2040))
      life_table_delay_year_3   <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=Delay_onset_by_years_smr(smr_matrix,years,delay_years=3) ,years_range=as.character(1960:2040))
      life_table_delay_year_5   <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=Delay_onset_by_years_smr(smr_matrix,years,delay_years=5) ,years_range=as.character(1960:2040))
      life_table_delay_year_8   <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=Delay_onset_by_years_smr(smr_matrix,years,delay_years=8) ,years_range=as.character(1960:2040))
      life_table_delay_year_13   <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=Delay_onset_by_years_smr(smr_matrix,years,delay_years=13) ,years_range=as.character(1960:2040))
    }

    if(switches$run_site_data )
    {
      # Site data
      lever_year_range        <- as.character(min(years):max(years))  #  only apply lever to years after lever_change_at

      # calculating years gained  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
      Smr_ratio_adjusting_site <- function(smr_matrix_n,smr_value,lever_year_range)
      {
        ratio        <- smr_value/( apply(smr_matrix_n[,,as.character(20:50),drop=FALSE], c(1, 2), sum)/length(20:50))
        ratio_matrix <- array(ratio, dim = c(dim(ratio), 100))
        # ratio_matrix[ratio_matrix>=1] <- 1
        smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
        # smr_matrix_adjusted[smr_matrix_adjusted<1] <- smr_matrix_n[smr_matrix_adjusted<1]
        smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1
        smr_matrix_n_new <- smr_matrix_n
        smr_matrix_n_new[,as.character(lever_year_range),] <- smr_matrix_adjusted[,as.character(lever_year_range),]
        # smr_matrix_m_new[as.character(lever_year_range),] <- smr_matrix_adjusted[as.character(lever_year_range),]
        return(smr_matrix_n_new)
      }
      smr_matrix_n_strips_low <- Smr_ratio_adjusting_site (smr_matrix_n,smr_value=mean(5.3,6.3),lever_year_range)
      smr_matrix_n_strips_hig <- Smr_ratio_adjusting_site (smr_matrix_n,smr_value=mean(3.1,3.7),lever_year_range)
      smr_matrix_n_sensor_low <- Smr_ratio_adjusting_site (smr_matrix_n,smr_value=mean(4.0,4.7),lever_year_range)
      smr_matrix_n_sensor_hig <- Smr_ratio_adjusting_site (smr_matrix_n,smr_value=mean(2.2,2.8),lever_year_range)

      qT1D_percent_n_new <- qT1D_percent_n
      qT1D_percent_n_new [,as.character(lever_year_range),] <- 1

      smr_matrix         <- smr_matrix_m * (1- qT1D_percent_n_new)  +  smr_matrix_n_strips_low *     qT1D_percent_n_new
      life_table_strips_low <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=smr_matrix ,years_range=as.character(1960:2040))

      smr_matrix         <- smr_matrix_m * (1- qT1D_percent_n_new)  +  smr_matrix_n_strips_hig *     qT1D_percent_n_new
      life_table_strips_hig <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=smr_matrix ,years_range=as.character(1960:2040))

      smr_matrix         <- smr_matrix_m * (1- qT1D_percent_n_new)  +  smr_matrix_n_sensor_low *     qT1D_percent_n_new
      life_table_sensor_low <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=smr_matrix ,years_range=as.character(1960:2040))

      smr_matrix         <- smr_matrix_m * (1- qT1D_percent_n_new)  +  smr_matrix_n_sensor_hig *     qT1D_percent_n_new
      life_table_sensor_hig <- calculate_ex_lifetime_years_lost_matrix(qB=qB, smr_matrix=smr_matrix ,years_range=as.character(1960:2040))
    }

    #  quantify the HLYs actually delivered by research?   1970 -----------------------------------------------------------------------------------------------------------------------------------
    if(TRUE)
    {

      Smr_ratio_adjusting_site_stagnant <- function(smr_matrix_n,smr_value)
      {
        ratio        <- smr_value/( apply(smr_matrix_n[,,as.character(20:50),drop=FALSE], c(1, 2), sum)/length(20:50))
        ratio_matrix <- array(ratio, dim = c(dim(ratio), 100))
        ratio_matrix[ratio_matrix < 1] <- 1
        smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
        smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1
        smr_matrix_n_new    <- smr_matrix_adjusted
        # smr_matrix_n_new[,as.character(lever_year_range),] <- smr_matrix_adjusted[,as.character(lever_year_range),]
        return(smr_matrix_n_new)
      }

      # ratio        <- smr_min/(rowSums(smr_matrix_m[,as.character(20:50)])/length(20:50))
      # ratio_matrix <- t(matrix(ratio, nrow=100, ncol=length(ratio), byrow=TRUE))
      # ratio_matrix[ratio_matrix < 1] <- 1
      # smr_matrix_adjusted <- smr_matrix_m * ratio_matrix
      # smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1
      # smr_matrix_m <- smr_matrix_adjusted


      smr_stagnant      <- exp(10.5 *  0.3545 - 1.5274)
      smr_matrix_n_stagnant <-       Smr_ratio_adjusting_site_stagnant (smr_matrix_n,smr_stagnant)
      smr_matrix_m_stagnant <-       Smr_ratio_adjusting_site_stagnant (smr_matrix_m,smr_stagnant)


      diab_odds    <- qB_full / (1 - qB_full) * smr_matrix_n_stagnant
      qT1D_n_stagnant  <- diab_odds / (1 + diab_odds)
      diab_odds    <- qB_full / (1 - qB_full) * smr_matrix_m_stagnant
      qT1D_m_stagnant  <- diab_odds / (1 + diab_odds)

      prev_stagnant_lever1970  <- calculate_prevalence_tensor(i_full,      qB, qT1D_n_stagnant,       qT1D_m_stagnant,           qT1D_percent_n,      dDx     )

    }

    if(switches$run_stages_1_2)
    {
      # prev stage 1 scenario -----------------------------------
      i <- i_full; qB <- qB_full; qT1D_percent_n <- qT1D_percent_n_full; dDx <- dDx_full

      i   <- (i/(1-dDx))
      dDx[] <- 0

      years_from_s1_to_s2  <- 8  # shift i matrix backward 8 years
      i[,as.character( 1860 :(2040-years_from_s1_to_s2)), ]    <- i[,as.character((1860+years_from_s1_to_s2):2040),,drop=FALSE]
      # fill most recent years same as 2040
      for( years_ in (2040 - years_from_s1_to_s2):2040 )
      {
        i[,as.character(years_),]     <- i[,as.character(2040),]
      }
      # i[i< i_full] <- i_full[i< i_full]
      i[i< i_full/(1-dDx_full)] <- (i_full/(1-dDx_full))[i< i_full/(1-dDx_full)]


      prev_stage1    <- calculate_prevalence_tensor(i,      qB, qB,       qB,           qT1D_percent_n,      dDx     )

      # prev stage 2 scenario -----------------------------------
      i <- i_full; qB <- qB_full; qT1D_percent_n <- qT1D_percent_n_full; dDx <- dDx_full

      i   <- (i/(1-dDx))
      dDx[] <- 0

      years_from_s1_to_s2  <- 3  # shift i matrix backward 3 years
      i[,as.character( 1860 :(2040-years_from_s1_to_s2)), ]    <- i[,as.character((1860+years_from_s1_to_s2):2040),,drop=FALSE]
      # fill most recent years same as 2040
      for( years_ in (2040 - years_from_s1_to_s2):2040 )
      {
        i[,as.character(years_),]     <- i[,as.character(2040),]
      }
      # i[i< i_full] <- i_full[i< i_full]
      i[i< i_full/(1-dDx_full)] <- (i_full/(1-dDx_full))[i< i_full/(1-dDx_full)]

      prev_stage2    <- calculate_prevalence_tensor(i,      qB, qB,       qB,           qT1D_percent_n,      dDx     )



    }


    if(FALSE)
    {
      # prev stage 1 Approach 2 ----------------------------------- ----------------------------------- ----------------------------------- -----------------------------------
      i <- i_full; qB <- qB_full; qT1D_percent_n <- qT1D_percent_n_full; dDx <- dDx_full

      i   <- (i/(1-dDx))
      dDx[] <- 0


      # incidence cases
      i_c    <- pop_full * i
      i_c_s1 <- i_c ; i_c_s1[] <- 0

      # i_c_extra_80  <- i_c[,"2040",,drop=FALSE][rep(1, times = 80), ]
      i_c_extra_80  <- abind(replicate(80, i_c[,"2040",,drop=FALSE], simplify = FALSE), along = 2)

      dimnames(i_c_extra_80)[[2]] <- as.character(2041:2120)
      # i_c <- rbind(i_c, i_c_extra_80)
      i_c  <- abind(i_c,i_c_extra_80 , along = 2)

      # calculate stage1 incidence cases by assuming 60% population turn to Stage 3 by age 19.
      for(j in seq(1861, 2036, 5) )

        { # j <- 1970
        i_c_s1[,as.character( (j+ (4*0)): (j+(4*1)) ),as.character(0:4) ] <-  (     i_c[,as.character( (j+ (4*0)): (j+(4*1)) ),as.character(0:4) ] *  +
                                                                                     i_c[,as.character( (j+ (4*1)): (j+(4*2)) ),as.character(5:9) ]+
                                                                                     i_c[,as.character( (j+ (4*2)): (j+(4*3)) ),as.character(10:14) ]+
                                                                                     i_c[,as.character( (j+ (4*3)): (j+(4*4)) ),as.character(15:19) ]+
                                                                                     i_c[,as.character( (j+ (4*4)): (j+(4*5)) ),as.character(20:24) ]+
                                                                                     i_c[,as.character( (j+ (4*5)): (j+(4*6)) ),as.character(25:29) ]+
                                                                                     i_c[,as.character( (j+ (4*6)): (j+(4*7)) ),as.character(30:34) ]+
                                                                                     i_c[,as.character( (j+ (4*7)): (j+(4*8)) ),as.character(35:39) ]+
                                                                                     i_c[,as.character( (j+ (4*8)): (j+(4*9)) ),as.character(40:44) ]+
                                                                                     i_c[,as.character( (j+ (4*9)): (j+(4*10)) ),as.character(45:49) ]+
                                                                                     i_c[,as.character( (j+ (4*10)): (j+(4*11)) ),as.character(50:54) ]+
                                                                                     i_c[,as.character( (j+ (4*11)): (j+(4*12)) ),as.character(55:59) ]+
                                                                                     i_c[,as.character( (j+ (4*12)): (j+(4*13)) ),as.character(60:64) ]+
                                                                                     i_c[,as.character( (j+ (4*13)): (j+(4*14)) ),as.character(65:69) ]+
                                                                                     i_c[,as.character( (j+ (4*14)): (j+(4*15)) ),as.character(70:74) ]+
                                                                                     i_c[,as.character( (j+ (4*15)): (j+(4*16)) ),as.character(75:79) ]
        )
      }
      i_c_s1[,"1860",] <-  i_c_s1[,"1861",]
      i_c_s1          <-  i_c_s1 / (pop_full+0.00001)

      # prev_stage1                     <- calculate_prevalence(i_c_s1          , qB, qB, qB,   qT1D_percent_n, dDx ,smr_matrix,  years)
      prev_stage1                  <- calculate_prevalence_tensor(i_c_s1,      qB, qB,       qB,           qT1D_percent_n,      dDx     )

      prev_stage2                     <- prev_stage1
    }

    # Calculate population =========================================================================================================================================================
    pop_scale_factor       <- pop_full / (prev$P + prev$S)   # exclude the dead from population
    P_level         <- pop_scale_factor * prev$P  #  sum(P_level["2021",]) ; sum( (pop_scale_factor_100d * prev_100d$P) ["2021",])

    PC_level        <- pop_scale_factor * adrop(prev$C_P["AtLeast1C",,,,drop=FALSE], drop = 1)

    S_level         <- pop_scale_factor * prev$S
    I_flow          <- pop_scale_factor * prev$I
    DBGP_flow       <- pop_scale_factor * prev$DBGP

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

    BD_flow         <- pop_full * qB_full

    # Calculate missing prevalence  and delta of levers. --------------------------------------------------------------------------------
    ghost_ddx_level   <- (prev_100d$P              - prev$P) * pop_scale_factor #
    ghost_hba1c_level <- (prev_100d_cure$P         - prev_100d$P) * pop_scale_factor #

    ghost_basic_care  <- (prev_100d_basic_care$P   - prev_100d$P) * pop_scale_factor #  delta
    ghost_best_care   <- (prev_100d_best_care$P    - prev_100d_basic_care$P) * pop_scale_factor #
    ghost_cure        <- (prev_100d_cure$P         - prev_100d_best_care$P) * pop_scale_factor #  sum(ghost_hba1c_level["2021",])
    ghost_level <- ghost_hba1c_level + ghost_ddx_level #   sum(ghost_level["2021",])

    if(switches$run_lever2023)
    {
      # Calculate missing prevalence  and delta of levers2023. --------------------------------------------------------------------------------
      ghost_ddx_level_lever2023   <- (prev_100d_lever2023$P              - prev$P) * pop_scale_factor #
      ghost_hba1c_level_lever2023 <- (prev_100d_cure_lever2023$P         - prev_100d_lever2023$P) * pop_scale_factor #
      ghost_basic_care_lever2023  <- (prev_100d_basic_care_lever2023$P   - prev_100d_lever2023$P) * pop_scale_factor #  delta
      ghost_best_care_lever2023   <- (prev_100d_best_care_lever2023$P    - prev_100d_basic_care_lever2023$P) * pop_scale_factor #
      ghost_cure_lever2023        <- (prev_100d_cure_lever2023$P         - prev_100d_best_care_lever2023$P) * pop_scale_factor #  sum(ghost_hba1c_level["2021",])
      # ghost population is the sum of these two
      ghost_level_lever2023 <- ghost_hba1c_level_lever2023 + ghost_ddx_level_lever2023 #   sum(ghost_level["2021",])
    }

    # ,sim_parameters=list(sim_start_year=2023,sim_min_diag_rates=0,sim_min_non_minimal_care_perc=0,sim_min_non_minimal_care_level="base_care")

    # Merge result to dataframe -------------------------------------------------------------------------------------------------------------------------------
    prev_merge <-       as.data.frame.table(pop_full[,as.character(1960:2040),,drop=FALSE],stringsAsFactors = F)[,1:3]
    colnames(prev_merge) <- c("loc_id","Year","Age")
    prev_merge$"Country"                               <- data_long$world_bank_name[data_long$year %in% as.character(1960:2040)]
    prev_merge$"dDx_full"                              <- as.data.frame.table(dDx_full[,as.character(1960:2040),,drop=FALSE],stringsAsFactors = F)[,4]
    prev_merge$"qT1D_percent_n_full"                   <- as.data.frame.table(qT1D_percent_n_full[,as.character(1960:2040),,drop=FALSE],stringsAsFactors = F)[,4]

    prev_merge$"sim_start_year"                        <- sim_parameters$sim_start_year
    prev_merge$"sim_min_diag_rates"                    <- sim_parameters$sim_min_diag_rates
    prev_merge$"sim_min_non_minimal_care_perc"         <- sim_parameters$sim_min_non_minimal_care_perc
    prev_merge$"sim_min_non_minimal_care_level"        <- sim_parameters$sim_min_non_minimal_care_level

    prev_merge$"Ann. background population"           <- as.data.frame.table(pop_full[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ann. background mortality"            <- as.data.frame.table(BD_flow[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Prevalence"                           <- as.data.frame.table(P_level[,as.character(1960:2040),,drop=FALSE])[,4]
    # prev_merge$"Prevalence with AtLeast1C"            <- as.data.frame.table(PC_level[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Incidence (1 base)"                   <- as.data.frame.table(I_flow[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Incidence (2 diagnosis)"              <- as.data.frame.table(I_flow_diagnosis[,as.character(1960:2040),,drop=FALSE])[,4]

    prev_merge$"Ghosts"                                        <- as.data.frame.table(ghost_level[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ghosts (onset death)"                          <- as.data.frame.table(ghost_ddx_level[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ghosts (early death)"                          <- as.data.frame.table(ghost_hba1c_level[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ghosts (delta basic care)"                     <- as.data.frame.table(ghost_basic_care[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ghosts (delta best care)"                      <- as.data.frame.table(ghost_best_care[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ghosts (delta cure)"                           <- as.data.frame.table(ghost_cure[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ann. onset deaths"                             <- as.data.frame.table(DDx_flow[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ann. early deaths"                             <- as.data.frame.table(DT1D_flow[,as.character(1960:2040),,drop=FALSE])[,4]

    prev_merge$"Ann. early deaths (background)"                               <- as.data.frame.table(DBGP_flow[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ann. early deaths (2 diagnosis)"                              <- as.data.frame.table(DT1D_flow_100d[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ann. early deaths (3 basic care)"                             <- as.data.frame.table(DT1D_flow_basic_care[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ann. early deaths (4 best care)"                              <- as.data.frame.table(DT1D_flow_best_care[,as.character(1960:2040),,drop=FALSE])[,4]
    prev_merge$"Ann. early deaths (5 cure)"                                   <- as.data.frame.table(DT1D_flow_cure[,as.character(1960:2040),,drop=FALSE])[,4]

    if(switches$run_lever2023)
    {
      prev_merge$"Ghosts sim_start_year"                                                      <- as.data.frame.table(ghost_level_lever2023[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Ghosts (onset death) sim_start_year"                                        <- as.data.frame.table(ghost_ddx_level_lever2023[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Ghosts (early death) sim_start_year"                                        <- as.data.frame.table(ghost_hba1c_level_lever2023[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Ghosts (delta basic care) sim_start_year"                                   <- as.data.frame.table(ghost_basic_care_lever2023[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Ghosts (delta best care) sim_start_year"                                    <- as.data.frame.table(ghost_best_care_lever2023[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Ghosts (delta cure) sim_start_year"                                         <- as.data.frame.table(ghost_cure_lever2023[,as.character(1960:2040),,drop=FALSE])$Freq

    }

    prev_merge$"Life expectency (1 background)"                                         <- as.data.frame.table(ex_background[,as.character(1960:2040),,drop=FALSE])$Freq

    prev_merge$"Life expectency (2 t1d base)"                                           <-  as.data.frame.table(life_table_base_scenario$ex[,as.character(1960:2040),,drop=FALSE])$Freq
    diagnosis_rate <- (prev_merge$`Incidence (1 base)`+0.001) / (prev_merge$`Ann. onset deaths`  + prev_merge$`Incidence (1 base)`+0.001)
    prev_merge$"Life expectency (2 t1d base)"                           <- prev_merge$"Life expectency (2 t1d base)" * diagnosis_rate +  (1-diagnosis_rate)*0.5
    prev_merge$"Lifetime years lost (2 t1d base) (complication)"        <- as.data.frame.table(life_table_base_scenario$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq *  prev_merge$"Life expectency (2 t1d base)"/ as.data.frame.table(life_table_base_scenario$ex[,as.character(1960:2040),])$Freq
    prev_merge$"Lifetime years lost (2 t1d base) (treatment)"           <- as.data.frame.table(life_table_base_scenario$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq    *  prev_merge$"Life expectency (2 t1d base)"/ as.data.frame.table(life_table_base_scenario$ex[,as.character(1960:2040),])$Freq



    prev_merge$"Life expectency (3 t1d diagnosis)"                                         <- as.data.frame.table(life_table_base_scenario$ex[,as.character(1960:2040),,drop=FALSE])$Freq
    prev_merge$"Lifetime years lost (3 t1d diagnosis) (complication)"                      <- as.data.frame.table(life_table_base_scenario$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
    prev_merge$"Lifetime years lost (3 t1d diagnosis) (treatment)"                         <- as.data.frame.table(life_table_base_scenario$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq

    prev_merge$"Life expectency (4 t1d basic care)"                                         <- as.data.frame.table(life_table_lever_2$ex[,as.character(1960:2040),,drop=FALSE])$Freq
    prev_merge$"Lifetime years lost (4 t1d basic care) (complication)"                      <- as.data.frame.table(life_table_lever_2$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
    prev_merge$"Lifetime years lost (4 t1d basic care) (treatment)"                         <- as.data.frame.table(life_table_lever_2$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq

    prev_merge$"Life expectency (5 t1d best care)"                                         <- as.data.frame.table(life_table_lever_3$ex[,as.character(1960:2040),,drop=FALSE])$Freq
    prev_merge$"Lifetime years lost (5 t1d best care) (complication)"                      <- as.data.frame.table(life_table_lever_3$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
    prev_merge$"Lifetime years lost (5 t1d best care) (treatment)"                         <- as.data.frame.table(life_table_lever_3$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq

    prev_merge$"Life expectency (6 t1d cure)"                                         <- as.data.frame.table(life_table_lever_4$ex[,as.character(1960:2040),,drop=FALSE])$Freq
    prev_merge$"Lifetime years lost (6 t1d cure) (complication)"                      <- as.data.frame.table(life_table_lever_4$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
    prev_merge$"Lifetime years lost (6 t1d cure) (treatment)"                         <- as.data.frame.table(life_table_lever_4$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq

    # Delay by 1,3,5,8,12 years -------------------------------------------------------------------------------------------------------------------------------------------

    prev_merge$"Prevalence stagnant 1970"                                                 <- as.data.frame.table((pop_scale_factor * prev_stagnant_lever1970$P)[,as.character(1960:2040),,drop=FALSE])$Freq

    if(switches$run_stages_1_2)
    {
      prev_merge$"Prevalence Stage 1 only"                                        <- as.data.frame.table((pop_scale_factor * (prev_stage1$P -  prev_stage2$P))[,as.character(1960:2040),,drop=FALSE] )$Freq
      prev_merge$"Prevalence Stage 2 only"                                        <- as.data.frame.table((pop_scale_factor * (prev_stage2$P -  prev_100d_cure$P))[,as.character(1960:2040),,drop=FALSE] )$Freq
    }

    if(switches$run_delay_onset)
    {
      prev_merge$"Prevalence delay onset by 1 years"                                        <- as.data.frame.table((pop_scale_factor * prev_delay_year_1$P)[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Prevalence delay onset by 3 years"                                        <- as.data.frame.table((pop_scale_factor * prev_delay_year_3$P)[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Prevalence delay onset by 5 years"                                        <- as.data.frame.table((pop_scale_factor * prev_delay_year_5$P)[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Prevalence delay onset by 8 years"                                        <- as.data.frame.table((pop_scale_factor * prev_delay_year_8$P)[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Prevalence delay onset by 13 years"                                       <- as.data.frame.table((pop_scale_factor * prev_delay_year_13$P)[,as.character(1960:2040),,drop=FALSE])$Freq

      prev_merge$"Lifetime years lost (delay onset year 1) (treatment)"      <- as.data.frame.table(life_table_delay_year_1$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (delay onset year 1) (complication)"   <- as.data.frame.table(life_table_delay_year_1$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Life expectency (delay onset year 1)"                      <- as.data.frame.table(life_table_delay_year_1$ex[,as.character(1960:2040),,drop=FALSE])$Freq

      prev_merge$"Lifetime years lost (delay onset year 3) (treatment)"      <- as.data.frame.table(life_table_delay_year_3$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (delay onset year 3) (complication)"   <- as.data.frame.table(life_table_delay_year_3$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Life expectency (delay onset year 3)"                      <- as.data.frame.table(life_table_delay_year_3$ex[,as.character(1960:2040),,drop=FALSE])$Freq

      prev_merge$"Lifetime years lost (delay onset year 5) (treatment)"      <- as.data.frame.table(life_table_delay_year_5$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (delay onset year 5) (complication)"   <- as.data.frame.table(life_table_delay_year_5$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Life expectency (delay onset year 5)"                      <- as.data.frame.table(life_table_delay_year_5$ex[,as.character(1960:2040),,drop=FALSE])$Freq

      prev_merge$"Lifetime years lost (delay onset year 8) (treatment)"      <- as.data.frame.table(life_table_delay_year_8$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (delay onset year 8) (complication)"   <- as.data.frame.table(life_table_delay_year_8$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Life expectency (delay onset year 8)"                      <- as.data.frame.table(life_table_delay_year_8$ex[,as.character(1960:2040),,drop=FALSE])$Freq

      prev_merge$"Lifetime years lost (delay onset year 13) (treatment)"      <- as.data.frame.table(life_table_delay_year_13$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (delay onset year 13) (complication)"   <- as.data.frame.table(life_table_delay_year_13$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Life expectency (delay onset year 13)"                      <- as.data.frame.table(life_table_delay_year_13$ex[,as.character(1960:2040),,drop=FALSE])$Freq

    }
    if(switches$run_site_data )
    {
      #  For site
      prev_merge$"Life expectency (strip low)"            <- as.data.frame.table(life_table_strips_low$ex[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (strip low)"        <- as.data.frame.table(life_table_strips_low$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq + as.data.frame.table(life_table_strips_low$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq

      prev_merge$"Life expectency (strip hig)"            <- as.data.frame.table(life_table_strips_hig$ex[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (strip hig)"        <- as.data.frame.table(life_table_strips_hig$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq + as.data.frame.table(life_table_strips_hig$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq

      prev_merge$"Life expectency (sensor low)"            <- as.data.frame.table(life_table_sensor_low$ex[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (sensor low)"        <- as.data.frame.table(life_table_sensor_low$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq + as.data.frame.table(life_table_sensor_low$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq

      prev_merge$"Life expectency (sensor hig)"            <- as.data.frame.table(life_table_sensor_hig$ex[,as.character(1960:2040),,drop=FALSE])$Freq
      prev_merge$"Lifetime years lost (sensor hig)"        <- as.data.frame.table(life_table_sensor_hig$ex_complication[,as.character(1960:2040),,drop=FALSE])$Freq + as.data.frame.table(life_table_sensor_hig$ex_treatment[,as.character(1960:2040),,drop=FALSE])$Freq
    }

    prev_merge$"% Odds living to"  <- prev_merge$Prevalence / (prev_merge$Prevalence + prev_merge$Ghosts) * 100
    prev_merge$"% Odds living to"[is.na(prev_merge$"% Odds living to")] <- 0

    # 1 in x families
    if(TRUE)
    {
      # 1 in ** families. country_indicator------------------------------------------------------------------------------------------------------------------------------------
      # country_indicator <- dplyr::select(readRDS("data_wb/country_indicator_imputed_avg_0.4.13.Rds"),country,year,fertility)
      country_indicator <- dplyr::select(readRDS("data_wb/country_indicator_imputed_avg_fertility_0.4.13.Rds"),country,year,fertility,loc_id_country_level=loc_id)
      country_indicator$year <- as.character(country_indicator$year)
      days_              <- setDT(prev_merge)[,list(background_population=sum(`Ann. background population`),`Total`= sum(Ghosts+Prevalence)),by=c("loc_id","Year")]

      days_$loc_id_country_level <- as.numeric(str_extract(days_$loc_id, "\\d+"))

      days_              <- days_ %>% dplyr::left_join(dplyr::select(country_indicator,loc_id_country_level,Year=year,fertility),by=c("loc_id_country_level","Year"))
      days_$family_size  <- days_$fertility^2 + days_$fertility +2
      days_$`1 in x families` <- 1/(1-((days_$background_population/days_$`Total`-1)/(days_$background_population/days_$`Total`) ) ^ days_$family_size)
      days_$`1 in x families`[is.nan(days_$`1 in x families`)] <- 0

      prev_merge   <- setDF(prev_merge) %>% dplyr::inner_join(dplyr::select(days_,loc_id,Year,`1 in x families` = `1 in x families`
                                                                               # ,`Lifetime years lost (2 t1d base) (treatment)`   =lifetime_years_lost_treatment
      ), by=c("loc_id","Year"))
    }

    # Final cleaning, Order and round and cut ,clean ---------------------------------------------------------

    prev_merge <- prev_merge %>% mutate_if(is.numeric, round, digits=2)
    prev_merge$Year    <- as.numeric(prev_merge$Year)
    prev_merge$Age     <- as.numeric(prev_merge$Age)
    prev_merge$loc_id  <- as.factor(prev_merge$loc_id)

    prev_merge <- prev_merge[prev_merge$Year>=1960,]
    prev_merge <- prev_merge[ with(prev_merge, order( loc_id, Year,Age)), ]   # order properly for transforming to matrix

    # prev_merge$"Life expectency (1 background)"      <- life_table_background$ex[life_table_background$year>=1960]

  # },
  # error = function(cond) {
  #   # Log('Error running model for %s:')
  #   sink("log.txt",append=TRUE);cat(paste0( "  ",  cond, " \n") );sink()
  #   NA
  # },
  # warning = function(cond) {
  #   # Log('Warning running model for %s:')
  #   NA
  # },
  # finally = {
  # })
  # write_parquet(prev_merge, paste0("../../temp/prev_merge_",input) )


  return(prev_merge)
}

if(FALSE)
{
  library(parallel)
  source("code_R/refresh_country_files_tensor_utils.R")
  num_thread <- 10
  clust      <- parallel::makeCluster(num_thread, type = "PSOCK")  # stopCluster(clust)
  clusterExport(cl=clust, varlist=c("calculate_prevalence_tensor","calculate_ex","calculate_ex_gpu","extract_for_purpose"))
  system.time({a <- clusterApply(clust, 1:10, refresh_country_files_tensor)})
  stopCluster(clust)
  Sys.time()

  system.time({refresh_country_files_tensor(1)})

}
