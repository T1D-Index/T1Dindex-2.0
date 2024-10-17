
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
prevalence_and_ghost_pop <- function(
    loc_id,
    # hba1c=make_age_function(7.5),
    start_year,
    end_year,
    data_dir,
    smr_scale_factor=1,
    incidence_scale_factor=1,
    data_long = NULL,
    matrices_list = NULL,
    log_timing=FALSE,
    full_run=TRUE,
    diagnosis_input=1 # 1 is t1d, 2 is ward , 3 is average
) {
  # smr_scale_factor=1 ;   incidence_scale_factor=1;matrices_list = NULL;log_timing=FALSE; start_year=1960; end_year=2040;full_run=TRUE ; data_long <- data_long_; diagnosis_input=1

  #----------- life expectancy t1ds -----------------------------------------------------------------------
  Get_life_expectancy_t1d <- function(qB, smr_matrix)
  {

    qB              <- qB[as.character(1960:2040),]
    smr_matrix      <- smr_matrix[as.character(1960:2040),]
    # qT1D_percent_n  <- qT1D_percent_n[as.character(1960:2040),]
    #
    # smr_matrix <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n
    smr_matrix_long    <- as.data.frame.table(smr_matrix,stringsAsFactors = F)%>% mutate_all(as.numeric) # when shiny run simulation, apply percent back to data_long
    colnames(smr_matrix_long) <- c("year",  "age",   "smr")

    t1d_mortality_matrix <- Get_qT1D_from_smr_matrix (qB, smr_matrix)

    t1d_mortality_long <- as.data.frame.table(t1d_mortality_matrix,stringsAsFactors = F)%>% mutate_all(as.numeric) # when shiny run simulation, apply percent back to data_long
    colnames(t1d_mortality_long) <- c("year",  "age",   "Value")
    t1d_mortality_long$smr <- smr_matrix_long$smr
    t1d_mortality_long <- t1d_mortality_long[t1d_mortality_long$year >= 1900,]
    t1d_mortality_long <- t1d_mortality_long[with(t1d_mortality_long, order( year,age)), ]

    t1d_mortality_long    <- dplyr::select(t1d_mortality_long, year,age,smr,qx =Value)
    if(FALSE)
    { # try speed up, did not work.
      t1d_mortality_long1 <- t1d_mortality_long; t1d_mortality_long1$run_no <- 1
      t1d_mortality_long2 <- t1d_mortality_long; t1d_mortality_long2$run_no <- 2
      life_table_t1d        <- calculate_ex_lifetime_years_lost(rbind(t1d_mortality_long2,t1d_mortality_long1 ))
    }

    life_table_t1d        <- calculate_ex_lifetime_years_lost(t1d_mortality_long)
    life_table_t1d
  }


  if(log_timing){ sink("log.txt",append=TRUE);cat(paste0(Sys.time()," prevalence_and_ghost_pop() \n") );sink()}

  ptm <- proc.time()
  # function to drop warm up periods
  dwu <- function(X) X[-(1:MAX_AGE),]
  # warm up period to seed prevalence: start model MAX_AGE years earlier
  years <- (start_year - MAX_AGE):end_year

  # Take draws from prevalence posterior distribution by simulating
  # country_wb_name <- "Morocco"
  # loc_id      <- get_loc_id(country_wb_name)
  # data_long   <- run_query_df (paste0("SELECT * FROM input_rates_combined  WHERE loc_id = '",loc_id,"'" ) )

  if(is.null(data_long))
    {
      loc_id_file_name <- as.numeric(str_extract(loc_id, "\\d+"))
      # data_long    <- setDF(readRDS(paste0(data_dir,"/",loc_id_file_name,".Rds")))
      data_long   <- setDF(arrow::read_feather(paste0(data_dir,"/",loc_id_file_name,".feather")))
      data_long   <- data_long[data_long$loc_id==loc_id,]
    }

  if(diagnosis_input==2)
  {
    data_long$mortality_undiagnosed_rate <- data_long$mortality_undiagnosed_rate_ward
  }
  if(diagnosis_input==3)
  {
    data_long$mortality_undiagnosed_rate <- data_long$mortality_undiagnosed_rate_average
  }
  # data_long_default_run   <- setDF(readRDS(paste0(strsplit(data_dir,"_mcmc")[[1]][1],"/",as.character(loc_id),".Rds")))

  country_wb_name <- unique(data_long$world_bank_name)

  if(is.null(matrices_list))
  {
    matrices_list <- data_long_2_matrices (data_long=data_long,data_long_default_run=data_long)
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
  lever_year_range        <- as.character(min(years):max(years))  #  only apply lever to years after lever_change_at


  # life expectency , background ---------------------------------------------------------------------
  life_table_background    <- dplyr::select(data_long, year,age,qx =background_mortality_rate)
  # life_table_background    <- calculate_ex(life_table_background)
  life_table_background    <- calculate_ex(life_table_background)
  # life_table_background[life_table_background$year==2022&life_table_background$age==10 ,]

  # make sure t1d mortality ( smr * backgroun_mortality  ) is always <=1 -------
  Get_qT1D_from_smr_matrix <- function(qB, smr_matrix)
  {
    diab_odds  <- qB / (1 - qB) * smr_matrix
    qT1D       <- diab_odds / (1 + diab_odds)
    qT1D
  }


  # life_table_base_scenario[life_table_base_scenario$year==2022&life_table_base_scenario$age==10 ,]
  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
  smr_matrix               <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n
  # prev base scenario  --------------------------------------------------------------------

  prev                     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years,run_complications=TRUE)
  # prev                     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years,run_complications=FALSE)

  life_table_base_scenario <- Get_life_expectancy_t1d (qB, smr_matrix)

  if(TRUE)
  {
    # prev stage 1 Approach 1 -----------------------------------
    i                 <- matrices_list$i   # mean(i[as.character(2010:2019),19]) *100000
    dDx               <- matrices_list$dDx
    qT1D_percent_n    <- matrices_list$qT1D_percent_n

    years_from_s1_to_s2  <- 8  # shift i matrix backward 8 years

    i   <- (i/(1-dDx))
    dDx[] <- 0

    i_cure <- i

    i[as.character( (1860) :(2040-years_from_s1_to_s2)), ]    <- i[as.character(( (1860)+years_from_s1_to_s2):2040),]
    i[as.character((2040 - years_from_s1_to_s2):2040),]       <-  matrix(i[as.character(2040),], nrow=years_from_s1_to_s2+1, ncol=ncol(i), byrow=TRUE)

    i[i< matrices_list$i/(1-matrices_list$dDx)] <- (matrices_list$i/(1-matrices_list$dDx))[i< matrices_list$i/(1-matrices_list$dDx)]


    # ratio_width    <-  ceiling(100/((20-years_from_s1_to_s2)/20))
    # i_extend_age   <-  i[,1:(ratio_width-99) ]
    # i_extend_age[] <- 0
    # dimnames(i_extend_age)[[2]] <- as.character(100:ratio_width)
    # i <- cbind(i, i_extend_age)
    # i[,  as.character(80:ratio_width)] <- i[,as.character(79)]

    # i_new  <- i ; i_new[] <- 0
    # for(j in 0:80)
    # {
    #   i_new[,as.character(j)] <- apply(i[,as.character((j+years_from_s1_to_s2):(j+years_from_s1_to_s2*2))], c(1), median) # lapply(i[,as.character(j:(j+years_from_s1_to_s2))])
    # }

    # i_new[i_new< matrices_list$i/(1-matrices_list$dDx)] <- (matrices_list$i/(1-matrices_list$dDx))[i_new< matrices_list$i/(1-matrices_list$dDx)]

    # i <- i_new

    if(FALSE)
    {
      data.frame(x= 0:99
                 ,y1=  i["2000",]
                 # ,y2=  i_new["2000",]
                 # y= pop["2000",]
      ) %>% echarts4r::e_chart(x) %>% echarts4r::e_line(y1)%>% echarts4r::e_line(y2) %>% e_tooltip(trigger="axis")

      data.frame(x= as.character(1960:2040)
                 , y1= i[as.character(1960:2040),"19"]
                 , y2= (matrices_list$i/(1-matrices_list$dDx)) [as.character(1960:2040),"19"]
      ) %>%
        echarts4r::e_chart(x) %>%
        echarts4r::e_line(y1) %>%
        echarts4r::e_line(y2) %>%
        e_tooltip(trigger="axis")

      data.frame(x= as.character(1960:2040)
                 , y1= dDx[as.character(1960:2040),"19"]
                 , y1= dDx[as.character(1960:2040),"19"]
      ) %>%
        echarts4r::e_chart(x) %>% echarts4r::e_line(y1) %>% e_tooltip(trigger="axis")
    }

    prev_stage1                     <- calculate_prevalence(i        , qB, qB, qB,   qT1D_percent_n, dDx ,smr_matrix,  years)

    # prev stage  2 scenario -----------------------------------
    # prev stage 1 Approach 1 -----------------------------------
    i                 <- matrices_list$i   # mean(i[as.character(2010:2019),19]) *100000
    dDx               <- matrices_list$dDx
    qT1D_percent_n    <- matrices_list$qT1D_percent_n

    years_from_s1_to_s2  <- 3  # shift i matrix backward 15 years

    i   <- (i/(1-dDx))
    dDx[] <- 0

    i_cure <- i

    i[as.character( (1860) :(2040-years_from_s1_to_s2)), ]    <- i[as.character(( (1860)+years_from_s1_to_s2):2040),]
    i[as.character((2040 - years_from_s1_to_s2):2040),]       <-  matrix(i[as.character(2040),], nrow=years_from_s1_to_s2+1, ncol=ncol(i), byrow=TRUE)

    i[i< matrices_list$i/(1-matrices_list$dDx)] <- (matrices_list$i/(1-matrices_list$dDx))[i< matrices_list$i/(1-matrices_list$dDx)]

    prev_stage2                     <- calculate_prevalence(i        , qB, qB, qB,   qT1D_percent_n, dDx ,smr_matrix,  years)

  }

  if(FALSE)
  {
    # prev stage 1 Approach 2 ----------------------------------- ----------------------------------- ----------------------------------- -----------------------------------
    i                 <- matrices_list$i   # mean(i[as.character(2010:2019),19]) *100000
    dDx               <- matrices_list$dDx
    qT1D_percent_n    <- matrices_list$qT1D_percent_n
    i   <- (i/(1-dDx))
    dDx[] <- 0
    years_from_s1_to_s3  <- 20  # shift i matrix backward 15 years
    # incidence cases
    i_c    <- pop * i
    i_c_s1 <- i_c ; i_c_s1[] <- 0

    i_c_extra_80  <- i_c["2040",,drop=FALSE][rep(1, times = 80), ]
    dimnames(i_c_extra_80)[[1]] <- as.character(2041:2120)
    i_c <- rbind(i_c, i_c_extra_80)

    # calculate stage1 incidence cases by assuming 60% population turn to Stage 3 by age 19.
    for(j in seq(1861, 2036, 5) )
    { # j <- 1970
      i_c_s1[as.character( (j+ (4*0)): (j+(4*1)) ),as.character(0:4) ] <-  (     i_c[as.character( (j+ (4*0)): (j+(4*1)) ),as.character(0:4) ] *  +
                                                                                   i_c[as.character( (j+ (4*1)): (j+(4*2)) ),as.character(5:9) ]+
                                                                                   i_c[as.character( (j+ (4*2)): (j+(4*3)) ),as.character(10:14) ]+
                                                                                   i_c[as.character( (j+ (4*3)): (j+(4*4)) ),as.character(15:19) ]+
                                                                                   i_c[as.character( (j+ (4*4)): (j+(4*5)) ),as.character(20:24) ]+
                                                                                   i_c[as.character( (j+ (4*5)): (j+(4*6)) ),as.character(25:29) ]+
                                                                                   i_c[as.character( (j+ (4*6)): (j+(4*7)) ),as.character(30:34) ]+
                                                                                   i_c[as.character( (j+ (4*7)): (j+(4*8)) ),as.character(35:39) ]+
                                                                                   i_c[as.character( (j+ (4*8)): (j+(4*9)) ),as.character(40:44) ]+
                                                                                   i_c[as.character( (j+ (4*9)): (j+(4*10)) ),as.character(45:49) ]+
                                                                                   i_c[as.character( (j+ (4*10)): (j+(4*11)) ),as.character(50:54) ]+
                                                                                   i_c[as.character( (j+ (4*11)): (j+(4*12)) ),as.character(55:59) ]+
                                                                                   i_c[as.character( (j+ (4*12)): (j+(4*13)) ),as.character(60:64) ]+
                                                                                   i_c[as.character( (j+ (4*13)): (j+(4*14)) ),as.character(65:69) ]+
                                                                                   i_c[as.character( (j+ (4*14)): (j+(4*15)) ),as.character(70:74) ]+
                                                                                   i_c[as.character( (j+ (4*15)): (j+(4*16)) ),as.character(75:79) ]
      )
    }
    i_c_s1["1860",] <-  i_c_s1["1861",]
    i_c_s1          <-  i_c_s1 / (pop+0.00001)

    prev_stage1                     <- calculate_prevalence(i_c_s1        , qB, qB, qB,   qT1D_percent_n, dDx ,smr_matrix,  years)
    prev_stage2                     <- prev_stage1

  }



  if(FALSE)
  {
    data_output <- data.frame(i); data_output$Year <- row.names(data_output); data_output <- gather(data_output ,key = "Age", value = "Value",-Year);data_output$Age <- as.numeric(gsub("X","", data_output$Age))
    e1 <- generate_3D_matrix ( data_output, "Prevalence")
    e1
  }
  # Delay onset by 1,3,5,8,13 years.------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # if(FALSE)
  if(full_run)
  {
    Delay_onset_by_years_incidence <- function(matrices_list,years,delay_years=1)
    {
      i                 <- matrices_list$i   # mean(i[as.character(2010:2019),19]) *100000
      for(j in 99:delay_years)
      {# j<- 99 #  shift incidence matrix to the right by delay_years
        i[as.character(2023:max(years)),as.character(j)] <- i[as.character(2023:max(years)),as.character(j-delay_years)]
      }
      for(j in (delay_years-1):0)
      {# j<- 99 #  fill 0 to  incidence matrix from age 0 to age delay_years
        i[as.character(2023:max(years)),as.character(j)] <- 0
      }
      i[,as.character(80:99)] <- 0
      return(i)
    }

    i  <- Delay_onset_by_years_incidence (matrices_list, years ,delay_years=1)
    prev_delay_year_1                     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years)
    i  <- Delay_onset_by_years_incidence (matrices_list, years ,delay_years=3)
    prev_delay_year_3                     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years)
    i  <- Delay_onset_by_years_incidence (matrices_list, years ,delay_years=5)
    prev_delay_year_5                     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years)
    i  <- Delay_onset_by_years_incidence (matrices_list, years ,delay_years=8)
    prev_delay_year_8                     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years)
    i  <- Delay_onset_by_years_incidence (matrices_list, years ,delay_years=13)
    prev_delay_year_13                    <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years)

    Delay_onset_by_years_smr <- function(smr_matrix,years,delay_years=1)
    {
      smr_matrix[as.character(2023:max(years)),as.character(10:(10+delay_years-1))] <- 1
      return(smr_matrix)
    }

    life_table_delay_year_1 <- Get_life_expectancy_t1d (qB, Delay_onset_by_years_smr(smr_matrix,years,delay_years=1) )
    life_table_delay_year_3 <- Get_life_expectancy_t1d (qB, Delay_onset_by_years_smr(smr_matrix,years,delay_years=3) )
    life_table_delay_year_5 <- Get_life_expectancy_t1d (qB, Delay_onset_by_years_smr(smr_matrix,years,delay_years=5) )
    life_table_delay_year_8 <- Get_life_expectancy_t1d (qB, Delay_onset_by_years_smr(smr_matrix,years,delay_years=8) )
    life_table_delay_year_13 <- Get_life_expectancy_t1d (qB, Delay_onset_by_years_smr(smr_matrix,years,delay_years=13) )
  }



  # prev                    <- calculate_prevalence_gpu(i        , qB, qT1D_n, qT1D_m, qT1D_percent_n, dDx ,  years)
  # Lever 1  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Lever 1  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  matrices_list_lever1 <- Apply_levers_to_input_matrices (matrices_list, lever=1, lever_year_range=lever_year_range)
  i                 <- matrices_list_lever1$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list_lever1$dDx
  smr_matrix_n      <- matrices_list_lever1$smr_matrix_n
  smr_matrix_m      <- matrices_list_lever1$smr_matrix_m
  qT1D_percent_n    <- matrices_list_lever1$qT1D_percent_n

  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
  smr_matrix               <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

  prev_100d               <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years) # deaths on diagnosis are converted to higer incidence

  # run lever start at 2023 -----------------------------------------------------------------------------------------------------------------------------------

    {

    matrices_list_lever1 <- Apply_levers_to_input_matrices (matrices_list, lever=1, lever_year_range= as.character(2023:max(years)))
    i                 <- matrices_list_lever1$i   # mean(i[as.character(2010:2019),19]) *100000
    dDx               <- matrices_list_lever1$dDx
    smr_matrix_n      <- matrices_list_lever1$smr_matrix_n
    smr_matrix_m      <- matrices_list_lever1$smr_matrix_m
    qT1D_percent_n    <- matrices_list_lever1$qT1D_percent_n

    qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
    qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
    smr_matrix               <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

    prev_100d_lever2023               <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx,smr_matrix ,  years) # deaths on diagnosis are converted to higer incidence


    matrices_list_lever2 <- Apply_levers_to_input_matrices (matrices_list, lever=2, lever_year_range=as.character(2023:max(years)))
    i                 <- matrices_list_lever2$i   # mean(i[as.character(2010:2019),19]) *100000
    dDx               <- matrices_list_lever2$dDx
    smr_matrix_n      <- matrices_list_lever2$smr_matrix_n
    smr_matrix_m      <- matrices_list_lever2$smr_matrix_m
    qT1D_percent_n    <- matrices_list_lever2$qT1D_percent_n

    qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
    qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
    smr_matrix               <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

    prev_100d_basic_care_lever2023    <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx,smr_matrix ,  years) # smr all non-minimal


    matrices_list_lever3 <- Apply_levers_to_input_matrices (matrices_list, lever=3, lever_year_range=as.character(2023:max(years)))
    # matrices_list_lever3 <- Apply_levers_to_input_matrices (matrices_list, lever=2.5, lever_year_range=as.character(2023:max(years)))
    i                 <- matrices_list_lever3$i   # mean(i[as.character(2010:2019),19]) *100000
    dDx               <- matrices_list_lever3$dDx
    smr_matrix_n      <- matrices_list_lever3$smr_matrix_n
    smr_matrix_m      <- matrices_list_lever3$smr_matrix_m
    qT1D_percent_n    <- matrices_list_lever3$qT1D_percent_n
    qT1D_n            <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
    qT1D_m            <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
    smr_matrix               <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

    prev_100d_best_care_lever2023     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx,smr_matrix ,  years) # deaths on diagnosis are converted to higer incidence

    matrices_list_lever4 <- Apply_levers_to_input_matrices (matrices_list, lever=4, lever_year_range=as.character(2023:max(years)))
    i                 <- matrices_list_lever4$i   # mean(i[as.character(2010:2019),19]) *100000
    dDx               <- matrices_list_lever4$dDx
    smr_matrix_n      <- matrices_list_lever4$smr_matrix_n
    smr_matrix_m      <- matrices_list_lever4$smr_matrix_m
    qT1D_percent_n    <- matrices_list_lever4$qT1D_percent_n
    qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
    qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
    smr_matrix               <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

    prev_100d_cure_lever2023          <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx,smr_matrix ,  years)


  }

  #  quantify the HLYs actually delivered by research?   1970 -----------------------------------------------------------------------------------------------------------------------------------
  if(full_run)
  {

    matrices_list_lever_1970 <- matrices_list
    i                 <- matrices_list_lever_1970$i   # mean(i[as.character(2010:2019),19]) *100000
    dDx               <- matrices_list_lever_1970$dDx
    smr_matrix_n      <- matrices_list_lever_1970$smr_matrix_n
    smr_matrix_m      <- matrices_list_lever_1970$smr_matrix_m
    # for(j in 1970:2040)
    # {
    #   smr_matrix_n[as.character(j),] <- smr_matrix_n[as.character( 1970),]
    #   smr_matrix_m[as.character(j),] <- smr_matrix_m[as.character( 1970),]
    # }
    # make minimal smr same as srm_min
    smr_min    <- exp(10.5 *  0.3545 - 1.5274)

    ratio        <- smr_min/(rowSums(smr_matrix_n[,as.character(20:50)])/length(20:50))
    ratio_matrix <- t(matrix(ratio, nrow=100, ncol=length(ratio), byrow=TRUE))
    ratio_matrix[ratio_matrix < 1] <- 1
    smr_matrix_adjusted <- smr_matrix_n * ratio_matrix
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1
    smr_matrix_n <- smr_matrix_adjusted

    ratio        <- smr_min/(rowSums(smr_matrix_m[,as.character(20:50)])/length(20:50))
    ratio_matrix <- t(matrix(ratio, nrow=100, ncol=length(ratio), byrow=TRUE))
    ratio_matrix[ratio_matrix < 1] <- 1
    smr_matrix_adjusted <- smr_matrix_m * ratio_matrix
    smr_matrix_adjusted[smr_matrix_adjusted<1] <- 1
    smr_matrix_m <- smr_matrix_adjusted

    qT1D_percent_n    <- matrices_list_lever_1970$qT1D_percent_n
    qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
    qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)


    prev_stagnant_lever1970               <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx,smr_matrix ,  years) # deaths on diagnosis are converted to higer incidence

  }

  # Lever 2  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Lever 2  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # prev lever 2: basic care Insulin, strips and education SMR  median( c(3.7, 4.4)) , 4.05--------------------------------------------------------------------

  matrices_list_lever2 <- Apply_levers_to_input_matrices (matrices_list, lever=2, lever_year_range=lever_year_range)
  i                 <- matrices_list_lever2$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list_lever2$dDx
  smr_matrix_n      <- matrices_list_lever2$smr_matrix_n
  smr_matrix_m      <- matrices_list_lever2$smr_matrix_m
  qT1D_percent_n    <- matrices_list_lever2$qT1D_percent_n

  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
  smr_matrix               <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

  prev_100d_basic_care    <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx,smr_matrix ,  years) # smr all non-minimal

  life_table_lever_2      <- Get_life_expectancy_t1d (qB, smr_matrix)

  # Lever 3  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Lever 3  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # # prev best care Pumps and FGMs/CGMs SMR median( c(2.2, 2.6)) , median 2.4--------------------------------------------------------------------
  matrices_list_lever3 <- Apply_levers_to_input_matrices (matrices_list, lever=3, lever_year_range=lever_year_range)
  i                 <- matrices_list_lever3$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list_lever3$dDx
  smr_matrix_n      <- matrices_list_lever3$smr_matrix_n
  smr_matrix_m      <- matrices_list_lever3$smr_matrix_m
  qT1D_percent_n    <- matrices_list_lever3$qT1D_percent_n
  qT1D_n            <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m            <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
  smr_matrix         <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

  prev_100d_best_care     <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx,smr_matrix  ,  years) # deaths on diagnosis are converted to higer incidence


  life_table_lever_3 <- Get_life_expectancy_t1d (qB, smr_matrix)


  # counterfactual: ghost pop is diff between prevalence and counterfactual where deaths on diagnosis are converted to incidence, and hba1c is fully controlled (4.3)
  # prev cure care --------------------------------------------------------------------
  # Lever 4  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Lever 4  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  matrices_list_lever4 <- Apply_levers_to_input_matrices (matrices_list, lever=4, lever_year_range=lever_year_range)
  i                 <- matrices_list_lever4$i   # mean(i[as.character(2010:2019),19]) *100000
  dDx               <- matrices_list_lever4$dDx
  smr_matrix_n      <- matrices_list_lever4$smr_matrix_n
  smr_matrix_m      <- matrices_list_lever4$smr_matrix_m
  qT1D_percent_n    <- matrices_list_lever4$qT1D_percent_n
  qT1D_n       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_n)
  qT1D_m       <- Get_qT1D_from_smr_matrix (qB, smr_matrix_m)
  smr_matrix         <- smr_matrix_m * (1-qT1D_percent_n)  +  smr_matrix_n  *    qT1D_percent_n

  prev_100d_cure          <- calculate_prevalence(i        , qB, qT1D_n, qT1D_m,   qT1D_percent_n, dDx ,smr_matrix,  years)


  life_table_lever_4 <- Get_life_expectancy_t1d (qB, smr_matrix)


  if(full_run)
  {
    # calculating years gained  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # calculating years gained  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
    matrices_list_strips_low  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 5.3,6.3,lever_year_range)
    matrices_list_strips_hig <- Apply_smr_to_input_matrices(matrices_list, lever=1, 3.1,3.7,lever_year_range)
    matrices_list_sensor_low  <- Apply_smr_to_input_matrices(matrices_list, lever=1, 4.0,4.7,lever_year_range)
    matrices_list_sensor_hig <- Apply_smr_to_input_matrices(matrices_list, lever=1, 2.2,2.8,lever_year_range)

    smr_matrix         <- matrices_list_strips_low$smr_matrix_m * (1- matrices_list_strips_low$qT1D_percent_n)  +  matrices_list_strips_low$smr_matrix_n  *     matrices_list_strips_low$qT1D_percent_n
    life_table_strips_low <- Get_life_expectancy_t1d (qB, smr_matrix)

    smr_matrix         <- matrices_list_strips_hig$smr_matrix_m * (1- matrices_list_strips_hig$qT1D_percent_n)  +  matrices_list_strips_hig$smr_matrix_n  *     matrices_list_strips_hig$qT1D_percent_n
    life_table_strips_hig <- Get_life_expectancy_t1d (qB, smr_matrix)

    smr_matrix         <- matrices_list_sensor_low$smr_matrix_m * (1- matrices_list_sensor_low$qT1D_percent_n)  +  matrices_list_sensor_low$smr_matrix_n  *     matrices_list_sensor_low$qT1D_percent_n
    life_table_sensor_low <- Get_life_expectancy_t1d (qB, smr_matrix)

    smr_matrix         <- matrices_list_sensor_hig$smr_matrix_m * (1- matrices_list_sensor_hig$qT1D_percent_n)  +  matrices_list_sensor_hig$smr_matrix_n  *     matrices_list_sensor_hig$qT1D_percent_n
    life_table_sensor_hig <- Get_life_expectancy_t1d (qB, smr_matrix)

  }

  life_table_final_all <- data.frame()

  if(FALSE)
  {
    if(full_run & country_wb_name == "New Zealand")
    {
      # calculate hab1cs from 5 to 12
      hba1c_list  <- c(5,5.5,6,6.5,7,7.5,8,8.5,9,10,12)
      # hba1c_matrix <- (log(smr) +  1.5274 )/ 0.3545
      smr_list    <- exp(hba1c_list *  0.3545 - 1.5274)

      # i <- 1
      matrices_list_        <- Apply_smr_to_input_matrices(matrices_list, lever=1, smr_list[1],smr_list[1],lever_year_range)
      smr_matrix            <- matrices_list_$smr_matrix_m * (1- matrices_list_$qT1D_percent_n)  +  matrices_list_$smr_matrix_n  *     matrices_list_$qT1D_percent_n
      life_table_final_ <- Get_life_expectancy_t1d (qB, smr_matrix,Enable_DKA_Hypo=TRUE)

      life_table_final_ <- dplyr::select(life_table_final_, year,age,ex,ex_complication,ex_treatment,events_DKA,events_Hypo)
      life_table_final_all <- life_table_final_[life_table_final_$age==10,]
      life_table_final_all$hba1c <- hba1c_list[1]

      for(i in 2:length(smr_list))
      {
        # i <- 2
        matrices_list_        <- Apply_smr_to_input_matrices(matrices_list, lever=1, smr_list[i],smr_list[i],lever_year_range)
        smr_matrix            <- matrices_list_$smr_matrix_m * (1- matrices_list_$qT1D_percent_n)  +  matrices_list_$smr_matrix_n  *     matrices_list_$qT1D_percent_n
        life_table_final_ <- Get_life_expectancy_t1d (qB, smr_matrix,Enable_DKA_Hypo=TRUE)
        life_table_final_ <- dplyr::select(life_table_final_, year,age,ex,ex_complication,ex_treatment,events_DKA,events_Hypo)

        life_table_final <- life_table_final_[life_table_final_$age==10,]
        life_table_final$hba1c <- hba1c_list[i]
        life_table_final_all <- rbind(life_table_final_all, life_table_final)

      }
      life_table_final_all <- setDF(life_table_final_all) %>% mutate_if(is.numeric, round, digits=4)

      # write.csv(life_table_final_all,"test.csv")
      life_table_final_all$loc_id <- loc_id

    }
  }


  # ----------------------------------------------------------------------------------------------------------------------------------------
  # translate proportions to population levels
  # P + S + D = pop_scale_factor
  # start simulation from 1800 ----------------------------------------------------

  pop_scale_factor       <- pop / (prev$P + prev$S)   # exclude the dead from population

  P_level                             <- pop_scale_factor * prev$P  #  sum(P_level["2021",]) ; sum( (pop_scale_factor_100d * prev_100d$P) ["2021",])
  PC_level                            <- pop_scale_factor * prev$C_P["AtLeast1C",,]  #  sum(P_level["2021",]) ; sum( (pop_scale_factor_100d * prev_100d$P) ["2021",])

  S_level         <- pop_scale_factor * prev$S
  I_flow          <- pop_scale_factor * prev$I

  I_flow_diagnosis           <- pop_scale_factor * prev_100d$I
  I_flow_basic_care          <- pop_scale_factor * prev_100d_basic_care$I
  I_flow_best_care           <- pop_scale_factor * prev_100d_best_care$I
  I_flow_cure                <- pop_scale_factor * prev_100d_cure$I

  DDx_flow                   <- pop_scale_factor * prev$DDx   # sum(DDx_flow["2021",])

  DT1D_flow                  <- pop_scale_factor * prev$DT1D
  DBGP_flow                  <- pop_scale_factor * prev$DBGP
  DT1D_flow_100d             <- pop_scale_factor * prev_100d$DT1D
  DT1D_flow_basic_care       <- pop_scale_factor * prev_100d_basic_care$DT1D
  DT1D_flow_best_care        <- pop_scale_factor * prev_100d_best_care$DT1D
  DT1D_flow_cure             <- pop_scale_factor * prev_100d_cure$DT1D

  BD_flow         <- pop * qB

  # Calculate missing prevalence  and delta of levers. --------------------------------------------------------------------------------
  # Calculate missing prevalence  and delta of levers. --------------------------------------------------------------------------------
  ghost_ddx_level   <- (prev_100d$P              - prev$P) * pop_scale_factor #
  ghost_hba1c_level <- (prev_100d_cure$P         - prev_100d$P) * pop_scale_factor #
  ghost_basic_care  <- (prev_100d_basic_care$P   - prev_100d$P) * pop_scale_factor #  delta
  ghost_best_care   <- (prev_100d_best_care$P    - prev_100d_basic_care$P) * pop_scale_factor #
  ghost_cure        <- (prev_100d_cure$P         - prev_100d_best_care$P) * pop_scale_factor #  sum(ghost_hba1c_level["2021",])
  ghost_stage1        <- (prev_stage1$P         - prev_stage2$P) * pop_scale_factor #  sum(ghost_hba1c_level["2021",])
  ghost_stage2        <- (prev_stage2$P         - prev_100d_cure$P) * pop_scale_factor #  sum(ghost_hba1c_level["2021",])
  # ghost population is the sum of these two
  ghost_level <- ghost_hba1c_level + ghost_ddx_level #   sum(ghost_level["2021",])

  print(paste0("ghost onset death 2021: ",round(sum(ghost_ddx_level["2021",]),2),"; ghost early death 2021: ",round(sum(ghost_hba1c_level["2021",]),2) ) )


  # Calculate missing prevalence  and delta of levers2023. --------------------------------------------------------------------------------
  # Calculate missing prevalence  and delta of levers2023. --------------------------------------------------------------------------------
  ghost_ddx_level_lever2023   <- (prev_100d_lever2023$P              - prev$P) * pop_scale_factor #
  ghost_hba1c_level_lever2023 <- (prev_100d_cure_lever2023$P         - prev_100d_lever2023$P) * pop_scale_factor #
  ghost_basic_care_lever2023  <- (prev_100d_basic_care_lever2023$P   - prev_100d_lever2023$P) * pop_scale_factor #  delta
  ghost_best_care_lever2023   <- (prev_100d_best_care_lever2023$P    - prev_100d_basic_care_lever2023$P) * pop_scale_factor #
  ghost_cure_lever2023        <- (prev_100d_cure_lever2023$P         - prev_100d_best_care_lever2023$P) * pop_scale_factor #  sum(ghost_hba1c_level["2021",])
  # ghost population is the sum of these two
  ghost_level_lever2023 <- ghost_hba1c_level_lever2023 + ghost_ddx_level_lever2023 #   sum(ghost_level["2021",])

  print(paste0("ghost onset death: ",round(sum(ghost_ddx_level_lever2023["2021",]),2),"; ghost early death: ",round(sum(ghost_hba1c_level_lever2023["2021",]),2) ) )



  if(FALSE)  # Plot checking
  {
    df <- data.frame(x=rownames(ghost_level)
                     ,diagnosed           = (rowSums(P_level))
                     ,ghost_undiagnosed   = (rowSums(ghost_ddx_level))
                     ,ghost_diagnosed     = (rowSums(ghost_hba1c_level))
                     ,ghost_basic_care    = (rowSums(ghost_basic_care))
                     ,ghost_best_care     = (rowSums(ghost_best_care))
                     ,ghost_cure          = (rowSums(ghost_cure))
                     ,ghost_stage1          = (rowSums(ghost_stage1))
                     ,ghost_stage2          = (rowSums(ghost_stage2))
    )
    df <- setDF(df) %>% mutate_if(is.numeric, round, digits=2)

    df[df$x>=1960 & df$x<=2040,] %>%
      e_charts(x) %>%
      e_bar(diagnosed, stack = "grp",name="Prevalence") %>%
      e_bar(ghost_stage1, stack = "grp",name="Prevalence (Stage1)") %>%
      e_bar(ghost_stage2, stack = "grp",name="Prevalence (Stage2)") %>%
      e_bar(ghost_undiagnosed, stack = "grp",name="Missing Prevalence (Undiagnosis)") %>%
      e_bar(ghost_basic_care, stack = "grp",name="Missing Prevalence (Basic Care)") %>%
      e_bar(ghost_best_care, stack = "grp",name="Missing Prevalence (Best Care)") %>%
      e_bar(ghost_cure, stack = "grp",name="Missing Prevalence (Cure)") %>%
      e_title(paste0("T1D population - ",country_wb_name) )%>%
      e_legend(left=80,top=50,orient='vertical')%>%      e_tooltip(trigger = 'axis')

  }

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
  prev_merge_new <-       as.data.frame.table(pop,stringsAsFactors = F)[,1:2]
  colnames(prev_merge_new) <- c("Year","Age")
  prev_merge_new      <- cbind(Country= country_wb_name,loc_id=loc_id,prev_merge_new)

  prev_merge_new$"dDx_full"                              <- as.data.frame.table(matrices_list$dDx)$Freq
  prev_merge_new$"qT1D_percent_n_full"                   <- as.data.frame.table(matrices_list$qT1D_percent_n)$Freq

  prev_merge_new$"sim_start_year"                        <- 0
  prev_merge_new$"sim_min_diag_rates"                    <- 0
  prev_merge_new$"sim_min_non_minimal_care_perc"         <- 0
  prev_merge_new$"sim_min_non_minimal_care_level"        <- 0


  prev_merge_new$"Ann. background population"           <- as.data.frame.table(pop)$Freq
  prev_merge_new$"Ann. background mortality"            <- as.data.frame.table(BD_flow)$Freq
  prev_merge_new$"Prevalence"                           <- as.data.frame.table(P_level)$Freq
  prev_merge_new$"Prevalence with AtLeast1C"            <- as.data.frame.table(PC_level)$Freq
  prev_merge_new$"Incidence (1 base)"                   <- as.data.frame.table(I_flow)$Freq
  prev_merge_new$"Incidence (2 diagnosis)"              <- as.data.frame.table(I_flow_diagnosis)$Freq

  prev_merge_new$"Ghosts"                                        <- as.data.frame.table(ghost_level)$Freq
  prev_merge_new$"Ghosts (onset death)"                          <- as.data.frame.table(ghost_ddx_level)$Freq
  prev_merge_new$"Ghosts (early death)"                          <- as.data.frame.table(ghost_hba1c_level)$Freq
  prev_merge_new$"Ghosts (delta basic care)"                     <- as.data.frame.table(ghost_basic_care)$Freq
  prev_merge_new$"Ghosts (delta best care)"                      <- as.data.frame.table(ghost_best_care)$Freq
  prev_merge_new$"Ghosts (delta cure)"                           <- as.data.frame.table(ghost_cure)$Freq
  prev_merge_new$"Ann. onset deaths"                             <- as.data.frame.table(DDx_flow)$Freq
  prev_merge_new$"Ann. early deaths"                             <- as.data.frame.table(DT1D_flow)$Freq

  prev_merge_new$"Ann. early deaths (background)"                               <- as.data.frame.table(DBGP_flow)$Freq
  prev_merge_new$"Ann. early deaths (2 diagnosis)"                              <- as.data.frame.table(DT1D_flow_100d)$Freq
  prev_merge_new$"Ann. early deaths (3 basic care)"                             <- as.data.frame.table(DT1D_flow_basic_care)$Freq
  prev_merge_new$"Ann. early deaths (4 best care)"                              <- as.data.frame.table(DT1D_flow_best_care)$Freq
  prev_merge_new$"Ann. early deaths (5 cure)"                                   <- as.data.frame.table(DT1D_flow_cure)$Freq

  prev_merge_new$"Ghosts lever2023"                                                      <- as.data.frame.table(ghost_level_lever2023)$Freq
  prev_merge_new$"Ghosts (onset death) lever2023"                                        <- as.data.frame.table(ghost_ddx_level_lever2023)$Freq
  prev_merge_new$"Ghosts (early death) lever2023"                                        <- as.data.frame.table(ghost_hba1c_level_lever2023)$Freq
  prev_merge_new$"Ghosts (delta basic care) lever2023"                                   <- as.data.frame.table(ghost_basic_care_lever2023)$Freq
  prev_merge_new$"Ghosts (delta best care) lever2023"                                    <- as.data.frame.table(ghost_best_care_lever2023)$Freq
  prev_merge_new$"Ghosts (delta cure) lever2023"                                         <- as.data.frame.table(ghost_cure_lever2023)$Freq


  prev_merge_new$"Prevalence Stage 1 only"                                        <- as.data.frame.table(pop_scale_factor * (prev_stage1$P -  prev_stage2$P))$Freq
  prev_merge_new$"Prevalence Stage 2 only"                                        <- as.data.frame.table(pop_scale_factor * (prev_stage2$P -  prev_100d_cure$P))$Freq



  if(full_run)
  {
    prev_merge_new$"Prevalence stagnant 1970"                                        <- as.data.frame.table(pop_scale_factor * prev_stagnant_lever1970$P)$Freq
    prev_merge_new$"Prevalence delay onset by 1 years"                                        <- as.data.frame.table(pop_scale_factor * prev_delay_year_1$P)$Freq
    prev_merge_new$"Prevalence delay onset by 3 years"                                        <- as.data.frame.table(pop_scale_factor * prev_delay_year_3$P)$Freq
    prev_merge_new$"Prevalence delay onset by 5 years"                                        <- as.data.frame.table(pop_scale_factor * prev_delay_year_5$P)$Freq
    prev_merge_new$"Prevalence delay onset by 8 years"                                        <- as.data.frame.table(pop_scale_factor * prev_delay_year_8$P)$Freq
    prev_merge_new$"Prevalence delay onset by 13 years"                                        <- as.data.frame.table(pop_scale_factor * prev_delay_year_13$P)$Freq

  }

  # --------------------------
  prev_merge_new$Year <- as.numeric(prev_merge_new$Year)
  prev_merge_new$Age  <- as.numeric(prev_merge_new$Age)
  prev_merge_new <- prev_merge_new[prev_merge_new$Year>=1960,]
  prev_merge_new <- prev_merge_new[ with(prev_merge_new, order(  Year,Age)), ]   # order properly for transforming to matrix


  prev_merge_wide          <- prev_merge_new
  prev_merge_wide$loc_id   <- as.factor(prev_merge_wide$loc_id)
  prev_merge_wide$`Life expectency (1 background)`                         <- life_table_background$ex[life_table_background$year>=1960]

  # for t1d base life expectancy, include diagnosis rate effect in the begining
  prev_merge_wide$`Life expectency (2 t1d base)`                           <- life_table_base_scenario$ex[life_table_base_scenario$year>=1960]
  diagnosis_rate <- (prev_merge_wide$`Incidence (1 base)`+0.001) / (prev_merge_wide$`Ann. onset deaths`  + prev_merge_wide$`Incidence (1 base)`+0.001)
  prev_merge_wide$`Life expectency (2 t1d base)`                           <- prev_merge_wide$`Life expectency (2 t1d base)` * diagnosis_rate +  (1-diagnosis_rate)*0.5
  prev_merge_wide$`Lifetime years lost (2 t1d base) (complication)`        <- life_table_base_scenario$ex_complication[life_table_base_scenario$year>=1960] *  prev_merge_wide$`Life expectency (2 t1d base)`/ life_table_base_scenario$ex[life_table_base_scenario$year>=1960]
  prev_merge_wide$`Lifetime years lost (2 t1d base) (treatment)`           <- life_table_base_scenario$ex_treatment[life_table_base_scenario$year>=1960] *  prev_merge_wide$`Life expectency (2 t1d base)`/ life_table_base_scenario$ex[life_table_base_scenario$year>=1960]

  prev_merge_wide$`Life expectency (3 t1d diagnosis)`                      <- life_table_base_scenario$ex[life_table_base_scenario$year>=1960]
  prev_merge_wide$`Lifetime years lost (3 t1d diagnosis) (complication)`   <- life_table_base_scenario$ex_complication[life_table_base_scenario$year>=1960]
  prev_merge_wide$`Lifetime years lost (3 t1d diagnosis) (treatment)`      <- life_table_base_scenario$ex_treatment[life_table_base_scenario$year>=1960]

  prev_merge_wide$`Life expectency (4 t1d basic care)`                     <- life_table_lever_2$ex[life_table_lever_2$year>=1960]
  prev_merge_wide$`Lifetime years lost (4 t1d basic care) (complication)`  <- life_table_lever_2$ex_complication[life_table_lever_2$year>=1960]
  prev_merge_wide$`Lifetime years lost (4 t1d basic care) (treatment)`     <- life_table_lever_2$ex_treatment[life_table_lever_2$year>=1960]

  prev_merge_wide$`Life expectency (5 t1d best care)`                      <- life_table_lever_3$ex[life_table_lever_3$year>=1960]
  prev_merge_wide$`Lifetime years lost (5 t1d best care) (complication)`   <- life_table_lever_3$ex_complication[life_table_lever_3$year>=1960]
  prev_merge_wide$`Lifetime years lost (5 t1d best care) (treatment)`      <- life_table_lever_3$ex_treatment[life_table_lever_3$year>=1960]

  prev_merge_wide$`Life expectency (6 t1d cure)`                           <- life_table_lever_4$ex[life_table_lever_4$year>=1960]
  prev_merge_wide$`Lifetime years lost (6 t1d cure) (complication)`        <- life_table_lever_4$ex_complication[life_table_lever_4$year>=1960]
  prev_merge_wide$`Lifetime years lost (6 t1d cure) (treatment)`           <- life_table_lever_4$ex_treatment[life_table_lever_4$year>=1960]


  if(full_run)
  {

    prev_merge_wide$`Lifetime years lost (delay onset year 1) (treatment)`     <- life_table_delay_year_1$ex_treatment[life_table_delay_year_1$year>=1960]
    prev_merge_wide$`Lifetime years lost (delay onset year 1) (complication)`  <- life_table_delay_year_1$ex_complication[life_table_delay_year_1$year>=1960]
    prev_merge_wide$`Life expectency (delay onset year 1)`                     <- life_table_delay_year_1$ex[life_table_delay_year_1$year>=1960]

    prev_merge_wide$`Lifetime years lost (delay onset year 3) (treatment)`     <- life_table_delay_year_3$ex_treatment[life_table_delay_year_3$year>=1960]
    prev_merge_wide$`Lifetime years lost (delay onset year 3) (complication)`  <- life_table_delay_year_3$ex_complication[life_table_delay_year_3$year>=1960]
    prev_merge_wide$`Life expectency (delay onset year 3)`                     <- life_table_delay_year_3$ex[life_table_delay_year_3$year>=1960]

    prev_merge_wide$`Lifetime years lost (delay onset year 5) (treatment)`     <- life_table_delay_year_5$ex_treatment[life_table_delay_year_5$year>=1960]
    prev_merge_wide$`Lifetime years lost (delay onset year 5) (complication)`  <- life_table_delay_year_5$ex_complication[life_table_delay_year_5$year>=1960]
    prev_merge_wide$`Life expectency (delay onset year 5)`                     <- life_table_delay_year_5$ex[life_table_delay_year_5$year>=1960]

    prev_merge_wide$`Lifetime years lost (delay onset year 8) (treatment)`     <- life_table_delay_year_8$ex_treatment[life_table_delay_year_8$year>=1960]
    prev_merge_wide$`Lifetime years lost (delay onset year 8) (complication)`  <- life_table_delay_year_8$ex_complication[life_table_delay_year_8$year>=1960]
    prev_merge_wide$`Life expectency (delay onset year 8)`                     <- life_table_delay_year_8$ex[life_table_delay_year_8$year>=1960]

    prev_merge_wide$`Lifetime years lost (delay onset year 13) (treatment)`     <- life_table_delay_year_13$ex_treatment[life_table_delay_year_13$year>=1960]
    prev_merge_wide$`Lifetime years lost (delay onset year 13) (complication)`  <- life_table_delay_year_13$ex_complication[life_table_delay_year_13$year>=1960]
    prev_merge_wide$`Life expectency (delay onset year 13)`                     <- life_table_delay_year_13$ex[life_table_delay_year_13$year>=1960]

    # for site data.
    prev_merge_wide$`Life expectency (strip low)`            <- life_table_strips_low$ex[life_table_strips_low$year>=1960]
    prev_merge_wide$`Lifetime years lost (strip low)`        <- (life_table_strips_low$ex_complication + life_table_strips_low$ex_treatment)  [life_table_strips_low$year>=1960]

    prev_merge_wide$`Life expectency (strip hig)`            <- life_table_strips_hig$ex[life_table_strips_hig$year>=1960]
    prev_merge_wide$`Lifetime years lost (strip hig)`        <- (life_table_strips_hig$ex_complication + life_table_strips_hig$ex_treatment)  [life_table_strips_hig$year>=1960]

    prev_merge_wide$`Life expectency (sensor low)`            <- life_table_sensor_low$ex[life_table_sensor_low$year>=1960]
    prev_merge_wide$`Lifetime years lost (sensor low)`        <- (life_table_sensor_low$ex_complication + life_table_sensor_low$ex_treatment)  [life_table_sensor_low$year>=1960]

    prev_merge_wide$`Life expectency (sensor hig)`            <- life_table_sensor_hig$ex[life_table_sensor_hig$year>=1960]
    prev_merge_wide$`Lifetime years lost (sensor hig)`        <- (life_table_sensor_hig$ex_complication + life_table_sensor_hig$ex_treatment)  [life_table_sensor_hig$year>=1960]
  }

  prev_merge_wide$`% Odds living to`  <- prev_merge_wide$Prevalence / (prev_merge_wide$Prevalence + prev_merge_wide$Ghosts) * 100
  prev_merge_wide$`% Odds living to`[is.na(prev_merge_wide$`% Odds living to`)] <- 0

  if(TRUE)
  {
    # 1 in ** families. country_indicator------------------------------------------------------------------------------------------------------------------------------------
    # country_indicator <- dplyr::select(readRDS("data_wb/country_indicator_imputed_avg_0.4.13.Rds"),country,year,fertility)
    country_indicator <- dplyr::select(readRDS("data_wb/country_indicator_imputed_avg_fertility_0.4.13.Rds"),country,year,fertility)

    days_              <- setDT(prev_merge_wide)[,list(background_population=sum(`Ann. background population`),`Total`= sum(Ghosts+Prevalence)),by=c("Country","Year")]
    days_              <- days_ %>% dplyr::left_join(dplyr::select(country_indicator,Country=country,Year=year,fertility),by=c("Country","Year"))
    days_$Country      <- as.factor(days_$Country)
    days_$family_size  <- days_$fertility^2 + days_$fertility +2

    days_$`1 in x families` <- 1/(1-((days_$background_population/days_$`Total`-1)/(days_$background_population/days_$`Total`) ) ^ days_$family_size)

    days_$`1 in x families`[is.nan(days_$`1 in x families`)] <- 0

    prev_merge_wide   <- prev_merge_wide %>% dplyr::inner_join(dplyr::select(days_,Country,Year,`1 in x families` = `1 in x families`
                                                                             # ,`Lifetime years lost (2 t1d base)`=lifetime_years_lost_complication_and_treatment
                                                                             # ,`Lifetime years lost (2 t1d base) (complication)`=lifetime_years_lost_complication
                                                                             # ,`Lifetime years lost (2 t1d base) (treatment)`   =lifetime_years_lost_treatment
    ), by=c("Country","Year"))
  }

  prev_merge_wide <- setDF(prev_merge_wide) %>% mutate_if(is.numeric, round, digits=2)

  # if(FALSE)
  # {
  #   prev_merge_wide_old  <- read_parquet( "C:/DropboxT1D/reruns_scenario_2/rerun_0.4.15_0/AFG.binary")
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
      prev_merge_wide=prev_merge_wide,
      life_table_final_all=life_table_final_all
    )),
    class='prevalence')
}
