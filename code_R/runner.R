
refresh_one_country_file <- function(
    partition_id,data_dir,cache_dir, start_year=1960, end_year=2040,use_3d_array=TRUE,run_type="national"
) {  # start_year=1960; end_year=2040 ; run_sim=TRUE ; partition_id <- "partition_1"
  source('code_R/prevalence_utils.R')
  source('code_R/001_main_util.R')

  # Load Library ---------------------------------------------------------------------------------------
  library(arrow)
  # library(data.table)
  library(dplyr)
  library(abind)
  library(stringr)
  library(parallel)

  # Load Dependency ---------------------------------------------------------------------------------------

  source("code_R/refresh_country_files_tensor.R")
  source("code_R/refresh_country_files_tensor_utils.R")
  source("code_R/utils.R")
  source("code_R/data.R")
  source("code_R/runner.R")


  # country_code <- countries$world_bank_code[countries$world_bank_name==country_wb_name]
  # cache_file   <- file.path(cache_dir, paste0(country_code, ".binary"))
  # log_timing <- grepl("mcmc_1_0\\/AFG\\.binary",cache_file )


  # country_name <- country_name_list[1]
  # library(RSQLite)
  # library(dplyr)
  # library(tidyr)
  # library(arrow)
  Log=function(fmt, ...) { cat(sprintf(paste0(fmt, '\n'), ...)) }

  Log('Starting model for %s...', partition_id)

  tryCatch({

   # run_type="sim"
   # run_type="sim_baseline"
   # run_type="national"

    cache_file                <- ( paste0(cache_dir,"/",partition_id))
    cache_file_lifetime       <- ( paste0(cache_dir,"/",partition_id, "_lifetime"))
    # sink("log.txt",append=TRUE);cat(paste0( "test shit" , "  ", " \n") );sink()
    prev_merge_wide      <- data.frame()
    life_table_final     <- data.frame()

    data_long                   <- setDF(arrow::read_feather(paste0(data_dir,"/",partition_id)))

    if(run_type=="sim")
    {
      # simulator 2500 runs
      for(sim_start_year in c(2000,2010,2020,2024,2025,2030) ) #   important to have 2023 2024 2025
        # for(sim_start_year in c(2000,2010,2020,2030) ) #   important to have 2023 2024 2025
        {
          prev_merge_list <- list()
          for(sim_min_diag_rates in c(0,20,40,60,80,100) )  #  not meaningful to have negatives.  we assume this set minimal diangosis rate to diag_rates
            # for(sim_min_diag_rates in c(0,30,100) )  #  not meaningful to have negatives.  we assume this set minimal diangosis rate to diag_rates
            {
              for(sim_min_non_minimal_care_perc in c(0,20,40,60,80,100) )  #  not meaningful to have negatives.  we assume this set minimal diangosis rate to diag_rates
                # for(sim_min_non_minimal_care_perc in c(0,30,100) )  #  not meaningful to have negatives.  we assume this set minimal diangosis rate to diag_rates
                {
                for(sim_min_non_minimal_care_level in c("default_care","basic_care","best_care", "cure") )  #  this is more meaning full than percentage.
                  # for(sim_min_non_minimal_care_level in c("default_care","minimal_care","basic_care","best_care", "cure") )  #  this is more meaning full than percentage.
                    # for(sim_min_non_minimal_care_level in c("base_care","basic_care") )  #  this is more meaning full than percentage.
                    {
              #  sim_start_year=2000;sim_min_diag_rates=0;sim_min_non_minimal_care_perc=60;sim_min_non_minimal_care_level="default_care"
                # if( !(sim_diag_rates_min==0 & sim_care_level_min== "default_care") )
                # {
                  prev_merge        <- data.frame()

                  sim_parameters <- list(  sim_start_year   = sim_start_year
                                          ,sim_min_diag_rates = sim_min_diag_rates
                                          ,sim_min_non_minimal_care_perc = sim_min_non_minimal_care_perc
                                          ,sim_min_non_minimal_care_level = sim_min_non_minimal_care_level
                  )
                  print(unlist(sim_parameters))
                  system.time( prev_merge <- refresh_country_files_tensor( data_long
                                                                          ,sim_enable=TRUE
                                                                          ,sim_parameters=sim_parameters
                                                                          ,switches =list(run_site_data=FALSE ,run_delay_onset=FALSE,run_lever2023=FALSE,run_stages_1_2=FALSE)
                                                                          ) )

                  prev_merge_list <- append(prev_merge_list,list(prev_merge))

                # }

                    }

              }


          }
          length(prev_merge_list)
          prev_merge_list <- rbindlist(prev_merge_list)
          prev_merge_list$diagnosis_input  <- "gregory"
          country_data_merge               <- merge_agg(country_data_merge = prev_merge_list,all_data_points=TRUE,all_age_brackets=FALSE)
          country_data_merge_lifetime      <- merge_agg_lifetime(country_data_merge_age_10= prev_merge_list)


          country_data_merge          <- country_data_merge[country_data_merge$Year >=2000,]
          country_data_merge_lifetime <- country_data_merge_lifetime[country_data_merge_lifetime$Year >=2000,]

          arrow::write_feather(country_data_merge         , paste0(cache_file,"_",sim_start_year,".binary") )
          arrow::write_feather(country_data_merge_lifetime, paste0(cache_file_lifetime,"_",sim_start_year,".binary"))


        }




    }
    if(run_type=="sim_baseline")
    {

      system.time( prev_merge_list <- refresh_country_files_tensor( data_long
                                                               ,sim_enable=FALSE
                                                               ,switches =list(run_site_data=TRUE ,run_delay_onset=FALSE,run_lever2023=FALSE,run_stages_1_2=FALSE)
      ) )
      # length(prev_merge_list)
      # prev_merge_list <- rbindlist(prev_merge_list)
      prev_merge_list$diagnosis_input  <- "gregory"
      country_data_merge               <- merge_agg(country_data_merge = prev_merge_list,all_data_points=TRUE,all_age_brackets=FALSE)
      country_data_merge_lifetime      <- merge_agg_lifetime(country_data_merge_age_10= prev_merge_list)

      arrow::write_feather(country_data_merge         , cache_file)
      arrow::write_feather(country_data_merge_lifetime, cache_file_lifetime)

    }

    if(run_type=="national")
    {
      if(use_3d_array)
      {
        data_long$loc_id    <- as.factor( as.character(data_long$loc_id))
        loc_id_list         <- unique(data_long$loc_id)

        print( paste0("debug: nrow(data_long) ",nrow(data_long)) )
        # which(loc_id_list=="484_fema_urban_15")

        # data_long   <- data_long[data_long$loc_id %in% loc_id_list[1:20],] ;data_long$loc_id    <- as.factor( as.character(data_long$loc_id));  loc_id_list         <- unique(data_long$loc_id)   # for quick testing

        # # # not working due to torch not working under parellel cpus.

        # system.time({
        prev_merge_1        <- refresh_country_files_tensor(data_long,sim_enable=FALSE)

        prev_merge_1$diagnosis_input  <- "gregory"
        prev_merge_wide               <- prev_merge_1
      }else
      {
        # print(i)
        loc_id_list         <- unique(data_long$loc_id)
        for(i in 1:length(loc_id_list))
        {
          loc_id       <- loc_id_list[i]
          # loc_id <- '100_fema_rural'
          data_long_   <- data_long[data_long$loc_id==loc_id,]
          # write.csv(data_long_,"temp/Input(India).csv",row.names = FALSE)

          system.time({
            prev1 <- prevalence_and_ghost_pop(
              # loc_id=loc_id_,
              loc_id=loc_id,
              # hba1c=hba1c_function(country_name),
              start_year=start_year,
              end_year=end_year,
              data_dir=data_dir,
              log_timing=FALSE,
              data_long=data_long_,
              # full_run=FALSE,
              diagnosis_input = 1
            )
          })

          if(FALSE)
          { # check consistency between 3D and 2D ------ ------ ------ ------ ------ ------ ------ ------ ------
            prev1$prev_merge_wide$diagnosis_input  <- "gregory"
            prev_merge_ <- prev_merge_1[prev_merge_1$loc_id == loc_id,]
            sum_check <- colSums(  prev1$prev_merge_wide[,colnames(prev_merge_)] == prev_merge_)
            sum_check[sum_check!=nrow(prev_merge_)]

            # asdf1<- prev_merge_           [prev_merge_$`Prevalence Stage 2 only`!= prev1$prev_merge_wide$`Prevalence Stage 2 only`,]
            # asdf2<- prev1$prev_merge_wide [prev_merge_$`Prevalence Stage 2 only`!= prev1$prev_merge_wide$`Prevalence Stage 2 only`,]
            # asdf1$`Prevalence Stage 2 only` - asdf2$`Prevalence Stage 2 only`
            colnames(prev1$prev_merge_wide)[!colnames(prev1$prev_merge_wide) %in% colnames(prev_merge_)]
          }

          prev1$prev_merge_wide$diagnosis_input  <- "gregory"
          if(nrow(prev1$life_table_final)){ prev1$life_table_final$diagnosis_input <- "gregory"}
          prev_merge_wide   <- rbind(prev_merge_wide , prev1$prev_merge_wide)
          life_table_final  <- rbind(life_table_final, prev1$life_table_final)
        }

      }

      country_data_merge               <- merge_agg(country_data_merge = prev_merge_wide,all_data_points=TRUE,all_age_brackets=FALSE)
      country_data_merge_lifetime      <- merge_agg_lifetime(prev_merge_wide)

      arrow::write_feather(prev_merge_wide            , paste0(cache_file,"_raw") )
      arrow::write_feather(country_data_merge         , cache_file)
      arrow::write_feather(country_data_merge_lifetime, cache_file_lifetime)

    }
    if(run_type=="ward")
    {
      data_long$loc_id    <- as.factor( as.character(data_long$loc_id))
      loc_id_list         <- unique(data_long$loc_id)

      print( paste0("debug: nrow(data_long) ",nrow(data_long)) )
      # which(loc_id_list=="484_fema_urban_15")

      # # # not working due to torch not working under parellel cpus.

      # system.time({

      data_long$mortality_undiagnosed_rate <- data_long$mortality_undiagnosed_rate_ward
      prev_merge_1        <- refresh_country_files_tensor(data_long,sim_enable=FALSE)
      prev_merge_1$diagnosis_input  <- "ward"

      prev_merge_wide               <- prev_merge_1

      country_data_merge               <- merge_agg(country_data_merge = prev_merge_wide,all_data_points=TRUE,all_age_brackets=FALSE)
      country_data_merge_lifetime      <- merge_agg_lifetime(prev_merge_wide)

      arrow::write_feather(prev_merge_wide            , paste0(cache_file,"_raw") )
      arrow::write_feather(country_data_merge         , cache_file)
      arrow::write_feather(country_data_merge_lifetime, cache_file_lifetime)
    }
    if(run_type=="blend")
    {
      data_long$loc_id    <- as.factor( as.character(data_long$loc_id))
      loc_id_list         <- unique(data_long$loc_id)

      print( paste0("debug: nrow(data_long) ",nrow(data_long)) )
      # which(loc_id_list=="484_fema_urban_15")

      # # # not working due to torch not working under parellel cpus.

      # system.time({

      data_long$mortality_undiagnosed_rate <- data_long$mortality_undiagnosed_rate_average
      prev_merge_1        <- refresh_country_files_tensor(data_long,sim_enable=FALSE)
      prev_merge_1$diagnosis_input  <- "blend"

      prev_merge_wide               <- prev_merge_1

      country_data_merge               <- merge_agg(country_data_merge = prev_merge_wide,all_data_points=TRUE,all_age_brackets=FALSE)
      country_data_merge_lifetime      <- merge_agg_lifetime(prev_merge_wide)

      arrow::write_feather(prev_merge_wide            , paste0(cache_file,"_raw") )
      arrow::write_feather(country_data_merge         , cache_file)
      arrow::write_feather(country_data_merge_lifetime, cache_file_lifetime)
    }



  },
  error = function(cond) {
    Log('Error running model for %s:', partition_id)
    sink("log.txt",append=TRUE);cat(paste0( partition_id , "  ",  cond, " \n") );sink()
    NA
  },
  warning = function(cond) {
    Log('Warning running model for %s:', paste0(partition_id," : ",cond) )
  },
  finally = {
  })

}


generate_3D_matrix <- function(data_output, title_string)
{
  data_output <- data_output[as.numeric(data_output$Year) %%2 ==0,]
  data_output <- data_output[as.numeric(data_output$Age ) %%5 ==0,]
  data_output %>%
    e_charts(Age) %>%
    e_scatter_3d(Year, Value) %>%
    e_x_axis_3d(axisLine = list(),name="Age")%>%
    e_y_axis_3d(axisLine = list(),name="Year",min=min(data_output$Year))%>%
    e_z_axis_3d(axisLine = list(),name=title_string)%>%
    e_title(title_string)
}

