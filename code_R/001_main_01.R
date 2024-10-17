#!/usr/bin/env Rscript
# setwd( "E:/GithubCode/t1d_global_model/t1dGlobalModel")
main_ <- function(version_no,run_per_country=TRUE,run_merge_country=TRUE,run_few_countries= c() ,scenario=2)
{
  # version_no <- "0.4.12_mcmc_1"
  # version_no_input <- "0.4.15"; version_no <- "0.4.15"; scenario <- 2 ; run_per_country=TRUE; run_merge_country=TRUE; run_few_countries= c()
  # version_no_input <- "1.0.1.wpp2022"; version_no <- "1.0.1.wpp2022"; scenario <- 2 ; run_per_country=TRUE; run_merge_country=TRUE; run_few_countries= c()
  # version_no_input <- "1.0.1.wpp2022.gender"; version_no <- "1.0.1.wpp2022.gender"; scenario <- 2 ; run_per_country=TRUE; run_merge_country=TRUE; run_few_countries= c()
  host_name <<- "localhost"
  # library(RSQLite)
  library(RPostgreSQL)
  library(dplyr)
  library(purrr)
  library(testthat)
  library(tidyverse)
  library(doParallel)
  library(futile.logger)
  library(data.table)
  library("writexl")
  library(zip)
  library(fst)
  library(arrow)
  library(echarts4r)
  library(abind)

  source('code_R_data_prep/DATA_CONSTANTS.R')
  source('code_shiny/QUERY_DATA_POINTS.R')
  source('code_R/runner.R')
  source('code_R/prevalence.R')
  source('code_R/prevalence_utils.R')
  source('code_R/001_main_util.R')
  source('code_R/data.R')
  source('code_R/utils.R')
  source('code_R/complications.R')
  source('code_R/burden.R')
  # source('code_quick_job_scripts/extract_for_purpose.R')
  source("code_quick_job_scripts/extract_for_purpose_v2.R")

  source("code_R/refresh_country_files_tensor.R")
  source("code_R/refresh_country_files_tensor_utils.R")

  time_start <- Sys.time()
  if (!dir.exists('cache')) {
    dir.create('cache')
  }

  source('code_R/000_parameters.R')

  # version_no <- "0.4.12"; scenario <- 2; config$run_days_lost <- TRUE;  config$run_days_lost_lever <- TRUE

  if(scenario==1)
  {
    config$run_projection      <- FALSE
  }else
  {
    config$run_projection      <- TRUE
  }
  config$version_no_input <- version_no_input
  print(config)
  # config$run_projection      <- FALSE
  # version_no <- "0.4.13"
  # for(i in 1:15)  {
  for(i in 1) {
    i <- 0  # default run
    # i <- 14  # diagnosis_rate_left run
    # i <- 16  # diagnosis_rate_right run
    print(i)
    scenarios[] <- FALSE
    ifelse(i==0,scenarios[i] <- FALSE,scenarios[i] <- TRUE )
    # scenarios[] <- FALSE , i <- cha
    # Sys.setenv(T1D_DATA_FILE = paste0("data_outputs/data.db") )
    data_dir  <- paste0( "../../data_outputs/data_",version_no_input)
    # cache_dir <- paste0('reruns/rerun_',version_no)

    if(scenarios$adult_incidence_all_studies)
    { data_dir  <- paste0( "../../data_outputs/data_",version_no,".adult_incidence_all_studies")   }
    if(scenarios$iot_no_growth)
    {  data_dir  <- paste0( "../../data_outputs/data_",version_no,".iot_no_growth")}
    if(scenarios$iot_global_curve_for_all)
    { data_dir  <- paste0( "../../data_outputs/data_",version_no,".iot_global_curve_for_all")}
    if(scenarios$smr_age_cruve_flat_average)
    { data_dir  <- paste0( "../../data_outputs/data_",version_no,".smr_age_cruve_flat_average")}
    if(scenarios$adult_onset_zero)
    { data_dir  <- paste0( "../../data_outputs/data_",version_no,".adult_onset_zero")}

    # Create SQLite database  data.db ------------------------------------------------------------------------------
    # unlink(cache_dir, recursive=TRUE)

    cat('Global summary script. Output is to calc.log.\n')
    flog.appender(appender.file('calc.log'), name='ROOT')
    flog.info('Starting cache refresh')

    # Parallell
    Log=function(fmt, ...) { cat(sprintf(paste0(fmt, '\n'), ...)) }

    countries <- Data_Run_query_return_df ("SELECT * FROM index_parameters.country;")

    # write.csv(dplyr::select(countries,-world_bank_classification ), "temp/countries.csv",row.names = F)
    # country_name_list <- c("France","Germany")
    run_level <- "subnational"
    run_level <- "national"
    run_level <- "national_v_1_0"
    run_level <- "national_v_1_1"
    run_level <- "national_v_2_0"; data_dir  <- paste0("data_prep_temp/inputs_",run_level)
    # use_3d_array=FALSE ; run_type=FALSE
    # use_3d_array=TRUE  ; run_type="sim"
    # use_3d_array=TRUE  ; run_type="sim_baseline"
    # use_3d_array=TRUE  ; run_type="national"
    # use_3d_array=TRUE  ; run_type="ward"
    # use_3d_array=TRUE  ; run_type="blend"

    cache_dir <- paste0(data_dir,"_output_3D",use_3d_array,"_",run_type);


    Sys.time()
    if(FALSE)
    { # test
      list.files(data_dir)
      # partition_id <-  "partition_16"
      partition_id <-  "528"
      refresh_one_country_file(partition_id=partition_id,data_dir=data_dir,cache_dir=cache_dir,use_3d_array=use_3d_array,run_type=run_type)
    }


    # Index Run ------------------------------------------------------------------------------
    Sys.time()
    if(FALSE)
    {
      # Exclude existing output loc_ids.
      # unlink(cache_dir,recursive = T)
      loc_id_input_list <-list.files(data_dir)

      dir.create(cache_dir)
      # input_list_exclude <- gsub("_life_table.binary","",list.files(cache_dir))
      # input_list_exclude <- unique(gsub(".binary","",input_list_exclude))
      # loc_id_input_list  <- loc_id_input_list[!loc_id_input_list %in% input_list_exclude]
      # loc_id_input_list <- loc_id_input_list[1]
      num_thread <- 15
      clust <- parallel::makeCluster(num_thread, setup_strategy = "sequential", outfile="makeCluster_Log.txt")
      clusterExport(cl=clust, varlist=c('Log','setDF','rbindlist','config' , 'scenarios','write.fst','setDT','str_extract','read_feather'
                                        ,'cache_dir','data_dir','run_query_df',"adrop"
                                        ,'countries'
                                        ,'abind','refresh_country_files_tensor','calculate_ex_matrix','calculate_ex_lifetime_years_lost_matrix','calculate_prevalence_tensor'
                                        ,'get_loc_id','get_database_connection'
                                        ,'is_testing','get_prevalence_reference'
                                        ,'MAX_AGE'
                                        ,'get_incidence_curve' ,'pivot_wider'
                                        ,'get_incidence_growth_rates','AGES'
                                        ,'spread'
                                        ,'matrix_from_function','make_age_function'
                                        ,'calculate_prevalence'
                                        ,'complication_prevalence'
                                        ,'get_complication_parameters','get_hba1c_assumption'
                                        ,'weib_survival','calculate_dalys'
                                        ,'get_disease_weights','prevalence_and_ghost_pop','calculate_ex','calculate_ex_lifetime_years_lost','data_long_2_matrices','Data_Run_query_return_df','host_name','assert'))


      # system.time({a    <- parLapply(clust, loc_id_input_list, refresh_one_country_file,data_dir=data_dir,cache_dir=cache_dir)})
      # system.time({a <- clusterApplyLB(clust, loc_id_input_list, refresh_one_country_file,data_dir=data_dir,cache_dir=cache_dir,run_type=run_type)})
      system.time({a <- clusterApplyLB(clust, loc_id_input_list, refresh_one_country_file,data_dir=data_dir,cache_dir=cache_dir,use_3d_array=use_3d_array,run_type=run_type)})

      stopCluster(clust)
      Sys.time()
    }


  # Read run results and write to DB
    if(FALSE)
    {

      # cache_dir  <- paste0("data_prep_temp/inputs_national_v_2_0_output")
      # Read run results ------------------------------------------------------------------------------------
      # loc_id_input_list <- loc_id_input_list[loc_id_input_list$world_bank_name %in% c("Finland","Denmark","New Zealand") ,]
      run_results <- read_run_results (cache_dir)

      country_data_merge2 <- run_results$country_data_merge2
      country_data_2      <- run_results$country_data_2
      run_sequence <- setDT(country_data_merge2)[,list(.N),by=c( 'sim_start_year','sim_min_diag_rates','sim_min_non_minimal_care_perc', 'sim_min_non_minimal_care_level')]



      # Export sim runs to CF dashboard -------------------------------------------------------------------------------------
      if(FALSE)
      {
        detach("package:RPostgreSQL", unload = TRUE) # it will collide with sqldf
        library(sqldf)

        country <- countries

        country_parquet_list <- data.frame()
        country_parquet_list <- rbind(country_parquet_list, data.frame(Country="GLOBAL", file_name=paste0("GLOBAL") ))
        country_parquet_list <- rbind(country_parquet_list, data.frame(Country=unique(country$wd_income_category), file_name=paste0(unique(country$wd_income_category), "")) )
        country_parquet_list <- rbind(country_parquet_list, data.frame(Country=unique(country$wd_region), file_name=paste0(unique(country$wd_region), "")))
        country_parquet_list$loc_id <- 1001:(1000 + nrow(country_parquet_list) )
        country_parquet_list <- rbind(country_parquet_list, dplyr::select(country, Country=world_bank_name,file_name=world_bank_code,loc_id))


        country_data_merge3 <- country_data_merge2 %>% dplyr::left_join(country_data_2,by=c("sim_start_year","sim_min_diag_rates","sim_min_non_minimal_care_perc","sim_min_non_minimal_care_level","diagnosis_input","loc_id","Country","Year"))
        country_data_merge3 <- country_data_merge3 %>% dplyr::left_join(dplyr::select(countries, Country=world_bank_name, wd_region, wd_income_category),by=c("Country"))


        country_data_merge3$diagnosis_input <- NULL


        group_by_set <- paste0("
                                 sim_start_year
                               , sim_min_diag_rates
                               , sim_min_non_minimal_care_perc
                               , sim_min_non_minimal_care_level
                               , age_bracket
                               , Year
                               ")

        # result1_2 <- setDT(country_data_merge3)[,list("Prevalence+Missing"=round(sum(Prevalence+Ghosts) ,0)  ),by=c('sim_start_year'
        #                                                , 'sim_min_diag_rates'
        #                                                , 'sim_min_non_minimal_care_perc'
        #                                                , 'sim_min_non_minimal_care_level'
        #                                                , 'age_bracket'
        #                                                , 'Year', 'Country' ,'wd_region' , 'wd_income_category')]

        result1 <- sqldf(paste0('SELECT ',group_by_set,', Country           as Country, wd_region    as wd_region, wd_income_category as wd_income_category, '
                                ,query_data_points,query_data_points_index_inputs,'
                               FROM country_data_merge3  as t1
                               GROUP BY ',group_by_set ,', Country ,wd_region , wd_income_category
                               ') )

        result2 <- sqldf(paste0('SELECT ',group_by_set,', wd_region         as Country, wd_region    as wd_region, \'\'                as wd_income_category, '
                                ,query_data_points,query_data_points_index_inputs,'
                               FROM country_data_merge3  as t1
                               GROUP BY ',group_by_set ,', wd_region
                               ') )
        result3 <- sqldf(paste0('SELECT ',group_by_set,', wd_income_category as Country,  \'\'      as wd_region, wd_income_category   as wd_income_category, '
                                ,query_data_points,query_data_points_index_inputs,'
                               FROM country_data_merge3  as t1
                               GROUP BY ',group_by_set ,', wd_income_category
                               ') )

        result4 <- sqldf(paste0('SELECT ',group_by_set,', \'GLOBAL\'         as Country, \'\'       as wd_region, \'\'                 as wd_income_category, '
                                ,query_data_points,query_data_points_index_inputs,'
                               FROM country_data_merge3  as t1
                               GROUP BY ',group_by_set ,'
                               ') )

        result <- rbind(result4,result3,result2,result1)


        # split in to loc_ids.
        # result  <- read_parquet("t1dindex_simulator_448_runs.parquet")
        result  <- result %>% dplyr::left_join(country_parquet_list,by=c("Country"))

        input_rates_combined_path <-  paste0(cache_dir,"_by_country")

        dir.create(input_rates_combined_path)


        loc_id_list <- unique(result$loc_id)
        for( i in 1: length(loc_id_list) )
        { # i <- 1
          result_t <-  result[result$loc_id==loc_id_list[i],]
          # write.fst(input_rates_combined_t, paste0(paste0("data_outputs/data_",version_no,"/",loc_id_list[i],".binary") ),compress = 100)
          # saveRDS(input_rates_combined_t,paste0(paste0(input_rates_combined_path,loc_id_list[i],".Rds") ) )
          arrow::write_feather(result_t,paste0(paste0(input_rates_combined_path,"/",loc_id_list[i],"") ) )
        }


        # arrow::write_parquet(result,"t1dindex_simulator_448_runs_V3.parquet")
        arrow::write_parquet(result,"t1dindex_simulator_256_runs_V3.parquet")
        arrow::write_parquet(result,"t1dindex_simulator_runs_V4_as_minimum.parquet")
        arrow::write_parquet(result,"t1dindex_simulator_v1_as_sealevel.parquet")
        arrow::write_parquet(result,"t1dindex_simulator_v1_as_increment.parquet")
        arrow::write_parquet(result,"data_export/CleverFranke/t1dindex_v2_simulator_as_value_2.parquet")  # fix care level to replace , instead of sealevel
        result <-  arrow::read_parquet(result,"data_export/CleverFranke/t1dindex_v2_simulator_as_value_2.parquet")  # fix care level to replace , instead of sealevel


      }


      # Start writing to DB --------------------------------------------------------------------------------------------------------------------------
      diagnosis_input_values <- unique(country_data_merge2$diagnosis_input)
      version_db  <- paste0(run_level,"_3d",tolower(use_3d_array),"_",tolower(run_type))
      table_name_ <- paste0("main_",version_db)



      host_name   <<- "localhost"
      run_tracker <- read_parquet("data_prep_temp/input_rates_combined_national_run_tracker.parquet")

      i <- 1
      print(i)
      print(diagnosis_input_values[i])
      print(paste0("writing to DB \n ", host_name))
      table_name <- paste0(table_name_ )
      country_data_merge2_ <- country_data_merge2[country_data_merge2$diagnosis_input==diagnosis_input_values[i],]
      country_data_2_      <- country_data_2[country_data_2$diagnosis_input==diagnosis_input_values[i],]
      # life_table_all_      <- life_table_all[life_table_all$diagnosis_input==diagnosis_input_values[i],]
      # write_fst(life_table_all_, paste0("data_outputs/life_table_",table_name,".binary") )

      # country_data_merge2_ <- subset(country_data_merge2, diagnosis_input == diagnosis_input_values[i])
      # Generate Site data ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
      # query <- paste0( 'SELECT * FROM ',table_name,' WHERE 1=1 AND age_bracket = \'00_99\' ORDER BY  "Country" ASC ' )
      # index_data  <- Data_Run_query_return_df(query)
      index_data  <- country_data_merge2_[country_data_merge2_$age_bracket=="00_99",]
      index_data  <- cbind(index_data,country_data_2_[,-1:-3])


      data_export_final_all_year <- data.frame()
      if(FALSE)
      {
        for(year_at in 1960:2040)
        { # year_at <- 2022
          data_export_final <-  extract_for_purpose_v2 (year_at = year_at,index_data=index_data,include_gbd=FALSE)
          # data_export_final[,26:71] <- NULL # do not save growth of T1D , population
          data_export_final_all_year <- rbind(data_export_final_all_year,data_export_final )
        }
      }else
      {
        num_thread <- 10
        clust <- parallel::makeCluster(num_thread, setup_strategy = "sequential")
        clusterExport(cl=clust, varlist=c("extract_for_purpose_v2","Data_Run_query_return_df","dbConnect","host_name","dbGetQuery","dbDisconnect","setDT","helathy_years_restored"))
        # system.time({a    <- parLapply(clust, loc_id_input_list, refresh_one_country_file,data_dir=data_dir,cache_dir=cache_dir)})
        system.time({a <- clusterApplyLB(clust, 1960:2040, extract_for_purpose_v2,index_data=index_data,include_gbd=FALSE)})
        stopCluster(clust)
        data_export_final <- rbindlist(a)
        # data_export_final[,26:71] <- NULL # do not save growth of T1D , population
        data_export_final_all_year <- data_export_final
      }

      #-------------------------------------------------------------------------------------
      country_data_2_$year_at <- country_data_2_$Year
      country_data_2_$diagnosis_input <- NULL
      country_data_2_$Year            <- NULL
      country_data_2_$Age             <- NULL
      country_data_2_$`Healthy years restored with device uptake`             <- NULL
      country_data_2_$`Healthy years restored with insulin and strips`            <- NULL
      country_data_2_$`Healthy years restored with onset diagnosis`            <- NULL

      data_export_final_all_year_ <- data_export_final_all_year %>% left_join(country_data_2_, by=c("loc_id","year_at","Country") )

      input_rates_combined_final_age_specific <- readRDS( "data_prep_temp/input_rates_combined_final_age_specific.binary")
      input_rates_combined_final_age_specific <- input_rates_combined_final_age_specific[input_rates_combined_final_age_specific$year %in% data_export_final_all_year_$year_at,]
      input_rates_combined_final_age_specific$year_at <- input_rates_combined_final_age_specific$year
      colnames(input_rates_combined_final_age_specific)[colnames(input_rates_combined_final_age_specific)=="% of Nonmininal Care 1"] <- "% of Nonmininal Care"
      input_rates_combined_final_age_specific$year <- NULL

      data_export_final_all_year_ <- data_export_final_all_year_ %>% left_join(input_rates_combined_final_age_specific, by=c("loc_id","year_at") )

      # country_indicator_imputed_avg <- readRDS("data_wb/country_indicator_imputed_avg_1_5_subnational_70.Rds")
      # country_indicator_imputed_avg <- readRDS("data_wb/country_indicator_imputed_avg_1_5_subnational_70_2.Rds")
      country_indicator_imputed_avg <- readRDS("data_prep_temp/country_indicator_imputed_avg_1_5_subnational_72_2.Rds")

      data_export_final_all_year_ <- data_export_final_all_year_ %>% left_join(dplyr::select(country_indicator_imputed_avg
                                                                                             ,loc_id,year_at=year
                                                                                             ,gdp_pc_2010,gni_pc
                                                                                             ,imr,mr_u5,doctors_per_capital
                                                                                             ,pop_urban), by=c("loc_id","year_at") )
      # -----------------------------------------------------------------------------------------------------
      print(paste0(Sys.time()," | ","DB Write Tables ..."))# ------------------------------------------------------------------

      Data_dump_data_frame(df=country_data_merge2_      ,schema_name="public",table_name=paste0(table_name) )
      Data_dump_data_frame(data_export_final_all_year_  ,schema_name="public",table_name=paste0(table_name,"_site") )
      Data_dump_data_frame(run_tracker                  ,schema_name="public",table_name=paste0(table_name,"_run_tracker") )

      # arrow::write_feather(country_data_merge2_,"temp/test.feather" )

      # Data_Run_query_return_df( query_get_create_partition_table( paste0(table_name) , '"Year"' ) )
      # partition_list <- unique(country_data_merge2_$Year)
      # for(i in 1:length(partition_list))
      # {
      #   Data_Run_query_return_df( paste0(" CREATE TABLE ",table_name,"_",partition_list[i]," PARTITION OF ",table_name," FOR VALUES IN ('",partition_list[i],"');") )
      # }
      # Data_Run_query_return_df( paste0("INSERT INTO ",table_name," SELECT * FROM ",table_name,"_unpartitioned ;") )
      # Data_Run_query_return_df(paste0('CREATE INDEX idx_country_year_age_',table_name,' ON ',table_name,' ("Country", "Year","age_bracket");'))
      # Data_Run_query_return_df(paste0('CREATE INDEX idx_country_year_'    ,table_name,' ON ',table_name,' ("Country", "Year");'))
      # Data_Run_query_return_df(paste0('CREATE INDEX idx_country_'         ,table_name,' ON ',table_name,' ("Country");'))
      #------------------------------------------------------------------------------------------------------
      print(paste0(Sys.time()," | ","DB Creating Index..."))# ------------------------------------------------------------------

      Data_Run_query_return_df(paste0('CREATE INDEX idx_age_bracket_'                   ,table_name,' ON ',table_name,' ("age_bracket");'))
      Data_Run_query_return_df(paste0('CREATE INDEX idx_loc_id_year_'                   ,table_name,' ON ',table_name,' ("loc_id","Year");'))
      Data_Run_query_return_df(paste0('CREATE INDEX idx_year_'                          ,table_name,' ON ',table_name,' ("Year");'))
      Data_Run_query_return_df(paste0('CLUSTER  ',table_name,' USING idx_year_'            ,table_name,' ;'))

      Data_Run_query_return_df(paste0('CREATE INDEX idx_loc_id_year_'            ,paste0(table_name,"_site"),' ON ',paste0(table_name,"_site"),' ("loc_id",year_at);')) # Essential for table joining
      Data_Run_query_return_df(paste0('CREATE INDEX idx_year_'            ,paste0(table_name,"_site"),' ON ',paste0(table_name,"_site"),' (year_at);')) # Essential for table joining
      Data_Run_query_return_df(paste0('CLUSTER  ',table_name,'_site USING idx_loc_id_year_',table_name,'_site ;'))

      Data_Run_query_return_df(paste0('CREATE INDEX idx_loc_id'                   ,table_name,'_run_tracker ON ',table_name,'_run_tracker ("loc_id");'))
      print(paste0(Sys.time()," | ","DB Write Tables ..."))# ------------------------------------------------------------------

      # Data_dump_data_frame(run_tracker ,schema_name="public",table_name=paste0("main_1_1_60_ward_run_tracker") )
      # asdf <- read_fst("test.binary")
      # Generate Site data ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------


      # update files for Shiny
      version_list <- paste0(version_db )
      # for(i in 1:length(version_list))
      # {  #  i  <- 1
      country_all <-  get_country_stats_query_df (in_version=version_list
                                                  , filter_clause='  AND age_bracket in(\'00_99\') and "Year" in (2022 )  '
                                                  ,query_data_points_use=' string_agg(t1.loc_id::varchar(255), \',\') as loc_id_list ,round(sum("Prevalence")::numeric,0) "Prevalence"'  )

      country_all$world_bank_name_display <- world_bank_name_convert_site(country_all$world_bank_name)
      saveRDS(country_all,paste0("data_temp/country_all_",version_list))

        # table_name <- paste0("main_",gsub("\\.","_",version_list[i]),"_name_loc_mapping")
        # # country_data_merge2 <- country_data_merge2 %
        # Data_dump_data_frame(country_all ,schema_name="public",table_name=table_name )
      # }

    }

  }

  gc()

}

if(FALSE)
{
  library(data.table)
  library(dplyr)
  library(echarts4r)
  library(parallel)
  library(arrow)
  library(RPostgreSQL)

  rm(list=ls())
  source('code_R_data_prep/000_0_build_database.R')
  source('code_R/001_main.R')

  sink("log.txt",append=TRUE);cat(paste0(" Begin \n") );sink()
  # main_(version_no="0.4.12",run_per_country=TRUE,run_merge_country=TRUE)  # default run
  main_mcmc <- function(i)
  {  # i <- 1
    run_name <- paste0("0.4.12_mcmc_",i)
    # build_database(version_no=run_name,random_data_ci=TRUE)
    # main_(run_name,run_per_country=TRUE,run_merge_country=FALSE,run_few_countries=c("India","United States","Netherlands","Democratic Republic of the Congo"))
    # main_(run_name,run_per_country=TRUE,run_merge_country=FALSE,run_few_countries=c("Afghanistan"))
    # main_(version_no=run_name,run_per_country=TRUE,run_merge_country=FALSE,scenario = 1)
    main_(version_no=run_name,run_per_country=TRUE,run_merge_country=FALSE,scenario = 2)

    # sink("log.txt",append=TRUE);cat(paste0(Sys.time()," main_mcmc, process id: ", i ," finished \n") );sink()
    # closeAllConnections()
    gc()
  }
  if(FALSE)
  {# check number
    data_new <- read_parquet("../../reruns_scenario_2_launch/rerun_0.4.12_mcmc_490_0/AFG.binary")
    data_old <- read_parquet("../../reruns_scenario_2/rerun_0.4.12_mcmc_490_0/AFG.binary")
    sum(data_old$Prevalence==data_new$Prevalence)
    sum(data_old$Ghosts==data_new$Ghosts)
    sum(data_old$`Life expectency (2 t1d base)`==data_new$`Life expectency (2 t1d base)`)
  }

  sink("log.txt",append=TRUE);cat(paste0(Sys.time()," main_mcmc, Start \n") );sink()

  num_thread <- 2
  clust      <- parallel::makeCluster(num_thread, type = "PSOCK")  # stopCluster(clust)
  # clusterExport(cl=clust, varlist=c("build_database","main_"))
  clusterExport(cl=clust, varlist=c("main_"))
  system.time({a <- clusterApply(clust, 76:78, main_mcmc)})
  stopCluster(clust)
  Sys.time()

  # if(FALSE)
  # {
  #   #  SparkR install guild: https://phoenixnap.com/kb/install-spark-on-windows-10
  #   library(SparkR, lib.loc = c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib")))
  #   sparkR.session(master = "local[10]", sparkConfig = list(
  #     spark.executor.memory='1G'
  #     ,spark.driver.memory='2G'
  #     ,spark.r.backendConnectionTimeout	=60000
  #     # ,spark.executor.rpc.netty.dispatcher.numThreads = 10
  #     ,spark.serializer="org.apache.spark.serializer.KryoSerializer"
  #     # ,spark.sql.execution.arrow.sparkr.enabled    =TRUE
  #   ))
  #   # sparkR.session.stop()
  #   system.time({SparkR::spark.lapply(571:1000,main_mcmc)})
  #   # sparkR.session.stop()
  #   detach("package:SparkR", unload = TRUE)
  # }



  # merge to GLOBAL, REGION, ICOME LEVEL
  num_thread <- 40
  clust      <- parallel::makeCluster(num_thread, type = "PSOCK")  # stopCluster(clust)
  clusterExport(cl=clust, varlist=c("build_database","main_"))
  system.time({a <- clusterApply(clust, 1:1000, function(i) {
    sink("log.txt",append=TRUE);cat(paste0(Sys.time()," Merge coutry to GLOBAL , process id: ", i ," \n") );sink()
    main_( version_no=paste0("0.4.12_mcmc_",i),run_per_country=FALSE,run_merge_country=TRUE,scenario = 1)
    # main_( version_no=paste0("0.4.12_mcmc_",i),run_per_country=FALSE,run_merge_country=TRUE,scenario = 2)
  })})
  stopCluster(clust)
  Sys.time()



  # Save to AWS RDS DB Cis for all ----------------------------------------------------------------------------------------------------------------------------
  con <-  dbConnect(RSQLite::SQLite(), paste0("../../data_outputs/data_0.4.12/","data.db") , read_only=read_only)
  country <- dbGetQuery( con, "SELECT *   FROM country ")
  dbDisconnect(con)
  country_parquet_list <- data.frame()
  country_parquet_list <- rbind(country_parquet_list, data.frame(country_name="GLOBAL", file_name=paste0("GLOBAL") ))
  country_parquet_list <- rbind(country_parquet_list, data.frame(country_name=unique(country$world_bank_classification), file_name=paste0(unique(country$world_bank_classification), "")) )
  country_parquet_list <- rbind(country_parquet_list, data.frame(country_name=unique(country$wd_region), file_name=paste0(unique(country$wd_region), "")))
  country_parquet_list <- rbind(country_parquet_list, dplyr::select(country, country_name=world_bank_name,file_name=world_bank_code))


  scenario  <- 2
  schema    <- paste0("main_ci_partition_scenario_",scenario)


  file_list_country_list <- list.files(paste0("../../reruns_scenario_",scenario,"/rerun_0.4.12_mcmc_1_0"))
  file_list_country_list <- file_list_country_list[grepl(".binary",file_list_country_list)]
  # file_list_country_list <- "GLOBAL.binary"
  run_query_return_df_aws <- function(query)
  {
    # host <- "database-t1d-dev-2.ccasxjjwcbmy.ap-southeast-2.rds.amazonaws.com"
    host <- "localhost"
    connec <- dbConnect(RPostgres::Postgres(),  dbname = "t1d", host = host,
                        port = "5432",  user = "postgres",   password = "postgrest1d")
    data_put <- RPostgres::dbSendQuery(connec, query )
    dbDisconnect(connec)
    return(data_put)
  }

  get_ci_per_country <-  function(i,schema,scenario)
  {
    # i <- 2
    library(data.table)
    library(arrow)
    library(RPostgreSQL)
    library(dplyr)
    # host <- "database-t1d-dev-2.ccasxjjwcbmy.ap-southeast-2.rds.amazonaws.com"
    host <- "localhost"
    connec <- dbConnect(RPostgres::Postgres(),  dbname = "t1d", host = host,
                        port = "5432",  user = "postgres",   password = "postgrest1d")
    file_list_runs <- list.files(paste0("../../reruns_scenario_",scenario),full.names = TRUE)
    file_list_runs <- file_list_runs[grepl("0.4.12_mcmc_",file_list_runs)]
    file_list_runs <- paste0(file_list_runs, "/",file_list_country_list[i])
    country_data_merge      <- rbindlist( lapply(1:1000,function(x,file_list_runs){df <- read_parquet(file_list_runs[x]);df$run<- x; df},file_list_runs ),use.names=TRUE)
    country_data_merge      <- dplyr::select(country_data_merge,run, Country,Year,Age,`Ann. background population`
                                             ,Prevalence
                                             ,Ghosts
                                             ,`Ghosts (early death)`
                                             ,`Ghosts (onset death)`
                                             , `Incidence (1 base)`
                                             , `Life expectency (1 background)`
                                             ,`Life expectency (2 t1d base)`
                                             ,`Ann. early deaths`
                                             ,`Ann. onset deaths`
    )
    # for(y in c(2000:2022,2040) )
    for(y in c(2025,2030,2035))
    {  # y <- 2000
       print(y)
      country_data_merge_ <- country_data_merge[country_data_merge$Year %in% c(y),]
      RPostgreSQL::dbWriteTable(connec, DBI::SQL(paste0(schema,".main_ci_",y))    , setDF(country_data_merge_), overwrite=FALSE,row.names=FALSE,append=TRUE)

    }

    dbDisconnect(connec)

    sink("log.txt",append=TRUE);cat(paste0(Sys.time()," get_ci_per_country, process id: ", i ," \n") );sink()

  }

  # create tables
  get_query_create <- function(year,schema){ paste0('

          CREATE TABLE  ',schema,'.main_ci_',year,'
          (
              run integer,
              "Country" text COLLATE pg_catalog."default",
              "Year" double precision,
              "Age" double precision,
              "Ann. background population" double precision,
              "Prevalence" double precision,
              "Ghosts" double precision,
              "Ghosts (early death)" double precision,
              "Ghosts (onset death)" double precision,
              "Incidence (1 base)" double precision,
              "Life expectency (1 background)" double precision,
              "Life expectency (2 t1d base)" double precision,
              "Ann. early deaths" double precision,
              "Ann. onset deaths" double precision
          )

          PARTITION BY LIST ("Country");
  ')}



  # for(year in c(2000:2022,2040))
    for(year in c(2025,2030,2035))
    {   # year <- 2000
    print(year)
    run_query_return_df_aws(paste0('DROP TABLE IF EXISTS ',schema,'.main_ci_',year,''))
    run_query_return_df_aws(query=get_query_create(year,schema))
    for(i in 1: nrow(country_parquet_list))
    {   # i <- 191
      # print(i)
      query_partition <-  paste0("CREATE TABLE \" ",scenario,".main_ci_",year,"main_ci_",year,"_",tolower(country_parquet_list$file_name[i]),"\"
             PARTITION OF ",schema,".main_ci_",year," FOR VALUES IN  ('",gsub("'","''",(country_parquet_list$country_name[i])),"');")
      run_query_return_df_aws(query_partition)

    }

    }
  # system.time({get_ci_per_country(1,file_list_country_list)})
  num_thread <- 6
  clust      <- parallel::makeCluster(num_thread, type = "PSOCK")  # stopCluster(clust)
  clusterExport(cl=clust, varlist=c("file_list_country_list"))
  system.time({a <- clusterApply(clust, 1:213, get_ci_per_country,schema,scenario)})
  stopCluster(clust)
  Sys.time()



  # check global prevalence , missing, le diff -------------------------------
  country    <- read.csv('data_internal/country.csv',stringsAsFactors = F,encoding='UTF-8')
  # USA, GLOBAL, India, demo of congo.  HTML files, csv . upper and lower. 2000-2040. ------- ------- ------- ------- ------- -------
  country_name <- "GLOBAL"
  # country_name <- "USA"
  # country_name <- "IND"
  # country_name <- "COD"
  # country_name <- "NLD"

  file_list <- list.files(paste0("../../reruns_scenario_2"),full.names = TRUE)
  # file_list <- paste0( "../t1dGlobalModel_data/reruns/rerun_0.4.12_mcmc_",1:500,"_0")
  file_list <- file_list[grepl("0.4.12_mcmc_",file_list)]
  file_list <- paste0(file_list, "/",country_name,".binary")
  country_name <- country$world_bank_name[country$world_bank_code==country_name]
  output_file <- "t1dGlobalModel_data/mcmc_plot_paper"
  country_data_merge      <- rbindlist( lapply(file_list,function(x){df <- read_parquet(x);df$run<- x; df} ))

  data_plot <- setDT(country_data_merge)[,list(Value=sum(Prevalence[Age<100])),by=c("Country","Year","run")]
  data_plot <- data_plot[data_plot$Year>=2000,]
  # data_plot <- setDT(country_data_merge)[,list(Value=`Life expectency (2 t1d base)`[Age==10]),by=c("Country","Year","run")]
  data_plot <- data_plot[,list(value_median= median(Value)
                               ,lower= quantile(Value,probs=c(.025))
                               ,upper= quantile(Value,probs=c(.975))
                               ,n_run=.N),by=c("Year")]






  write.csv(data_plot,paste0("../",output_file,"/",country_name,"_prevalence.csv") )

  e1 <- data_plot %>%
    e_charts(Year) %>%
    e_line(value_median,name="Value") %>%
    e_band(lower,upper,
           stack = "confidence-band",
           symbol = c("none", "none"),)%>%e_x_axis(name="Year")%>%e_y_axis(name="Prevalence") %>%
    e_tooltip(trigger = "axis") %>% e_legend(bottom = 0)

  data_plot <- setDT(country_data_merge)[,list(Value=sum(Ghosts)),by=c("Country","Year","run")]
  data_plot <- data_plot[data_plot$Year>=2000,]
  # data_plot <- setDT(country_data_merge)[,list(Value=`Life expectency (2 t1d base)`[Age==10]),by=c("Country","Year","run")]
  data_plot <- data_plot[,list(value_median= median(Value)
                               ,lower= quantile(Value,probs=c(.025))
                               ,upper= quantile(Value,probs=c(.975))),by=c("Year")]
  write.csv(data_plot,paste0("../",output_file,"/",country_name,"_missing_prevalence.csv") )

  e2 <- data_plot %>%
    e_charts(Year) %>%
    e_line(value_median,name="Value") %>%
    e_band(lower,upper,
           stack = "confidence-band",
           symbol = c("none", "none"),)%>%e_x_axis(name="Year")%>%e_y_axis(name="Missing Prevalence") %>%
    e_tooltip(trigger = "axis") %>% e_legend(bottom = 0)


  data_plot <- setDT(country_data_merge)[,list(Value=sum(`Incidence (1 base)`)),by=c("Country","Year","run")]
  data_plot <- data_plot[data_plot$Year>=2000,]
  # data_plot <- data_plot[data_plot$Year>=2000 & data_plot$Year<=2020,]
  # data_plot <- setDT(country_data_merge)[,list(Value=`Life expectency (2 t1d base)`[Age==10]),by=c("Country","Year","run")]
  data_plot <- data_plot[,list(value_median= median(Value)
                               ,lower= quantile(Value,probs=c(.025))
                               ,upper= quantile(Value,probs=c(.975))),by=c("Year")]

  write.csv(data_plot,paste0("../",output_file,"/",country_name,"_incidence.csv") )

  e3 <- data_plot %>%
    e_charts(Year) %>%
    e_line(value_median,name="Value") %>%
    e_band(lower,upper,
           stack = "confidence-band",
           symbol = c("none", "none"),)%>%e_x_axis(name="Year")%>%e_y_axis(name="Incidence") %>%
    e_tooltip(trigger = "axis") %>% e_legend(bottom = 0)


  # data_plot <- setDT(country_data_merge)[,list(Value=sum(`Incidence (1 base)`)),by=c("Country","Year","run")]
  data_plot <- setDT(country_data_merge)[,list(Value=`Life expectency (2 t1d base)`[Age==10]),by=c("Country","Year","run")]
  data_plot <- data_plot[data_plot$Year>=2000,]
  data_plot <- data_plot[,list(value_median= median(Value)
                               ,lower= quantile(Value,probs=c(.025))
                               ,upper= quantile(Value,probs=c(.975))),by=c("Year")]
  write.csv(data_plot,paste0("../",output_file,"/",country_name,"_life_expectency_at10.csv") )

  e4 <- data_plot %>%
    e_charts(Year) %>%
    e_line(value_median,name="Value") %>%
    e_band(lower,upper,
           stack = "confidence-band",
           symbol = c("none", "none"),)%>%e_x_axis(name="Year")%>%e_y_axis(name="Life Expectency") %>%
    e_tooltip(trigger = "axis") %>% e_legend(bottom = 0)


  e_all <-  e_arrange( e1,e2,e3,e4,   cols = 2)





  # get ci for all -----------------------------



}



