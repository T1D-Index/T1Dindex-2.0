#!/usr/bin/env Rscript
# setwd( "E:/GithubCode/t1d_global_model/t1dGlobalModel")
library(RSQLite)
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
# library(arrow)

source('code_R/runner.R')
source('code_R/prevalence.R')
source('code_R/data.R')
source('code_R/utils.R')
source('code_R/incidence.R')
source('code_R/mortality.R')
source('code_R/complications.R')
source('code_R/burden.R')
source('code_quick_job_scripts/extract_for_purpose.R')

time_start <- Sys.time()
if (!dir.exists('cache')) {
    dir.create('cache')
}

source('code_R/000_parameters.R')
# config$run_projection      <- FALSE

for(i in 8:length(scenarios))
{
  # i <- 0  # default run
  # i <- 14  # diagnosis_rate_left run
  # i <- 16  # diagnosis_rate_right run
    print(i)
    scenarios[] <- FALSE
    ifelse(i==0,scenarios[i] <- FALSE,scenarios[i] <- TRUE )

    # scenarios[] <- FALSE , i <- cha

    # version_no <- "0.4.2"
    # version_no <- "0.4.9"
    # version_no <- "0.4.11"
    version_no <- "0.4.13"

    Sys.setenv(T1D_DATA_FILE = paste0("data_outputs/data_",version_no,"/","data.db") )
    data_dir  <- paste0( "data_outputs/data_",version_no)
    # cache_dir <- paste0('reruns/rerun_',version_no)
    cache_dir <- paste0('../reruns/rerun_',version_no,"_",i)
    # cache_dir <- paste0('reruns/rerun_',version_no,"_c",)


    if(scenarios$adult_incidence_all_studies)
    {Sys.setenv(T1D_DATA_FILE = paste0("data_outputs/data_",version_no,".adult_incidence_all_studies/","data.db") )
        data_dir  <- paste0( "data_outputs/data_",version_no,".adult_incidence_all_studies")   }
    if(scenarios$iot_no_growth)
    {Sys.setenv(T1D_DATA_FILE = paste0("data_outputs/data_",version_no,".iot_no_growth/","data.db") )
        data_dir  <- paste0( "data_outputs/data_",version_no,".iot_no_growth")}
    if(scenarios$iot_global_curve_for_all)
    {Sys.setenv(T1D_DATA_FILE = paste0("data_outputs/data_",version_no,".iot_global_curve_for_all/","data.db") )
        data_dir  <- paste0( "data_outputs/data_",version_no,".iot_global_curve_for_all")}
    if(scenarios$smr_age_cruve_flat_average)
    {Sys.setenv(T1D_DATA_FILE = paste0("data_outputs/data_",version_no,".smr_age_cruve_flat_average/","data.db") )
        data_dir  <- paste0( "data_outputs/data_",version_no,".smr_age_cruve_flat_average")}
    if(scenarios$adult_onset_zero)
    {Sys.setenv(T1D_DATA_FILE = paste0("data_outputs/data_",version_no,".adult_onset_zero/","data.db") )
        data_dir  <- paste0( "data_outputs/data_",version_no,".adult_onset_zero")}

    # Create SQLite database  data.db ------------------------------------------------------------------------------
    # unlink(cache_dir, recursive=TRUE)
    dir.create(cache_dir)

    cat('Global summary script. Output is to calc.log.\n')
    flog.appender(appender.file('calc.log'), name='ROOT')
    flog.info('Starting cache refresh')

    # Parallell
    Log=function(fmt, ...) { cat(sprintf(paste0(fmt, '\n'), ...)) }
    countries         <- get_countries_and_regions()
    countries$wd_income_category[countries$wd_income_category==''] <- countries$world_bank_classification[countries$wd_income_category=='' ]
    countries <- countries[countries$idf_country_name!='',]
    # write.csv(dplyr::select(countries,-world_bank_classification ), "temp/countries.csv",row.names = F)
    country_name_list <- countries$world_bank_name
    # country_name_list <- c("France","Germany")
    Sys.time()
    num_thread <- 5
    clust <- parallel::makeCluster(num_thread, setup_strategy = "sequential")
    clusterExport(cl=clust, varlist=c('Log','config' , 'scenarios','write.fst','extract_for_purpose','setDT'
                                      ,'cache_dir','data_dir','run_query_df'
                                      ,'countries'
                                      ,'get_loc_id','get_database_connection'
                                      ,'is_testing','get_prevalence_reference'
                                      ,'MAX_AGE','incidence_function'
                                      ,'get_incidence_curve' ,'pivot_wider'
                                      ,'get_incidence_growth_rates','AGES'
                                       ,'spread'
                                      ,'matrix_from_function','make_age_function'
                                      ,'matrix_from_function'  ,'background_mortality_function','t1d_mortality_function'
                                      ,'calculate_prevalence'
                                      ,'complication_prevalence'
                                      ,'get_complication_parameters','get_hba1c_assumption'
                                      ,'weib_survival','calculate_dalys'
                                      ,'get_disease_weights','prevalence_and_ghost_pop','calculate_ex','data_long_2_matrices'))

    # num_loops <- ceiling( length(country_name_list)/num_thread)
    # for(p in 1:num_loops)
    # {   # p <- 51
    #     # incProgress(1/num_loops, detail = paste("Doing part ", p, "/",num_loops)) # shiny run
    #     start <- 1+(p-1)*num_thread
    #     end   <- min(p*num_thread, length(country_name_list))
    #     country_name_list_t <-  country_name_list[start:end]
    #     print(country_name_list_t) ; print(Sys.time())
    #     system.time({a <- parLapply(clust, country_name_list_t, refresh_one_country_file)})
    # }
    system.time({a <- parLapply(clust, country_name_list, refresh_one_country_file)})
    # refresh_one_country_file(country_wb_name="Rwanda")
    # refresh_one_country_file(country_wb_name="United States")


    stopCluster(clust)
    Sys.time()

    # Merge and write to xlsx ----------------------------------------------------------------------------------------------
    flog.info('Generate summary_break_down_age.xlsx..')
    # combine data
    files_db <- list.files(cache_dir,full.names=TRUE)
    files_db <- files_db[files_db %in% paste0(cache_dir,"/",countries$world_bank_code,".binary")]
    # files_db <- files_db[!grepl(".parquet", files_db)]
    country_data_merge      <- rbindlist( lapply(files_db, read_fst))

    country_data_merge      <- gather(country_data_merge ,key = "Value_type", value = "Value",-Country,-Type,-Year,-Age)


    country_data_merge_cname      <- colnames(country_data_merge)
    country_data_merge            <- dplyr::select(countries,Country=world_bank_name,wd_region,wd_income_category  ) %>% dplyr::inner_join(country_data_merge,by="Country")

    # calculate regions and income class and global  -----------------------------
    index_caculate_mean <- grepl("Life expectency",country_data_merge$Value_type)|country_data_merge$Value_type=="1 in x families"
    country_data_global1 <- as.data.frame( setDT(country_data_merge)[index_caculate_mean,
                                                                     list(Country= 'GLOBAL',Type  = 'GLOBAL',Value=mean(Value) ),by=c('Value_type','Year','Age')] )

    country_data_global2 <- as.data.frame( setDT(country_data_merge)[!index_caculate_mean,
                                                                     list(Country= 'GLOBAL',Type  = 'GLOBAL',Value=sum(Value)),by=c('Value_type','Year','Age')] )
    country_data_global <- rbind(country_data_global1,country_data_global2)

    # country_data_global$Value[country_data_global$Value_type=="Incidence (1 base)" & country_data_global$Year== 2021]/country_data_global$Value[country_data_global$Value_type=="Ann. background population" & country_data_global$Year== 2021]*100000

    country_data_region1 <- as.data.frame( setDT(country_data_merge)[wd_region!='' & index_caculate_mean,
                                                                     list(Country= wd_region,Type   = 'wd_region',Value=mean(Value )),by=c('wd_region','Value_type','Year','Age')] )
    country_data_region2 <- as.data.frame( setDT(country_data_merge)[wd_region!='' & (!index_caculate_mean),
                                                                     list(Country= wd_region,Type   = 'wd_region',Value=sum(Value)),by=c('wd_region','Value_type','Year','Age')] )
    country_data_region <- rbind(country_data_region1,country_data_region2)

    # country_data_region$Value[country_data_region$Value_type=="Incidence (1 base)" & country_data_region$Year== 2021]/country_data_region$Value[country_data_region$Value_type=="Ann. background population" & country_data_region$Year== 2021]*100000


    country_data_region_broad1 <- as.data.frame( setDT(country_data_merge)[wd_region!=''  & index_caculate_mean,
                                                                          list( Country= wd_income_category,Type   = 'wd_income_category',Value=mean(Value )),by=c('wd_income_category','Value_type','Year','Age')] )
    country_data_region_broad2 <- as.data.frame( setDT(country_data_merge)[wd_region!=''  & (!index_caculate_mean),
                                                                          list( Country= wd_income_category,Type   = 'wd_income_category',Value=sum(Value)),by=c('wd_income_category','Value_type','Year','Age')] )

    country_data_region_broad <- rbind(country_data_region_broad1,country_data_region_broad2)


    setDF(country_data_merge)

    countrys_and_regions <-  data.frame(country_name = c(unique(country_data_global$Country),
                                                         unique(country_data_region_broad$Country),
                                                         unique(country_data_region$Country),
                                                         unique(country_data_merge$Country)
                                                                ))
    countrys_and_regions$file_name <- 1:nrow(countrys_and_regions)

    country_data_merge <- rbind(country_data_merge[,country_data_merge_cname]
                                  ,country_data_global[,country_data_merge_cname]
                                  ,country_data_region[,country_data_merge_cname]
                                  ,country_data_region_broad[,country_data_merge_cname] )

    dir.create(paste0(cache_dir,"/prev_merge_wide"))
    dir.create(paste0(cache_dir,"/aggregates"))

    country_stats_over_year_wide <- spread(country_data_merge, Value_type, Value )
    # over-write life-expectency calculation for Global, region, and Income level ------- ------- ------- ------- ------- -------

    country_stats_over_year_wide_with_region <- dplyr::select(countries,Country=world_bank_name,wd_region,wd_income_category  )  %>% dplyr::inner_join(country_stats_over_year_wide, by = "Country")

    life_expec_  <- setDT(country_stats_over_year_wide_with_region)[,list(`Life expectency (1 background)` = sum(`Life expectency (1 background)` *`Ann. background population` )/sum(`Ann. background population`) ),by= c('Year','Age')]
    country_stats_over_year_wide$`Life expectency (1 background)`[country_stats_over_year_wide$Country=="GLOBAL"] <- life_expec_$`Life expectency (1 background)`

    life_expec_  <- setDT(country_stats_over_year_wide_with_region)[,list(
      `Life expectency (2 t1d base)`       = sum(`Life expectency (2 t1d base)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (3 t1d diagnosis)`  = sum(`Life expectency (3 t1d diagnosis)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (4 t1d basic care)` = sum(`Life expectency (4 t1d basic care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (5 t1d best care)`  = sum(`Life expectency (5 t1d best care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (6 t1d cure)`       = sum(`Life expectency (6 t1d cure)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),


      `Lifetime years lost (2 t1d base)`       = sum(`Lifetime years lost (2 t1d base)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (3 t1d diagnosis)`  = sum(`Lifetime years lost (3 t1d diagnosis)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (4 t1d basic care)` = sum(`Lifetime years lost (4 t1d basic care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (5 t1d best care)`  = sum(`Lifetime years lost (5 t1d best care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (6 t1d cure)`       = sum(`Lifetime years lost (6 t1d cure)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts)

    ),by= c('Year','Age')]

    country_stats_over_year_wide$`Life expectency (2 t1d base)`[country_stats_over_year_wide$Country=="GLOBAL"]         <- life_expec_$`Life expectency (2 t1d base)`
    country_stats_over_year_wide$`Life expectency (3 t1d diagnosis)`[country_stats_over_year_wide$Country=="GLOBAL"]    <- life_expec_$`Life expectency (3 t1d diagnosis)`
    country_stats_over_year_wide$`Life expectency (4 t1d basic care)`[country_stats_over_year_wide$Country=="GLOBAL"]   <- life_expec_$`Life expectency (4 t1d basic care)`
    country_stats_over_year_wide$`Life expectency (5 t1d best care)`[country_stats_over_year_wide$Country=="GLOBAL"]    <- life_expec_$`Life expectency (5 t1d best care)`
    country_stats_over_year_wide$`Life expectency (6 t1d cure)`[country_stats_over_year_wide$Country=="GLOBAL"]         <- life_expec_$`Life expectency (6 t1d cure)`

    country_stats_over_year_wide$`Lifetime years lost (2 t1d base)`[country_stats_over_year_wide$Country=="GLOBAL"]         <- life_expec_$`Lifetime years lost (2 t1d base)`
    country_stats_over_year_wide$`Lifetime years lost (3 t1d diagnosis)`[country_stats_over_year_wide$Country=="GLOBAL"]    <- life_expec_$`Lifetime years lost (3 t1d diagnosis)`
    country_stats_over_year_wide$`Lifetime years lost (4 t1d basic care)`[country_stats_over_year_wide$Country=="GLOBAL"]   <- life_expec_$`Lifetime years lost (4 t1d basic care)`
    country_stats_over_year_wide$`Lifetime years lost (5 t1d best care)`[country_stats_over_year_wide$Country=="GLOBAL"]    <- life_expec_$`Lifetime years lost (5 t1d best care)`
    country_stats_over_year_wide$`Lifetime years lost (6 t1d cure)`[country_stats_over_year_wide$Country=="GLOBAL"]         <- life_expec_$`Lifetime years lost (6 t1d cure)`



    life_expec_wd_region  <- setDT(country_stats_over_year_wide_with_region)[,list(`Life expectency (1 background)` = sum(`Life expectency (1 background)` *`Ann. background population` )/sum(`Ann. background population`) ),by= c('wd_region','Year','Age')]
    life_expec_wd_region$Country <- life_expec_wd_region$wd_region
    life_expec_wd_region_ <- country_stats_over_year_wide[country_stats_over_year_wide$Type=="wd_region",]
    life_expec_wd_region_2 <- life_expec_wd_region_ %>% dplyr::inner_join(life_expec_wd_region, by = c('Country','Year','Age'))
    country_stats_over_year_wide$`Life expectency (1 background)`[country_stats_over_year_wide$Type=="wd_region"] <- life_expec_wd_region_2$`Life expectency (1 background).y`


    life_expec_wd_region  <- setDT(country_stats_over_year_wide_with_region)[,list(
      `Life expectency (2 t1d base)`       = sum(`Life expectency (2 t1d base)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (3 t1d diagnosis)`  = sum(`Life expectency (3 t1d diagnosis)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (4 t1d basic care)` = sum(`Life expectency (4 t1d basic care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (5 t1d best care)`  = sum(`Life expectency (5 t1d best care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (6 t1d cure)`       = sum(`Life expectency (6 t1d cure)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),

      `Lifetime years lost (2 t1d base)`       = sum(`Lifetime years lost (2 t1d base)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (3 t1d diagnosis)`  = sum(`Lifetime years lost (3 t1d diagnosis)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (4 t1d basic care)` = sum(`Lifetime years lost (4 t1d basic care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (5 t1d best care)`  = sum(`Lifetime years lost (5 t1d best care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (6 t1d cure)`       = sum(`Lifetime years lost (6 t1d cure)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts)
      ),by= c('wd_region','Year','Age')]
    life_expec_wd_region$Country <- life_expec_wd_region$wd_region
    life_expec_wd_region_ <- country_stats_over_year_wide[country_stats_over_year_wide$Type=="wd_region",]
    life_expec_wd_region_2 <- life_expec_wd_region_ %>% dplyr::inner_join(life_expec_wd_region, by = c('Country','Year','Age'))
    country_stats_over_year_wide$`Life expectency (2 t1d base)`[country_stats_over_year_wide$Type=="wd_region"]         <- life_expec_wd_region_2$`Life expectency (2 t1d base).y`
    country_stats_over_year_wide$`Life expectency (3 t1d diagnosis)`[country_stats_over_year_wide$Type=="wd_region"]    <- life_expec_wd_region_2$`Life expectency (3 t1d diagnosis).y`
    country_stats_over_year_wide$`Life expectency (4 t1d basic care)`[country_stats_over_year_wide$Type=="wd_region"]   <- life_expec_wd_region_2$`Life expectency (4 t1d basic care).y`
    country_stats_over_year_wide$`Life expectency (5 t1d best care)`[country_stats_over_year_wide$Type=="wd_region"]    <- life_expec_wd_region_2$`Life expectency (5 t1d best care).y`
    country_stats_over_year_wide$`Life expectency (6 t1d cure)`[country_stats_over_year_wide$Type=="wd_region"]         <- life_expec_wd_region_2$`Life expectency (6 t1d cure).y`

    country_stats_over_year_wide$`Lifetime years lost (2 t1d base)`[country_stats_over_year_wide$Type=="wd_region"]         <- life_expec_wd_region_2$`Lifetime years lost (2 t1d base).y`
    country_stats_over_year_wide$`Lifetime years lost (3 t1d diagnosis)`[country_stats_over_year_wide$Type=="wd_region"]    <- life_expec_wd_region_2$`Lifetime years lost (3 t1d diagnosis).y`
    country_stats_over_year_wide$`Lifetime years lost (4 t1d basic care)`[country_stats_over_year_wide$Type=="wd_region"]   <- life_expec_wd_region_2$`Lifetime years lost (4 t1d basic care).y`
    country_stats_over_year_wide$`Lifetime years lost (5 t1d best care)`[country_stats_over_year_wide$Type=="wd_region"]    <- life_expec_wd_region_2$`Lifetime years lost (5 t1d best care).y`
    country_stats_over_year_wide$`Lifetime years lost (6 t1d cure)`[country_stats_over_year_wide$Type=="wd_region"]         <- life_expec_wd_region_2$`Lifetime years lost (6 t1d cure).y`


    life_expec_wd_region  <- setDT(country_stats_over_year_wide_with_region)[,list(`Life expectency (1 background)` = sum(`Life expectency (1 background)` *`Ann. background population` )/sum(`Ann. background population`) ),by= c('wd_income_category','Year','Age')]
    life_expec_wd_region$Country <- life_expec_wd_region$wd_income_category
    life_expec_wd_region_ <- country_stats_over_year_wide[country_stats_over_year_wide$Type=="wd_income_category",]
    life_expec_wd_region_2 <- life_expec_wd_region_ %>% dplyr::inner_join(life_expec_wd_region, by = c('Country','Year','Age'))
    country_stats_over_year_wide$`Life expectency (1 background)`[country_stats_over_year_wide$Type=="wd_income_category"] <- life_expec_wd_region_2$`Life expectency (1 background).y`


    life_expec_wd_region  <- setDT(country_stats_over_year_wide_with_region)[,list(
      `Life expectency (2 t1d base)`       = sum(`Life expectency (2 t1d base)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (3 t1d diagnosis)`  = sum(`Life expectency (3 t1d diagnosis)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (4 t1d basic care)` = sum(`Life expectency (4 t1d basic care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (5 t1d best care)`  = sum(`Life expectency (5 t1d best care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Life expectency (6 t1d cure)`       = sum(`Life expectency (6 t1d cure)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),

      `Lifetime years lost (2 t1d base)`       = sum(`Lifetime years lost (2 t1d base)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (3 t1d diagnosis)`  = sum(`Lifetime years lost (3 t1d diagnosis)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (4 t1d basic care)` = sum(`Lifetime years lost (4 t1d basic care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (5 t1d best care)`  = sum(`Lifetime years lost (5 t1d best care)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts),
      `Lifetime years lost (6 t1d cure)`       = sum(`Lifetime years lost (6 t1d cure)` *(Prevalence+Ghosts) )/sum(Prevalence+Ghosts)
      ),by= c('wd_income_category','Year','Age')]
    life_expec_wd_region$Country <- life_expec_wd_region$wd_income_category
    life_expec_wd_region_ <- country_stats_over_year_wide[country_stats_over_year_wide$Type=="wd_income_category",]
    life_expec_wd_region_2 <- life_expec_wd_region_ %>% dplyr::inner_join(life_expec_wd_region, by = c('Country','Year','Age'))

    country_stats_over_year_wide$`Life expectency (2 t1d base)`[country_stats_over_year_wide$Type=="wd_income_category"]         <- life_expec_wd_region_2$`Life expectency (2 t1d base).y`
    country_stats_over_year_wide$`Life expectency (3 t1d diagnosis)`[country_stats_over_year_wide$Type=="wd_income_category"]    <- life_expec_wd_region_2$`Life expectency (3 t1d diagnosis).y`
    country_stats_over_year_wide$`Life expectency (4 t1d basic care)`[country_stats_over_year_wide$Type=="wd_income_category"]   <- life_expec_wd_region_2$`Life expectency (4 t1d basic care).y`
    country_stats_over_year_wide$`Life expectency (5 t1d best care)`[country_stats_over_year_wide$Type=="wd_income_category"]    <- life_expec_wd_region_2$`Life expectency (5 t1d best care).y`
    country_stats_over_year_wide$`Life expectency (6 t1d cure)`[country_stats_over_year_wide$Type=="wd_income_category"]         <- life_expec_wd_region_2$`Life expectency (6 t1d cure).y`

    country_stats_over_year_wide$`Lifetime years lost (2 t1d base)`[country_stats_over_year_wide$Type=="wd_income_category"]         <- life_expec_wd_region_2$`Lifetime years lost (2 t1d base).y`
    country_stats_over_year_wide$`Lifetime years lost (3 t1d diagnosis)`[country_stats_over_year_wide$Type=="wd_income_category"]    <- life_expec_wd_region_2$`Lifetime years lost (3 t1d diagnosis).y`
    country_stats_over_year_wide$`Lifetime years lost (4 t1d basic care)`[country_stats_over_year_wide$Type=="wd_income_category"]   <- life_expec_wd_region_2$`Lifetime years lost (4 t1d basic care).y`
    country_stats_over_year_wide$`Lifetime years lost (5 t1d best care)`[country_stats_over_year_wide$Type=="wd_income_category"]    <- life_expec_wd_region_2$`Lifetime years lost (5 t1d best care).y`
    country_stats_over_year_wide$`Lifetime years lost (6 t1d cure)`[country_stats_over_year_wide$Type=="wd_income_category"]         <- life_expec_wd_region_2$`Lifetime years lost (6 t1d cure).y`


    country_parquet_list <- data.frame()
    country_parquet_list <- rbind(country_parquet_list, data.frame(country_name="GLOBAL", file_name=paste0("GLOBAL") ))
    country_parquet_list <- rbind(country_parquet_list, data.frame(country_name=unique(countries$world_bank_classification), file_name=paste0(unique(countries$world_bank_classification), "")) )
    country_parquet_list <- rbind(country_parquet_list, data.frame(country_name=unique(countries$wd_region), file_name=paste0(unique(countries$wd_region), "")))
    # country_parquet_list <- rbind(country_parquet_list, dplyr::select(country, country_name=world_bank_name,file_name=world_bank_code))

     # save parquet for each aggregates , including regions,
    for(p in 1:nrow(country_parquet_list) )
    { print(paste0("save prev_wide parquet for ",p) )  # p  <- 1
      prev_merge_wide_t <-  country_stats_over_year_wide[country_stats_over_year_wide$Country==country_parquet_list$country_name[p],]
      # write_parquet(prev_merge_wide_t, paste0(cache_dir,"/",country_parquet_list$file_name[p],".binary"))
      # saveRDS(prev_merge_wide_t, paste0(cache_dir,"/prev_merge_wide/",countrys_and_regions$file_name[p],".test"))
      write_fst(prev_merge_wide_t, paste0(cache_dir,"/",country_parquet_list$file_name[p],".binary"),compress = 90)
    }

    # arrow::write_parquet(country_stats_over_year_wide,paste0(cache_dir,"/aggregates/country_data_merge.parquet"),compression ="uncompressed")
    # write to SQlite database ----------------------------------------------------------------------------------------------------------------------
    db_file_path <-  paste0("data_outputs_final/",version_no,"_",i,".db")
    file.remove(db_file_path)
    con <- dbConnect(RSQLite::SQLite(),db_file_path)
    # dbExecute(con,"VACUUM;")
    dbWriteTable(con, "main"       , country_stats_over_year_wide, overwrite=TRUE)
    dbExecute(con,"CREATE INDEX idx_country_year ON main (Country, Year);")
    dbExecute(con,"CREATE INDEX idx_country ON main (Country);")
    dbExecute(con,"CREATE INDEX idx_year ON main (Year);")
    dbDisconnect(con)

    # run_query_df_path("select  * from main limit 10",db_file_path)


    # exports to excel , purpose  -----------------------------  -----------------------------  -----------------------------  -----------------------------
    if(FALSE)
    {
      # paper info, get proportion of prevalence from extrapolated countries idf ------
      # what's the proportion of stats when a country is using extrapolated incidence number from IDF ? 2021
      data <- country_stats_over_year[country_stats_over_year$Year==2021,]
      data <- data[data$world_bank_name %in% countries$world_bank_name,]

      data$is_extrapolated <- data$idf_reference_country!= 'N/A'
      data <- setDT(data)[,list(
        number_of_countries= .N,
        HIC= sum(world_bank_classification=="HIC"),
        LIC= sum(world_bank_classification=="LIC"),
        LMIC= sum(world_bank_classification=="LMIC"),
        UMIC = sum(world_bank_classification=="UMIC"),
        `Ann. background population`=sum(`Ann. background population`),
        prevalence=sum(Prevalence),
        Ghosts=sum(Ghosts)
      ),by=is_extrapolated]

    }

    if(FALSE)
    {
        setDF(country_data_merge)
        country_data_merge_wide  <- round(country_stats_over_year_wide[,c(-1,-2,-3,-4)],2)
        country_data_merge_wide  <- cbind(country_stats_over_year_wide[,c(1,2,3,4)], country_data_merge_wide)

        # fwrite(country_data_merge_wide[,], "summary_break_down_age.csv")
        # saveRDS(country_data_merge_wide, paste0(cache_dir,"/aggregates/t1dGlobalModel.Rds") )
        fwrite(country_data_merge_wide[1:1000000,],  paste0(cache_dir,"/aggregates/summary_break_down_age_1-1000000.csv") )
        fwrite(country_data_merge_wide[1000001:nrow(country_data_merge_wide),],  paste0(cache_dir,"/aggregates/summary_break_down_age_1000000-end.csv"))


        # age bin 5 years , eg: Age=0 denotes  0 <= age <=4 ------------------------------------------------
        prev_merge_long      <- gather(country_data_merge_wide ,key = "Value_type", value = "Value",-Country,-Type,-Year,-Age)
        prev_merge_long$Year <- as.numeric(prev_merge_long$Year)
        prev_merge_long$Age  <- as.numeric(prev_merge_long$Age)
        prev_merge_long <- prev_merge_long[,colnames(country_data_merge)]


        country_data_merge_bin       <- prev_merge_long
        country_data_merge_bin$Age   <- (country_data_merge_bin$Age) -  (country_data_merge_bin$Age) %% 5

        country_data_merge_bin_le       <- setDT(country_data_merge_bin)[ grepl("Life expectency",Value_type),list(Value=Value[1] ),by=c( 'Country', 'Type','Value_type','Year','Age')]
        country_data_merge_bin_non_le   <- setDT(country_data_merge_bin)[!grepl("Life expectency",Value_type),list(Value=sum(Value) ),by=c( 'Country', 'Type','Value_type','Year','Age')]

        country_data_merge_bin       <- rbind(country_data_merge_bin_le,country_data_merge_bin_non_le )
        country_data_merge_bin       <- spread(country_data_merge_bin, Value_type, Value )
        fwrite(country_data_merge_bin,  paste0(cache_dir,"/aggregates/summary_break_down_age_bin.csv" ) )

        print(Sys.time() - time_start)
        # Read the 2 CSV file names from working directory
        Zip_Files <- list.files( paste0(cache_dir,"/aggregates/"),pattern = "*.csv")
        # Zip the files and place the zipped file in working directory
        zip::zip(zipfile  =  paste0("t1dGlobalModel_",version_no,".zip"),
                 files    =  Zip_Files ,
                 root     =  paste0(cache_dir,"/aggregates/"))

        unlink(paste0(cache_dir,"/aggregates/",Zip_Files))
    }

}

if(FALSE)
{

  # sensitivity analysis

  compare_to_published_estimates <- function(country_data_merge)
  {
    library(sqldf)
    published_estimates <- read.csv("data_sensitivity_analysis/published estimates.csv", stringsAsFactors = FALSE)
    published_estimates <- published_estimates[1:27,]
    published_estimates <- published_estimates[!is.na(published_estimates$Year),]
    published_estimates$Prevalence <- as.integer(trimws(gsub("\\,","", published_estimates$Prevalence )))
    published_estimates$row_no <- 1:nrow(published_estimates)
    published_estimates$Country[published_estimates$Country=="Taiwan"] <- "China, Taiwan Province of China"
    published_estimates <-  sqldf("select t1.row_no,t1.Country,t1.Year,t1.Source,t1.Low, t1.High,t1.Prevalence, round(sum(t2.Prevalence),0) as Prevalence_T1D_index_estimates
                  FROM published_estimates as  t1 left join  country_data_merge t2
                  WHERE t1.Year=t2.Year
                  AND t1.Country= t2.Country
                  AND (t2.Age >= t1.Low AND t2.Age <= t1.High)
                  GROUP BY t1.row_no, t1.Country,t1.Year,t1.Source,t1.Low, t1.High,t1.Prevalence
                  ORDER BY t1.row_no asc
                  --limit 100")
    published_estimates$Country[published_estimates$Country=="China, Taiwan Province of China"] <- "Taiwan"

    published_estimates$Variance <- round((published_estimates$Prevalence_T1D_index_estimates - published_estimates$Prevalence)/  published_estimates$Prevalence_T1D_index_estimates * 100,0)

    published_estimates$"Var(abs)" <- abs(published_estimates$Variance)

    return(published_estimates)
  }

  get_data_merge <- function( file_path)
  {
    filenames <- list.files(paste0(file_path), pattern="*.binary", full.names=TRUE)
    ldf <- lapply(filenames, read_parquet)
    country_data_merge   <- rbindlist(ldf)

  }

  # sensitivity summary
  # country_data_merge <-  arrow::read_parquet( paste0("reruns/rerun_0.4.10_0/aggregates/country_data_merge.parquet") )
  country_data_merge <-  get_data_merge( file_path= paste0("reruns/rerun_0.4.12_0/"))

  # sum(country_data_merge[ country_data_merge$Year==2021  & country_data_merge$Country=="Mali"& country_data_merge$Age <= 25 ])

  # sum(country_data_merge$Prevalence[ country_data_merge$Year==2021  & country_data_merge$Country=="Mali"& country_data_merge$Age <= 25 ])

  prevalance_base                  <- sum(country_data_merge$Prevalence[ country_data_merge$Year==2021  & country_data_merge$Country=="GLOBAL" ])
  Ghosts_base                           <- sum(country_data_merge$Ghosts    [ country_data_merge$Year==2021  & country_data_merge$Country=="GLOBAL" ])
  published_estimates <- compare_to_published_estimates(country_data_merge)
  Var_base            <- round( median( published_estimates$`Var(abs)`) ,0)

  scenarios_df <- as.data.frame( t(as.data.frame(scenarios)) )
  scenarios_df$Scenario  <- rownames(scenarios_df)
  scenarios_df$Prevalance <- NA
  scenarios_df$Ghosts     <- NA
  scenarios_df$"Var(abs)"     <- NA

  for(i in 1:length(scenarios))
  { # i <- 13
    print(i)
    country_data_merge <-  get_data_merge( file_path= paste0("reruns/rerun_0.4.12_",i,"/"))

    scenarios_df$Prevalance[i] <- sum(country_data_merge$Prevalence[ country_data_merge$Year==2021  & country_data_merge$Country=="GLOBAL" ])
    scenarios_df$Ghosts[i]    <- sum(country_data_merge$Ghosts     [ country_data_merge$Year==2021  & country_data_merge$Country=="GLOBAL" ])
    published_estimates <- compare_to_published_estimates(country_data_merge)
    scenarios_df$`Var(abs)` [i]    <-         round( median( published_estimates$`Var(abs)`) ,0)

  }

  scenarios_df_base <- scenarios_df[1,]
  scenarios_df_base$Scenario <- "base_scenario"
  scenarios_df_base$Prevalance <- prevalance_base
  scenarios_df_base$Ghosts <- Ghosts_base
  scenarios_df_base$`Var(abs)`  <- Var_base
  scenarios_df <- rbind(scenarios_df_base,scenarios_df)

  scenarios_df$`Change in global prevalence (2021)`         <- round( (scenarios_df$Prevalance - scenarios_df$Prevalance[1]) / scenarios_df$Prevalance[1] *100,0)
  scenarios_df$`Change in global missing prevalence (2021)` <- round( (scenarios_df$Ghosts - scenarios_df$Ghosts[1]) / scenarios_df$Ghosts[1] *100,0)
  scenarios_df$`Change in ‘fit’ with observed prevalence(registries & studies)` <- round( (scenarios_df$`Var(abs)` - scenarios_df$`Var(abs)`[1]) / scenarios_df$`Var(abs)`[1] *100,0)

  scenarios_df$V1 <- NULL

  # write.csv(scenarios_df, "data_sensitivity_analysis/sensivitity_analysis_output_v2.csv",row.names = FALSE)
  write.csv(scenarios_df, "data_sensitivity_analysis/sensivitity_analysis_output_v0.4.12.csv",row.names = FALSE)




}
