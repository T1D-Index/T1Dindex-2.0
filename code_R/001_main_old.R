{
  
  
  # for ( i in 0:15) { print(i); cache_dir <- paste0('../../reruns_scenario_',scenario,'/rerun_',version_no,"_",i)
  # Generate incomelevel/region/global and write to output -------------------------------------------------------------------------
  if(run_merge_country)
  {
    countries <- Data_Run_query_return_df ("SELECT * FROM index_parameters.country;")
    
    # Merge and write to xlsx ----------------------------------------------------------------------------------------------
    flog.info('Generate summary_break_down_age.xlsx..')
    # combine data
    files_db <- list.files(cache_dir,full.names=TRUE)
    files_db <- files_db[files_db %in% paste0(cache_dir,"/",countries$loc_id,".binary")]
    # files_db <- files_db[!grepl(".parquet", files_db)]
    country_data_merge      <- rbindlist( lapply(files_db, read_parquet))
    
    colnames_saved          <- colnames(country_data_merge)
    country_data_merge            <- dplyr::select(countries,Country=world_bank_name,wd_region,wd_income_category  ) %>% dplyr::inner_join(country_data_merge,by="Country")
    
    country_data_merge_agg <- Agg_country(country_data_merge ,key_list = c("Year","Age"))
    country_data_merge_agg$Country <-  'GLOBAL'
    country_data_merge_agg$Type    <-  'GLOBAL'
    country_data_merge_global <- setDF(country_data_merge_agg)[,colnames_saved]
    
    country_data_merge_agg <- Agg_country(country_data_merge ,key_list = c("Year","Age","wd_region") )
    country_data_merge_agg$Country <-  country_data_merge_agg$wd_region
    country_data_merge_agg$Type    <-  'wd_region'
    country_data_merge_wd_region <- setDF(country_data_merge_agg)[,colnames_saved]
    
    
    country_data_merge_agg <- Agg_country(country_data_merge ,key_list = c("Year","Age","wd_income_category") )
    country_data_merge_agg$Country <-  country_data_merge_agg$wd_income_category
    country_data_merge_agg$Type    <-  'wd_income_category'
    country_data_merge_wd_income_category <- setDF(country_data_merge_agg)[,colnames_saved]
    
    
    country_data_merge_all <- rbind(setDF(country_data_merge)[,colnames_saved]
                                    ,country_data_merge_global
                                    ,country_data_merge_wd_region
                                    ,country_data_merge_wd_income_category)
    
    country_data_merge_all <- country_data_merge_all %>% mutate_if(is.numeric, round, digits=2)
    
    country_parquet_list <- data.frame()
    country_parquet_list <- rbind(country_parquet_list, data.frame(country_name="GLOBAL", file_name=paste0("GLOBAL") ))
    country_parquet_list <- rbind(country_parquet_list, data.frame(country_name=unique(countries$world_bank_classification), file_name=paste0(unique(countries$world_bank_classification), "")) )
    country_parquet_list <- rbind(country_parquet_list, data.frame(country_name=unique(countries$wd_region), file_name=paste0(unique(countries$wd_region), "")))
    # country_parquet_list <- rbind(country_parquet_list, dplyr::select(country, country_name=world_bank_name,file_name=world_bank_code))
    
    
    
    # save parquet for each aggregates , including regions,
    for(p in 1:nrow(country_parquet_list) )
    { print(paste0("save prev_wide parquet for ",p) )  # p  <- 1
      prev_merge_wide_t <-  country_data_merge_all[country_data_merge_all$Country==country_parquet_list$country_name[p],colnames_saved]
      # write_parquet(prev_merge_wide_t, paste0(cache_dir,"/",country_parquet_list$file_name[p],".binary"))
      # saveRDS(prev_merge_wide_t, paste0(cache_dir,"/prev_merge_wide/",countrys_and_regions$file_name[p],".test"))
      write_parquet(prev_merge_wide_t, paste0(cache_dir,"/",country_parquet_list$file_name[p],".binary"))
    }
    
  }
  
  # arrow::write_parquet(country_stats_over_year_wide,paste0(cache_dir,"/aggregates/country_data_merge.parquet"),compression ="uncompressed")
  if(FALSE)
  {
    files_db <- list.files(cache_dir,full.names=TRUE)
    files_db <- files_db[files_db %in% paste0(cache_dir,"/",countries$loc_id,".binary")]
    files_db <- files_db[grepl(".binary", files_db)]
    country_data_merge      <- rbindlist( lapply(files_db, read_parquet))
    
    if(TRUE)
    {
      merge_agg_age_bracket <- function( country_data_merge, start,end)
      { # start = 0 ; end=99;
        country_data_merge_agg <- setDT(country_data_merge)[Age<=end & Age >= start,list(
          Age= 10
          ,`Ann. background mortality`=sum(`Ann. background mortality`)
          ,`Ann. background population`=sum(`Ann. background population`)
          
          ,`Ann. background population (age 10)`=sum(`Ann. background population`[Age==10])
          ,`prevalence plus missing (age 10)`=sum(`Ghosts`[Age==10]) + sum(`Prevalence`[Age==10])
          
          ,`Ann. days lost (1 base)`=sum(`Ann. days lost (1 base)`)
          ,`Ann. days lost (2 diagnosis)`=sum(`Ann. days lost (2 diagnosis)`)
          ,`Ann. days lost (3 basic care)`=sum(`Ann. days lost (3 basic care)`)
          ,`Ann. days lost (4 best care)`=sum(`Ann. days lost (4 best care)`)
          ,`Ann. days lost (5 cure)`=sum(`Ann. days lost (5 cure)`)
          ,`Ann. early deaths`=sum(`Ann. early deaths`)
          ,`Ann. early deaths (2 diagnosis)`=sum(`Ann. early deaths (2 diagnosis)`)
          ,`Ann. early deaths (3 basic care)`=sum(`Ann. early deaths (3 basic care)`)
          ,`Ann. early deaths (4 best care)`=sum(`Ann. early deaths (4 best care)`)
          ,`Ann. early deaths (5 cure)`=sum(`Ann. early deaths (5 cure)`)
          ,`Ann. onset deaths`=sum(`Ann. onset deaths`)
          # ,`diagnosis rate`=sum(`Ann. onset deaths`[Age==10])/ (sum(`Ann. onset deaths`[Age==10]) + `Incidence (1 base)`[Age==10])
          ,`Ghosts`=sum(`Ghosts`)
          ,`Ghosts (delta basic care)`=sum(`Ghosts (delta basic care)`)
          ,`Ghosts (delta best care)`=sum(`Ghosts (delta best care)`)
          ,`Ghosts (delta cure)`=sum(`Ghosts (delta cure)`)
          ,`Ghosts (early death)`=sum(`Ghosts (early death)`)
          ,`Ghosts (onset death)`=sum(`Ghosts (onset death)`)
          
          ,`Ghosts lever2023`=sum(`Ghosts lever2023`)
          ,`Ghosts (delta basic care) lever2023`=sum(`Ghosts (delta basic care) lever2023`)
          ,`Ghosts (delta best care) lever2023`=sum(`Ghosts (delta best care) lever2023`)
          ,`Ghosts (delta cure) lever2023`=sum(`Ghosts (delta cure) lever2023`)
          ,`Ghosts (early death) lever2023`=sum(`Ghosts (early death) lever2023`)
          ,`Ghosts (onset death) lever2023`=sum(`Ghosts (onset death) lever2023`)
          
          
          ,`Incidence (1 base)`=sum(`Incidence (1 base)`)
          ,`Incidence (2 diagnosis)`=sum(`Incidence (2 diagnosis)`)
          ,`Prevalence`=sum(`Prevalence`)
          ,`Life expectency (1 background)`=sum(`Life expectency (1 background)`[Age==10])
          ,`Life expectency (2 t1d base)`=sum(`Life expectency (2 t1d base)`[Age==10])
          ,`Life expectency (3 t1d diagnosis)`=sum(`Life expectency (3 t1d diagnosis)`[Age==10])
          ,`Life expectency (4 t1d basic care)`=sum(`Life expectency (4 t1d basic care)`[Age==10])
          ,`Life expectency (5 t1d best care)`=sum(`Life expectency (5 t1d best care)`[Age==10])
          ,`Life expectency (6 t1d cure)`=sum(`Life expectency (6 t1d cure)`[Age==10])
          ,`Lifetime years lost (2 t1d base) (complication)`=sum(`Lifetime years lost (2 t1d base) (complication)`[Age==10])
          ,`Lifetime years lost (3 t1d diagnosis) (complication)`=sum(`Lifetime years lost (3 t1d diagnosis) (complication)`[Age==10])
          ,`Lifetime years lost (4 t1d basic care) (complication)`=sum(`Lifetime years lost (4 t1d basic care) (complication)`[Age==10])
          ,`Lifetime years lost (5 t1d best care) (complication)`=sum(`Lifetime years lost (5 t1d best care) (complication)`[Age==10])
          ,`Lifetime years lost (6 t1d cure) (complication)`=sum(`Lifetime years lost (6 t1d cure) (complication)`[Age==10])
          ,`Lifetime years lost (2 t1d base) (treatment)`=sum(`Lifetime years lost (2 t1d base) (treatment)`[Age==10])
          ,`Lifetime years lost (3 t1d diagnosis) (treatment)`=sum(`Lifetime years lost (3 t1d diagnosis) (treatment)`[Age==10])
          ,`Lifetime years lost (4 t1d basic care) (treatment)`=sum(`Lifetime years lost (4 t1d basic care) (treatment)`[Age==10])
          ,`Lifetime years lost (5 t1d best care) (treatment)`=sum(`Lifetime years lost (5 t1d best care) (treatment)`[Age==10])
          ,`Lifetime years lost (6 t1d cure) (treatment)`=sum(`Lifetime years lost (6 t1d cure) (treatment)`[Age==10])
          ,`1 in x families`=sum(`1 in x families`[Age==10])
          ,`% Odds living to 55`= sum(`Prevalence`[Age==55])/ (sum(`Ghosts`[Age==55]) + sum(`Prevalence`[Age==55])) * 100
          ,`% Odds living to 60`= sum(`Prevalence`[Age==60])/ (sum(`Ghosts`[Age==60]) + sum(`Prevalence`[Age==60])) * 100
          ,`% Odds living to 65`= sum(`Prevalence`[Age==65])/ (sum(`Ghosts`[Age==65]) + sum(`Prevalence`[Age==65])) * 100
        ),by=c("Country","Type","Year")]
        country_data_merge_agg
      }
      
      country_data_merge_agg2_00_99 <-  merge_agg_age_bracket(country_data_merge,00,99)
      country_data_merge_agg2_00_99$age_bracket <- "00_99"
      country_data_merge_agg2_20_99 <-  merge_agg_age_bracket(country_data_merge,20,99)
      country_data_merge_agg2_20_99$age_bracket <- "20_99"
      country_data_merge_agg2_20_59 <-  merge_agg_age_bracket(country_data_merge,20,59)
      country_data_merge_agg2_20_59$age_bracket <- "20_59"
      country_data_merge_agg2_60_99 <-  merge_agg_age_bracket(country_data_merge,60,99)
      country_data_merge_agg2_60_99$age_bracket <- "60_99"
      country_data_merge_agg2_00_19 <-  merge_agg_age_bracket(country_data_merge,00,19)
      country_data_merge_agg2_00_19$age_bracket <- "00_19"
      country_data_merge_agg2_00_24 <-  merge_agg_age_bracket(country_data_merge,00,24)
      country_data_merge_agg2_00_24$age_bracket <- "00_24"
      country_data_merge <- rbind( country_data_merge_agg2_00_99
                                   ,country_data_merge_agg2_20_99
                                   ,country_data_merge_agg2_20_59
                                   ,country_data_merge_agg2_60_99
                                   ,country_data_merge_agg2_00_19
                                   ,country_data_merge_agg2_00_24
      )
      country_data_merge
    }
    # table_name <- "main_0_4_15_lifetime_low_high"
    table_name <- paste0("main_",gsub("\\.","_",version_no))
    
    Data_dump_data_frame(country_data_merge ,schema_name="public",table_name=table_name )
    
    Data_Run_query_return_df(paste0('CREATE INDEX idx_country_year_',table_name,' ON ',table_name,' ("Country", "Year");'))
    Data_Run_query_return_df(paste0('CREATE INDEX idx_country_'     ,table_name,' ON ',table_name,' ("Country");'))
    Data_Run_query_return_df(paste0('CREATE INDEX idx_year_'        ,table_name,' ON ',table_name,' ("Year");'))
    
    
    
  }
}   # old code













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