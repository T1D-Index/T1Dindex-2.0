Agg_country <- function(country_data_merge ,key_list)
{
  agg_life_expectancy_median <- function(year,frequency)
  {
    median <- sum(year * (frequency+0.00001) )/ sum(frequency+0.00001)
    median
  }

  country_data <- setDT(country_data_merge)[,list(
    `Ann. background mortality` = sum( `Ann. background mortality` )
    ,`Ann. background population (age 10)` = sum( `Ann. background population (age 10)` )
    ,`Ann. background population` = sum( `Ann. background population` )
    ,`prevalence plus missing (age 10)`= sum(`prevalence plus missing (age 10)`)

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


    ,`Incidence (1 base)`  = sum(`Incidence (1 base)`)
    ,`Incidence (2 diagnosis)`  = sum(`Incidence (2 diagnosis)`)
    ,`Prevalence`  = sum(`Prevalence`)
    ,`Prevalence (minimal care)`  = sum(`Prevalence (minimal care)`)

    ,`1 in x families`  = sum(`1 in x families` *`Ann. background population` )/sum(`Ann. background population`)

  ),by=key_list]


  # country_data_age10 <- setDT(country_data_merge)[Age==10,list(
  country_data_age10 <- setDT(country_data_merge)[,list(

    `Life expectency (1 background)`      =  agg_life_expectancy_median(`Life expectency (1 background)`      , `Ann. background population (age 10)`  )
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

  ),by=key_list]

  # country_data_age10$`Lifetime years lost (2 t1d base)`       <- country_data_age10$`Lifetime years lost (2 t1d base) (complication)` + country_data_age10$`Lifetime years lost (2 t1d base) (treatment)`
  # country_data_age10$`Lifetime years lost (3 t1d diagnosis)`  <- country_data_age10$`Lifetime years lost (3 t1d diagnosis) (complication)` + country_data_age10$`Lifetime years lost (3 t1d diagnosis) (treatment)`
  # country_data_age10$`Lifetime years lost (4 t1d basic care)` <- country_data_age10$`Lifetime years lost (4 t1d basic care) (complication)` + country_data_age10$`Lifetime years lost (4 t1d basic care) (treatment)`
  # country_data_age10$`Lifetime years lost (5 t1d best care)`  <- country_data_age10$`Lifetime years lost (5 t1d best care) (complication)` + country_data_age10$`Lifetime years lost (5 t1d best care) (treatment)`
  # country_data_age10$`Lifetime years lost (6 t1d cure)`       <- country_data_age10$`Lifetime years lost (6 t1d cure) (complication)` + country_data_age10$`Lifetime years lost (6 t1d cure) (treatment)`

  country_data_merge_global <- country_data %>% left_join(country_data_age10, by=key_list)
  country_data_merge_global <- country_data_merge_global %>% mutate_if(is.numeric, round, digits=2)

  country_data_merge_global
}



median_per_country <- function(i)
{   # i <- 1

  library(data.table)
  library(dplyr)
  library(arrow)
  country    <- read.csv('data_internal/country.csv',stringsAsFactors = F,encoding='UTF-8')

  country_code <- country$world_bank_code[i]
  loc_id       <- country$loc_id[i]
  # country_code <- "AFG"
  # country_name <- "IND"
  # country_name <- "COD"
  # country_name <- "NLD"

  file_list <- list.files(paste0("C:/DropboxT1D/reruns_scenario_2_launch/"),full.names = TRUE)
  # file_list <- paste0( "../t1dGlobalModel_data/reruns/rerun_0.4.12_mcmc_",1:50,"_0")
  file_list <- file_list[grepl("0.4.12_mcmc_",file_list)]
  file_list <- paste0(file_list, "/",country_code,".binary")

  # file_list <- file_list[1:5]
  country_data_merge_temp      <- rbindlist( lapply(file_list,function(x){ # x <- file_list[1]
    df <- read_parquet(x);
    df$run<- x;
    df <- df[df$Year >=1960,]
    df} ))

  index_input <- readRDS(paste0("data_outputs/index_input_1_0_0/",loc_id,".Rds") )
  country_data_merge_temp <- country_data_merge_temp %>% dplyr::left_join(dplyr::select(index_input, Year=year,Age=age, value_percent_non_minimal_care), by = c("Year","Age") )
  country_data_merge_temp$`Prevalence (minimal care)` = country_data_merge_temp$Prevalence * (1-country_data_merge_temp$value_percent_non_minimal_care)
  # country_data_merge_temp$`Prevalence (non-minimal care)` = country_data_merge_temp$Prevalence * (country_data_merge_temp$value_percent_non_minimal_care)
  # cat(paste0(paste0("`",colnames(country_data_merge_temp),"`=sum(`", colnames(country_data_merge_temp),"`)" ),collapse = "\n,"))

  agg_1000_run_age_bracket <- function(country_data_merge_temp,start,end)
  { # start <- 0 ; end <- 99
    country_data_merge_agg <- setDT(country_data_merge_temp)[Age<=end & Age >= start,list(
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
      ,`Incidence (1 base)`=sum(`Incidence (1 base)`)
      ,`Incidence (2 diagnosis)`=sum(`Incidence (2 diagnosis)`)
      ,`Prevalence`=sum(`Prevalence`)
      ,`Prevalence (minimal care)`=sum(`Prevalence (minimal care)`)
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
    ),by=c("Country","Type","run","Year")]


    # cat(paste0(paste0("`",colnames(country_data_merge_temp),"`=median(`", colnames(country_data_merge_temp),"`)" ),collapse = "\n,"))

    country_data_merge_agg2 <- setDT(country_data_merge_agg)[,list(
      age_bracket= paste0(sprintf("%02d", start), "_",sprintf("%02d", end))

      ,`Ann. background mortality`=median(`Ann. background mortality`)
      ,`Ann. background population`=median(`Ann. background population`)
      ,`Ann. background population (age 10)`=median(`Ann. background population (age 10)`)
      ,`prevalence plus missing (age 10)`=median(`prevalence plus missing (age 10)`)

      ,`Ann. days lost (1 base)`=median(`Ann. days lost (1 base)`)
      ,`Ann. days lost (2 diagnosis)`=median(`Ann. days lost (2 diagnosis)`)
      ,`Ann. days lost (3 basic care)`=median(`Ann. days lost (3 basic care)`)
      ,`Ann. days lost (4 best care)`=median(`Ann. days lost (4 best care)`)
      ,`Ann. days lost (5 cure)`=median(`Ann. days lost (5 cure)`)
      ,`Ann. early deaths`=median(`Ann. early deaths`)
      ,`Ann. early deaths (2 diagnosis)`=median(`Ann. early deaths (2 diagnosis)`)
      ,`Ann. early deaths (3 basic care)`=median(`Ann. early deaths (3 basic care)`)
      ,`Ann. early deaths (4 best care)`=median(`Ann. early deaths (4 best care)`)
      ,`Ann. early deaths (5 cure)`=median(`Ann. early deaths (5 cure)`)
      ,`Ann. onset deaths`=median(`Ann. onset deaths`)
      ,`Ghosts`=median(`Ghosts`)
      ,`Ghosts (delta basic care)`=median(`Ghosts (delta basic care)`)
      ,`Ghosts (delta best care)`=median(`Ghosts (delta best care)`)
      ,`Ghosts (delta cure)`=median(`Ghosts (delta cure)`)
      ,`Ghosts (early death)`=median(`Ghosts (early death)`)
      ,`Ghosts (onset death)`=median(`Ghosts (onset death)`)
      ,`Incidence (1 base)`=median(`Incidence (1 base)`)
      ,`Incidence (2 diagnosis)`=median(`Incidence (2 diagnosis)`)
      ,`Prevalence`=median(`Prevalence`)
      ,`Prevalence (minimal care)`=median(`Prevalence (minimal care)`)
      ,`Life expectency (1 background)`=median(`Life expectency (1 background)`)
      ,`Life expectency (2 t1d base)`=median(`Life expectency (2 t1d base)`)
      ,`Life expectency (3 t1d diagnosis)`=median(`Life expectency (3 t1d diagnosis)`)
      ,`Life expectency (4 t1d basic care)`=median(`Life expectency (4 t1d basic care)`)
      ,`Life expectency (5 t1d best care)`=median(`Life expectency (5 t1d best care)`)
      ,`Life expectency (6 t1d cure)`=median(`Life expectency (6 t1d cure)`)

      ,`Lifetime years lost (2 t1d base) (complication)`=median(`Lifetime years lost (2 t1d base) (complication)`)
      ,`Lifetime years lost (3 t1d diagnosis) (complication)`=median(`Lifetime years lost (3 t1d diagnosis) (complication)`)
      ,`Lifetime years lost (4 t1d basic care) (complication)`=median(`Lifetime years lost (4 t1d basic care) (complication)`)
      ,`Lifetime years lost (5 t1d best care) (complication)`=median(`Lifetime years lost (5 t1d best care) (complication)`)
      ,`Lifetime years lost (6 t1d cure) (complication)`=median(`Lifetime years lost (6 t1d cure) (complication)`)

      ,`Lifetime years lost (2 t1d base) (treatment)`=median(`Lifetime years lost (2 t1d base) (treatment)`)
      ,`Lifetime years lost (3 t1d diagnosis) (treatment)`=median(`Lifetime years lost (3 t1d diagnosis) (treatment)`)
      ,`Lifetime years lost (4 t1d basic care) (treatment)`=median(`Lifetime years lost (4 t1d basic care) (treatment)`)
      ,`Lifetime years lost (5 t1d best care) (treatment)`=median(`Lifetime years lost (5 t1d best care) (treatment)`)
      ,`Lifetime years lost (6 t1d cure) (treatment)`=median(`Lifetime years lost (6 t1d cure) (treatment)`)

      ,`1 in x families`=median(`1 in x families`)
      ,`% Odds living to 55`=median(`% Odds living to 55`)
    ),by=c("Country","Type","Year","Age")]

    country_data_merge_agg2
  }


  country_data_merge_agg2 <- data.frame()

  age_specific <- FALSE
  # age_specific <- TRUE
  if(age_specific)
  {
    # Export for Graham,  Ethisopia, individual age 0 to 30
    country_data_merge_temp <- country_data_merge_temp[country_data_merge_temp$Year==2019,]
    for(i in 65:99)
    {
      # print(i)
      country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,i,i); country_data_$age_bracket <- i ;country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
    }
    country_data_merge_agg2 <- dplyr::select(country_data_merge_agg2,Country,Year,age_bracket, Prevalence,`Ann. background population`)
    # write.csv(country_data_merge_agg2,"temp/Ethiopia_2022_individual_age_0_30.csv")
  }else
  {

  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,00,99); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,20,99); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,20,59); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,20,25); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,60,99); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,00,19); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,00,14); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,00,11); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)

  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,00,24); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)

  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,20,39); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,40,59); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,60,79); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,80,99); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,85,99); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)


  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,20,79); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)

  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,0,4); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,5,9); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,10,14); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,15,19); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,20,24); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,25,29); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,30,34); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,35,39); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,40,44); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,45,49); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,50,54); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,55,59); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,60,64); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,65,69); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,70,74); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,75,79); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,80,84); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,85,89); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,90,94); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,95,99); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)

  country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,30,99); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)

  }


  if(FALSE)
  {
    # canada export requisted by Sarah,
    country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,60,69); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
    country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,70,79); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
    country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,80,89); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)
    country_data_ <-  agg_1000_run_age_bracket(country_data_merge_temp,90,99); country_data_merge_agg2 <- rbind(country_data_merge_agg2,country_data_)

    export <- dplyr::select(country_data_merge_agg2[country_data_merge_agg2$Year %in% c(2021,2022)], Country,Year,age_bracket,Prevalence,Incidence= `Incidence (1 base)`)
    export <- export[order(export$Year),]
    write.csv(export,"temp/export_canada_2021_2022.csv",row.names = FALSE)
  }

  country_data_merge_agg2


}


if(FALSE)
{

  host_name <<- "localhost"

  source("code_R_data_prep/DATA_CONSTANTS.R")
  # source("code_R/prevalence_utils.R")
  library(data.table)
  library(dplyr)
  library(echarts4r)
  library(parallel)
  library(RPostgreSQL)
  library(arrow)


  num_thread <- 10
  clust      <- parallel::makeCluster(num_thread, type = "PSOCK")  # stopCluster(clust)
  # clusterExport(cl=clust, varlist=c("build_database","main_"))
  # clusterExport(cl=clust, varlist=c("main_"))
  # system.time({a <- clusterApply(clust, 1:10, median_per_country)})
  # system.time({a <- clusterApply(clust, 96, median_per_country)})  # keynya only
  system.time({a <- clusterApply(clust, 1:201, median_per_country)})

  stopCluster(clust)
  Sys.time()

  country_data_merge <- rbindlist(a)

  # country_data_merge$Prevalence_rate <- round(country_data_merge$Prevalence/ country_data_merge$`Ann. background population` * 100000 ,0)
  country_data_merge <- setDT(country_data_merge)[,Background_population_world_sum := sum(`Ann. background population`) , by=c("age_bracket","Year")]



  current <- country_data_merge[country_data_merge$Country=="Malaysia" & country_data_merge$age_bracket %in% c("00_11", "00_24", "00_99" ) & country_data_merge$Year=="2022",c("Country","Year","Prevalence","age_bracket")]
  current$assume <- 3000
  current$assume <- current$Prevalence * current$assume/583.705
  # write_parquet(country_data_merge,"../t1dGlobalModel_temp/country_data_merge_all.binary")
  # write_parquet(country_data_merge,"../t1dGlobalModel_temp/country_data_merge_all_graham.binary")
  # write_parquet(country_data_merge,"../t1dGlobalModel_temp/country_data_merge_all_dianna.binary")
  write.csv(country_data_merge,"temp/all_country_65_99_age_specific.csv")
  # write_parquet(country_data_merge,"../t1dGlobalModel_temp/country_data_merge_all_for_finland_registry.binary")
  # country_data_merge <- read_parquet("../t1dGlobalModel_temp/country_data_merge_all_for_finland_registry.binary")

  # merge with lever2023
  main_0_4_15_lever2023 <- Data_Run_query_return_df ("SELECT * FROM public.main_0_4_15_lever2023;")
  merge_agg_age_bracket <- function( df_input, start,end)
  { # start = 0 ; end=99;
    country_data_merge_agg <- setDT(df_input)[,list(
      Age= 10
      , age_bracket= paste0(sprintf("%02d", start), "_",sprintf("%02d", end))
      ,`Ghosts`=sum(`Ghosts`[Age<=end & Age >= start])
      ,`Ghosts (delta basic care)`=sum(`Ghosts (delta basic care)`[Age<=end & Age >= start])
      ,`Ghosts (delta best care)`=sum(`Ghosts (delta best care)`[Age<=end & Age >= start])
      ,`Ghosts (delta cure)`=sum(`Ghosts (delta cure)`[Age<=end & Age >= start])
      ,`Ghosts (early death)`=sum(`Ghosts (early death)`[Age<=end & Age >= start])
      ,`Ghosts (onset death)`=sum(`Ghosts (onset death)`[Age<=end & Age >= start])

    ),by=c("Country","Type","Year")]
    country_data_merge_agg
  }

  country_data_merge2 <- data.frame()
  age_bracket_list <- unique(country_data_merge$age_bracket)
  for( i in 1: length(age_bracket_list))
  {
    print(age_bracket_list[i])
    start <- as.numeric(strsplit(age_bracket_list[i],"_")[[1]][1])
    end   <- as.numeric(strsplit(age_bracket_list[i],"_")[[1]][2])
    country_data_ <-  merge_agg_age_bracket(main_0_4_15_lever2023,start,end)
    country_data_merge2 <- rbind(country_data_merge2,country_data_)
  }


  main_0_4_15_lever2023_ <- dplyr::select(country_data_merge2 , "Country" ,"Year","Age","age_bracket"
                                         ,"Ghosts lever2023"= "Ghosts"
                                         ,"Ghosts (delta basic care) lever2023"="Ghosts (delta basic care)"
                                         ,"Ghosts (delta best care) lever2023"="Ghosts (delta best care)"
                                         ,"Ghosts (delta cure) lever2023"="Ghosts (delta cure)"
                                         ,"Ghosts (early death) lever2023"="Ghosts (early death)"
                                         ,"Ghosts (onset death) lever2023"="Ghosts (onset death)")

  country_data_merge            <- country_data_merge %>% dplyr::inner_join(main_0_4_15_lever2023_,by=c("Country" ,"Year","Age","age_bracket"))



  # Global, HiC etc... -----------------------------------------------------------------------------
  countries <- Data_Run_query_return_df ("SELECT * FROM index_parameters.country;")
  # colnames_saved                <- colnames(country_data_merge)
  country_data_merge            <- dplyr::select(countries,Country=world_bank_name,wd_region,wd_income_category  ) %>% dplyr::inner_join(country_data_merge,by="Country")

  country_data_merge_agg <- Agg_country(country_data_merge ,key_list = c("Year","Age","age_bracket"))
  country_data_merge_agg$Country <-  'GLOBAL'
  country_data_merge_agg$Type    <-  'GLOBAL'
  colnames_saved                <- colnames(country_data_merge_agg)


  country_data_merge_global <- setDF(country_data_merge_agg)[,colnames_saved]

  country_data_merge_agg <- Agg_country(country_data_merge ,key_list = c("Year","Age","age_bracket","wd_region") )
  country_data_merge_agg$Country <-  country_data_merge_agg$wd_region
  country_data_merge_agg$Type    <-  'wd_region'
  country_data_merge_wd_region <- setDF(country_data_merge_agg)[,colnames_saved]


  country_data_merge_agg <- Agg_country(country_data_merge ,key_list = c("Year","Age","age_bracket","wd_income_category") )
  country_data_merge_agg$Country <-  country_data_merge_agg$wd_income_category
  country_data_merge_agg$Type    <-  'wd_income_category'
  country_data_merge_wd_income_category <- setDF(country_data_merge_agg)[,colnames_saved]


  country_data_merge_all <- rbind(setDF(country_data_merge)[,colnames_saved]
                                  ,country_data_merge_global
                                  ,country_data_merge_wd_region
                                  ,country_data_merge_wd_income_category)

  country_data_merge_all <- country_data_merge_all %>% mutate_if(is.numeric, round, digits=2)


  # country_data_merge <- read_parquet("temp/country_data_merge_all.binary")
  # country_data_merge_all <- setDF(country_data_merge)


  # merge with single run ---------------------------------------------------------
  query <- paste0( 'SELECT *
                      -- FROM main_0_4_15_paper_median
                       FROM main_0_4_15
                       WHERE 1=1 AND "Age" = 10
                      --  AND age_bracket = \'00_99\'
                      --  AND "Year"=  2021
                       ORDER BY  "Country" ASC , "Year", "Age"  '
  )
  main_0_4_15 <- Data_Run_query_return_df(query)
  main_0_4_15 <- main_0_4_15[main_0_4_15$Country %in% country_data_merge_all$Country,]

  query <- paste0( 'SELECT "Country","Year", "Age", "% Odds living to"
                       FROM main_0_4_15  WHERE 1=1
                       AND "Age" = 55
                       ORDER BY  "Country" ASC , "Year", "Age"  '
  )
  main_0_4_15_odds_living_to <- Data_Run_query_return_df(query)
  main_0_4_15_odds_living_to <- main_0_4_15_odds_living_to[main_0_4_15_odds_living_to$Country %in% country_data_merge_all$Country,]

  main_0_4_15$`% Odds living to 55` <- main_0_4_15_odds_living_to$`% Odds living to`

  query <- paste0( 'SELECT "Country","Year", "Age", "% Odds living to"
                       FROM main_0_4_15  WHERE 1=1
                       AND "Age" = 60
                       ORDER BY  "Country" ASC , "Year", "Age"  '
  )
  main_0_4_15_odds_living_to <- Data_Run_query_return_df(query)
  main_0_4_15_odds_living_to <- main_0_4_15_odds_living_to[main_0_4_15_odds_living_to$Country %in% country_data_merge_all$Country,]

  main_0_4_15$`% Odds living to 60` <- main_0_4_15_odds_living_to$`% Odds living to`

  query <- paste0( 'SELECT "Country","Year", "Age", "% Odds living to"
                       FROM main_0_4_15  WHERE 1=1
                       AND "Age" = 65
                       ORDER BY  "Country" ASC , "Year", "Age"  '
  )
  main_0_4_15_odds_living_to <- Data_Run_query_return_df(query)
  main_0_4_15_odds_living_to <- main_0_4_15_odds_living_to[main_0_4_15_odds_living_to$Country %in% country_data_merge_all$Country,]

  main_0_4_15$`% Odds living to 65` <- main_0_4_15_odds_living_to$`% Odds living to`



  main_0_4_15$true_incidence_age_10 <- main_0_4_15$`Incidence (1 base)`
  main_0_4_15$`Life expectency (2 t1d base) (Default run)` <- main_0_4_15$`Life expectency (2 t1d base)`

  name_list_replace <-  c("true_incidence_age_10","Life expectency (2 t1d base) (Default run)","Life expectency (3 t1d diagnosis)",                     "Life expectency (4 t1d basic care)"
                         ,"Life expectency (5 t1d best care)",                     "Life expectency (6 t1d cure)",                          "Lifetime years lost (2 t1d base) (complication)"
                          ,"Lifetime years lost (3 t1d diagnosis) (complication)",  "Lifetime years lost (4 t1d basic care) (complication)", "Lifetime years lost (5 t1d best care) (complication)"
                          ,"Lifetime years lost (6 t1d cure) (complication)",       "Lifetime years lost (2 t1d base) (treatment)",          "Lifetime years lost (3 t1d diagnosis) (treatment)"
                          ,"Lifetime years lost (4 t1d basic care) (treatment)",    "Lifetime years lost (5 t1d best care) (treatment)",     "Lifetime years lost (6 t1d cure) (treatment)"
                          ,"Life expectency (strip low)",                          "Life expectency (strip hig)",                          "Lifetime years lost (strip low)"
                          ,"Lifetime years lost (strip hig)",                      "Life expectency (sensor low)",                         "Life expectency (sensor hig)"
                          ,"Lifetime years lost (sensor low)",                     "Lifetime years lost (sensor hig)",                     "1 in x families"
                         , "% Odds living to 55" , "% Odds living to 60" , "% Odds living to 65"
                         ,"Ann. days lost (1 base)", "Ann. days lost (2 diagnosis)", "Ann. days lost (3 basic care)" ,"Ann. days lost (4 best care)", "Ann. days lost (5 cure)")

  country_data_merge_all <- country_data_merge_all[, setdiff(colnames(country_data_merge_all), name_list_replace)]


  country_data_merge_all_conbine <- country_data_merge_all %>% dplyr::inner_join(main_0_4_15[,c("Country" ,"Year","Age",name_list_replace) ],by=c("Country" ,"Year","Age") )


  # add loc id ----------------------------
  country  <- Data_Run_query_return_df ("SELECT * FROM country order by world_bank_name ;")
  country_parquet_list <- data.frame()
  country_parquet_list <- rbind(country_parquet_list, data.frame(country_name="GLOBAL", file_name=paste0("GLOBAL") ))
  country_parquet_list <- rbind(country_parquet_list, data.frame(country_name=unique(country$world_bank_classification), file_name=paste0(unique(country$world_bank_classification), "")) )
  country_parquet_list <- rbind(country_parquet_list, data.frame(country_name=unique(country$wd_region), file_name=paste0(unique(country$wd_region), "")))
  country_parquet_list$loc_id <- 1001:(1000 + nrow(country_parquet_list) )
  country_parquet_list <- rbind(country_parquet_list, dplyr::select(country, country_name=world_bank_name,file_name=world_bank_code,loc_id))
  country_data_merge_all_conbine  <-country_data_merge_all_conbine %>% left_join(select(country_parquet_list, Country=country_name,  loc_id))


#   country_data_merge_all_conbine_read <-     Data_Run_query_return_df("select * from main_0_4_15_paper_median")
  # Dump to DB ---------------------------------------------------------------------------------------
  # table_name <- "main_0_4_15_paper_median"
  table_name <- "main_1_0_0"  # 1.0 is the same as main_0_4_15_paper_median, reorganised to work on new dashboard
  Data_dump_data_frame(country_data_merge_all_conbine ,schema_name="public",table_name=table_name )

  #------------------------------------------------------------------------------------------------------
  print(paste0(Sys.time()," | ","DB Creating Index..."))

  Data_Run_query_return_df(paste0('CREATE INDEX idx_age_bracket_'                   ,table_name,' ON ',table_name,' ("age_bracket");'))
  Data_Run_query_return_df(paste0('CREATE INDEX idx_loc_id_year_'                   ,table_name,' ON ',table_name,' ("loc_id","Year");'))
  Data_Run_query_return_df(paste0('CREATE INDEX idx_year_'                          ,table_name,' ON ',table_name,' ("Year");'))
  Data_Run_query_return_df(paste0('CLUSTER  ',table_name,' USING idx_year_'            ,table_name,' ;'))

  Data_Run_query_return_df(paste0('CREATE INDEX idx_loc_id_year_'            ,paste0(table_name,"_site"),' ON ',paste0(table_name,"_site"),' ("loc_id",year_at);')) # Essential for table joining
  Data_Run_query_return_df(paste0('CREATE INDEX idx_year_'            ,paste0(table_name,"_site"),' ON ',paste0(table_name,"_site"),' (year_at);')) # Essential for table joining
  Data_Run_query_return_df(paste0('CLUSTER  ',table_name,'_site USING idx_loc_id_year_',table_name,'_site ;'))

  Data_Run_query_return_df(paste0('CREATE INDEX idx_loc_id'                   ,table_name,'_run_tracker ON ',table_name,'_run_tracker ("loc_id");'))

  # write_parquet(country_data_merge,"temp/test.binary")
  # country_data_merge_all <- arrow::read_parquet(("temp/test.binary"))
}
