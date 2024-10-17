library(data.table)
library(dplyr)
library(RPostgreSQL)
library(arrow)

source("code_R/data.R")
source("code_shiny/QUERY_DATA_POINTS.R")
source("code_R_data_prep/DATA_CONSTANTS.R")
host_name   <<- "localhost"

# in_version <- "1_1_72_gregory_national" ;
# in_version <- "1_0_0" ;
in_version <- "1_1_72_gregory_subnational"
# Load 1.0 and reformat same as 2.0
main_national            <- Data_Run_query_return_df ( paste0("SELECT * FROM  public.main_1_0_0;") )

prevalence_plus_missing  <- main_national[main_national$age_bracket=="00_99" ,c("Country","Year","prevalence plus missing (age 10)" )]
colnames(prevalence_plus_missing)[colnames(prevalence_plus_missing)=="Year"] <- "year_at"


main_national_site   <- Data_Run_query_return_df ( paste0("SELECT * FROM  public.main_1_0_0_site;") )
main_national_site   <- main_national_site[main_national_site$Country %in% main_national$Country,]


main_national_site1  <- (main_national[main_national$age_bracket=="00_99",c(1,30:49) ])
colnames(main_national_site1)[colnames(main_national_site1)=="Year"] <- "year_at"

main_national_site <- main_national_site %>% left_join (main_national_site1, by=c("Country" , "year_at") )
main_national_site <- main_national_site %>% left_join (prevalence_plus_missing, by=c("Country" , "year_at") )


# Load subnational
main_sub_national <- Data_Run_query_return_df ( paste0("SELECT * FROM  main_1_1_72_gregory_subnational_original;") )
main_sub_national_site <- Data_Run_query_return_df ( paste0("SELECT * FROM  main_1_1_72_gregory_subnational_site_original;") )
# Data_dump_data_frame (main_sub_national     ,"public","public.main_1_1_72_gregory_subnational_original")
# Data_dump_data_frame (main_sub_national_site,"public","main_1_1_72_gregory_subnational_site_original")
main_sub_national_new      <- main_sub_national
main_sub_national_site_new <- main_sub_national_site

main_national <- main_national[main_national$age_bracket %in% main_sub_national$age_bracket,]
main_national <- main_national[main_national$Country %in% main_sub_national$Country,]

country_list <- unique(main_sub_national$Country)
age_bra_list <- unique(main_sub_national$age_bracket)

aligh_number_by_name <- function(main_sub_national,main_national , column_name )
{
  #  column_name <- "Ann. background population"
  main_national_temp <- main_national[,c( c("age_bracket","Country"  ,"Year", column_name ))]
  colnames(main_national_temp)[4] <- "column_name"
  main_sub_national[,"column_name_sub"] <-  main_sub_national[,column_name]
  main_sub_national <- main_sub_national  %>% left_join (main_national_temp, by=c("Country" ,"age_bracket", "Year") )
  main_sub_national <- setDT(main_sub_national)[,column_name_ratio := column_name/sum(column_name_sub ) ,by=c("Country" ,"age_bracket", "Year") ]
  setDF(main_sub_national)
  main_sub_national[,column_name] <- main_sub_national[,column_name] * main_sub_national[,"column_name_ratio"]


  main_sub_national$column_name_sub     <- NULL
  main_sub_national$column_name        <- NULL
  main_sub_national$column_name_ratio  <- NULL

  return( main_sub_national)

}
name_list <- c("Ann. background population"
               , "Ann. early deaths"
               , "Ann. onset deaths"
               , "Incidence (1 base)"
               , "Prevalence"
               , "Ann. early deaths (2 diagnosis)"
               , "Ann. early deaths (3 basic care)"
               , "Ann. early deaths (4 best care)"
               , "Ann. early deaths (5 cure)"
               , "Ghosts"
               , "Ghosts (delta basic care)"
               , "Ghosts (delta best care)"
               , "Ghosts (delta cure)"
               , "Ghosts (early death)"
               , "Ghosts (onset death)"
               , "Incidence (2 diagnosis)"
              )
for(i in 1:length(name_list))
{
  print(i)
  main_sub_national_new <- aligh_number_by_name (main_sub_national_new,main_national , column_name=name_list[i] )
}


colnames_list  <- colnames(main_sub_national)
colnames_list  <- colnames(main_sub_national_site)
for(i in 1:length(colnames_list))
{
  contains <- grepl(colnames_list[i], query_data_points, fixed = TRUE)
  if(contains)
  {
    print(colnames_list[i])
  }
}


# adjust "prevalence plus missing (age 10)"
main_sub_national_site_new$age_bracket <- 1
main_national_site$age_bracket         <- 1
main_sub_national_site_new$Year <- main_sub_national_site_new$year_at
main_national_site$Year         <- main_national_site$year_at
main_sub_national_site_new <- aligh_number_by_name (main_sub_national_site_new,main_national_site , column_name="prevalence plus missing (age 10)" )
main_sub_national_site_new$age_bracket <- NULL
main_national_site$age_bracket         <- NULL
main_sub_national_site_new$Year <- NULL
main_national_site$Year         <- NULL


aligh_number_by_name_site <- function(main_sub_national_site,main_national_site , column_name )
{
  #  column_name <- "Lifetime years lost (2 t1d base) (complication)"
  main_national_temp <- main_national_site[,c( c("Country"  ,"year_at", column_name ))]
  colnames(main_national_temp)[3] <- "column_name"
  main_sub_national_site[,"column_name_sub"] <-  main_sub_national_site[,column_name]
  main_sub_national_site <- main_sub_national_site  %>% left_join (main_national_temp, by=c("Country" , "year_at") )
  main_sub_national_site <- setDT(main_sub_national_site)[,column_name_ratio := column_name/ (sum(column_name_sub *`prevalence plus missing (age 10)` ) /sum(`prevalence plus missing (age 10)`) )  ,by=c("Country" , "year_at") ]
  setDF(main_sub_national_site)
  main_sub_national_site[,column_name] <- main_sub_national_site[,column_name] * main_sub_national_site[,"column_name_ratio"]

  # view2 <- main_sub_national_site[main_sub_national_site$Country=="Australia" & main_sub_national_site$year_at==2020 ,c("Country","loc_id","year_at","prevalence plus missing (age 10)",column_name )]
  # view1 <- main_national_site[main_national_site$Country=="Australia" & main_national_site$year_at==2020,c("Country","loc_id","year_at",column_name )]


  main_sub_national_site$column_name_sub     <- NULL
  main_sub_national_site$column_name        <- NULL
  main_sub_national_site$column_name_ratio  <- NULL
  return( main_sub_national_site)
}

name_list_site <- c(
                 "Life expectency (1 background)"
               , "Life expectency (2 t1d base)"
               , "Life expectency (3 t1d diagnosis)"
               , "Life expectency (4 t1d basic care)"
               , "Life expectency (5 t1d best care)"
               , "Lifetime years lost (2 t1d base) (complication)"
               , "Lifetime years lost (3 t1d diagnosis) (complication)"
               , "Lifetime years lost (4 t1d basic care) (complication)"
               , "Lifetime years lost (5 t1d best care) (complication)"
               , "Lifetime years lost (2 t1d base) (treatment)"
               , "Lifetime years lost (3 t1d diagnosis) (treatment)"
               , "Lifetime years lost (4 t1d basic care) (treatment)"
               , "Lifetime years lost (5 t1d best care) (treatment)"
               , "Healthy years restored with onset diagnosis"
               , "Healthy years restored with insulin and strips"
               , "Healthy years restored with device uptake"
)

# replace "prevalence plus missing (age 10)" with new one



for(i in 1:length(name_list_site))
{
  print(i)
  main_sub_national_site_new <- aligh_number_by_name_site (main_sub_national_site_new,main_national_site , column_name=name_list_site[i] )
}

# Data_dump_data_frame (main_sub_national_new     ,"public","main_1_1_72_gregory_subnational")
# Data_dump_data_frame (main_sub_national_site_new,"public","main_1_1_72_gregory_subnational_site")







