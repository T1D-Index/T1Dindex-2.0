#!/usr/bin/env Rscript
# draft paper : https://docs.google.com/document/d/1JDVJVxQ-4XFap5CTxPGUEdG5r0SKQm2QEsxCYYcBRcU/edit
library(RSQLite)
library(dplyr)
library(purrr)
library(testthat)
library(tidyverse)
library(data.table)

source('code_R/runner.R')
source('code_R/prevalence.R')
source('code_R/data.R')
source('code_R/population.R')
source('code_R/utils.R')
source('code_R/incidence.R')
source('code_R/onset_death.R')
source('code_R/mortality.R')
source('code_R/hba1c.R')
source('code_R/complications.R')
source('code_R/burden.R')



# Sys.setenv(T1D_DATA_FILE = "/Users/feiwang/data.db")
# Sys.setenv(T1D_DATA_FILE = "data_outputs/data.db")


con <- get_database_connection(read_only=TRUE)

dbListTables(con)

country <- run_query_df ('SELECT * FROM country;')


hba1c <- dbGetQuery(con, 'SELECT * FROM hba1c;')
classical_model_smr <- dbGetQuery(con, 'SELECT * FROM classical_model_smr;')

incidence_curve <- dbGetQuery(con, 'SELECT * FROM incidence_curve;')

#  saveRDS(incidence_curve,"article_markdown/incidence_curve.Rds")

incidence_growth_rate <- dbGetQuery(con, 'SELECT * FROM incidence_growth_rate;')

onset_death_rate <- dbGetQuery(con, 'SELECT * FROM onset_death_rate;')
onset_death_rate_single_year <- dbGetQuery(con, 'SELECT * FROM onset_death_rate_single_year;')

life_table <- dbGetQuery(con, 'SELECT * FROM life_table;')
population <- dbGetQuery(con, 'SELECT * FROM population;')
population_single_year <- dbGetQuery(con, 'SELECT * FROM population_single_year;')



prevalence_reference <- dbGetQuery(con, 'SELECT  * FROM prevalence_reference;')
weibull_complications <- dbGetQuery(con, 'SELECT * FROM weibull_complications;')
model_constants <- dbGetQuery(con, 'SELECT * FROM model_constants;')
disease_weights <- dbGetQuery(con, 'SELECT * FROM disease_weights;')
write.csv(disease_weights,"disease_weights.csv")
constant_complications <- dbGetQuery(con, 'SELECT * FROM constant_complications;')



countries <- get_countries_and_regions()


country_code = countries$world_bank_code[countries$world_bank_name=="Australia"]
loc_id       = countries$loc_id         [countries$world_bank_name=="Australia"]
country_name = countries$world_bank_name[countries$world_bank_name=="Australia"]
country_wb_name <- country_name
Log=function(fmt, ...) { cat(sprintf(paste0(fmt, '\n'), ...)) }
cache_file <- file.path(paste0(country_code, ".Rda"))



# refresh_one_country_file(country_name, 2000, 2040, Log, cache_file)

start_year  <-  2000
end_year    <-  2040
result <- run_model(
  country_name,
  start_year=start_year,
  end_year=end_year
)


value<- result$dalys_years
Year <- as.numeric(substring(names(value),1,4))
Age  <- as.numeric(substring(names(value),5,10))
deaths <- data.frame(Country=result$country, Type="Country", Year=Year,Age=Age, Value= value)
dalys_years <- sum(deaths$Value[deaths$Year==2021])
dalys_years


value<- result$pop_scale_factor
Year <- as.numeric(substring(names(value),1,4))
Age  <- as.numeric(substring(names(value),5,10))
background_population <- data.frame(Country=result$country, Type="Country", Year=Year,Age=Age, Value= value)
population <- sum(background_population$Value[background_population$Year==2021])
population

value<- result$pop
Year <- as.numeric(substring(names(value),1,4))
Age  <- as.numeric(substring(names(value),5,10))
background_population <- data.frame(Country=result$country, Type="Country", Year=Year,Age=Age, Value= value)
population2 <- sum(background_population$Value[background_population$Year==2021])
population2

value<- result$BD_flow
Year <- as.numeric(substring(names(value),1,4))
Age  <- as.numeric(substring(names(value),5,10))
deaths <- data.frame(Country=result$country, Type="Country", Year=Year,Age=Age, Value= value)
deaths1 <- sum(deaths$Value[deaths$Year==2021])
deaths1

value<- result$dt1d_years
Year <- as.numeric(substring(names(value),1,4))
Age  <- as.numeric(substring(names(value),5,10))
deaths <- data.frame(Country=result$country, Type="Country", Year=Year,Age=Age, Value= value)
dt1d_years <- sum(deaths$Value[deaths$Year==2021])
dt1d_years

value<- result$prev_years
Year <- as.numeric(substring(names(value),1,4))
Age  <- as.numeric(substring(names(value),5,10))
deaths <- data.frame(Country=result$country, Type="Country", Year=Year,Age=Age, Value= value)
prev_years <- sum(deaths$Value[deaths$Year==2021])
prev_years

value<- result$ddx_years
Year <- as.numeric(substring(names(value),1,4))
Age  <- as.numeric(substring(names(value),5,10))
deaths <- data.frame(Country=result$country, Type="Country", Year=Year,Age=Age, Value= value)
ddx_years <- sum(deaths$Value[deaths$Year==2021])
ddx_years








variables <- tribble(
  ~variable,                    ~Name,
  "prev_years",             "Prevalence",
  "prev_lt20_years",        "Prevalence <20 yr",
  "ghost_years",            "Ghosts",
  "ghost_years_ddx",        "Ghosts (onset death)",
  "ghost_years_dt1d",       "Ghosts (early death)",
  "incidence_years",        "Incidence",
  "ddx_years",              "Ann. onset deaths",
  "dt1d_years",             "Ann. early deaths",
  "dalys_years",            "Ann. days lost",
  "qB",                     "Ann. background deaths",
  "pop",                    "Ann. background population"
)













