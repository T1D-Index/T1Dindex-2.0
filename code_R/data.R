# Functions for loading data files
#
# Confidential
# Copyright (c) JDRF 2020, All rights reserved
#
# This file contains functions for loading data from the filesystem. These
# functions are used for loading custom parameters, as well as by the data prep
# scripts.
#

#' Converts an incidence curve tibble to an array
#'
#' Also divides by 1e5 so that incidence is given as an annual probability
#'
#' @param tbl curves tibble loaded from incidence.xlsx
#' @return Nx100 matrix of incidence rates, with year as rownames
#' @noRd
#' @importFrom dplyr select
age_matrix <- function(tbl) {
  m <- select(tbl, .data$`0-4`:.data$`95-99`) %>%
    as.matrix(ncol = 20)
  UNBIN_IDXS <- sort(rep(1:20, 5))
  m100 <- matrix(m[, UNBIN_IDXS], ncol = 100)  # convert to 1-year age bins
  rownames(m100) <- tbl$year
  m100 * 1e-5
}

#' Obtains a (possibly read-only) database connection
#'
#' If the database is not available in ~/data.db then it uses the built-in test database.
#' Override the default location by setting the environment variable T1D_DATA_FILE.
#'
#' @param read_only if TRUE create r/o connection
#' @param use_test if TRUE use test db (always TRUE when under test)
#' @return DBI connection
#' @importFrom testthat is_testing
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbConnect
get_database_connection <- function(read_only=TRUE, use_test=FALSE) {
  db_file <- path.expand(Sys.getenv('T1D_DATA_FILE', unset = '~/data.db'))
  if (use_test | is_testing() | !file.exists(db_file)) {
    db_file <- system.file('data/test.db', package='T1DModel', mustWork = TRUE)
  }
  if (!file.exists(db_file)) {
    stop(sprintf('Database %s does not exist', db_file))
  }
  dbConnect(RSQLite::SQLite(), db_file, read_only=read_only)
}




run_query_df <- function(query)
{

  con <- get_database_connection(read_only=TRUE, use_test=FALSE)
  df  <- dbGetQuery(con,  query)
  dbDisconnect(con)

  return(df)
}


run_query_df_path <- function(query,db_file_path)
{
  # host_name <- "localhost"  # local
  con <-      connec <- dbConnect(RPostgres::Postgres(),  dbname = "t1d"
                                  , host = host_name
                                  , port = "5432",  user = "postgres",   password = "postgrest1d")
  df  <- dbGetQuery(con,  query)
  dbDisconnect(con)

  return(df)
}



#' Obtain location ID for country with given WB name
#'
#' @param country_wb_name World Bank name of country
#' @param use_test force use of test database
#' @return integer location ID
#' @export
#' @importFrom DBI dbDisconnect dbGetQuery

get_loc_id <- function(country_wb_name, use_test=FALSE) {
  con <- get_database_connection(read_only=TRUE, use_test=use_test)
  loc_id_df <- dbGetQuery(
    con,
    'SELECT loc_id FROM country WHERE world_bank_name = :wbname;',
    list(wbname=country_wb_name))
  dbDisconnect(con)
  if (nrow(loc_id_df) == 0) {
    stop(sprintf('Location ID for country %s not found', country_wb_name))
  }
  loc_id_df$loc_id
}

#' Obtain country name associated with given loc_id
#'
#' @param loc_id location ID for country
#' @param use_test force use of test database
#' @return integer location ID
#' @export
#' @importFrom DBI dbDisconnect
get_country_name <- function(loc_id, use_test=FALSE) {
  ctry <- get_country_details(loc_id, use_test=use_test)
  ctry$world_bank_name
}

#' Obtain details for a country by loc_id
#'
#' @param loc_id Location ID of country
#' @param use_test force use of test database
#' @return one-row data frame
#' @export
#' @importFrom DBI dbDisconnect dbGetQuery
get_country_details <- function(loc_id, use_test=FALSE) {
  con <- get_database_connection(read_only=TRUE, use_test=use_test)
  details <- dbGetQuery(
    con,
    'SELECT
         loc_id, world_bank_name, world_bank_code,
         jdrf_region_name, jdrf_broad_region_name,
         broad_region_id, od_region_id
    FROM country
    WHERE loc_id = :loc_id',
    list(loc_id=loc_id))
  dbDisconnect(con)
  if (nrow(details) == 0) {
    stop(sprintf('Details for location ID %d not found', loc_id))
  }
  details
}

#' Obtain list of countries
#'
#' @importFrom DBI dbGetQuery dbDisconnect
#' @export
get_countries <- function() {
  con <- get_database_connection(read_only=TRUE)
  details <- dbGetQuery(
    con, 'SELECT loc_id, world_bank_name, world_bank_code FROM country;')
  dbDisconnect(con)
  details
}

#' Obtain list of countries and regions
#'
#' @importFrom DBI dbGetQuery dbDisconnect
#' @export
get_countries_and_regions <- function() {
  con <- get_database_connection()
  details <- dbGetQuery(
    con,
    "SELECT
      loc_id,
      idf_country_name,
      world_bank_name,
      world_bank_code,
      jdrf_region_name,
      wd_region,
      wd_income_category,
      world_bank_classification,
      jdrf_broad_region_name,
      idf_reference_country
    FROM country
    --WHERE drop_country='FALSE';
    ")
  dbDisconnect(con)
  details
}

#' Get incidence curve data for location
#'
#' @param loc_id numeric location ID
#' @return data frame
#' @importFrom DBI dbGetQuery dbDisconnect
#' @export
get_incidence_curve <- function(loc_id) {
  con <- get_database_connection()
  curve_details <- dbGetQuery(
    con,
    'SELECT year, age_from, age_to, incidence/1e5 as incidence
    FROM incidence_curve
    WHERE loc_id = :loc_id
    ORDER BY year, age_from;',
    list(loc_id=loc_id)
  )
  dbDisconnect(con)
  curve_details %>%
    pivot_wider(id_cols='year', values_from='incidence', names_from='age_from') %>%
    as('matrix')
}

#' Get incidence growth rate for location
#'
#' @param loc_id numeric location ID
#' @return data frame
#' @importFrom DBI dbGetQuery dbDisconnect
#' @export
get_incidence_growth_rates <- function(loc_id) {
  con <- get_database_connection()
  gr <- dbGetQuery(
    con,
    'SELECT past_growth_rate, future_growth_rate
    FROM incidence_growth_rate
    WHERE loc_id = :loc_id',
    list(loc_id=loc_id)
  )
  dbDisconnect(con)
  list(
    past=gr$past_growth_rate,
    future=gr$future_growth_rate)
}

#' Get complication parameters
#'
#' @return list of weibull and constant parameters
#' @export
get_complication_parameters <- function() {
  # con <- get_database_connection()
  # const <- dbGetQuery(con,
  #   'SELECT complication, hba1c_lt_8, hba1c_8_to_9, hba1c_gteq_9
  #   FROM constant_complications;'
  # )
  # weib <- dbGetQuery(
  #   con,
  #   'SELECT abbrev, intercept, slope, scale, name, graph_title
  #   FROM weibull_complications;'
  # )
  # dbDisconnect(con)

  const <- read.csv("data_internal/constant_complications.csv")
  weib  <- read.csv("data_internal/weibull_complications.csv")
  list(constant=const, weibull=weib)
}

#' Get disease weights
#'
#' @return data frame
#' @export
#' @importFrom DBI dbDisconnect dbGetQuery
get_disease_weights <- function() {
  # con <- get_database_connection()
  # dw <- dbGetQuery(con, 'SELECT complication, DALYs FROM disease_weights;')
  # dbDisconnect(con)

  dw  <- read.csv("data_internal/disease_weights.csv")

  wts <- dw$DALYs
  names(wts) <- dw$complication



  as.list(wts)
}

#' Get assumed mean HbA1cs
#'
#' @param loc_id location ID
#' @return data frame of year and assumed HbA1c
#' @export
#' @importFrom DBI dbDisconnect
get_hba1c_assumption <- function(loc_id) {
  con <- get_database_connection()
  hb <- dbGetQuery(
    con,
    'SELECT year, hba1c
    FROM hba1c
    WHERE loc_id = :loc_id
    ORDER BY year;',
    list(loc_id=loc_id)
  )
  dbDisconnect(con)
  hb
}

#' Get prevalence reference point, if any
#'
#' @param loc_id country location ID
#' @return a structure containing year and prevalence point estimate,
#'     if available, or NULL otherwise
#' @export
get_prevalence_reference <- function(loc_id) {
  con <- get_database_connection()
  refs <- dbGetQuery(
    con,
    'SELECT year, prev
    FROM prevalence_reference
    WHERE loc_id = :loc_id;',
    list(loc_id=loc_id)
  )
  dbDisconnect(con)
  if (nrow(refs) > 0) {
    list(year=refs$year[1], prev=refs$prev[1])
  } else {
    NULL
  }
}
