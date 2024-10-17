# Incidence calculation functions
#
# Confidential
# Copyright (c) JDRF and affiliates 2019-2020. All rights reserved.
#
# Our incidence curves are constructed as follows:
#   - Incidence up to age 20 is taken from the IDF atlas
#   - Adult incidence profiles are constructed using the Diaz-Valencia paper
#   - Growth rates based on income are applied to incidence rates, expressed
#     as odds ratios
#
# See data-raw/10_incidence_curve_data.R and data-raw/70_country_parameters.R
# for data sources.
#
# A current limitation is that only data from a single study is used.

#' Construct incidence curve by age
#'
#' This is a higher-order function that returns a function of year. The
#' incidence function returns a curve as a numeric array of length 100,
#' corresponding to incidence rates for ages 0-99.
#'
#' \code{country_wb_name} should be the world bank name for that country.
#' \code{past_growth_rate} and \code{future_growth_rate} are given as a ratio,
#' e.g. 0.05 for 5% annual growth. If the growth rate is not specified, we take
#' the default country parameter for the given country.
#'
#' Note: odds ratio grows at the given rate, not the incidence rate itself.
#'
#' @param country_wb_name (scalar character) World Bank name for country
#' @return function of year that returns vector of length 100, representing
#'   incidence rates for ages 0-99
#'
#' @examples
#' # Obtain incidence curves for India
#' # ind_inc is a function of year that returns incidence for ages 0-99.
#' ind_inc <- incidence_function('India')
#'
#' # Plot 2020 incidence curve for India
#' plot(0:99, 1e5*ind_inc(2020), type='l',
#'      main='India 2020 Incidence by age',
#'      ylab='Annual cases per 100,000')
#'
#' # Plot age 5 incidence curve for India, 1990 - 2020
#' inc_age_5 <- Vectorize(function(y) 1e5*ind_inc(y)[5+1])
#' plot(inc_age_5, from=1990, to=2020,
#'      main='Incidence for children aged 5, India',
#'      ylab='Annual cases per 100,000')
#' @importFrom dplyr filter
#' @importFrom stats approxfun
#' @export
incidence_function <- function(country_wb_name,years) {

  loc_id   <- run_query_df (paste0("SELECT loc_id FROM country WHERE world_bank_name = '",country_wb_name,"'" ) )$loc_id
  inc_data <- run_query_df (paste0("SELECT year, age, incidence_rate/1e5 as incidence
                                    FROM incidence_curve
                                    WHERE loc_id = ",loc_id,"
                                    ORDER BY year;" ) )
  inc_data <- inc_data %>%
              pivot_wider(id_cols='year', values_from='incidence', names_from='age') %>%
              as('matrix')


  growth_rates   <- run_query_df (paste0("SELECT past_growth_rate, future_growth_rate
    FROM incidence_growth_rate
    WHERE loc_id =" , loc_id ) )

  growth_rates <- list(
    past=growth_rates$past_growth_rate,
    future=growth_rates$future_growth_rate)



#
#   loc_id <- get_loc_id(country_wb_name)
#   inc_data <- get_incidence_curve(loc_id)
#   growth_rates <- get_incidence_growth_rates(loc_id)




  sample_range <- range(inc_data[,1])
  # incidence is passed to the returned function as a list of incidence functions,
  # one for each age
  if (sample_range[1] == sample_range[2]) {
    # single curve available - constant function for only inc rate available
    inc_fs <- lapply(AGES, function(age) {function(years) inc_data[floor(age/5)+2]})
  } else {
    # multiple curves available - interpolate incidence for each age
    inc_fs <- lapply(AGES, function(age) {
      approxfun(inc_data[, 1], inc_data[, floor(age/5)+2])
    })
  }
  # function of year & draw that either interpolates or applies growth rate,
  # returning a vector of length 100 for ages 0..99.
  draw <- 1
  stopifnot(is.numeric(years))

  result <- matrix(NA, nrow=length(years), ncol=MAX_AGE,
                   dimnames=list(years, AGES))


  for (i in seq_along(years)) {
    yr <- years[i]
    if (sample_range[1] < yr && yr < sample_range[2]) {
      # strictly within sample: interpolate
      result[i,] <- as.vector(sapply(seq_len(MAX_AGE), function(i) inc_fs[[i]](yr)))
    } else {
      # before/after sample - apply growth rate
      gr <- ifelse(yr <= sample_range[1], growth_rates$past, growth_rates$future)
      # print(gr)
      base_year <- sample_range[1+(yr > sample_range[1])]
      # print(base_year)

      curve <- as.vector(inc_data[-1])
      # print(base_year)

      result[i,] <- curve * (1 + gr) ^ (yr - base_year)
    }
  }
  result

}
