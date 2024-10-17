# Functions for plotting model outputs
#
# Confidential
# Copyright (c) JDRF 2020, All rights reserved
#
# These plots are used by various model outputs, including the paper and shiny
# dashboard. They will eventually be superseded by the web site.

color_order <- c(
  7,  9, 3, 13, 2,
  18, 6, 15, 10, 5, 12, 14,
  8, 16, 17, 4, 11, 1, 19)
light_colors <- c(
  '#EF9A9A',
  '#F48FB1',
  '#CE93D8',
  '#B39DDB',
  '#9FA8DA',
  '#90CAF9',
  '#81D4FA',
  '#80DEEA',
  '#80CBC4',
  '#80DEEA',
  '#80CBC4',
  '#A5D6A7',
  '#C5E1A5',
  '#E6EE9C',
  '#FFF59D',
  '#FFE082',
  '#FFCC80',
  '#FFAB91',
  '#BCAAA4',
  '#EEEEEE',
  '#B0BEC5'
)
mid_colors <- c(
  '#F44336', # 1. red
  '#E91E63', # 2. pink
  '#9C27B0', # 3. violet
  '#673AB7', # 4. purple
  '#3F51B5', # 5. navy
  '#2196F3', # 6. blue
  '#03A9F4', # 7. light blue
  '#00BCD4', # 8. teal
  '#009688', # 9. dark green
  '#4CAF50', # 10. sea green
  '#8BC34A', # 11. apple
  '#CDDC39', # 12. lemon ginger
  '#FFC107', # 13. sea buckthorn
  '#FF9800', # 14. burnt orange
  '#FF5722', # 15. orange
  '#FF5722', # 16. red-orange
  '#795548', # 17. brown
  '#9E9E9E', # 18. gray
  '#607D8B'  # 19. bluemetal
)

dark_colors <- c(
  '#C62828', # 1. red
  '#AD1457', # 2. pink
  '#6A1B9A', # 3. violet
  '#4527A0', # 4. purple
  '#283593', # 5. navy
  '#1565C0', # 6. blue
  '#0277BD', # 7. light blue
  '#00838F', # 8. teal
  '#00695C', # 9. dark green
  '#2E7D32', # 10. sea green
  '#558B2F', # 11. apple
  '#9E9D24', # 12. lemon ginger
  '#F9A825', # 13. sea buckthorn
  '#FF8F00', # 14. burnt orange
  '#EF6C00', # 15. orange
  '#D84315', # 16. red-orange
  '#4E342E', # 17. brown
  '#424242', # 18. gray
  '#37474F'  # 19. bluemetal
)

#' JDRF ggplot outline colors
#'
#' @importFrom ggplot2 scale_color_manual
#' @export
theme_color <- function() {
  scale_color_manual(values=mid_colors[color_order])
}

#' ggplot fill colors
#'
#' @importFrom ggplot2 scale_fill_manual
#' @export
theme_fill <- function() {
  set.seed(1)
  scale_fill_manual(values=mid_colors[color_order])
}

#' Simple theme to use for plotting outputs
#'
#' @param size relative size of text in theme
#' @importFrom ggplot2 theme theme_minimal element_blank element_line unit element_text
#' @export
plot_theme <- function(size=15) {
  theme_minimal(size) +
    theme(legend.position = 'bottom',
          panel.grid.major.y = element_line(color = "#DDDDDD", size=0.5),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_text(vjust=-1),
          axis.ticks.length.x.bottom = unit(3, 'pt'),
          axis.ticks = element_line(color = '#DDDDDD', size=0.5))
}

#' Simple theme to use for plotting horizontal graphs
#'
#' @param size relative size of text in theme
#' @importFrom ggplot2 theme theme_minimal element_blank element_line unit element_text
#' @export
horizontal_plot_theme <- function(size=15) {
  theme_minimal(size) +
    theme(legend.position = 'bottom',
          panel.grid.major.x = element_line(color = "#DDDDDD", size=0.5),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
    )
}

#' Override default line to make much thicker
#'
#' @noRd
#' @param size relative line thickness
#' @param ... other parameters to pass to `ggplot2::geom_line()`
geom_line <- function(size=1.5, ...) {
  ggplot2::geom_line(size=size, ...)
}

#' Single lines should be colored
#'
#' @param size relative line thickness
#' @param ... other parameters to pass to `ggplot2::geom_line()`
#' @export
geom_single_line <- function(size=1.5, ...) {
  ggplot2::geom_line(size=size,
                     color=mid_colors[color_order[1]], ...)
}

#' Single columns should be colored
#'
#' @param ... other parameters to pass to `ggplot2::geom_col()`
#' @export
geom_single_col <- function(...) {
  ggplot2::geom_col(fill=mid_colors[color_order[1]], ...)
}

#' ggplot scales that include units on top line
#'
#' This is a specialization of scale_y_units
#'
#' @param units description of units
#' @param fmt sprintf format to apply to labels
#' @param sep separator between value and units on top line
#' @param ... parameters to pass to scale_y_units
#' @export
scale_y_units <- function(units, fmt=NULL, sep='\n', ...) {
  labelfun <- function(breaks) {
    non_na_breaks <- breaks[!is.na(breaks)]
    idx <- match(non_na_breaks[length(non_na_breaks)], breaks)
    if (is.null(fmt)) {
      fmt_breaks <- as.character(breaks)
    }
    else {
      fmt_breaks <- sprintf(fmt, breaks)
    }
    fmt_breaks[idx] <- paste0(fmt_breaks[idx], sep, units)
    fmt_breaks
  }
  scale_y_continuous(labels = labelfun, ...)
}

#' Summarize incidence function
#'
#' Calculates a simple average for broad age groups, 0-20, 20-40, and 40+ years
#'
#' @param prev Prevalence, output of \link{prevalence_and_ghost_pop}()
#' @param inc_fn Prevalence function
#' @return ggplot object
#' @importFrom dplyr bind_cols group_by summarize mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line ylab labs
#' @export
plot_incidence_timeseries <- function(prev, inc_fn) {
  inc_matrix <- 1e5 * t(mapply(inc_fn, prev$years))
  colnames(inc_matrix) <- AGES
  to_plot <-
    bind_cols(year=prev$years, as_tibble(inc_matrix)) %>%
    pivot_longer(-.data$year, names_to='age', values_to='inc') %>%
    mutate(age_group = cut(as.integer(.data$age), right = FALSE,
                           breaks=c(0, 20, 60, Inf),
                           labels=c('< 20 years', '20-40 years', '60+ years'))) %>%
    group_by(.data$year, .data$age_group) %>%
    summarize(inc = mean(.data$inc))
  ggplot(to_plot, aes(x=.data$year, y=.data$inc, color=.data$age_group)) +
    geom_line() +
    theme_color() +
    scale_y_units('cases per\n100,000')+
    labs(title=paste0(prev$country, ' - T1D incidence rates by age'),
         color='Age', x=NULL, y=NULL) +
    plot_theme()
}

#' Plot incidence cross section
#'
#' @param prev Prevalence, output of \link{prevalence_and_ghost_pop}()
#' @param inc_fn Prevalence function
#' @param view_year data will be subset to this year
#' @return ggplot2 object
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_line ylab scale_y_continuous labs
#' @export
plot_incidence_xsection <- function(prev, inc_fn, view_year) {
  to_plot <- tibble(
    age = seq(MAX_AGE)-1,
    incidence = 1e5 * inc_fn(view_year)
  )

  ggplot(to_plot, aes(x=.data$age+0.5, y=.data$incidence)) +
    geom_single_line() +
    scale_y_units('cases per\n100,000', limits = c(0, NA)) +
    labs(title=paste0(prev$country, ' - T1D incidence rates in ', view_year),
         color='Age', y=NULL, x='Age') +
    plot_theme()
}

#' Plot incidence multi cross section
#'
#' @param prev prevalence object, output of \link{prevalence_and_ghost_pop}()
#' @param inc_fn incidence callback, a function of year
#' @param years vector of years to show cross-sections for
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_line ylab scale_y_continuous labs
#' @importFrom dplyr bind_rows
#' @export
plot_incidence_multi_xsection <- function(prev, inc_fn, years) {
  to_plot <- NULL
  for (year in years) {
    to_plot <- bind_rows(
      to_plot,
      tibble(
        age = AGES,
        incidence = 1e5 * inc_fn(year),
        year = year
      )
    )
  }
  ggplot(to_plot, aes(x=.data$age + 0.5, y=.data$incidence, color=factor(.data$year))) +
    geom_line() +
    scale_y_units('cases per\n100,000', limits = c(0, NA)) +
    labs(title=paste0(prev$country, ' - T1D incidence rates by age'),
         color=NULL, y=NULL, x='Age') +
    theme_color() +
    plot_theme()
}

#' Plot death on diagnosis rates
#'
#' @param prev prevalence object, output of \link{prevalence_and_ghost_pop}()
#' @param ddx_fn death at onset function
#' @return a ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_line ylab labs
#' @importFrom tibble tibble
#' @export
plot_ddx_timeseries <- function(prev, ddx_fn) {
  ddx <- 100 * apply(t(mapply(ddx_fn, prev$years)), 1, mean)
  to_plot <- tibble(year=prev$years, ddx=ddx)
  ggplot(to_plot, aes(x=.data$year, y=.data$ddx)) +
    geom_single_line() +
    scale_y_units('%') +
    labs(
      title=paste0(prev$country, ' - Death on Incidence Rates'),
      x=NULL, y=NULL) +
    plot_theme()
}

#' Plot average hba1c levels over time
#'
#' @param prev prevalence object, output of \link{prevalence_and_ghost_pop}()
#' @param hba1c_fn hba1c callback, a function of year
#' @return a ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_line ylab labs
#' @importFrom tibble tibble
#' @export
plot_hba1c_timeseries <- function(prev, hba1c_fn) {
  hba1c <- apply(t(mapply(hba1c_fn, prev$years)), 1, mean)
  to_plot <- tibble(year=prev$years, hba1c=hba1c)
  ggplot(to_plot, aes(x=.data$year, y=.data$hba1c)) +
    geom_single_line() +
    scale_y_units('%') +
    labs(
      title=paste0(prev$country, ' - Mean Glycated Hemoglobin Concentration (HbA1c)'),
      y=NULL) +
    theme_color() +
    plot_theme()
}

#' Plot prevalence cross section for a given year
#'
#' @param prevdata output of \link{prevalence_and_ghost_pop}()
#' @param view_year year of cross section to plot
#' @param binwidth width of bins to plot
#' @return ggplot object
#' @export
#' @examples
#' prev <- prevalence_and_ghost_pop('India')
#' plot_cross_section(prev, 2020, binwidth = 1)
#' @importFrom dplyr filter mutate group_by summarize_at select recode vars ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_col ggtitle xlab ylab aes
#' @importFrom tidyselect ends_with
plot_cross_section <- function(prevdata, view_year, binwidth=5) {
  longdata <- with(prevdata, matrixes_to_long_format(
    P_level=P_level, ghost_ddx_level=ghost_ddx_level, ghost_hba1c_level=ghost_hba1c_level,
    pop=pop, years=years))
  d <- filter(longdata, .data$year == view_year)

  ctry_name <- prevdata$country
  prev <- with(d, sum(P_level))
  pop <- with(d, sum(pop))
  title <- sprintf('%s - T1D Prevalence in %d', ctry_name, view_year)
  caption <- sprintf('Estimated total prevalence is %.0fk (%.0f per 100k)',
                      prev/1e3, prev/pop*1e5)
  breaks <- seq(0, 100, binwidth)
  labels <- paste0(breaks[-length(breaks)],'-',breaks[-1])

  if (binwidth == 1) {
    toplot <- d
  } else {
    toplot <- d %>%
      mutate(age = cut(.data$age, breaks=breaks, labels=labels, right=FALSE)) %>%
      group_by(.data$age, .data$year) %>%
      summarize_at(vars(ends_with('level')), ~sum(.)) %>%
      ungroup()
  }

  toplot %>%
    select(.data$age, .data$year, .data$P_level, .data$ghost_ddx_level, .data$ghost_hba1c_level) %>%
    pivot_longer(.data$P_level:.data$ghost_hba1c_level, names_to='group', values_to='prev') %>%
    mutate(group = recode(.data$group,
                          P_level='Living',
                          ghost_ddx_level='Ghost (death on onset)',
                          ghost_hba1c_level='Ghost (inadequate care)')) %>%
    ggplot(aes(x=.data$age, y=.data$prev*1e-3, fill=.data$group)) +
    geom_col() +
    scale_y_units("'000\ncases") +
    labs(
      title=title,
      caption=caption,
      fill=NULL, x='Age', y=NULL) +
    theme_fill() +
    plot_theme()
}

#' Multiple year cross-section
#'
#' @param prevdata prevalence data, output of a call to
#'   \link{prevalence_and_ghost_pop}()
#' @param view_years years of data to show
#' @return ggplot object
#' @export
#' @examples
#' prevalence_and_ghost_pop('India') %>%
#'   plot_multi_cross_section(view_years=c(2000, 2010, 2020))
#' @importFrom dplyr filter group_by summarize mutate inner_join select
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_col facet_wrap aes
plot_multi_cross_section <- function(prevdata, view_years=c(2000,2010,2020,2030)) {
  title <- sprintf('%s - T1D Prevalence', prevdata$country)
  longdata <- with(prevdata, matrixes_to_long_format(
    P_level=P_level, ghost_level=ghost_level,
    pop=pop, years=years))

  labels <- longdata %>%
    filter(.data$year %in% view_years) %>%
    group_by(.data$year) %>%
    summarize(pop=sum(.data$pop), prev=sum(.data$P_level)) %>%
    mutate(
      prev_100k=round(.data$prev/.data$pop*1e5),
      label = sprintf('%d (prevalence %.0fk, %.0f per 100k)', .data$year, .data$prev*1e-3, .data$prev_100k))

  longdata %>%
    inner_join(labels, by='year') %>%
    select(.data$year, .data$age, .data$P_level, .data$ghost_level, .data$label) %>%
    pivot_longer(.data$P_level:.data$ghost_level, names_to='group', values_to='prevalence') %>%
    mutate(group = recode(.data$group,
                          P_level='T1D prevalence',
                          ghost_level='Ghost population')) %>%
    ggplot(aes(x=.data$age, y=.data$prevalence*1e-3, fill=.data$group)) +
    geom_col() +
    facet_wrap(~.data$label, ncol=1) +
    scale_y_units("'000\ncases") +
    labs(title=title, fill=NULL, x='Age', y=NULL) +
    theme_fill() +
    plot_theme()
}

#' Plot prevalence by age at onset
#'
#' This graph overplots a single line per age at onset
#'
#' @param prevdata output of \link{prevalence_and_ghost_pop}()
#' @param xsec_year year of cross section to plot
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab labs scale_color_gradient
#' @export
plot_prevalence_by_age_at_onset_and_age <- function(prevdata, xsec_year) {
  data <- NULL
  for (aoo in seq_len(MAX_AGE)) {
    data <- bind_rows(
      data,
      tibble(
        age=aoo:MAX_AGE,
        age_of_onset=rep(aoo, MAX_AGE - aoo + 1),
        cohort_prev = prevdata$P_cohorts_level['2020',aoo:MAX_AGE,aoo])
    )
  }
  ggplot(data, aes(x=.data$age, y=.data$cohort_prev, color=.data$age_of_onset, group=.data$age_of_onset)) +
    geom_line() +
    scale_y_units('persons') +
    labs(
      title=sprintf('%s - T1DM prevalence - %d', prevdata$country, xsec_year),
      subtitle=sprintf('Contributions to prevalence by age at onset'),
      color='Age at T1DM onset',
      group='Age at T1DM onset',
      y=NULL, x=sprintf('Age in %d', xsec_year)) +
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    plot_theme()
}

#' Prevalence cross-section by age at onset
#'
#' This graph looks a little bit like the incidence curve,
#' even though each point is the sum of many years of
#' incidence.
#'
#' @param prevdata prevalence data
#' @param xsec_year cross-section year
#' @return a ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes labs
#' @importFrom tibble tibble
plot_prevalence_by_age_at_onset <- function(prevdata, xsec_year) {
  age_cohort <- prevdata$P_cohorts_level[as.character(xsec_year),,]
  data <- tibble(
    age=1:MAX_AGE,
    inc=apply(age_cohort, 2, sum)
  )
  ggplot(data, aes(x=.data$age, y=1e-3*.data$inc)) +
    geom_single_col() +
    scale_y_units("'000\ncases") +
    labs(
      title=sprintf('%s - T1DM prevalence in %d',
                    prevdata$country, xsec_year),
      x='Age at T1DM onset', y=NULL) +
    plot_theme()
}

#' Plot prevalence timeseries
#'
#' @param prevdata prevalence data, output of a call to
#'   \link{prevalence_and_ghost_pop}()
#' @param age_breaks where to divide the age distribution; use Inf to specify
#'   an unbounded upper interval
#' @return ggplot object
#' @examples
#' prevalence_and_ghost_pop('Australia') %>%
#'   plot_prevalence_timeseries
#' @importFrom dplyr mutate group_by summarize
#' @importFrom ggplot2 ggplot geom_col ggtitle ylab aes position_stack labs
#' @importFrom stringr str_replace
#' @export
plot_prevalence_timeseries <- function(prevdata, age_breaks=c(0,20,40,60,Inf)) {
  labels <- paste0(age_breaks[-length(age_breaks)],'-',age_breaks[-1]) %>%
    str_replace('-Inf', '+')
  longdata <- with(prevdata, matrixes_to_long_format(
    P_level=P_level, ghost_level=ghost_level,
    pop=pop, years=years))
  title <- sprintf('%s - T1D Prevalence', prevdata$country)
  longdata %>%
    mutate(age_group = cut(.data$age, breaks=age_breaks, labels=labels, right=FALSE)) %>%
    group_by(.data$year, .data$age_group) %>%
    summarize(
      prevalence = 1e-3*sum(.data$P_level),
      max_age = max(.data$age)
    ) %>%
    ggplot(aes(x=.data$year, y=.data$prevalence)) +
    geom_col(aes(fill=.data$age_group), position=position_stack(reverse=TRUE)) +
    scale_y_units("'000\ncases") +
    labs(
      title=title,
      fill='Age',
      x=NULL, y=NULL) +
    theme_fill() +
    plot_theme()
}

#' Complication prevalence as stacked bar chart
#'
#' @param prev prevalence object, output of \link{prevalence_and_ghost_pop}()
#' @param comp complications object, output of \link{complication_prevalence}()
#' @param long_names if TRUE use long names not abbreviations in legend
#' @return ggplot2 object
#' @importFrom dplyr bind_cols filter
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_col labs ylab guide_legend guides
#' @export
plot_complication_prevalence_timeseries <- function(prev, comp, long_names=FALSE) {
  to_plot <- bind_cols(
      year=as.integer(rownames(comp$year_comp)),
      as_tibble(comp$year_comp)
    ) %>%
    pivot_longer(
      -.data$year,
      names_to='complication',
      values_to='incidence') %>%
    filter(!(.data$complication %in% c('RF_dialysis', 'RF_transplant')))

  name_map <- tibble::tribble(
    ~complication,                              ~name, ~const_haz,
    "DSP",      "Distal symm. polyneuropathy",  TRUE,
    "PR",        "Proliferative retinopathy",  TRUE,
    "BL",                        "Blindness",  TRUE,
    "ON",                "Overt nephropathy",  TRUE,
    "RF",                    "Renal failure",  TRUE,
    "UoA",         "Foot ulcer or amputation",  TRUE,
    "HoM", "Hypertension or microalbuminuria",  TRUE,
    "nfMI",                     "Non-fatal MI",  TRUE,
    "nfCBVD",                 "Non-fatal stroke",  TRUE,
    "Hypo",                    "Hypoglycaemia",  FALSE,
    "DKA",            "Diabetic Ketoacidosis",  FALSE
  )
  lg <- guide_legend(ncol = 11)

  if (long_names) {
    to_plot <- to_plot %>%
      inner_join(name_map, by='complication') %>%
      mutate(complication = .data$name)
    lg <- guide_legend(ncol = 3)
  }

  title <- paste0(prev$country, ' - prevalence of T1D complications*')

  ggplot(to_plot, aes(x=.data$year, y=.data$incidence*1e-3, fill=.data$complication)) +
      geom_col() +
      scale_y_units("'000") +
      labs(title=title,
           caption='* Individuals may be counted multiple times',
           fill=NULL, x=NULL, y=NULL) +
      guides(fill=lg) +
      theme_fill() +
      plot_theme()
}

#' Complication prevalence by age
#'
#' @param prev prevalence object, output of \link{prevalence_and_ghost_pop}()
#' @param comp complications object, output of \link{complication_prevalence}()
#' @param view_year data will be subset to this year
#' @return ggplot2 object
#' @export
plot_complication_prevalence_xsection <- function(
    prev,
    comp,
    view_year
  ) {
  year_subset_m <- comp$year_comp_age[as.character(view_year),,]
  to_plot <- bind_cols(
      age=seq(MAX_AGE),
      as_tibble(t(year_subset_m))
    ) %>%
    pivot_longer(
      -.data$age,
      names_to='complication',
      values_to='incidence') %>%
    filter(!(.data$complication %in% c('RF_dialysis', 'RF_transplant')))

  title <- paste0(prev$country, ' - T1D Complications in ',view_year)

  ggplot(to_plot, aes(x=.data$age, y=.data$incidence*1e-3)) +
    geom_single_col() +
    facet_wrap(~.data$complication) +
    scale_y_units("'000") +
    labs(title=title, y=NULL) +
    plot_theme()
}

#' Compare complication prevalence levels between two years
#'
#' @param start_year integer like 2000, start of comparison
#' @param end_year integer like 2040, end of comparison
#' @param comp complication estimates, output of \link{complication_prevalence}()
#' @param country character name of country like 'India'
#' @export
#' @importFrom tibble tribble
#' @importFrom dplyr filter group_by summarize_if inner_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_col labs ylab xlab coord_flip
#' @importFrom forcats fct_rev
plot_complication_comparison <- function(start_year, end_year, comp, country) {
  name_map <- tribble(
    ~complication,                              ~name, ~const_haz,
    "DSP",      "Distal symmetric\npolyneuropathy",  TRUE,
    "PR",        "Proliferative\nretinopathy",  TRUE,
    "BL",                        "Blindness",  TRUE,
    "ON",                "Overt nephropathy",  TRUE,
    "RF",                    "Renal failure",  TRUE,
    "UoA",         "Foot ulcer or\namputation",  TRUE,
    "HoM", "Hypertension or\nmicroalbuminuria",  TRUE,
    "nfMI",                     "Non-fatal MI",  TRUE,
    "nfCBVD",                 "Non-fatal stroke",  TRUE,
    "AN",                 "Autonomic neuropathy",  TRUE,
    "Hypo",                    "Hypoglycaemia",  FALSE,
    "DKA",            "Diabetic Ketoacidosis",  FALSE
  )

  tidy_complications(comp) %>%
    filter(.data$year %in% c(start_year, end_year)) %>%
    group_by(.data$year) %>%
    summarize_if(is.numeric, ~round(1e-3 * sum(.), 2)) %>%
    pivot_longer(-.data$year, names_to='complication') %>%
    filter(.data$complication != 'age') %>%
    inner_join(name_map, by='complication') %>%
    mutate(year = fct_rev(ordered(.data$year))) %>%
    ggplot(aes(x=.data$name, y=.data$value, fill=.data$year)) +
    geom_col(position='dodge') +
    labs(
      title=sprintf('%s - T1DM complication prevalence', country),
      fill=NULL, x=NULL, y="Prevalence ('000)") +
    coord_flip() +
    theme_fill() +
    horizontal_plot_theme()
}

#' Plot disease burden timeseries
#'
#' @param prev Prevalence, results of \link{prevalence_and_ghost_pop}()
#' @param dalys Burden, results of \link{calculate_dalys}()
#' @return ggplot object
#' @export
#' @importFrom dplyr group_by summarize
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_col labs ylab xlab scale_fill_manual
plot_burden_timeseries <- function(prev, dalys) {
  plot_levels <- c("Due to disablity", "Death at onset", "Excess T1D deaths")
  to_plot <- tidy_dalys(prev, dalys) %>%
    group_by(.data$year) %>%
    summarize(
      `Due to disablity`=sum(.data$disability),
      `Death at onset`=sum(.data$ddx),
      `Excess T1D deaths`=sum(.data$excess)
    ) %>%
    pivot_longer(-.data$year, names_to='Cause', values_to='dalys') %>%
    mutate(Cause = ordered(.data$Cause, levels=plot_levels))

  ggplot(to_plot, aes(x=.data$year+.5, y=.data$dalys*1e-3, fill=.data$Cause)) +
    geom_col() +
    scale_y_units("'000\nDALYs") +
    labs(title=sprintf('%s - T1D Disease Burden', prev$country),
         subtitle="Disability-adjusted life years",
         fill=NULL, y=NULL, x=NULL) +
    theme_fill() +
    plot_theme()
}

#' Plot disease burden cross-section
#'
#' @param prev Prevalence, results of \link{prevalence_and_ghost_pop}()
#' @param dalys Burden, results of \link{calculate_dalys}()
#' @param show_year which year cross-section to show
#' @return ggplot object
#' @export
#' @importFrom dplyr group_by summarize
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_col labs ylab xlab scale_fill_manual
#' @importFrom rlang .data
plot_burden_cross_section <- function(prev, dalys, show_year) {
  plot_levels <- c("Due to disablity", "Death at onset", "Excess T1D deaths")
  to_plot <- tidy_dalys(prev, dalys) %>%
    filter(.data$year == show_year) %>%
    transmute(
      .data$age,
      `Due to disablity`=.data$disability,
      `Death at onset`=.data$ddx,
      `Excess T1D deaths`=.data$excess
    ) %>%
    pivot_longer(-.data$age, names_to='Cause', values_to='dalys') %>%
    mutate(Cause = ordered(.data$Cause, levels=plot_levels))

  ggplot(to_plot, aes(x=.data$age+.5, y=1e-3*.data$dalys, fill=.data$Cause)) +
    geom_col() +
    scale_y_units("'000\nDALYs") +
    labs(title=sprintf('%s - T1D Disease Burden in %d', prev$country, show_year),
         subtitle="Disability-adjusted life years",
         fill=NULL, y=NULL, x='Age') +
    theme_fill() +
    plot_theme()
}
