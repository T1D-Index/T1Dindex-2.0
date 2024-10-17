library(echarts4r)
for(i in 1:10)
{
  country_wb_name <- "Mali"
  prevalence_year   <-  2021
  prevalence_number <-  1000

  intervention_year_start <- 2005
  intervention_type <- "diagnosis_rate"

  print(Sys.time())

  data_long   <- readRDS(paste0(data_dir,"/",get_loc_id(country_wb_name),".Rds"))

  data_plot  <- data_long[data_long$age==24 & data_long$year <= prevalence_year& data_long$year >= intervention_year_start ,]
  data_plot  <- data.frame(x= data_plot$year,y= 1-data_plot$mortality_undiagnosed_rate)
  data_plot$rate  <- NA
  data_plot$y_new <- NA
  for(i in 2:nrow(data_plot))
  {
    data_plot$rate[i] <- (data_plot$y[i] - data_plot$y[i-1])/data_plot$y[i-1]
  }
  data_plot$rate_new <- data_plot$rate * 3
  data_plot$y_new    <- data_plot$y
  for(i in 2:nrow(data_plot))
  {
    data_plot$y_new[i] <- data_plot$y_new[i-1] * (1 +  data_plot$rate_new[i])
  }


  data_plot %>%
    e_chart(x)%>%
    e_line( y, name='diagnosis_rate_base')%>%
    e_line( y_new, name='diagnosis_rate_new')%>%
    e_x_axis(min=intervention_year_start, max= prevalence_year)%>%
    e_y_axis(min=0.4, max= 1)


  for(i in 0:24)
  { data_long[data_long$age==i & data_long$year <= prevalence_year& data_long$year >= intervention_year_start ,]$mortality_undiagnosed_rate <- (1-  data_plot$y_new) }


  prev        <- prevalence_and_ghost_pop(country_wb_name="Mali",start_year=1960,end_year=2040,data_long=data_long)
  prev        <- prev$prev_merge_long



  data_plot  <- setDT(prev[prev$Value_type=="Prevalence",])[,list(y = sum(Value) ),by="Year"]

  data_plot  <- data.frame(x= data_plot$Year,y=data_plot$y)
  data_plot %>%
    e_chart(x)%>%
    e_line(  y, name='Prevalence')%>%
    e_x_axis(min=intervention_year_start, max= prevalence_year)





}
