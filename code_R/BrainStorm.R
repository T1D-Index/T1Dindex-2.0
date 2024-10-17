
# how to map incidence to stage 1 incidence 

library(dplyr)
library(echarts4r)


data.frame(
  x= ( sequence(11,0,10) ) , 
  y= c(1,1,20,40,10,10,10,10,10,10,10),
  x2= ( 100-(sqrt(100-sequence(11,0,10)))*10 ), 
  y2= c(1,1,20,40,10,10,10,10,10,10,10)
  ) %>%
  echarts4r::e_chart(x) %>% echarts4r::e_line(y) %>%
  echarts4r::e_chart(x2) %>% echarts4r::e_line(y2) %>%
  
  e_x_axis(name="Age") %>% e_y_axis(name="T1D Incidence (Stage 3 at 2024)")


data.frame(
  x= ( 100-(sqrt(100-sequence(11,0,10)))*10 ) 
  , y= c(1,1,20,40,10,10,10,10,10,10,10)
  ) %>% 
  echarts4r::e_chart(x) %>% echarts4r::e_line(y) %>%
  e_x_axis(name="Age") %>% e_y_axis(name="T1D Incidence (Stage 1 at 2014)")
