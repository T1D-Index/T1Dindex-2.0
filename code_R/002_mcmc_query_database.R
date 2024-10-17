library(data.table)
library(arrow)
library(RPostgreSQL)
library(dplyr)

connec <- dbConnect(dbDriver("PostgreSQL"),  dbname = "t1d", host = "database-t1d-dev-2.ccasxjjwcbmy.ap-southeast-2.rds.amazonaws.com",
                    port = "5432",  user = "postgres",   password = "postgres")


query_compose <-function(year,age_max,country)
{  
  query <- 
          paste0('SELECT "Country", 
          ROUND( (percentile_cont(0.025) WITHIN GROUP (ORDER BY "Prevalence"))::numeric,0) as Prevalence_lower ,
          ROUND( (percentile_cont(0.5)   WITHIN GROUP (ORDER BY "Prevalence"))::numeric,0) as Prevalence_median ,
          ROUND( (percentile_cont(0.975) WITHIN GROUP (ORDER BY "Prevalence"))::numeric,0) as Prevalence_upper 
          FROM
          (
          SELECT  "Country" , "run",SUM("Prevalence") as "Prevalence"
          FROM public.main_ci_',year,'
          WHERE "Age" <= ',age_max,'
              and "Country" = \'',country,'\'
          group by "Country","run"
          ) AS apple
          GROUP BY "Country" ')

return(query)
}

country    <- read.csv('data_internal/country.csv',stringsAsFactors = F,encoding='UTF-8')


table1 <- read.csv("paper_submission/plots_paper_revision/Appendix_tables/appendix6-table1.csv")




for(i in 1:nrow(table1))
{  # i <- 1
  print(i)
  if(table1$Country[i]!="")
  {
    country_name <- table1$Country[i]
    country_name <- gsub("Taiwan","China, Taiwan Province of China",country_name)
    ci <- RPostgreSQL::dbGetQuery(connec, query_compose(year=table1$Year.of.data[i],age_max = table1$X[i], country=country_name) )
    
    table1$Estimated.prevalence.ci.median[i] <- ci$prevalence_median
    table1$Estimated.prevalence.ci.interval[i] <- paste0("(95% CI ",ci$prevalence_lower,"-",ci$prevalence_upper, ")")
    
  }
  
}


table1$Variance.ci <- paste0(round((table1$Estimated.prevalence.ci.median- as.numeric(gsub(",","",table1$Published.prevalence)) )/ table1$Estimated.prevalence.ci.median *100,0),"%")
table1$Variance.absolute.ci <- paste0(abs(round((table1$Estimated.prevalence.ci.median- as.numeric(gsub(",","",table1$Published.prevalence)) )/ table1$Estimated.prevalence.ci.median *100,0)),"%")


write.csv(table1, "paper_submission/plots_paper_revision/Appendix_tables/appendix6-table1-ci.csv")












dbDisconnect(connec)
