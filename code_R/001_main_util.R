

merge_agg <- function(country_data_merge,all_data_points,all_age_brackets=TRUE)
{

  # country_data_merge <-  data_final_urban_rural
  merge_agg_age_bracket <- function( country_data_merge, start,end)
  { # start = 10 ; end=99;

    # setDF(country_data_merge)

    if(all_data_points)
    {
      country_data_merge_agg_1 <- setDT(country_data_merge)[Age<=end & Age >= start,list(
        # Basic stats of 1 run

        `Ann. background mortality`=sum(`Ann. background mortality`)
        ,`Ann. background population`=sum(`Ann. background population`)
        ,`Ann. early deaths`=sum(`Ann. early deaths`)
        ,`Ann. early deaths (background)`=sum(`Ann. early deaths (background)`)
        ,`Ann. onset deaths`=sum(`Ann. onset deaths`)
        ,`Incidence (1 base)`=sum(`Incidence (1 base)`)
        ,`Prevalence`=sum(`Prevalence`)
        ,`Prevalence with AtLeast1C`= if('Prevalence with AtLeast1C' %in% colnames(country_data_merge)) sum(`Prevalence with AtLeast1C`) else  NA
        # Advanced Stats
        ,`Ann. early deaths (2 diagnosis)`=sum(`Ann. early deaths (2 diagnosis)`)
        ,`Ann. early deaths (3 basic care)`=sum(`Ann. early deaths (3 basic care)`)
        ,`Ann. early deaths (4 best care)`=sum(`Ann. early deaths (4 best care)`)
        ,`Ann. early deaths (5 cure)`=sum(`Ann. early deaths (5 cure)`)
        # ,`diagnosis rate`=sum(`Ann. onset deaths`[Age==10])/ (sum(`Ann. onset deaths`[Age==10]) + `Incidence (1 base)`[Age==10])
        ,`Ghosts`=sum(`Ghosts`)
        ,`Ghosts (delta basic care)`=sum(`Ghosts (delta basic care)`)
        ,`Ghosts (delta best care)`=sum(`Ghosts (delta best care)`)
        ,`Ghosts (delta cure)`=sum(`Ghosts (delta cure)`)
        ,`Ghosts (early death)`=sum(`Ghosts (early death)`)
        ,`Ghosts (onset death)`=sum(`Ghosts (onset death)`)

        ,`Ghosts sim_start_year`= if('Ghosts sim_start_year' %in% colnames(country_data_merge))  sum(`Ghosts sim_start_year`) else NA
        ,`Ghosts (delta basic care) sim_start_year`= if('Ghosts (delta basic care) sim_start_year' %in% colnames(country_data_merge))sum(`Ghosts (delta basic care) sim_start_year`) else  NA
        ,`Ghosts (delta best care) sim_start_year`= if('Ghosts (delta best care) sim_start_year' %in% colnames(country_data_merge))sum(`Ghosts (delta best care) sim_start_year`) else  NA
        ,`Ghosts (delta cure) sim_start_year`= if('Ghosts (delta cure) sim_start_year' %in% colnames(country_data_merge)) sum(`Ghosts (delta cure) sim_start_year`) else  NA
        ,`Ghosts (early death) sim_start_year`= if('Ghosts (early death) sim_start_year' %in% colnames(country_data_merge))sum(`Ghosts (early death) sim_start_year`) else  NA
        ,`Ghosts (onset death) sim_start_year`=if('Ghosts (onset death) sim_start_year' %in% colnames(country_data_merge))sum(`Ghosts (onset death) sim_start_year`) else  NA



        ,`Incidence (2 diagnosis)`=sum(`Incidence (2 diagnosis)`)

        ,`Prevalence stagnant 1970`=sum(`Prevalence stagnant 1970`)

        ,`Prevalence delay onset by 1 years`  = if('Prevalence delay onset by 1 years' %in% colnames(country_data_merge)) sum(`Prevalence delay onset by 1 years`) else  NA
        ,`Prevalence delay onset by 3 years`  = if('Prevalence delay onset by 3 years' %in% colnames(country_data_merge)) sum(`Prevalence delay onset by 3 years`) else  NA
        ,`Prevalence delay onset by 5 years`  = if('Prevalence delay onset by 5 years' %in% colnames(country_data_merge)) sum(`Prevalence delay onset by 5 years`) else  NA
        ,`Prevalence delay onset by 8 years`  = if('Prevalence delay onset by 8 years' %in% colnames(country_data_merge)) sum(`Prevalence delay onset by 8 years`) else  NA
        ,`Prevalence delay onset by 13 years` =if('Prevalence delay onset by 13 years' %in% colnames(country_data_merge)) sum(`Prevalence delay onset by 13 years`) else  NA

        ,`Prevalence Stage 1 only` = if('Prevalence Stage 1 only' %in% colnames(country_data_merge)) sum(`Prevalence Stage 1 only`) else  NA
        ,`Prevalence Stage 2 only` = if('Prevalence Stage 2 only' %in% colnames(country_data_merge)) sum(`Prevalence Stage 2 only`) else  NA


        # ,`Ghosts delay onset by 1 years` = sum(`Ghosts delay onset by 1 years`)
        # ,`Ghosts delay onset by 3 years` = sum(`Ghosts delay onset by 3 years`)
        # ,`Ghosts delay onset by 5 years` = sum(`Ghosts delay onset by 5 years`)
        # ,`Ghosts delay onset by 8 years` = sum(`Ghosts delay onset by 8 years`)
        # ,`Ghosts delay onset by 13 years` = sum(`Ghosts delay onset by 13 years`)

        # # Lifetime years data   Age 10
        # , true_incidence_age_10                 = sum(`Incidence (1 base)`[Age==10])
        # ,`Ann. background population (age 10)`  = sum(`Ann. background population`[Age==10])
        # ,`prevalence plus missing (age 10)`     = sum(`Ghosts`[Age==10]) + sum(`Prevalence`[Age==10])
        # ,`Life expectency (1 background)`       = sum(`Life expectency (1 background)`[Age==10])
        # ,`Life expectency (2 t1d base)`         = sum(`Life expectency (2 t1d base)`[Age==10])
        # ,`Life expectency (3 t1d diagnosis)`    = sum(`Life expectency (3 t1d diagnosis)`[Age==10])
        # ,`Life expectency (4 t1d basic care)`   = sum(`Life expectency (4 t1d basic care)`[Age==10])
        # ,`Life expectency (5 t1d best care)`    = sum(`Life expectency (5 t1d best care)`[Age==10])
        # ,`Life expectency (6 t1d cure)`         = sum(`Life expectency (6 t1d cure)`[Age==10])
        # ,`Lifetime years lost (2 t1d base) (complication)`      = sum(`Lifetime years lost (2 t1d base) (complication)`[Age==10])
        # ,`Lifetime years lost (3 t1d diagnosis) (complication)` = sum(`Lifetime years lost (3 t1d diagnosis) (complication)`[Age==10])
        # ,`Lifetime years lost (4 t1d basic care) (complication)`= sum(`Lifetime years lost (4 t1d basic care) (complication)`[Age==10])
        # ,`Lifetime years lost (5 t1d best care) (complication)` = sum(`Lifetime years lost (5 t1d best care) (complication)`[Age==10])
        # ,`Lifetime years lost (6 t1d cure) (complication)`      = sum(`Lifetime years lost (6 t1d cure) (complication)`[Age==10])
        # ,`Lifetime years lost (2 t1d base) (treatment)`=sum(`Lifetime years lost (2 t1d base) (treatment)`[Age==10])
        # ,`Lifetime years lost (3 t1d diagnosis) (treatment)`=sum(`Lifetime years lost (3 t1d diagnosis) (treatment)`[Age==10])
        # ,`Lifetime years lost (4 t1d basic care) (treatment)`=sum(`Lifetime years lost (4 t1d basic care) (treatment)`[Age==10])
        # ,`Lifetime years lost (5 t1d best care) (treatment)`=sum(`Lifetime years lost (5 t1d best care) (treatment)`[Age==10])
        # ,`Lifetime years lost (6 t1d cure) (treatment)`=sum(`Lifetime years lost (6 t1d cure) (treatment)`[Age==10])
        #
        # ,`Life expectency (strip low)`   =   if('Life expectency (strip low)' %in% colnames(country_data_merge))   sum(`Life expectency (strip low)` [Age==10]) else  NA
        # ,`Life expectency (strip hig)`   =   if('Life expectency (strip hig)' %in% colnames(country_data_merge))   sum(`Life expectency (strip hig)` [Age==10]) else  NA
        # ,`Lifetime years lost (strip low)`   =   if('Lifetime years lost (strip low)' %in% colnames(country_data_merge))   sum(`Lifetime years lost (strip low)`  [Age==10]) else  NA
        # ,`Lifetime years lost (strip hig)`   =   if('Lifetime years lost (strip hig)' %in% colnames(country_data_merge))   sum(`Lifetime years lost (strip hig)` [Age==10]) else  NA
        # ,`Life expectency (sensor low)`   =   if('Life expectency (sensor low)' %in% colnames(country_data_merge))   sum(`Life expectency (sensor low)`  [Age==10]) else  NA
        # ,`Life expectency (sensor hig)`   =   if('Life expectency (sensor hig)' %in% colnames(country_data_merge))   sum(`Life expectency (sensor hig)`  [Age==10]) else  NA
        # ,`Lifetime years lost (sensor low)`   =   if('Lifetime years lost (sensor low)' %in% colnames(country_data_merge))   sum(`Lifetime years lost (sensor low)`  [Age==10]) else  NA
        # ,`Lifetime years lost (sensor hig)`   =   if('Lifetime years lost (sensor hig)' %in% colnames(country_data_merge))   sum(`Lifetime years lost (sensor hig)`  [Age==10]) else  NA
        #
        # ,`Lifetime years lost (delay onset year 1) (treatment)`     = if('Lifetime years lost (delay onset year 1) (treatment)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 1) (treatment)`  [Age==10]) else  NA
        # ,`Lifetime years lost (delay onset year 3) (treatment)`     = if('Lifetime years lost (delay onset year 3) (treatment)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 3) (treatment)`  [Age==10]) else  NA
        # ,`Lifetime years lost (delay onset year 5) (treatment)`     = if('Lifetime years lost (delay onset year 5) (treatment)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 5) (treatment)`  [Age==10]) else  NA
        # ,`Lifetime years lost (delay onset year 8) (treatment)`     = if('Lifetime years lost (delay onset year 8) (treatment)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 8) (treatment)`  [Age==10]) else  NA
        # ,`Lifetime years lost (delay onset year 13) (treatment)`    = if('Lifetime years lost (delay onset year 13) (treatment)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 13) (treatment)`  [Age==10]) else  NA
        #
        # ,`Lifetime years lost (delay onset year 1) (complication)`  = if('Lifetime years lost (delay onset year 1) (complication)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 1) (complication)`  [Age==10]) else  NA
        # ,`Lifetime years lost (delay onset year 3) (complication)`  = if('Lifetime years lost (delay onset year 3) (complication)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 3) (complication)`  [Age==10]) else  NA
        # ,`Lifetime years lost (delay onset year 5) (complication)`  = if('Lifetime years lost (delay onset year 5) (complication)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 5) (complication)`  [Age==10]) else  NA
        # ,`Lifetime years lost (delay onset year 8) (complication)`  = if('Lifetime years lost (delay onset year 8) (complication)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 8) (complication)`  [Age==10]) else  NA
        # ,`Lifetime years lost (delay onset year 13) (complication)` = if('Lifetime years lost (delay onset year 13) (complication)' %in% colnames(country_data_merge)) sum(`Lifetime years lost (delay onset year 13) (complication)`  [Age==10]) else  NA
        #
        # ,`Life expectency (delay onset year 1)`                     = if('Life expectency (delay onset year 1)' %in% colnames(country_data_merge)) sum(`Life expectency (delay onset year 1)`  [Age==10]) else  NA
        # ,`Life expectency (delay onset year 3)`                     = if('Life expectency (delay onset year 3)' %in% colnames(country_data_merge)) sum(`Life expectency (delay onset year 3)`  [Age==10]) else  NA
        # ,`Life expectency (delay onset year 5)`                     = if('Life expectency (delay onset year 5)' %in% colnames(country_data_merge)) sum(`Life expectency (delay onset year 5)`  [Age==10]) else  NA
        # ,`Life expectency (delay onset year 8)`                     = if('Life expectency (delay onset year 8)' %in% colnames(country_data_merge)) sum(`Life expectency (delay onset year 8)`  [Age==10]) else  NA
        # ,`Life expectency (delay onset year 13)`                    = if('Life expectency (delay onset year 13)' %in% colnames(country_data_merge)) sum(`Life expectency (delay onset year 13)`  [Age==10]) else  NA
        #
        # ,`1 in x families`=sum(`1 in x families`[Age==10])
        # ,`% Odds living to 55`= (sum(`Prevalence`[Age==55])+ 0.0000001)/ (sum(`Ghosts`[Age==55]) + sum(`Prevalence`[Age==55])+ 0.0000001) * 100
        # ,`% Odds living to 60`= (sum(`Prevalence`[Age==60])+ 0.0000001)/ (sum(`Ghosts`[Age==60]) + sum(`Prevalence`[Age==60])+ 0.0000001) * 100
        # ,`% Odds living to 65`= (sum(`Prevalence`[Age==65])+ 0.0000001)/ (sum(`Ghosts`[Age==65]) + sum(`Prevalence`[Age==65])+ 0.0000001) * 100

      ),by=c("diagnosis_input","loc_id","Country","Year","sim_start_year","sim_min_diag_rates","sim_min_non_minimal_care_perc","sim_min_non_minimal_care_level")]

      country_data_merge_agg_1 <- setDF(country_data_merge_agg_1)[,colSums(is.na(country_data_merge_agg_1)) != nrow(country_data_merge_agg_1)]

    }else
    {
      country_data_merge_agg_1 <- setDT(country_data_merge)[Age<=end & Age >= start,list(
        # Basic stats of 1 run
        `Ann. background mortality`=sum(`Ann. background mortality`)
        ,`Ann. background population`=sum(`Ann. background population`)
        ,`Ann. early deaths`=sum(`Ann. early deaths`)
        ,`Ann. early deaths (background)`=sum(`Ann. early deaths (background)`)
        ,`Ann. onset deaths`=sum(`Ann. onset deaths`)
        ,`Incidence (1 base)`=sum(`Incidence (1 base)`)
        ,`Prevalence`=sum(`Prevalence`)

      ),by=c("loc_id","Year")]
    }
    country_data_merge_agg <- country_data_merge_agg_1
    country_data_merge_agg <- cbind(age_bracket=paste0(sprintf("%02d", start), "_",sprintf("%02d", end)), country_data_merge_agg)
    # country_data_merge_agg
    #  aggregate complications ---
    colnames_comp <- colnames(country_data_merge)[grepl("Comp ",  colnames(country_data_merge))]
    if(length(colnames_comp))
    {
      code_string <- ""
      for(i in 1:length(colnames_comp))
      { # i <- 1
        code_string <- paste0(code_string, ",`",colnames_comp[i],"`=sum(`",colnames_comp[i],"`[Age<=end & Age >= start]) \n")
      }
      eval(parse(text = paste0("
                country_data_merge_agg_comp <- setDT(country_data_merge)[,list(
                 dotN =.N
                 ",code_string,"
                ),by=c('diagnosis_input','loc_id','Year')]"
      )))
      country_data_merge_agg_comp$dotN <- NULL
      country_data_merge_agg <- cbind(country_data_merge_agg, country_data_merge_agg_comp[,-(1:4)])
    }

    country_data_merge_agg
  }

country_data_merge2 <- data.frame()
country_data_ <-  merge_agg_age_bracket(country_data_merge,00,99); country_data_merge2 <- rbind(country_data_merge2,country_data_)
country_data_ <-  merge_agg_age_bracket(country_data_merge,00,19); country_data_merge2 <- rbind(country_data_merge2,country_data_)
country_data_ <-  merge_agg_age_bracket(country_data_merge,20,59); country_data_merge2 <- rbind(country_data_merge2,country_data_)
country_data_ <-  merge_agg_age_bracket(country_data_merge,60,99); country_data_merge2 <- rbind(country_data_merge2,country_data_)

country_data_ <-  merge_agg_age_bracket(country_data_merge,00,14); country_data_merge2 <- rbind(country_data_merge2,country_data_) # ofen used age bracket from studies.
country_data_ <-  merge_agg_age_bracket(country_data_merge,00,24); country_data_merge2 <- rbind(country_data_merge2,country_data_) # this is for factsheet, defined as Young People
country_data_ <-  merge_agg_age_bracket(country_data_merge,20,99); country_data_merge2 <- rbind(country_data_merge2,country_data_) # for young and adult seperation



if(all_age_brackets)
{
  country_data_ <-  merge_agg_age_bracket(country_data_merge,20,99); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,00,24); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,00, 9); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,65,99); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,00,29); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,85,99); country_data_merge2 <- rbind(country_data_merge2,country_data_)

  country_data_ <-  merge_agg_age_bracket(country_data_merge, 0, 4); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge, 5, 9); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,10,14); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,15,19); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,20,24); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,25,29); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,30,34); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,35,39); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,40,44); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,45,49); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,50,54); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,55,59); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,60,64); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,65,69); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,70,74); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,75,79); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,80,84); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,85,89); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,90,94); country_data_merge2 <- rbind(country_data_merge2,country_data_)
  country_data_ <-  merge_agg_age_bracket(country_data_merge,95,99); country_data_merge2 <- rbind(country_data_merge2,country_data_)

  country_data_ <-  merge_agg_age_bracket(country_data_merge,90,99); country_data_merge2 <- rbind(country_data_merge2,country_data_)

}

# country_data_ <-  merge_agg_age_bracket(country_data_merge,20,39); country_data_$age_bracket <- "20_39" ;  country_data_merge2 <- rbind(country_data_merge2,country_data_)
# country_data_ <-  merge_agg_age_bracket(country_data_merge,40,59); country_data_$age_bracket <- "40_59" ;  country_data_merge2 <- rbind(country_data_merge2,country_data_)
# country_data_ <-  merge_agg_age_bracket(country_data_merge,60,79); country_data_$age_bracket <- "60_79" ;  country_data_merge2 <- rbind(country_data_merge2,country_data_)
# country_data_ <-  merge_agg_age_bracket(country_data_merge,80,99); country_data_$age_bracket <- "80_99" ;  country_data_merge2 <- rbind(country_data_merge2,country_data_)
country_data_merge2
}


merge_agg_lifetime <- function(country_data_merge_age_10)
{
  # country_data_merge_age_10 <- country_data_raw_global
  country_data_2 <- setDT(country_data_merge_age_10)[,list(
    Age= 10
    , true_incidence_age_10                 = sum(`Incidence (1 base)`[Age==10])
    ,`Ann. background population (age 10)`  = sum(`Ann. background population`[Age==10])
    ,`prevalence plus missing (age 10)`     = sum(`Ghosts`[Age==10]) + sum(`Prevalence`[Age==10])
    ,`dDx_full`                             = sum(`dDx_full`[Age==10])
    ,`qT1D_percent_n_full`                  = sum(`qT1D_percent_n_full`[Age==10])
    ,`Life expectency (1 background)`       = sum(`Life expectency (1 background)`[Age==10])
    ,`Life expectency (2 t1d base)`         = sum(`Life expectency (2 t1d base)`[Age==10])
    ,`Life expectency (3 t1d diagnosis)`    = sum(`Life expectency (3 t1d diagnosis)`[Age==10])
    ,`Life expectency (4 t1d basic care)`   = sum(`Life expectency (4 t1d basic care)`[Age==10])
    ,`Life expectency (5 t1d best care)`    = sum(`Life expectency (5 t1d best care)`[Age==10])
    ,`Life expectency (6 t1d cure)`         = sum(`Life expectency (6 t1d cure)`[Age==10])
    ,`Lifetime years lost (2 t1d base) (complication)`      = sum(`Lifetime years lost (2 t1d base) (complication)`[Age==10])
    ,`Lifetime years lost (3 t1d diagnosis) (complication)` = sum(`Lifetime years lost (3 t1d diagnosis) (complication)`[Age==10])
    ,`Lifetime years lost (4 t1d basic care) (complication)`= sum(`Lifetime years lost (4 t1d basic care) (complication)`[Age==10])
    ,`Lifetime years lost (5 t1d best care) (complication)` = sum(`Lifetime years lost (5 t1d best care) (complication)`[Age==10])
    ,`Lifetime years lost (6 t1d cure) (complication)`      = sum(`Lifetime years lost (6 t1d cure) (complication)`[Age==10])
    ,`Lifetime years lost (2 t1d base) (treatment)`=sum(`Lifetime years lost (2 t1d base) (treatment)`[Age==10])
    ,`Lifetime years lost (3 t1d diagnosis) (treatment)`=sum(`Lifetime years lost (3 t1d diagnosis) (treatment)`[Age==10])
    ,`Lifetime years lost (4 t1d basic care) (treatment)`=sum(`Lifetime years lost (4 t1d basic care) (treatment)`[Age==10])
    ,`Lifetime years lost (5 t1d best care) (treatment)`=sum(`Lifetime years lost (5 t1d best care) (treatment)`[Age==10])
    ,`Lifetime years lost (6 t1d cure) (treatment)`=sum(`Lifetime years lost (6 t1d cure) (treatment)`[Age==10])

    ,`Healthy years restored with onset diagnosis` = round( (((`Life expectency (3 t1d diagnosis)`- `Lifetime years lost (3 t1d diagnosis) (complication)`- `Lifetime years lost (3 t1d diagnosis) (treatment)`  )[ Age==10] )
                                                                - ((`Life expectency (2 t1d base)`- `Lifetime years lost (2 t1d base) (complication)`- `Lifetime years lost (2 t1d base) (treatment)`  )[ Age==10] ))
                                                               ,1)

    ,`Healthy years restored with insulin and strips` = round( (((`Life expectency (4 t1d basic care)`- `Lifetime years lost (4 t1d basic care) (complication)`- `Lifetime years lost (4 t1d basic care) (treatment)`  )[ Age==10] )
                                                                - ((`Life expectency (3 t1d diagnosis)`- `Lifetime years lost (3 t1d diagnosis) (complication)`- `Lifetime years lost (3 t1d diagnosis) (treatment)`  )[ Age==10] ))
                                                               ,1)

    ,`Healthy years restored with device uptake` = round( (((`Life expectency (5 t1d best care)`- `Lifetime years lost (5 t1d best care) (complication)`- `Lifetime years lost (5 t1d best care) (treatment)`  )[ Age==10] )
                                                           - ((`Life expectency (4 t1d basic care)`- `Lifetime years lost (4 t1d basic care) (complication)`- `Lifetime years lost (4 t1d basic care) (treatment)`  )[ Age==10] ))
                                                          ,1)

    ,`Life expectency (strip low)`   =   if('Life expectency (strip low)' %in% colnames(country_data_merge_age_10))   sum(`Life expectency (strip low)` [Age==10]) else  NA
    ,`Life expectency (strip hig)`   =   if('Life expectency (strip hig)' %in% colnames(country_data_merge_age_10))   sum(`Life expectency (strip hig)` [Age==10]) else  NA
    ,`Lifetime years lost (strip low)`   =   if('Lifetime years lost (strip low)' %in% colnames(country_data_merge_age_10))   sum(`Lifetime years lost (strip low)`  [Age==10]) else  NA
    ,`Lifetime years lost (strip hig)`   =   if('Lifetime years lost (strip hig)' %in% colnames(country_data_merge_age_10))   sum(`Lifetime years lost (strip hig)` [Age==10]) else  NA
    ,`Life expectency (sensor low)`   =   if('Life expectency (sensor low)' %in% colnames(country_data_merge_age_10))   sum(`Life expectency (sensor low)`  [Age==10]) else  NA
    ,`Life expectency (sensor hig)`   =   if('Life expectency (sensor hig)' %in% colnames(country_data_merge_age_10))   sum(`Life expectency (sensor hig)`  [Age==10]) else  NA
    ,`Lifetime years lost (sensor low)`   =   if('Lifetime years lost (sensor low)' %in% colnames(country_data_merge_age_10))   sum(`Lifetime years lost (sensor low)`  [Age==10]) else  NA
    ,`Lifetime years lost (sensor hig)`   =   if('Lifetime years lost (sensor hig)' %in% colnames(country_data_merge_age_10))   sum(`Lifetime years lost (sensor hig)`  [Age==10]) else  NA

    ,`Lifetime years lost (delay onset year 1) (treatment)`     = if('Lifetime years lost (delay onset year 1) (treatment)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 1) (treatment)`  [Age==10]) else  NA
    ,`Lifetime years lost (delay onset year 3) (treatment)`     = if('Lifetime years lost (delay onset year 3) (treatment)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 3) (treatment)`  [Age==10]) else  NA
    ,`Lifetime years lost (delay onset year 5) (treatment)`     = if('Lifetime years lost (delay onset year 5) (treatment)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 5) (treatment)`  [Age==10]) else  NA
    ,`Lifetime years lost (delay onset year 8) (treatment)`     = if('Lifetime years lost (delay onset year 8) (treatment)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 8) (treatment)`  [Age==10]) else  NA
    ,`Lifetime years lost (delay onset year 13) (treatment)`    = if('Lifetime years lost (delay onset year 13) (treatment)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 13) (treatment)`  [Age==10]) else  NA

    ,`Lifetime years lost (delay onset year 1) (complication)`  = if('Lifetime years lost (delay onset year 1) (complication)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 1) (complication)`  [Age==10]) else  NA
    ,`Lifetime years lost (delay onset year 3) (complication)`  = if('Lifetime years lost (delay onset year 3) (complication)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 3) (complication)`  [Age==10]) else  NA
    ,`Lifetime years lost (delay onset year 5) (complication)`  = if('Lifetime years lost (delay onset year 5) (complication)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 5) (complication)`  [Age==10]) else  NA
    ,`Lifetime years lost (delay onset year 8) (complication)`  = if('Lifetime years lost (delay onset year 8) (complication)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 8) (complication)`  [Age==10]) else  NA
    ,`Lifetime years lost (delay onset year 13) (complication)` = if('Lifetime years lost (delay onset year 13) (complication)' %in% colnames(country_data_merge_age_10)) sum(`Lifetime years lost (delay onset year 13) (complication)`  [Age==10]) else  NA

    ,`Life expectency (delay onset year 1)`                     = if('Life expectency (delay onset year 1)' %in% colnames(country_data_merge_age_10)) sum(`Life expectency (delay onset year 1)`  [Age==10]) else  NA
    ,`Life expectency (delay onset year 3)`                     = if('Life expectency (delay onset year 3)' %in% colnames(country_data_merge_age_10)) sum(`Life expectency (delay onset year 3)`  [Age==10]) else  NA
    ,`Life expectency (delay onset year 5)`                     = if('Life expectency (delay onset year 5)' %in% colnames(country_data_merge_age_10)) sum(`Life expectency (delay onset year 5)`  [Age==10]) else  NA
    ,`Life expectency (delay onset year 8)`                     = if('Life expectency (delay onset year 8)' %in% colnames(country_data_merge_age_10)) sum(`Life expectency (delay onset year 8)`  [Age==10]) else  NA
    ,`Life expectency (delay onset year 13)`                    = if('Life expectency (delay onset year 13)' %in% colnames(country_data_merge_age_10)) sum(`Life expectency (delay onset year 13)`  [Age==10]) else  NA

    ,`1 in x families`=sum(`1 in x families`[Age==10])
    ,`% Odds living to 55`= (sum(`Prevalence`[Age==55])+ 0.0000001)/ (sum(`Ghosts`[Age==55]) + sum(`Prevalence`[Age==55])+ 0.0000001) * 100
    ,`% Odds living to 60`= (sum(`Prevalence`[Age==60])+ 0.0000001)/ (sum(`Ghosts`[Age==60]) + sum(`Prevalence`[Age==60])+ 0.0000001) * 100
    ,`% Odds living to 65`= (sum(`Prevalence`[Age==65])+ 0.0000001)/ (sum(`Ghosts`[Age==65]) + sum(`Prevalence`[Age==65])+ 0.0000001) * 100

  ),by=c("diagnosis_input","loc_id","Country","Year","sim_start_year","sim_min_diag_rates","sim_min_non_minimal_care_perc","sim_min_non_minimal_care_level")]

  country_data_2 <- setDF(country_data_2)[,colSums(is.na(country_data_2)) != nrow(country_data_2)]


  country_data_2

}


read_run_results <- function(cache_dir)
{
  files_db <- list.files(cache_dir,full.names=TRUE)
  files_db <- files_db[!grepl("_lifetime", files_db) ]
  files_db <- files_db[!grepl("_raw", files_db) ]
  diagnosis_input <- "gregory"# diagnosis_input <- "ward"# diagnosis_input <- "blend"
  #----- ----------------------- Read all rows -------------------- --------------------
  read_rds_file <- function(file_path,diagnosis_input) {
    # file_path <- files_db[6]
    # country_data_merge_temp <- readRDS(file_path)
    country_data_merge_temp <- arrow::read_feather(file_path)
    # country_data_merge_temp <- country_data_merge_temp[country_data_merge_temp$diagnosis_input==diagnosis_input,]
    # rm(country_data_merge_temp)
    return(country_data_merge_temp)
  }
  country_data_merge2      <- rbindlist(lapply(files_db, read_rds_file,diagnosis_input) )

  run_sequence <- setDT(country_data_merge2)[,list(.N),by=c( 'sim_start_year','sim_min_diag_rates','sim_min_non_minimal_care_perc', 'sim_min_non_minimal_care_level')]
  print(paste0("number of NAs in country_data_merge: ", sum(is.na(country_data_merge2))))

  if(FALSE)
  { # data check
    country_data_merge2[country_data_merge2$Year==2040,list(`Ghosts (delta best care) lever2023`= sum(`Ghosts (delta best care) lever2023`)), by=c("age_bracket")]
  }

  #----- ----------------------- Read _lifetime  -------------------- --------------------
  files_db <- list.files(cache_dir,full.names=TRUE)
  files_db <- files_db[ grepl("_lifetime", files_db) ]
  country_data_2      <- rbindlist(lapply(files_db, read_rds_file,diagnosis_input) )

  run_sequence <- setDT(country_data_2)[,list(.N),by=c( 'sim_start_year','sim_min_diag_rates','sim_min_non_minimal_care_perc', 'sim_min_non_minimal_care_level')]

  #----- ----------------------- Read _raw  -------------------- --------------------
  files_db <- list.files(cache_dir,full.names=TRUE)
  files_db <- files_db[ grepl("_raw", files_db) ]
  country_data_raw      <- rbindlist(lapply(files_db, read_rds_file,diagnosis_input) )




  return(list(country_data_merge2=country_data_merge2 , country_data_2=country_data_2,country_data_raw=country_data_raw))

}
