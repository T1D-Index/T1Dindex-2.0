
calculate_ex_gpu <- function(life_table)
{

  # life_table <-  life_table_background_


  library(torch)
  device <- torch_device("cuda")

  life_table <- cbind(life_table, lx=c(NA) )
  life_table <- cbind(life_table, dx=c(NA) )
  life_table <- cbind(life_table, Lx=c(NA) )
  life_table <- cbind(life_table, Tx=c(NA) )
  life_table <- cbind(life_table, ex=c(NA) )
  life_table_gpu <- torch_tensor(life_table,  device = device)
  # life_table <- life_table_t1d
  # life_table
  # life_table$lx <- NA
  # life_table$dx <- NA

  life_table_gpu[,5][life_table_gpu[,3]==0] <- 100000
  life_table_gpu[,6][life_table_gpu[,3]==0] <- 100000 * life_table_gpu[,4][life_table_gpu[,3]==0]
  for(i in 1:max(life_table[,"age"]))
  {# i <- 1
    life_table_gpu[,5][life_table_gpu[,3]==i] <- life_table_gpu[,5][life_table_gpu[,3]==(i-1)] -  life_table_gpu[,6][life_table_gpu[,3]==(i-1)]
    life_table_gpu[,6][life_table_gpu[,3]==i] <- life_table_gpu[,5][life_table_gpu[,3]==i] * life_table_gpu[,4][life_table_gpu[,3]==i]
  }

  life_table_gpu[,7] <- torch_tensor(c(as_array(life_table_gpu$cpu())[,5][-1],0),device = device)
  life_table_gpu[,7][life_table_gpu[,3]==max(life_table_gpu[,3])] <- 0
  life_table_gpu[,7] <- life_table_gpu[,7] + 0.5 * life_table_gpu[,6]

  life_table_gpu[,8] <- 0
  life_table_gpu[,8][life_table_gpu[,3]==max(life_table_gpu[,3])] <- life_table_gpu[,8][life_table_gpu[,3]==max(life_table_gpu[,3])]

  for(i in (max(life_table[,3])-1):0)
  { # i <- 1
    life_table_gpu[,8][life_table_gpu[,3]==i] <- life_table_gpu[,8][life_table_gpu[,3]==(i+1) ]+ life_table_gpu[,7][life_table_gpu[,3]==i ]
  }


  life_table_gpu <- as_array(life_table_gpu$cpu())
  life_table_gpu[,9] <- ifelse(life_table_gpu[,8]==0 , 0 ,life_table_gpu[,8] / life_table_gpu[,5] )
  life_table_gpu <- as.data.frame(life_table_gpu)
  colnames(life_table_gpu) <- colnames(life_table)
  life_table_gpu <- life_table_gpu[ with(life_table_gpu, order(age,year, loc_id)), ]   # order properly for transforming to matrix

  return(life_table_gpu)
}



# # calculate ex  life expectency at age x ---------------
# calculate_ex <- function(life_table)
# {
#   # life_table <- dplyr::select(data_long, loc_id, year,age,qx =background_mortality_rate)
#   # life_table
#   life_table$lx <- NA
#   life_table$dx <- NA
#   life_table$lx[life_table$age==0] <- 100000
#   life_table$dx[life_table$age==0] <- 100000 * life_table$qx[life_table$age==0]
#   for(i in 1:max(life_table$age))
#   {# i <- 1
#     life_table$lx[life_table$age==i] <- life_table$lx[life_table$age==(i-1)] -  life_table$dx[life_table$age==(i-1)]
#     life_table$dx[life_table$age==i] <- life_table$lx[life_table$age==i] * life_table$qx[life_table$age==i]
#   }
#
#   life_table$Lx <- c(life_table$lx[-1],0)
#   life_table$Lx[life_table$age==max(life_table$age)] <- 0
#   life_table$Lx <- life_table$Lx + 0.5 * life_table$dx
#
#   life_table$Tx <- 0
#   life_table$Tx[life_table$age==max(life_table$age)] <- life_table$Lx[life_table$age==max(life_table$age)]
#
#   for(i in (max(life_table$age)-1):0)
#   { # i <- 1
#     life_table$Tx[life_table$age==i] <- life_table$Tx[life_table$age==(i+1) ]+ life_table$Lx[life_table$age==i ]
#   }
#
#   life_table$ex <- ifelse(life_table$Tx==0 , 0 ,life_table$Tx / life_table$lx )
#
#   return(life_table)
# }



calculate_prevalence_tensor <- function(i, qB,  qT1D_n, qT1D_m, qT1D_percent_n, dDx,smr_matrix,track_cohort=FALSE,run_complications=FALSE) {
  # run_complications=FALSE ; run_complications=TRUE
  # track_cohort=FALSE
  weibull <- data.frame()
  if(run_complications)  # not run complications
  {

    hba1c_full   <- (log(smr_matrix) +  1.5274 )/ 0.3545
    weibull      <- get_complication_parameters()$weibull
    # weibull <- weibull[1:2,]

    weib_survival_props <- function(hba1c_year,age_pair, weibull) {

      tail <- as.numeric(dimnames(age_pair)[[3]])
      head <- round(age_pair,0)
      props_full <- array(0,dim=c(nrow(weibull)+1,dim(hba1c_year)[[1]], length(age_pair)), dimnames = list(comp=c(weibull$abbrev,"AtLeast1C"),dimnames(hba1c_year)[[1]], age=1:length(age_pair)))

      for(a in 1:length( age_pair) )
      { # a <- 2
        if(tail[a] >  head[a])
        {
          hba1c <- hba1c_year
          hba1c <- hba1c[,,head[a]:tail[a] ,drop=FALSE]
          # weibull=complication_parameters$weibull  ;  hba1c <- c(10,10,10,10,10,10)
          props <- array(0,dim=c(nrow(weibull)+1,dim(hba1c)[[1]], dim(hba1c)[[3]] ), dimnames = list(comp=c(weibull$abbrev,"AtLeast1C"), loc_id=dimnames(hba1c)[[1]], age=dimnames(hba1c)[[3]]))
          for(i_ in 1:nrow(weibull))
          {
            # i_ <- 1
            intercept <- weibull$intercept[i_]
            slope     <- weibull$slope[i_]
            scale     <- weibull$scale[i_]
            hba1c_age_cut_10 <- hba1c
            T <- dim(hba1c_age_cut_10)[[3]]-1 # total time length (years) incl. T=0
            hs_odd <- 0.5 * (hba1c_age_cut_10[,,-1,drop=FALSE] + hba1c_age_cut_10[,,-(T+1),drop=FALSE]) # interpolate at half-year points
            ts <- pmin(0:T, 30)  # constant hazard after 30 years
            ts <- array(ts,dim=c(length(ts), dim(hba1c_age_cut_10)[[2]], dim(hba1c_age_cut_10)[[1]] ))
            ts <- aperm(ts,c(3,2,1))
            # Simpson's rule
            l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c_age_cut_10)/scale) / scale
            l_odd  <- (ts[,,-1,drop=FALSE] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale
            l_even_cumsum <- abind( array(0, dim=c(dim(hba1c_age_cut_10)[[1]], dim(hba1c_age_cut_10)[[2]], 2 )), l_even[,,-c(1, T+1),drop=FALSE],along=3)
            l_odd_cumsum  <- abind( array(0, dim=c(dim(hba1c_age_cut_10)[[1]], dim(hba1c_age_cut_10)[[2]], 1 )),  l_odd[,,-(T+1),drop=FALSE]    ,along=3)

            cumsum_matrix <-  2 * l_even_cumsum + 4 * l_odd_cumsum
            for(j in dim(hba1c_age_cut_10)[[3]]:2)
            { # j <- 100
              cumsum_matrix[,,j] <-    rowSums(cumsum_matrix[,,1:j,drop=FALSE], dims = 2)
            }
            Lambda <- 1/6 * (l_even + cumsum_matrix)
            props[weibull$abbrev[i_],,] <- 1- exp(-Lambda)

          }
          # This is Graham's and Gabriel's 'journeys' adjustments. Reduces ON and PR by survivor shares of RF and BL, respectively,
          # effectively making this a multistate model (rather than a pure time-to-event model)
          props['ON',,] <- props['ON',,] * (1 - props['RF',,])
          props['PR',,] <- props['PR',,] * (1 - props['BL',,])
          props['AtLeast1C',,] <- 1-(apply((1-props), 3, prod))
          # props[,length(hba1c),drop=FALSE]
          props_full[,,a] <-  props[,,dim(hba1c)[[3]],drop=FALSE]
        }
      }
      props_full
    }
  }


  loc_ids <- dimnames(i)[[1]]
  years   <- (1960 - 100):2040
  AGES    <- 0: 99
  MAX_AGE <- 100

  # Model compartments, all annual cohorts as a proportion of 1.
  S <-  array(NA, dim=list(length(loc_ids), length(years), MAX_AGE), dimnames=list(loc_ids,years, AGES))# S - susceptible (non-T1D)
  P <-  array(NA, dim=list(length(loc_ids), length(years), MAX_AGE), dimnames=list(loc_ids,years, AGES))# S - susceptible (non-T1D)
  D <-  array(NA, dim=list(length(loc_ids), length(years), MAX_AGE), dimnames=list(loc_ids,years, AGES))# S - susceptible (non-T1D)


  O    <-  array(NA, dim=list(                 length(loc_ids),length(years), MAX_AGE), dimnames=list(                              loc_ids,years, AGES))     # O , mean age of Onset
  C_I  <-  array(0,  dim=list(nrow(weibull)+1, length(loc_ids),length(years), MAX_AGE), dimnames=list(c(weibull$abbrev,"AtLeast1C"),loc_ids,years, AGES))  #  complication incidence  probability
  C_P  <-  array(NA, dim=list(nrow(weibull)+1, length(loc_ids),length(years), MAX_AGE), dimnames=list(c(weibull$abbrev,"AtLeast1C"),loc_ids,years, AGES)) #  complication prevalence cases



  # incidence flows are reflected the following year, for individuals one year older
  i_shift <- array(NA, dim=list(length(loc_ids), length(years), MAX_AGE), dimnames=list(loc_ids,years, AGES))
  i_shift[,, 1] <- 0
  i_shift[,,-1] <- i[,,-MAX_AGE]
  i_all <- i / (1-dDx)              # includes deaths at onset
  i_all_shift <- i_shift / (1-dDx)  # shifted & includes deaths at onset

  # Track P cohorts and in/outflows in 3D array: {year}x{age}x{onset age}
  Pcohorts <- NULL
  if(track_cohort)
  {
    Pcohorts  <- array(NA, dim=list(length(loc_ids),length(years), MAX_AGE, MAX_AGE), dimnames=list(loc_ids,years, AGES, AGES))
    Icohorts  <- array(0, dim=list(length(loc_ids),length(years), MAX_AGE, MAX_AGE) , dimnames=list(loc_ids,years, AGES, AGES))
    PDcohorts <- array(0, dim=list(length(loc_ids),length(years), MAX_AGE, MAX_AGE) , dimnames=list(loc_ids,years, AGES, AGES))
    Pcohorts[,1,1,1] <- 0

  }

  # The model proceeds yearwise, populating successive matrix rows.
  # Equation references below refer to the model summary.
  S[,1,1] <- 1
  P[,1,1] <- 0
  D[,1,1] <- 0

  O[,,1] <- 0
  C_P[,1,1,1] <- 0

  dim_cohorts <- 'if'(length(loc_ids)==1, c(MAX_AGE, MAX_AGE),c(length(loc_ids),MAX_AGE, MAX_AGE))
  diagnal_1s_100_99 <-  array(diag(1,100), c(100,99,length(loc_ids)))
  diagnal_1s_100_99 <-  aperm(diagnal_1s_100_99, c(3,1,2))
  diagnal_1s_100_100 <-  array(diag(1,100), c(100,100,length(loc_ids)))
  diagnal_1s_100_100 <-  aperm(diagnal_1s_100_100, c(3,1,2))


  for (t in seq_along(years[-1])) {
    # t <- 1
    # for (t in 1:10) {
    # bump previous period's data along one cell (ie one year older)
    Sshift <- abind( array(1, dim=c(length(loc_ids),1)), S[,t, -MAX_AGE,drop=FALSE],along=3); dimnames(Sshift)[[3]] <- as.character(0:99)
    Pshift <- abind( array(0, dim=c(length(loc_ids),1)), P[,t, -MAX_AGE,drop=FALSE],along=3); dimnames(Pshift)[[3]] <- as.character(0:99)
    Dshift <- abind( array(0, dim=c(length(loc_ids),1)), D[,t, -MAX_AGE,drop=FALSE],along=3); dimnames(Dshift)[[3]] <- as.character(0:99)

    if(track_cohort)
    {
      PCshift <- array(0, dim=c(length(loc_ids),MAX_AGE, MAX_AGE))
      PCshift[,1,] <- 0
      PCshift[,-1,] <- Pcohorts[,t,-MAX_AGE,]
    }

    #
    # # equation (3) - susceptible compartment, S
    S[,t+1,] <- Sshift[,1,] * (1 - i_all_shift[,t,]) * (1 - qB[,t,])

    # # equation (4) - prevalence compartment P
    P[,t+1,] <-    qT1D_percent_n[,t,]  *  Pshift[,1,]  * (1 - qT1D_n[,t,]) + i_shift[,t,] * Sshift[,1,] +
      (1-qT1D_percent_n[,t,]) *  Pshift[,1,]  * (1 - qT1D_m[,t,])
    #
    if(run_complications)  # not run complications
    {
      # Oshift <- c(0, O[,t, -MAX_AGE]) ; names(Oshift) <- as.character(0:99) # mean age of onset plus 1 with shifting
      Oshift <- abind( array(0, dim=c(length(loc_ids),1)), O[,t, -MAX_AGE,drop=FALSE],along=3) ; dimnames(Oshift)[[3]] <- as.character(0:99)

      # C_P_shift <- cbind(rep(0, nrow(weibull)+1), C_P[,t, -MAX_AGE]) ; colnames(C_P_shift) <- as.character(0:99) #
      C_P_shift <- abind(array(0, dim=c(nrow(weibull)+1,length(loc_ids),1)), C_P[,,t, -MAX_AGE,drop=FALSE],along=4) ; dimnames(C_P_shift)[[4]] <- as.character(0:99) #

      index_not_na <- !is.na(Pshift) ; index_not_na[1] <- FALSE

      if(sum(index_not_na))
      {
        O[,t+1,index_not_na] <- ( Oshift[index_not_na] * Pshift[index_not_na] +  c(0:99) [index_not_na]* i_shift[,t,index_not_na]* Sshift[index_not_na] + 1e-20 ) / ( Pshift[index_not_na] + i_shift[,t,index_not_na]* Sshift[index_not_na] + 1e-20 )

        C_I[,,t+1,index_not_na]  <- weib_survival_props (hba1c_year= hba1c_full[,t,,drop=FALSE], age_pair = O[,t+1,index_not_na,drop=FALSE], weibull)

        C_P[,,t+1,]  <-  qT1D_percent_n[,t,]  *  C_P_shift  * (1 - qT1D_n[,t,])  +
          (1-qT1D_percent_n[,t,]) *  C_P_shift  * (1 - qT1D_m[,t,])  +   C_I[,,t+1,,drop=FALSE] * (  abind(replicate(11, P[,t+1,,drop=FALSE], simplify = FALSE), along = 0) - C_P_shift  )
      }
    }



    if(track_cohort)
    {
      # # equation (4') - prevalence compartment P by cohort
      Ishift <-  array(NA, dim=c(length(loc_ids),MAX_AGE, MAX_AGE))

      Ishift_diag   <- i_shift[,t,] * Sshift[,1,]
      Icohorts_diag <- i[,t,] * S[,t,]
      # for(c in 1:length(loc_ids))
      # {  # c <- 1
      #   Ishift[c,,-MAX_AGE] <- diag(Ishift_diag[c,])[,-1] # Only want to shift age, not onset age
      #
      #   Icohorts[c,t,,] <- diag(Icohorts_diag[c,])
      #
      #   # PDcohorts[c,t,,] <- qT1D_percent_n[c,t,] * Pcohorts[c,t,,] * array(qT1D_n[c,t,], c(MAX_AGE, MAX_AGE)) +
      #   #   (1 - qT1D_percent_n[c,t,]) * Pcohorts[c,t,,] * array(qT1D_m[c,t,], c(MAX_AGE, MAX_AGE))
      #   # Pcohorts[c,t+1,,] <- qT1D_percent_n[c,t,] * PCshift[c,,] * array(1 - qT1D_n[c,t,], c(MAX_AGE, MAX_AGE)) + Ishift[c,,] +
      #   #   (1 - qT1D_percent_n[c,t,]) * PCshift[c,,] * array(1 - qT1D_m[c,t,], c(MAX_AGE, MAX_AGE))
      # }

      Ishift[,,-MAX_AGE] <- (Ishift_diag[,-1]) # Only want to shift age, not onset age
      Ishift[,,-MAX_AGE] <- Ishift[,,-MAX_AGE] *   diagnal_1s_100_99

      Icohorts[,t,,] <- Icohorts_diag[,]
      Icohorts[,t,,] <- Icohorts[,t,,] *   diagnal_1s_100_100


      Ishift[,,MAX_AGE]  <- 0  # NB: half-cycle adjustment populates any incidence for MAX_AGE

      PDcohorts[,t,,] <- array(qT1D_percent_n[,t,], dim_cohorts)  * Pcohorts[,t,,] * array(qT1D_n[,t,], dim_cohorts) +
        (1 - array(qT1D_percent_n[,t,], dim_cohorts)) * Pcohorts[,t,,] * array(qT1D_m[,t,], dim_cohorts)
      Pcohorts[,t+1,,] <- array(qT1D_percent_n[,t,], dim_cohorts) * PCshift[,,] * array(1 - qT1D_n[,t,], dim_cohorts) + Ishift[,,] +
        (1 - array(qT1D_percent_n[,t,], dim_cohorts)) * PCshift[,,] * array(1 - qT1D_m[,t,], dim_cohorts)


    }


    # # equation (5) - death compartment D
    D[,t+1,] <- (Dshift[,1,]
                + i_all_shift[,t,] * dDx[,t,] * Sshift[,1,]
                # + qT1D[t,] * Pshift
                + qT1D_n[,t,] * Pshift[,1,] *       qT1D_percent_n[,t,]
                + qT1D_m[,t,] * Pshift[,1,] *  (1 - qT1D_percent_n[,t,])
                + Sshift[,1,] * qB[,t,] * (1 - i_all_shift[,t,]))

  }
  # final year cohort flows - same as in the loop but for final period
  C_P[,,"1960","99"] <- 0 # Replace NA becasuse of 1 year lagging due to simulation start years

  t <- t + 1

  if(track_cohort)
  {
    for(c in 1:length(loc_ids))
    {
      Icohorts[c,t,,]   <- diag(i[c,t,] * S[c,t,])
      PDcohorts[c,t,,]  <- Pcohorts[c,t,,] * array(qT1D_n[c,t,], c(MAX_AGE, MAX_AGE)) *      qT1D_percent_n[c,t,]  +
        Pcohorts[c,t,,] * array(qT1D_m[c,t,], c(MAX_AGE, MAX_AGE)) * (1 - qT1D_percent_n[c,t,])
    }
    Pcohorts <- Pcohorts + 0.5 * Icohorts - 0.5 * PDcohorts

  }


  # # flows: based on unshifted versions of incidence
  I <- i * S                    # T1D incidence
  DDx <- i_all * dDx * S        # deaths at T1D onset
  # # DT1D <- (qT1D - qB) * P       # T1D-cause mortality
  DT1D <- (qT1D_n - qB) * P  *          qT1D_percent_n  +
          (qT1D_m - qB) * P  *     (1 - qT1D_percent_n)   # T1D-cause mortality
  DBGP <- P * qB                # background mortality for people w/ T1D
  DBGS <- S * qB * (1 - i_all)  # background mortality for susceptible pop'n
  # # DBGS <- S * qB  # background mortality for susceptible pop'n
  #
  # # Half-cycle adjustments - apply half of each flow in the reference
  # # year for that flow. We expect that on average, 1/2 of the flow
  # # occurs part-way through the year.
  S <- S - 0.5 * I - 0.5 * DDx - 0.5 * DBGS
  P <- P + 0.5 * I - 0.5 * (DT1D + DBGP)
  D <- D + 0.5 * (DDx + DBGS + DT1D + DBGP)
  list(S=S, P=P, Pcohorts=Pcohorts, D=D, I=I, DDx=DDx, DT1D=DT1D,DBGP=DBGP, DBGS=DBGS,C_P=C_P)
}


calculate_prevalence_cuda_test <- function(test) {

  device <- torch_device("cuda")
  i_gpu  <- torch_tensor(i, device = device)
  qB_gpu <- torch_tensor(qB, device = device)
  qT1D_n_gpu <- torch_tensor(qT1D_n, device = device)
  qT1D_m_gpu <- torch_tensor(qT1D_m, device = device)
  qT1D_percent_n_gpu <- torch_tensor(qT1D_percent_n, device = device)
  dDx_gpu <- torch_tensor(dDx, device = device)

  S_gpu <- torch_tensor( matrix(-1, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES)),device = device )     # S - susceptible (non-T1D)
  P_gpu <- torch_tensor( matrix(-1, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES)),device = device )     # S - susceptible (non-T1D)
  D_gpu <- torch_tensor( matrix(-1, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES)),device = device )     # S - susceptible (non-T1D)

  Pcohorts_gpu   <- torch_tensor( array(0, dim=list(length(years), MAX_AGE, MAX_AGE),dimnames=list(years, AGES, AGES)),device = device )
  Icohorts_gpu   <- torch_tensor( array(0, dim=list(length(years), MAX_AGE, MAX_AGE),dimnames=list(years, AGES, AGES)),device = device )
  PDcohorts_gpu  <- torch_tensor( array(0, dim=list(length(years), MAX_AGE, MAX_AGE),dimnames=list(years, AGES, AGES)),device = device )


  # S_gpu[1,1] <- 1
  # P_gpu[1,1] <- PDcohorts_gpu[1,1,1] <- D_gpu[1,1] <- 0
  # # Icohorts_gpu[,,] <- PDcohorts_gpu[,,] <- 0
  #
  # for (t in seq_along(years[-1])) {
  #  # t <- 1
  #   Sshift <-  S_gpu[t, 1:(MAX_AGE-1) ]
  #   Pshift <-  P_gpu[t, 1:(MAX_AGE-1) ]
  #   Dshift <-  D_gpu[t, 1:(MAX_AGE-1) ]
  #   #
  #   # PCshift <- torch_tensor(array(0, dim=c(MAX_AGE, MAX_AGE)),device = device )
  #   # PCshift[1,] <- 0
  #   # PCshift[2:MAX_AGE,] <- Pcohorts_gpu[t,-MAX_AGE,]
  #
  #   # Ishift <- torch_tensor(matrix(NA, nrow = MAX_AGE, ncol = MAX_AGE),device = device )
  #   # Ishift[,1:(MAX_AGE-1)] <- diag(i_shift[t,] * Sshift)[,-1] # Only want to shift age, not onset age
  #   # Ishift[,MAX_AGE] <- 0  # NB: half-cycle adjustment populates any incidence for MAX_AGE
  #
  # }


  # # Data matrixes: years on rows, ages on columns. Memory footprint of storing
  # # all this information really is tiny, a few MB.
  # model_matrix <- function() {
  #   matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))
  # }
  # # Similarly, a 3D array: {year}x{age}x{onset age}
  # model_array <- function() {
  #   array(NA, dim=list(length(years), MAX_AGE, MAX_AGE),
  #         dimnames=list(years, AGES, AGES))
  # }
  # # Model compartments, all annual cohorts as a proportion of 1.
  # S <- model_matrix()     # S - susceptible (non-T1D)
  # P <- model_matrix()     # P - prevalence (T1D)
  # D <- model_matrix()     # D - deaths (absorbing state)
  #
  # # incidence flows are reflected the following year, for individuals one year older
  # i_shift <- model_matrix()
  # i_shift[, 1] <- 0
  # i_shift[,-1] <- i[,-MAX_AGE]
  # i_all <- i / (1-dDx)              # includes deaths at onset
  # i_all_shift <- i_shift / (1-dDx)  # shifted & includes deaths at onset
  #
  # # Track P cohorts and in/outflows in 3D array: {year}x{age}x{onset age}
  # Pcohorts <- model_array()
  # Icohorts <- model_array()
  # PDcohorts <- model_array()
  #
  # # The model proceeds yearwise, populating successive matrix rows.
  # # Equation references below refer to the model summary.
  # S[1,1] <- 1
  # P[1,1] <- Pcohorts[1,1,1] <- D[1,1] <- 0
  # Icohorts[,,] <- PDcohorts[,,] <- 0
  # for (t in seq_along(years[-1])) {
  #   # for (t in 1:10) {
  #   # bump previous period's data along one cell (ie one year older)
  #   Sshift <- c(1, S[t, -MAX_AGE])
  #   Pshift <- c(0, P[t, -MAX_AGE])
  #   Dshift <- c(0, D[t, -MAX_AGE])
  #
  #   PCshift <- array(0, dim=c(MAX_AGE, MAX_AGE))
  #   PCshift[1,] <- 0
  #   PCshift[-1,] <- Pcohorts[t,-MAX_AGE,]
  #
  #   Ishift <- matrix(NA, nrow = MAX_AGE, ncol = MAX_AGE)
  #   Ishift[,-MAX_AGE] <- diag(i_shift[t,] * Sshift)[,-1] # Only want to shift age, not onset age
  #   Ishift[,MAX_AGE] <- 0  # NB: half-cycle adjustment populates any incidence for MAX_AGE
  #
  #   # equation (3) - susceptible compartment, S
  #   S[t+1,] <- Sshift * (1 - i_all_shift[t,]) * (1 - qB[t,])
  #
  #   # equation (4) - prevalence compartment P
  #   # P[t+1,] <- Pshift * (1 - qT1D[t,]) + i_shift[t,] * Sshift
  #
  #   P[t+1,] <-    qT1D_percent_n[t,]  *  Pshift  * (1 - qT1D_n[t,]) + i_shift[t,] * Sshift +
  #     (1-qT1D_percent_n[t,]) *  Pshift  * (1 - qT1D_m[t,])
  #
  #
  #
  #   # equation (4') - prevalence compartment P by cohort
  #   Icohorts[t,,] <- diag(i[t,] * S[t,])
  #   # PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D[t,], c(MAX_AGE, MAX_AGE))
  #   PDcohorts[t,,] <- qT1D_percent_n[t,] * Pcohorts[t,,] * array(qT1D_n[t,], c(MAX_AGE, MAX_AGE)) +
  #     (1 - qT1D_percent_n[t,]) * Pcohorts[t,,] * array(qT1D_m[t,], c(MAX_AGE, MAX_AGE))
  #
  #   # Pcohorts[t+1,,] <- PCshift * array(1 - qT1D[t,], c(MAX_AGE, MAX_AGE)) + Ishift
  #   Pcohorts[t+1,,] <- qT1D_percent_n[t,] * PCshift * array(1 - qT1D_n[t,], c(MAX_AGE, MAX_AGE)) + Ishift +
  #     (1 - qT1D_percent_n[t,]) * PCshift * array(1 - qT1D_m[t,], c(MAX_AGE, MAX_AGE))
  #
  #   # equation (5) - death compartment D
  #   D[t+1,] <- (Dshift
  #               + i_all_shift[t,] * dDx[t,] * Sshift
  #               # + qT1D[t,] * Pshift
  #               + qT1D_n[t,] * Pshift *       qT1D_percent_n[t,]
  #               + qT1D_m[t,] * Pshift *  (1 - qT1D_percent_n[t,])
  #               + Sshift * qB[t,] * (1 - i_all_shift[t,]))
  #
  # }
  # # final year cohort flows - same as in the loop but for final period
  # t <- length(years)
  # Icohorts[t,,] <- diag(i[t,] * S[t,])
  #
  # # PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D[t,], c(MAX_AGE, MAX_AGE))
  # PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D_n[t,], c(MAX_AGE, MAX_AGE)) *      qT1D_percent_n[t,]  +
  #   Pcohorts[t,,] * array(qT1D_m[t,], c(MAX_AGE, MAX_AGE)) * (1 - qT1D_percent_n[t,])
  #
  # # flows: based on unshifted versions of incidence
  # I <- i * S                    # T1D incidence
  # DDx <- i_all * dDx * S        # deaths at T1D onset
  # # DT1D <- (qT1D - qB) * P       # T1D-cause mortality
  # DT1D <- (qT1D_n - qB) * P  *          qT1D_percent_n[t,]  +
  #   (qT1D_m - qB) * P  *     (1 - qT1D_percent_n[t,])   # T1D-cause mortality
  # DBGP <- P * qB                # background mortality for people w/ T1D
  # DBGS <- S * qB * (1 - i_all)  # background mortality for susceptible pop'n
  # # DBGS <- S * qB  # background mortality for susceptible pop'n
  #
  # # Half-cycle adjustments - apply half of each flow in the reference
  # # year for that flow. We expect that on average, 1/2 of the flow
  # # occurs part-way through the year.
  # S <- S - 0.5 * I - 0.5 * DDx - 0.5 * DBGS
  # P <- P + 0.5 * I - 0.5 * (DT1D + DBGP)
  # D <- D + 0.5 * (DDx + DBGS + DT1D + DBGP)
  # Pcohorts <- Pcohorts + 0.5 * Icohorts - 0.5 * PDcohorts
  # # list(S=S, P=P, Pcohorts=Pcohorts, D=D, I=I, DDx=DDx, DT1D=DT1D,DBGP=DBGP, DBGS=DBGS)
}



complication_prevalence_tensor <- function(years,P_cohorts_level,smr_matrix_n, smr_matrix_m,qT1D_percent_n, hba1c_f=NULL,complication_parameters,disease_weights) {
  #   years=1960:2040 ; P_cohorts_level <-  P_cohorts_level[1,,,]; smr_matrix_n <- smr_matrix_n[1,,]; smr_matrix_m <- smr_matrix_m[1,,]; qT1D_percent_n <- qT1D_percent_n[1,,];hba1c_f=NULL
  # nyears          <- length(prev$years)

  start_year      <- min(years)
  end_year        <- max(years)

  data_start_year <- start_year-MAX_AGE+1
  all_years       <- seq(data_start_year, end_year)

  # if (is.null(hba1c_f)) {
  #   hba1c_f <- hba1c_function(prev$country)  # dimensions: year x age
  # }
  # hba1c_matrix <- hba1c_f(all_years)

  smr <- smr_matrix_n * qT1D_percent_n + smr_matrix_m * (1- qT1D_percent_n)
  hba1c_matrix <- (log(smr) +  1.5274 )/ 0.3545

  hba1c_matrix <- hba1c_matrix[-1,]


  comp_names              <- with(complication_parameters, c(weibull$abbrev, constant$complication))
  nweib                   <- nrow(complication_parameters$weibull)
  nconst                  <- nrow(complication_parameters$constant)

  # 4D comp_prev array tabulates complication prevalence by time, complication,
  # age, and age at diagnosis. Dimensions are `year` x `complication` x `age` x
  # `age at diagnosis`. Values are numbers of individuals in the reference
  # population. That is, a value of 1,000 would indicate 1,000 individuals in
  # the reference country, at that age and cohort.
  comp_prev <- array(NA,
                     dim=list(length(comp_names), length(years), MAX_AGE, MAX_AGE),
                dimnames=list(comp_names               , years , AGES   , AGES))

  # Weibull time-to-event complications
  risk <- array(0,
                dim=c(length(comp_names), length(years), MAX_AGE, MAX_AGE),
                dimnames = list(comp=comp_names, year=years, age=AGES, cohort=AGES))
  can_reuse_S <- any(('uniform_ages' %in% attributes(hba1c_f)) & c(attr(hba1c_f, 'uniform_ages'), FALSE))

  for (incidence_year in seq(end_year, data_start_year)) {
    # incidence_year <- 2020
    # print(incidence_year)

    for (incidence_age in seq(0, 99 - max(0, start_year - incidence_year))) {
      # incidence_age <- 20
      # print(incidence_age)

      max_age <- min(incidence_age + end_year - incidence_year, 99)
      # hba1c array indexes - potentially starts before start_year because
      # integration starts from age of incidence
      h_length <- max_age - incidence_age + 1
      h_age_idx <- seq(incidence_age + 1, max_age + 1)
      h_year_idx <- seq(incidence_year - data_start_year + 1, incidence_year - data_start_year + h_length)
      # survivor curve starting in incidence_year at age incidence_age
      hba1cs <- hba1c_matrix[cbind(h_year_idx, h_age_idx)]
      # risk indexes are potentially shorter, only starting at start_year
      first_risk_age <- incidence_age + max(start_year - incidence_year, 0)
      first_risk_year <- max(start_year, incidence_year)
      risk_length <- h_length - (first_risk_age - incidence_age)
      keep_risk_idx <- seq(1 + first_risk_year - incidence_year, first_risk_year - incidence_year + risk_length)
      cohort_idx <- incidence_age + 1
      risk_year_idx <- seq(max(0, incidence_year - start_year) + 1, max(0, incidence_year - start_year) + risk_length)
      risk_age_idx <- seq(first_risk_age + 1, first_risk_age + risk_length)
      # risk is the complement of survival
      if(FALSE)
      {
          hba1cs <- c(15,16,17)
          hba1cs <- c(15)

          hba1cs_rep <- rep(hba1cs,10)
          dim(hba1cs_rep) <- c(length(hba1cs_rep)/10, 10)
          hba1cs_rep <- t(hba1cs_rep)

          T <- length(hba1cs)-1 # total time length (years) incl. T=0
          hs_odd <- 0.5 * (hba1cs[-1] + hba1cs[-(T+1)]) # interpolate at half-year points
          ts <- pmin(0:T, 30)  # constant hazard after 30 years

          ts <- rep(ts,10)
          dim(ts) <- c(length(ts)/10, 10)
          ts <- t(ts)

          l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1cs_rep)/scale) / scale
          l_odd <- (ts[,-1] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale

          Lambda <- 1/6 * (l_even[1,] + cumsum(2 * c(0, 0, l_even[1,-c(1, T+1)]) + 4 * c(0, l_odd[1,-(T+1)])))

          exp(-Lambda) # survival function
      }

      inputs <- list()
      intercept <- complication_parameters$weibull$intercept
      slope <- complication_parameters$weibull$slope
      scale <- complication_parameters$weibull$scale
      for (comp_i in seq_len(nweib)) {
        # comp_i <- 1
        # print(comp_i)
        inputs[[comp_i]] <- list(hba1c=hba1cs,intercept=intercept[comp_i],slope=slope[comp_i],scale=scale[comp_i])
      }

      weib_survival <- function(inputs) {
        # inputs <- inputs[[1]]
        hba1c <- inputs$hba1c ; intercept<- inputs$intercept; slope<- inputs$slope; scale<- inputs$scale
        T <- length(hba1c)-1 # total time length (years) incl. T=0
        hs_odd <- 0.5 * (hba1c[-1] + hba1c[-(T+1)]) # interpolate at half-year points
        ts <- pmin(0:T, 30)  # constant hazard after 30 years
        l_even <- ts^((1 - scale)/scale) * exp(-(intercept + slope * hba1c)/scale) / scale
        l_odd <- (ts[-1] - 0.5)^((1 - scale)/scale) * exp(-(intercept + slope * hs_odd)/scale) / scale
        Lambda <- 1/6 * (l_even + cumsum(2 * c(0, 0, l_even[-c(1, T+1)]) + 4 * c(0, l_odd[-(T+1)])))
        exp(-Lambda) # survival function
      }
      S <-  lapply( inputs, weib_survival )
      # do.call(rbind, S)

      # hba1c <- inputs$hba1c ; intercept<- inputs$intercept; slope<- inputs$slope; scale<- inputs$scale
      # T <- length(hba1c)-1 # total time length (years) incl. T=0
      # hs_odd <- 0.5 * (hba1c[-1] + hba1c[-(T+1)])
     # pweibull(hba1c, shape = -slope[1], scale = scale[1])

      for (comp_i in seq_len(nweib)) {
        # comp_i <- 1
        # print(comp_i)
        risk[cbind(comp=comp_i, year=risk_year_idx, age=risk_age_idx, cohort=cohort_idx)] <- (1 - S[[comp_i]])[keep_risk_idx]
      }
    }
  }

  # This is Graham's and Gabriel's 'journeys' adjustments.
  # Reduces ON and PR by survivor shares of RF and BL, respectively,
  # effectively making this a multistate model (rather than a pure
  # time-to-event model)
  risk['ON',,,] <- risk['ON',,,] * (1 - risk['RF',,,])
  risk['PR',,,] <- risk['PR',,,] * (1 - risk['BL',,,])
  for (comp_i in seq_len(nweib)) {
    comp_prev[comp_i,,,] <- P_cohorts_level * risk[comp_i,,,]
  }

  # constant risk complications
  for (comp_i in seq_len(nconst)) {
    const_risk <- complication_parameters$constant[comp_i,]
    lt8 <- 1e-2 * const_risk$hba1c_lt_8
    gt8lt9 <- 1e-2 * const_risk$hba1c_8_to_9
    gt9 <- 1e-2 * const_risk$hba1c_gteq_9
    cr <- function(h) ifelse(h < lt8, lt8, ifelse(h < 9, gt8lt9, gt9))
    rsk <- apply(hba1c_matrix[-(1:99),], MARGIN=1:2, FUN=cr)
    for (cohort_i in seq_len(MAX_AGE)) {
      risk[nweib + comp_i,,,cohort_i] <- rsk
    }
    comp_prev[nweib + comp_i,,,] <- risk[nweib + comp_i,,,] * P_cohorts_level
  }
  # structure(list(
  #   year_comp_age_cohort=comp_prev,
  #   year_comp_age_cohort_prev_risk=risk
  # ), class='complications')
  #
  props <- risk

  disease_weights <- as.list(disease_weights)
  t1d <- disease_weights$T1D

  # Calculate disability weights to apply to prevalence. Since 100% of prevalent
  # cases have T1D, we apply the T1D disability weight directly. For all other
  # complications, we scale by the modeled complication probability.
  # Terminology note: disease weight == disability weight
  dsp <- props["DSP",,,] * disease_weights$DSP
  pr <- props["PR",,,] * disease_weights$PR + props["BL",,,] * disease_weights$BL
  on <- props["ON",,,] * disease_weights$ON
  #rf <- (props["RF_transplant",,,] * disease_weights$Transplant
  #  + props["RF_dialysis",,,] * disease_weights$Dialysis)
  rf <- props["RF",,,] * disease_weights$Dialysis
  uoa <- props["UoA",,,] * disease_weights$UoA
  hom <- props["HoM",,,] * disease_weights$HoM
  nfMI <- props["nfMI",,,] * disease_weights$nfMI
  nfCBVD <- props["nfCBVD",,,] * disease_weights$nfCBVD
  # disability weights are combined multiplicatively and applied to prevalence level
  disability_wts <- (
    1 - (1 - t1d) * (1 - dsp) * (1 - pr) * (1 - on - rf) * (1 - uoa) *
      (1 - hom) * (1 - nfMI) * (1 - nfCBVD))
  disability <- disability_wts * P_cohorts_level
  structure(
    list(
      year=apply(disability, 1, sum),
      year_age=apply(disability, c(1, 2), sum),
      year_cohort=apply(disability, c(1, 3), sum),
      year_age_cohort = disability
    ),
    class='DALYs')

}



