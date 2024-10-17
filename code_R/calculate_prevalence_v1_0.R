calculate_prevalence_v1_0 <- function(i, qB,  qT1D_n, qT1D_m, qT1D_percent_n, dDx, years) {

  # Data matrixes: years on rows, ages on columns. Memory footprint of storing
  # all this information really is tiny, a few MB. Similarly, a 3D array: {year}x{age}x{onset age}
  # Model compartments, all annual cohorts as a proportion of 1.
  S <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # S - susceptible (non-T1D)
  P <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # P - prevalence (T1D)
  D <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))     # D - deaths (absorbing state)

  # incidence flows are reflected the following year, for individuals one year older
  i_shift <- matrix(NA, nrow=length(years), ncol=MAX_AGE, dimnames=list(years, AGES))
  i_shift[, 1] <- 0
  i_shift[,-1] <- i[,-MAX_AGE]
  i_all <- i / (1-dDx)              # includes deaths at onset
  i_all_shift <- i_shift / (1-dDx)  # shifted & includes deaths at onset

  # Track P cohorts and in/outflows in 3D array: {year}x{age}x{onset age}
  Pcohorts  <- array(NA, dim=list(length(years), MAX_AGE, MAX_AGE), dimnames=list(years, AGES, AGES))
  Icohorts  <- array(0, dim=list(length(years), MAX_AGE, MAX_AGE), dimnames=list(years, AGES, AGES))
  PDcohorts <- array(0, dim=list(length(years), MAX_AGE, MAX_AGE), dimnames=list(years, AGES, AGES))

  # The model proceeds yearwise, populating successive matrix rows.
  # Equation references below refer to the model summary.
  S[1,1] <- 1
  P[1,1] <- Pcohorts[1,1,1] <- D[1,1] <- 0
  for (t in seq_along(years[-1])) {
    # for (t in 1:90) {
    # bump previous period's data along one cell (ie one year older)
    Sshift <- c(1, S[t, -MAX_AGE])
    Pshift <- c(0, P[t, -MAX_AGE])
    Dshift <- c(0, D[t, -MAX_AGE])

    # equation (3) - susceptible compartment, S
    S[t+1,] <- Sshift * (1 - i_all_shift[t,]) * (1 - qB[t,])

    # equation (4) - prevalence compartment P
    # P[t+1,] <- Pshift * (1 - qT1D[t,]) + i_shift[t,] * Sshift

    P[t+1,] <-    qT1D_percent_n[t,]  *  Pshift  * (1 - qT1D_n[t,]) + i_shift[t,] * Sshift +
      (1-qT1D_percent_n[t,]) *  Pshift  * (1 - qT1D_m[t,])


    if(config$run_days_lost)  # do not run levers for paper stats
    {
      PCshift <- array(0, dim=c(MAX_AGE, MAX_AGE))
      PCshift[1,] <- 0
      PCshift[-1,] <- Pcohorts[t,-MAX_AGE,]

      Ishift <- matrix(NA, nrow = MAX_AGE, ncol = MAX_AGE)
      Ishift[,-MAX_AGE] <- diag(i_shift[t,] * Sshift)[,-1] # Only want to shift age, not onset age
      Ishift[,MAX_AGE] <- 0  # NB: half-cycle adjustment populates any incidence for MAX_AGE

      # equation (4') - prevalence compartment P by cohort
      Icohorts[t,,] <- diag(i[t,] * S[t,])
      # PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D[t,], c(MAX_AGE, MAX_AGE))
      PDcohorts[t,,] <- qT1D_percent_n[t,] * Pcohorts[t,,] * array(qT1D_n[t,], c(MAX_AGE, MAX_AGE)) +
        (1 - qT1D_percent_n[t,]) * Pcohorts[t,,] * array(qT1D_m[t,], c(MAX_AGE, MAX_AGE))

      # Pcohorts[t+1,,] <- PCshift * array(1 - qT1D[t,], c(MAX_AGE, MAX_AGE)) + Ishift
      Pcohorts[t+1,,] <- qT1D_percent_n[t,] * PCshift * array(1 - qT1D_n[t,], c(MAX_AGE, MAX_AGE)) + Ishift +
        (1 - qT1D_percent_n[t,]) * PCshift * array(1 - qT1D_m[t,], c(MAX_AGE, MAX_AGE))
    }
    # equation (5) - death compartment D
    D[t+1,] <- (Dshift
                + i_all_shift[t,] * dDx[t,] * Sshift
                # + qT1D[t,] * Pshift
                + qT1D_n[t,] * Pshift *       qT1D_percent_n[t,]
                + qT1D_m[t,] * Pshift *  (1 - qT1D_percent_n[t,])
                + Sshift * qB[t,] * (1 - i_all_shift[t,]))

  }
  # final year cohort flows - same as in the loop but for final period
  t <- length(years)
  Icohorts[t,,] <- diag(i[t,] * S[t,])

  # PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D[t,], c(MAX_AGE, MAX_AGE))
  PDcohorts[t,,] <- Pcohorts[t,,] * array(qT1D_n[t,], c(MAX_AGE, MAX_AGE)) *      qT1D_percent_n[t,]  +
    Pcohorts[t,,] * array(qT1D_m[t,], c(MAX_AGE, MAX_AGE)) * (1 - qT1D_percent_n[t,])

  # flows: based on unshifted versions of incidence
  I <- i * S                    # T1D incidence
  DDx <- i_all * dDx * S        # deaths at T1D onset
  # DDx <- i_all * S - I        # deaths at T1D onset
  # DT1D <- (qT1D - qB) * P       # T1D-cause mortality
  if(config$lever_change_start_at==1860)
  {
    DT1D <- (qT1D_n - qB) * P  *          qT1D_percent_n[t,]  +
      (qT1D_m - qB) * P  *     (1 - qT1D_percent_n[t,])   # T1D-cause mortality
  }else
  {
    DT1D <- (qT1D_n - qB) * P  *          qT1D_percent_n  +
      (qT1D_m - qB) * P  *     (1 - qT1D_percent_n)   # T1D-cause mortality
  }

  DBGP <- P * qB                # background mortality for people w/ T1D
  DBGS <- S * qB * (1 - i_all)  # background mortality for susceptible pop'n
  # DBGS <- S * qB  # background mortality for susceptible pop'n

  # Half-cycle adjustments - apply half of each flow in the reference
  # year for that flow. We expect that on average, 1/2 of the flow
  # occurs part-way through the year.
  S <- S - 0.5 * I - 0.5 * DDx - 0.5 * DBGS
  P <- P + 0.5 * I - 0.5 * (DT1D + DBGP)
  D <- D + 0.5 * (DDx + DBGS + DT1D + DBGP)
  Pcohorts <- Pcohorts + 0.5 * Icohorts - 0.5 * PDcohorts
  list(S=S, P=P, Pcohorts=Pcohorts, D=D, I=I, DDx=DDx, DT1D=DT1D,DBGP=DBGP, DBGS=DBGS)
}
