# soft impute tools

#   SOftimpute -----------------------------------------------------------------
impute_matrix_completion <- function(Y,lambda=0.2)
{

  # Y <- country_indicator_wide_gdp

  Y <- as.matrix(Y)
  Y_scale <- Y
  tryCatch({
    Y_scale <- softImpute::biScale(Y, maxit = 100)
  },error=function(e){Y_scale <- Y})

  fit1 <- softImpute::softImpute(Y_scale, rank.max = min(dim(Y)) - 1,  lambda =lambda )

  Y_predict <- softImpute::complete(Y,fit1)

  Y_predict[Y_predict<0]<- 0
#
  return(Y_predict)

}



soft_cv <- function(Y,
                    k = 10,
                    lambda_grid = seq(0.1,1, 0.02),
                    print_update = FALSE
) {

  # svdY <- svd(Y)
  # d <- svdY$d
  # if(is.null(lambda_grid)) {
  #   lambda_grid <- seq(min(d), max(d), length = 20)
  # }
  k = 10
  lambda_grid = seq(0.1,1, 0.02)
  # lambda_grid = seq(1,15, 1)


  num_cell <- prod(dim(Y))
  kfold_index <- sample(1:k,num_cell,replace = TRUE)

  for(i in 1:length(lambda_grid))
  {
    # i <- 18
    print(i)
    lambda <- lambda_grid[i]

    Rsquared_per_lambda <- c()
    mean_error <- c()
    for(j in 1:k)
    {
      # j <- 1
      Y_train <- Y
      Y_train[kfold_index==j] <- NA


      Y_train <- softImpute::biScale(Y_train, maxit = 100)
      fit1 <- softImpute::softImpute(Y_train,
                                     rank.max = min(dim(Y)) - 1,
                                     lambda =lambda)
      Y_predict <- softImpute::complete(Y_train, fit1)

      # Y_predict <- impute_By_Age_Ratio_Naive(Y_train)



      index_predict <- is.na(Y_train)&(!is.na(Y))& (kfold_index==j)
      Rsquared <- cor(Y[index_predict], Y_predict[index_predict]) ^ 2
      Rsquared_per_lambda <- c(Rsquared_per_lambda,Rsquared)


      mean_error <- c(mean_error,mean(abs(Y[index_predict] - Y_predict[index_predict])))

    }

    print(mean(Rsquared_per_lambda,na.rm= TRUE))
    print(mean(mean_error))



  }


}
