CV_randomForest <- function(model.data,features_for_training, response_name="smr",Nfold = 10,nodesize = 5)
{

  # model.data <- d_trian; # model.data$smr <- log(model.data$smr+1)

  # CV_randomForest(d_trian, features_for_training,response= "smr",10,nodesize = 5)

  # fold_key="unique_key";response_name="smr";  model.data          <- model_data_global_studies  ;nodesize = 5 ,Nfold = 10
  model.data$response <- model.data[,response_name]
  # set.seed(11)
  # mtry <- 20
  # Nfold    <- floor(length(unique(model.data[,fold_key])) /7)
  Nfold_df <- unique(model.data[,"unique_key",drop=FALSE] )

  Nfold_df$Nfold <- (1:nrow(Nfold_df) %% Nfold)+1
  model.data <- model.data  %>% left_join(Nfold_df, by='unique_key')

  valid_predicted_all  <- data.frame()
  train_predicted_all  <- data.frame()
  R2list_train   <- c();    R2list_predict <- c()
  RMSElist_train   <- c();    RMSElist_predict <- c()
  MAPElist_predict <- c()

  for(i in 1:Nfold)
  {
    print(i)
    # i <- 4;  set.seed(11)
    # paste0("CV fold: ",print(i))

    train.data <- model.data[model.data$Nfold!=i,] ; train.data$Nfold <- NULL
    valid.data <- model.data[model.data$Nfold==i,] ; valid.data$Nfold <- NULL

    train.data_ <-  train.data[,features_for_training]
    rf <- randomForest(train.data_, train.data$response, ntree = 1000,nodesize = 2)

    # rf <- glm.nb(response ~ doctors_per_capita , data = train.data)

    valid.data$predicted <- predict(rf,valid.data[,features_for_training],type="response")
    valid_predicted_all  <- rbind(valid_predicted_all,valid.data)
    train.data$predicted <- predict(rf,setDF(train.data[,features_for_training]),type="response")
    train_predicted_all  <- rbind(train_predicted_all,train.data)

    RMSElist_train   <-  c(RMSElist_train    ,rmse(train.data$response, train.data$predicted))
    RMSElist_predict <-  c(RMSElist_predict  ,rmse(valid.data$response, valid.data$predicted))
    MAPElist_predict <-  c(MAPElist_predict  ,mean(abs(( (valid.data$response+1)- (valid.data$predicted+1) )/(valid.data$response+1) )) * 100)


    R2list_train     <-  c(R2list_train  ,cor(train.data$response, train.data$predicted) ^ 2*100)
    R2list_predict   <-  c(R2list_predict,cor(valid.data$response, valid.data$predicted) ^ 2*100)

  }

  # r squared
  R2_t <- cor(train_predicted_all$response, train_predicted_all$predicted) ^ 2*100
  R2_p <- cor(valid_predicted_all$response, valid_predicted_all$predicted) ^ 2*100

  index_non_zero  <- valid_predicted_all$response !=0
  list_x <- ( valid_predicted_all$predicted[index_non_zero]- valid_predicted_all$response[index_non_zero] )/ valid_predicted_all$response[index_non_zero]*100
  ci_pc_lower  <- median((list_x[list_x<0]) )
  ci_pc_upper  <- median(abs(list_x[list_x>0]) )
  # quantile(list_x,probs=c(.025,0.5,.975))
  print(paste0("ci_pc_lower: ",ci_pc_lower))
  print(paste0("ci_pc_upper: ",ci_pc_upper))
  #  rmse

  rmse(train_predicted_all$response, train_predicted_all$predicted)
  rmse(valid_predicted_all$response, valid_predicted_all$predicted)
#  variance

  median(1-abs(train_predicted_all$response - train_predicted_all$predicted) / train_predicted_all$predicted)
  median(1-abs(valid_predicted_all$response - valid_predicted_all$predicted) / valid_predicted_all$predicted)

  # randomForest::importance(rf)
  # print(paste0("CV R2 train   "  ,R2_t ))
  # print(paste0("CV R2 predict   ",R2_p ))

  # randomForest::importance(rf)
  # print(paste0("CV RMSE train   ",round(mean(RMSElist_train) ,2) ))
  # print(paste0("CV RMSE predict   ",round(mean(RMSElist_predict) ,2) ))

}
