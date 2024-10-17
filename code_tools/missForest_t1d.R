missForest_t1d <- function(impute_input, iteration = 6,node_size=1)
{
  # impute_input <- impute_input_train
  setDF(impute_input)
  response_name <- colnames(impute_input)[colSums(is.na(impute_input))>0]

  index_train <- !is.na(impute_input[,response_name])

  for(i in iteration)
  {
    rf <- randomForest( impute_input[index_train,setdiff(colnames(impute_input),response_name) ]   , (impute_input[index_train,response_name]), ntree = 1000,nodesize = node_size)
    impute_input[,paste0(response_name,"_imputed") ] <-  (predict(rf,data.frame(impute_input[,setdiff(colnames(impute_input),response_name) ]),type="response"))
    R2 <-  cor(impute_input$smr[index_train], impute_input[,paste0(response_name,"_imputed")][index_train]) ^ 2
    R2

    # rf <- randomForest( impute_input[index_train,setdiff(colnames(impute_input),response_name) ]   , (impute_input[index_train,response_name]), ntree = 1000,nodesize = node_size)
    # impute_input[,paste0(response_name,"_imputed") ] <-  (predict(rf,data.frame(impute_input[,setdiff(colnames(impute_input),response_name) ]),type="response"))
    # R2 <-  cor(impute_input$smr[index_train], impute_input[,paste0(response_name,"_imputed")][index_train]) ^ 2
    # R2
    #
    # rf <- randomForest( impute_input[index_train,setdiff(colnames(impute_input),response_name) ]   , (impute_input[index_train,response_name]), ntree = 1000,nodesize = node_size)
    # impute_input[,paste0(response_name,"_imputed") ] <-  (predict(rf,data.frame(impute_input[,setdiff(colnames(impute_input),response_name) ]),type="response"))
    # R2 <-  cor(impute_input$smr[index_train], impute_input[,paste0(response_name,"_imputed")][index_train]) ^ 2
    # R2
    #
    # rf <- randomForest( impute_input[index_train,setdiff(colnames(impute_input),response_name) ]   , (impute_input[index_train,response_name]), ntree = 1000,nodesize = node_size)
    # impute_input[,paste0(response_name,"_imputed") ] <-  (predict(rf,data.frame(impute_input[,setdiff(colnames(impute_input),response_name) ]),type="response"))
    # R2 <-  cor(impute_input$smr[index_train], impute_input[,paste0(response_name,"_imputed")][index_train]) ^ 2
    # R2

  }



  return(impute_input)

}
