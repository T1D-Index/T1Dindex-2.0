
impute_By_Age_Ratio_Naive <- function(matrix)
{
  # matrix = country_indicator_wide_gdp

  matrix_impute <- matrix
  for(row_i in 1:nrow(matrix))
  {
    # print(rownames(matrix) [row_i])
    # row_i <- 1
    for(col_j in 1:ncol(matrix))
    {
      # col_j <- 1
      if(is.na(matrix[row_i, col_j]))
      {
        rows_with_data <- which(!is.na(matrix[,col_j]))
        cols_with_data <- which(!is.na(matrix[row_i,]))

        col_vector     <- matrix[rows_with_data,col_j,drop=FALSE]
        row_vector     <- matrix[row_i,cols_with_data,drop=FALSE]

        matrix_temp    <- matrix[rows_with_data, cols_with_data,drop=FALSE]



        index_not_null <- colSums(is.na(matrix_temp),2)==0

        if(sum(index_not_null)==0)
        {
          # remove rows that are all empty
          index_not_null_row <- rowSums(is.na(matrix_temp),2)==0

          col_vector <- col_vector[ index_not_null_row ,,drop=FALSE ]

          matrix_temp    <- matrix_temp[index_not_null_row,,drop=FALSE]

        }

        index_not_null <- colSums(is.na(matrix_temp),2)==0


        matrix_temp    <- matrix_temp[,  index_not_null,drop=FALSE]

        row_vector     <- row_vector[ index_not_null]

        # matrix_impute[row_i,col_j] <- sum(row_vector) * (sum(col_vector)/sum(matrix_temp))

        matrix_impute[row_i,col_j] <- sum(row_vector) * mean(unlist(col_vector/rowSums(matrix_temp)))

        if(is.nan( matrix_impute[row_i,col_j] ))
        {
          matrix_impute[row_i,col_j] <- NA
        }



      }
    }
  }

  return(matrix_impute)
}
