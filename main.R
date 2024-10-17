
# Load Library ---------------------------------------------------------------------------------------
library(arrow)
library(data.table)
library(dplyr)
library(abind)
library(stringr)

# Load Dependency ---------------------------------------------------------------------------------------

source("code_R/refresh_country_files_tensor.R")
source("code_R/refresh_country_files_tensor_utils.R")
source("code_R/utils.R")
source("code_R/data.R")

# Read Index input ---------------------------------------------------------------------------------------

index_data_input  <- setDF(arrow::read_feather("data_prep_temp/inputs_national_v_2_0/partition_1"))


# Run Index get output  ---------------------------------------------------------------------------------------
index_data_output <- refresh_country_files_tensor(index_data_input,sim_enable=FALSE)
