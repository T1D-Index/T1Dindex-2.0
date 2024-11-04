
# Load Library ---------------------------------------------------------------------------------------
library(arrow)
library(data.table)
library(dplyr)
library(abind)
library(stringr)
library(parallel)

# Load Dependency ---------------------------------------------------------------------------------------

source("code_R/refresh_country_files_tensor.R")
source("code_R/refresh_country_files_tensor_utils.R")
source("code_R/utils.R")
source("code_R/data.R")
source("code_R/runner.R")

# Read Index input ---------------------------------------------------------------------------------------

index_data_input  <- setDF(arrow::read_feather("data_prep_temp/inputs_national_v_2_0/partition_1"))


# Run Index get output  ---------------------------------------------------------------------------------------
system.time(
  {
    index_data_output <- refresh_country_files_tensor(index_data_input,sim_enable=FALSE)
  }
)

# Parallel ---------------------------------------------------------------------------
data_dir <- "data_prep_temp/inputs_national_v_2_0/"
cache_dir <- "data_prep_temp/inputs_national_v_2_0_output/"
loc_id_input_list <-list.files(data_dir)

Log=function(fmt, ...) { cat(sprintf(paste0(fmt, '\n'), ...)) }

num_thread <- 15
clust <- parallel::makeCluster(num_thread, setup_strategy = "sequential", outfile="makeCluster_Log.txt")
clusterExport(cl=clust, varlist=c('Log','setDF','rbindlist','setDT','str_extract','read_feather'
                                  ,'cache_dir','data_dir','run_query_df',"adrop"
                                  ,'abind','refresh_country_files_tensor','calculate_ex_matrix','calculate_ex_lifetime_years_lost_matrix','calculate_prevalence_tensor'
                                  ,'get_loc_id','get_database_connection'
                                  ,'get_prevalence_reference'

                                  ,'get_incidence_curve'
                                  ,'get_incidence_growth_rates'

                                  ,'matrix_from_function','make_age_function'
                                  ,'get_complication_parameters','get_hba1c_assumption'

                                  ,'get_disease_weights','calculate_ex','calculate_ex_lifetime_years_lost','assert'))

system.time({a <- clusterApplyLB(clust, loc_id_input_list, refresh_one_country_file,data_dir=data_dir,cache_dir=cache_dir,use_3d_array=TRUE,run_type="national")})


stopCluster(clust)
Sys.time()
