library(FLOPART)
data_count = 5L
data_vec <- 1:data_count
weight_vec <- 1:data_count
penalty <- 0L
cost_mat <- numeric(data_count * 2)
end_vec <- integer(data_count)
mean_vec <- numeric(data_count)
intervals_mat <- integer(data_count*2)
label_starts <- 1L
label_ends <- 3L
label_types <- 0L
label_count <- 1L
result <- 
  .C("PeakSegFPOPLog_interface", data_vec = data_vec, weight_vec = weight_vec,
     data_count = data_count, penalty = penalty, cost_mat = cost_mat,
   end_vec = end_vec, mean_vec = mean_vec, intervals_mat = intervals_mat,
   label_starts = label_starts, label_ends = label_ends,
   label_types = label_types, label_count = label_count,
   PACKAGE="FLOPART")

costMatrix <- matrix(result[["cost_mat"]], data_count)


#PoissonLoss <- function(int meanVal, int dataVal){
  #meanVal - dataVal * log(meanVal)
  
#}