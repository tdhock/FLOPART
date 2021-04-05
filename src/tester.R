library(FLOPART)

data_to_test <- 1:10
weight_vec <- 1:10
data_count <- 10
penalty <- 10
cost_mat <- 1:10
end_vec <- 1:10
mean_vec <- 1:10
intervals_mat <- 1:10
label_starts <- as.integer(3)
label_ends <- as.integer(4)
label_types <- as.integer(0)
label_count <- 3
.C("PeakSegFPOPLog_interface", data_to_test, weight_vec, as.integer(data_count), penalty, cost_mat = cost_mat,
   end_vec, mean_vec, intervals_mat, label_starts, label_ends, label_types, label_count,
   PACKAGE="FLOPART")

costMatrix <- matrix(cost_mat, data_count, 2)

