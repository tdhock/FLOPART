library(FLOPART)
data_count = 5L
data_vec <- as.integer(c(1,2,3,4,5))
weight_vec <- as.double(rep(1, data_count))
penalty <- as.double(0)
cost_mat <- numeric(data_count * 2)
end_vec <- integer(data_count)
mean_vec <- numeric(data_count) 
intervals_mat <- integer(data_count*2)
label_starts <- 1L
label_ends <- 3L
label_count <- 1L

PoissonLoss <- function(meanVal, dataVal){
  meanVal - dataVal * log(meanVal)
  
}


test_that("costs are correct for noPeaks label",
          {
            label_types <- 0L
              .C("FLOPART_interface", data_vec = data_vec, weight_vec = weight_vec,
                 data_count = data_count, penalty = penalty, cost_mat = cost_mat,
                 end_vec = end_vec, mean_vec = mean_vec, intervals_mat = intervals_mat,
                 label_starts = label_starts, label_ends = label_ends,
                 label_types = label_types, label_count = label_count,
                 PACKAGE="FLOPART")
            
            costMatrix <- matrix(result[["cost_mat"]], nrow=2, ncol=data_count, 
                                 byrow=TRUE)
            
            expect_true(is.infinite(costMatrix[1,2]))
            expect_true(is.infinite(costMatrix[1,3]))
            expect_true(is.infinite(costMatrix[1,4]))
          })
test_that("first of peak start cost and last of peak start down costs are infinite",
          {
            label_types <- 1L
            result <- 
              .C("FLOPART_interface", data_vec = data_vec, weight_vec = weight_vec,
                 data_count = data_count, penalty = penalty, cost_mat = cost_mat,
                 end_vec = end_vec, mean_vec = mean_vec, intervals_mat = intervals_mat,
                 label_starts = label_starts, label_ends = label_ends,
                 label_types = label_types, label_count = label_count,
                 PACKAGE="FLOPART")
            
            costMatrix <- matrix(result[["cost_mat"]], nrow=2, ncol=data_count, 
                                 byrow=TRUE)
            
            expect_true(is.infinite(costMatrix[1,2]))
            expect_true(is.infinite(costMatrix[2,4]))
          })
test_that("last of peak end down cost and first of peak end up cost are infinite",
          {
            label_types <- -1L
            result <- 
              .C("FLOPART_interface", data_vec = data_vec, weight_vec = weight_vec,
                 data_count = data_count, penalty = penalty, cost_mat = cost_mat,
                 end_vec = end_vec, mean_vec = mean_vec, intervals_mat = intervals_mat,
                 label_starts = label_starts, label_ends = label_ends,
                 label_types = label_types, label_count = label_count,
                 PACKAGE="FLOPART")
            
            costMatrix <- matrix(result[["cost_mat"]], nrow=2, ncol=data_count, 
                                 byrow=TRUE)
            
            expect_true(is.infinite(costMatrix[2,2]))
            expect_true(is.infinite(costMatrix[1,4]))
          })

