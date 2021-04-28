library(testthat)
library(FLOPART)
data_vec <- 1:5
data_count <- length(data_vec)
weight_vec <- as.double(rep(1, data_count))
penalty <- as.double(0)
cost_mat <- numeric(data_count * 2)
end_vec <- integer(data_count)
mean_vec <- numeric(data_count)
intervals_mat <- integer(data_count*2)
label_starts <- 1L
label_ends <- 3L
label_count <- length(label_starts)
PoissonLoss <- function(meanVal, dataVal){
  meanVal - dataVal * log(meanVal)
}
seg <- function(mean, start, end){
  data.frame(mean, start, end)
}

label_types <- 0L
result <- FLOPART::FLOPART_interface(
  data_vec, weight_vec, penalty,
  label_types, label_starts, label_ends)
costMatrix <- t(result[["cost_mat"]])
test_that("correct for noPeaks label",{
  expected.segs <- rbind(
    seg(2.5, 1, 4),
    seg(5.0, 5, 5))
  expect_equal(result[["segments_df"]], expected.segs)
  last.up.cost <- mean(PoissonLoss(c(rep(2.5, 4), 5), data_vec))
  first.up.cost <- PoissonLoss(1, 1)
  expected.up.cost <- c(first.up.cost, rep(Inf, 3), last.up.cost)
  expect_equal(costMatrix[1,], expected.up.cost)
})

label_types <- 1L
result <- FLOPART::FLOPART_interface(
  data_vec, weight_vec, penalty,
  label_types, label_starts, label_ends)
costMatrix <- t(result[["cost_mat"]])
test_that("first of peak start cost and last of peak start down costs are infinite", {
  expect_true(is.infinite(costMatrix[1,2]))
  expect_true(is.infinite(costMatrix[2,4]))
})

label_types <- -1L
result <- FLOPART::FLOPART_interface(
  data_vec, weight_vec, penalty,
  label_types, label_starts, label_ends)
costMatrix <- t(result[["cost_mat"]])
test_that("last of peak end down cost and first of peak end up cost are infinite",{
  expect_true(is.infinite(costMatrix[2,2]))
  expect_true(is.infinite(costMatrix[1,4]))
})

