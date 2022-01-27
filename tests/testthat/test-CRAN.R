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
seg <- function(mean, firstRow, lastRow, state){
  data.frame(mean, firstRow, lastRow, state)
}

label_types <- 0L
result <- FLOPART::FLOPART_interface(
  data_vec, weight_vec, penalty,
  label_types, label_starts, label_ends)
costMatrix <- t(result[["cost_mat"]])
test_that("correct for noPeaks label",{
  expected.segs <- rbind(
    seg(2.5, 1, 4, 0),
    seg(5.0, 5, 5, 5))
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

label_code <- FLOPART::get_label_code()
test_that("get_label_code works", {
  expect_true(all(c("noPeaks", "peakStart", "peakEnd") %in% names(label_code)))
})

label <- function(type, start, end){
  data.frame(type=label_code[type], start, end, row.names=NULL)
}
two.labels <- rbind(
  label("noPeaks", 2, 3),
  label("peakEnd", 3, 4))
test_that("error for overlapping labels", {
  expect_error({
    with(two.labels, FLOPART::FLOPART_interface(
      data_vec, weight_vec, penalty,
      type, start, end))
  }, "label start should be greater than previous label end")
})

test_that("error for label start before coverage", {
  cov.df <- data.frame(chromStart=10, chromEnd=20, count=4)
  lab.df <- data.frame(chromStart=5, chromEnd=15, annotation="noPeaks")
  expect_error({
    FLOPART::FLOPART_data(cov.df, lab.df)
  }, "label starts must be on or after first coverage")
})

test_that("error for label end after coverage", {
  cov.df <- data.frame(chromStart=10, chromEnd=20, count=4)
  lab.df <- data.frame(chromStart=15, chromEnd=25, annotation="noPeaks")
  expect_error({
    FLOPART::FLOPART_data(cov.df, lab.df)
  }, "label ends must be on or before last coverage")
})

test_that("can start up", {
  seg.mean.vec <- c(11, 3)
  set.seed(1)
  count <- unlist(lapply(seg.mean.vec, function(m)rpois(10, m)))
  N <- length(count)
  cov.df <- data.frame(
    chromStart=seq(0, N-1),
    chromEnd=seq(1, N),
    count)
  result.list <- FLOPART::FLOPART(cov.df, penalty=5)
  segs <- result.list[["segments_dt"]]
  expect_identical(segs[["status"]], c("peak", "background"))
})

test_that("factor annotation ok", {
  seg.mean.vec <- c(11, 3)
  set.seed(1)
  count <- unlist(lapply(seg.mean.vec, function(m)rpois(10, m)))
  N <- length(count)
  cov.df <- data.frame(
    chromStart=seq(0, N-1),
    chromEnd=seq(1, N),
    count)
  label.df <- data.frame(
    chromStart=12,
    chromEnd=18,
    annotation=factor(
      "noPeaks", levels=c("peaks", "peakStart", "peakEnd", "noPeaks")))
  result.list <- FLOPART::FLOPART(cov.df, label.df, penalty=5)
  segs <- result.list[["segments_dt"]]
  expect_identical(segs[["status"]], c("peak", "background"))
})

test_that("error for non-df coverage", {
  expect_error({
    FLOPART::FLOPART(c(3, 4), penalty=5)
  }, "coverage must be a data frame")
})