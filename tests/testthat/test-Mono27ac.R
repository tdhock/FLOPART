data("Mono27ac", package="FLOPART")
library(ggplot2)
library(data.table)
Mono27ac$coverage[,index:=seq_along(chromStart)]
ggplot(data=Mono27ac$coverage)+
  geom_step(aes(index,count), color="grey50")+
  geom_rect(aes(xmin=150, xmax=250, ymin=-Inf, ymax=Inf),color="red",
            fill=NA, alpha=0.5)
library(FLOPART)
data_count = nrow(Mono27ac$coverage)
data_vec <- as.integer(Mono27ac$coverage$count)
weight_vec <- as.double(rep(1, data_count))
penalty <- as.double(0)
cost_mat <- numeric(data_count * 2)
end_vec <- integer(data_count)
mean_vec <- numeric(data_count)
intervals_mat <- integer(data_count*2)
label_starts <- 150L
label_ends <- 250L
label_types <- 1L
label_count <- 1L

end_vec
weight_vec
result <- 
  .C("PeakSegFPOPLog_interface", data_vec = data_vec, weight_vec = weight_vec,
     data_count = data_count, penalty = penalty, cost_mat = cost_mat,
     end_vec = end_vec, mean_vec = mean_vec, intervals_mat = intervals_mat,
     label_starts = label_starts, label_ends = label_ends,
     label_types = label_types, label_count = label_count,
     PACKAGE="FLOPART")

costMatrix <- matrix(result[["cost_mat"]], nrow=2, ncol=data_count, byrow=TRUE)

