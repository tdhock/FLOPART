counts <- data.table::fread("~/Downloads/counts.csv",drop=1)
labels <- data.table::fread("~/Downloads/labels.csv",drop=1)
fit <- FLOPART::FLOPART(counts,labels,penalty=5)
