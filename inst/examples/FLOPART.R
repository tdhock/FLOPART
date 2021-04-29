data("Mono27ac", package="FLOPART")
fit <- with(Mono27ac, FLOPART::FLOPART(coverage, labels, 1e6))
str(fit)
