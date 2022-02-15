label_colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c")

##' FLOPART needs at most one label per coverage data row, which may
##' not be the case for arbitrary coverage/labels.
##' @title Convert data for input to FLOPART
##' @param coverage data frame of coverage with columns chromStart, chromEnd, count
##' @param label data.frame of labels with with columns chromStart, chromEnd, annotation
##' @return named list: coverage and label
##' @author Toby Dylan Hocking
##' @example inst/examples/FLOPART_data.R
FLOPART_data <- function(coverage, label){
  i.chromStart <- i.chromEnd <- count <- chromStart <- chromEnd <-
    annotation <- labelStart <- labelEnd <- type <- firstRow <-
      lastRow <- . <- NULL
  if(missing(label)){
    label <- data.table(
      chromStart=integer(), chromEnd=integer(), annotation=character())
  }
  if(!is.data.frame(coverage)){
    stop("coverage must be a data frame")
  }
  if(any(max(coverage[["chromEnd"]]) < label[["chromEnd"]])){
    stop("label ends must be on or before last coverage")
  }
  if(any(label[["chromStart"]] < min(coverage[["chromStart"]]))){
    stop("label starts must be on or after first coverage")
  }
  label_code <- get_label_code()
  df.names <- list("label", "coverage")
  key.vec <- c("chromStart", "chromEnd")
  uniq.pos.list <- list()
  dt.list <- list()
  for(data.type in df.names){
    data.dt <- data.table(get(data.type), key=key.vec)
    for(k in key.vec){
      uniq.pos.list[[paste(data.type, k)]] <- data.dt[[k]]
    }
    dt.list[[data.type]] <- data.dt
  }
  uniq.pos <- unique(sort(unlist(uniq.pos.list)))
  uniq.dt <- data.table(
    chromStart=uniq.pos[-length(uniq.pos)],
    chromEnd=uniq.pos[-1])
  with.counts <- dt.list[["coverage"]][
    uniq.dt,
    .(chromStart=i.chromStart, chromEnd=i.chromEnd, count,
      weight = i.chromEnd - i.chromStart),
    on=.(chromStart < chromEnd, chromEnd > chromStart)]
  label.index.dt <- with.counts[
    dt.list[["label"]],
    .(firstRow=.I[1],
      lastRow=.I[.N],
      labelStart=i.chromStart,
      labelEnd=i.chromEnd,
      type=label_code[paste(annotation)],
      annotation),
    by=.EACHI,
    on=.(chromEnd > chromStart, chromStart < chromEnd)
  ][, .(
    chromStart=labelStart,
    chromEnd=labelEnd,
    annotation, type,
    firstRow, lastRow
  )]
  list(coverage_dt=with.counts, label_dt=label.index.dt)
}

##' Main function for computing optimal segmentation model with
##' Poisson loss, up-down constraints, and label constraints.
##' @title Functional Labeled Optimal Partitioning
##' @param coverage data frame of coverage
##' @param label data frame of labels
##' @param penalty non-negative penalty constant
##' @return augmented list of output from FLOPART_interface and
##'   FLOPART_data
##' @author Toby Dylan Hocking
##' @example inst/examples/FLOPART.R
FLOPART <- function(coverage, label, penalty){
  status <- state <- NULL
  data.list <- FLOPART_data(coverage, label)
  ##print(data.list)
  result <- with(data.list, FLOPART_interface(
    coverage_dt[["count"]],
    coverage_dt[["weight"]],
    penalty,
    label_dt[["type"]],
    label_dt[["firstRow"]]-1,
    label_dt[["lastRow"]]-1))
  segs <- as.data.table(result[["segments_df"]])
  name.vec <- c(firstRow="chromStart", lastRow="chromEnd")
  for(xRow in names(name.vec)){
    chromX <- name.vec[[xRow]]
    row.vec <- segs[[xRow]]
    pos.vec <- data.list[["coverage_dt"]][[chromX]]
    set(segs, j=chromX, value=pos.vec[row.vec])
  }
  segs[, status := ifelse(state==0, "background", "peak")]
  result[["segments_df"]] <- NULL
  result[["segments_dt"]] <- segs
  c(data.list, result)
}
