library(data.table)
library(ggplot2)
data(Mono27ac, package="PeakSegDisk")

future::plan("multisession")

pen.dt.list <- list()
for(bin.size in c(10, 50, 100, 500, 1000, 5000, 10000)){
  print(bin.size)
  bin.dt <- Mono27ac$coverage[, {
    edge.vec <- seq(chromStart[1], chromEnd[.N], by=bin.size)
    binStart <- edge.vec[-length(edge.vec)]
    data.table(
      binStart,
      binStart1=binStart+1L,
      binEnd=edge.vec[-1])
  }]
  setkey(bin.dt, binStart1, binEnd)
  n.folds <- 3
  bin.dt[, bin.i := 1:.N]
  bin.dt[, fold := bin.i %% n.folds + 1]
  Mono27ac$coverage[, chromStart1 := chromStart+1L]
  setkey(Mono27ac$coverage, chromStart1, chromEnd)
  cov.in.bins <- foverlaps(Mono27ac$coverage, bin.dt, nomatch=0L)
  cov.in.bins[, dataStart := ifelse(chromStart < binStart, binStart, chromStart)]
  cov.in.bins[, dataEnd := ifelse(chromEnd < binEnd, chromEnd, binEnd)]
  stopifnot(cov.in.bins[, dataStart[-1] == dataEnd[-.N]])
  ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    geom_point(aes(
      dataEnd, count),
      data=cov.in.bins)+
    geom_segment(aes(
      dataStart, count,
      xend=dataEnd, yend=count),
      data=cov.in.bins)+
    facet_grid(fold ~ ., labeller="label_both")
  for(validation.fold in 1:n.folds){
    cov.in.bins[, set := ifelse(fold==validation.fold, "validation", "train")]
    valid.bins <- bin.dt[fold==validation.fold]
    valid.bins[, binMid := as.integer((binStart+binEnd)/2)]
    train.dt <- rbind(
      valid.bins[, .(dataStart=binStart, dataEnd=binMid, count=NA)],
      valid.bins[, .(dataStart=binMid, dataEnd=binEnd, count=NA)],
      cov.in.bins[set=="train", .(dataStart, dataEnd, count)]
    )[order(dataStart)]
    valid.dt <- cov.in.bins[set=="validation"]
    end.i <- which(diff(is.na(c(train.dt[["count"]],0))) == -1)
    after.i <- ifelse(end.i==nrow(train.dt), end.i-2, end.i+1)
    after.end <- train.dt[after.i, count]
    train.dt[end.i, count := after.end]
    start.i <- end.i-1
    before.i <- ifelse(start.i==1, 3, start.i-1)
    before.start <- train.dt[before.i, count]
    train.dt[start.i, count := before.start]
    stopifnot(is.finite(train.dt$count))
    vdir <- file.path(
      "figure-search-valid-min-fill-two-constant", bin.size, validation.fold)
    dir.create(vdir, showWarnings=FALSE, recursive=TRUE)
    coverage.bedGraph <- file.path(vdir, "coverage.bedGraph")
    fwrite(
      train.dt[, .(
        chrom=cov.in.bins$chrom[1],
        as.integer(dataStart),
        as.integer(dataEnd),
        as.integer(count)
      )],
      coverage.bedGraph,
      sep="\t", col.names=FALSE, quote=FALSE)
    OnePenalty <- function(penalty){
      one.result <- PeakSegDisk::PeakSegFPOP_dir(vdir, penalty)
      seg.result <- one.result$segments[, .(
        segStart1=chromStart+1,
        segEnd=chromEnd,
        mean)]
      setkey(seg.result, segStart1, segEnd)
      over.dt <- foverlaps(
        valid.dt, seg.result, nomatch=0L)
      vloss <- over.dt[, PeakSegOptimal::PoissonLoss(
        count, mean, dataEnd-dataStart)]
      one.result$loss[, `:=`(
        valid.PoissonLoss = vloss,
        train.loss = total.loss/bases,
        valid.loss = vloss/valid.bases)]
      one.result
    }
    valid.bases <- valid.dt[, sum(dataEnd-dataStart)]
    initial.penalty.vec <- 10^seq(-6, 9, by=0.5)
    result.list <- list()
    next.penalty.vec <- initial.penalty.vec
    while(length(print(next.penalty.vec))){
      LAPPLY <- future.apply::future_lapply
      ##LAPPLY <- lapply
      result.list[paste(next.penalty.vec)] <-
        LAPPLY(next.penalty.vec, OnePenalty)
      loss.dt <- do.call(rbind, lapply(result.list, "[[", "loss"))
      loss.tall <- nc::capture_melt_single(
        loss.dt,
        set="train|valid",
        "[.]loss",
        value.name="mean.PoissonLoss")
      ggplot()+
        geom_label(aes(
          penalty, mean.PoissonLoss, color=set, label=peaks),
          data=loss.tall)+
        scale_x_log10()
      selected <- data.table(penaltyLearning::modelSelection(
        loss.dt, "total.loss", "peaks"))
      candidate.penalty.vec <- selected[
        which.min(valid.loss), c(min.lambda, max.lambda)]
      next.penalty.vec <- candidate.penalty.vec[
        !candidate.penalty.vec %in% names(result.list)]
    }#while(
    pen.dt.list[[paste(bin.size, validation.fold)]] <- data.table(
      bin.size, validation.fold, selected)
  }
}
pen.dt <- do.call(rbind, pen.dt.list)

ggplot()+
  geom_segment(aes(
    min.log.lambda, valid.PoissonLoss,
    xend=max.log.lambda, yend=valid.PoissonLoss),
    data=pen.dt)+
  facet_grid(validation.fold ~ bin.size, labeller=label_both)

mean.dt <- pen.dt[, {
  log.pen.vec <- unique(sort(c(min.log.lambda, max.log.lambda)))
  over.dt <- foverlaps(
    data.table(
      .SD, key=c("min.log.lambda","max.log.lambda")),
    data.table(
      min.log.penalty=log.pen.vec[-length(log.pen.vec)],
      max.log.penalty=log.pen.vec[-1],
      key=c("min.log.penalty","max.log.penalty")),
    nomatch=0L)[
      min.log.lambda < max.log.penalty & min.log.penalty < max.log.lambda]
  over.dt[, .(
    mean.valid.loss=mean(valid.loss), folds=.N
  ), by=.(min.log.penalty, max.log.penalty)]
}, by=bin.size]

gg <- ggplot()+
  geom_segment(aes(
    min.log.penalty, mean.valid.loss,
    xend=max.log.penalty, yend=mean.valid.loss),
    size=1,
    data=mean.dt)+
  facet_grid(. ~ bin.size, labeller=label_both)
png(
  "figure-search-valid-min-fill-two-constant.png",
  width=10, height=3, units="in", res=200)
print(gg)
dev.off()
