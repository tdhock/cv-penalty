library(data.table)
library(ggplot2)
data(Mono27ac, package="PeakSegDisk")

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
  for(shuffle.fold in 1:n.folds){
    bin.dt[, set := ifelse(fold==shuffle.fold, "to.shuffle", "to.keep")]
    cov.in.bins <- foverlaps(Mono27ac$coverage, bin.dt, nomatch=0L)
    cov.in.bins[, dataStart := ifelse(
      chromStart < binStart, binStart, chromStart)]
    cov.in.bins[, dataEnd := ifelse(
      chromEnd < binEnd, chromEnd, binEnd)]
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
    for(seed in 1:2){
      bin.perm <- data.table(bin.dt, key="bin.i")
      set.seed(seed)
      bin.perm[set=="to.shuffle", bin.i := sample(bin.i)]
      bin.join <- bin.perm[, .(
        newBinStart=binStart,
        bin.i,
        set=ifelse(set=="to.shuffle", "shuffled", "kept"))]
      ## train data has shuffled and kept.
      cov.permutation <- data.table(cov.in.bins)[bin.join, on="bin.i"]
      cov.permutation[, newDataStart := dataStart-binStart+newBinStart]
      cov.permutation[, newDataEnd := dataEnd-binStart+newBinStart]
      ggplot()+
        theme_bw()+
        theme(panel.spacing=grid::unit(0, "lines"))+
        geom_step(aes(
          dataEnd, count, color=set),
          data=data.table(cov.in.bins)[, set := "full"])+
        geom_point(aes(
          newDataEnd, count,
          color=set),
          data=cov.permutation)+
        geom_segment(aes(
          newDataStart, count,
          color=set,
          xend=newDataEnd, yend=count),
          data=cov.permutation)
      vdir <- file.path(
        "figure-search-valid-min-permutation", bin.size, shuffle.fold, seed)
      dir.create(vdir, showWarnings=FALSE, recursive=TRUE)
      coverage.bedGraph <- file.path(vdir, "coverage.bedGraph")
      fwrite(
        cov.permutation[, .(
          chrom,
          as.integer(newDataStart),
          as.integer(newDataEnd),
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
          cov.in.bins, seg.result, nomatch=0L)
        ## validation data has to.shuffle and to.keep.
        set.loss <- over.dt[, .(
          total.loss=PeakSegOptimal::PoissonLoss(
            count, mean, dataEnd-dataStart),
          bases=sum(dataEnd-dataStart)
        ), by=set][, mean.loss := total.loss/bases]
        ## TODO maybe also store mean valid loss?
        for(set.i in 1:nrow(set.loss)){
          set.loss[set.i, set(
            one.result$loss,
            j=paste0(set, ".loss"),
            value=mean.loss)]
        }
        one.result
      }
      initial.penalty.vec <- 10^seq(-6, 9, by=0.5)
      result.list <- list()
      next.penalty.vec <- initial.penalty.vec
      while(length(print(next.penalty.vec))){
        LAPPLY <- future.apply::future_lapply
        LAPPLY <- lapply
        result.list[paste(next.penalty.vec)] <-
          LAPPLY(next.penalty.vec, OnePenalty)
        loss.dt <- do.call(rbind, lapply(result.list, "[[", "loss"))
        loss.tall <- nc::capture_melt_single(
          loss.dt,
          set="to.shuffle|to.keep",
          "[.]loss",
          value.name="mean.PoissonLoss")
        ggplot()+
          geom_line(aes(
            penalty, mean.PoissonLoss),
            data=loss.tall)+
          geom_label(aes(
            penalty, mean.PoissonLoss, label=peaks),
            data=loss.tall)+
          scale_x_log10()+
          facet_grid(set~.)
        selected <- data.table(penaltyLearning::modelSelection(
          loss.dt, "total.loss", "peaks"))
        candidate.penalty.vec <- selected[
          which.min(to.shuffle.loss), c(min.lambda, max.lambda)]
        next.penalty.vec <- candidate.penalty.vec[
          !candidate.penalty.vec %in% names(result.list)]
      }#while(
      pen.dt.list[[paste(bin.size, shuffle.fold, seed)]] <- data.table(
        bin.size, shuffle.fold, seed, selected)
    }
  }
}
pen.dt <- do.call(rbind, pen.dt.list)

##TODO.
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
  "figure-search-valid-min-fill-noise.png",
  width=10, height=3, units="in", res=200)
print(gg)
dev.off()
