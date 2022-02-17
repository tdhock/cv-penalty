library(data.table)
library(ggplot2)
data(Mono27ac, package="PeakSegDisk")

pen.dt.list <- list()
for(bin.size in c(100, 500, 1000)){
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
  bin.dt[, bin.i := 1:.N]
  Mono27ac$coverage[, chromStart1 := chromStart+1L]
  setkey(Mono27ac$coverage, chromStart1, chromEnd)
  cov.in.bins <- foverlaps(Mono27ac$coverage, bin.dt, nomatch=0L)
  n.folds <- 3
  cov.in.bins[, fold := bin.i %% n.folds + 1]
  cov.in.bins[, dataStart := ifelse(chromStart < binStart, binStart, chromStart)]
  cov.in.bins[, dataEnd := ifelse(chromEnd < binEnd, chromEnd, binEnd)]
  stopifnot(cov.in.bins[, dataStart[-1] == dataEnd[-.N]])
  gg <- ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    geom_point(aes(
      dataEnd, count),
      data=cov.in.bins)+
    facet_grid(fold ~ .)
  for(validation.fold in 1:n.folds){
    cov.in.bins[, set := ifelse(fold==validation.fold, "validation", "train")]
    train.dt <- cov.in.bins[set=="train"]
    train.dt[, nextStart := c(dataStart[-1], NA)]
    train.dt[, nextCount := c(count[-1], NA)]
    ## impute via mean of count before and after gap.
    train.imp <- rbind(train.dt[dataEnd != nextStart, {
      data.table(
        chrom,
        dataStart=dataEnd,
        dataEnd=nextStart,
        count=as.integer((nextCount+count)/2))
    }], train.dt[, .(
          chrom,
          dataStart,
          dataEnd,
          count)])[order(dataStart)]
    stopifnot(train.imp[, dataStart[-1] == dataEnd[-.N]])
    show.dt <- rbind(
      train.imp[, .(
        dataStart, dataEnd, count,
        set="train", type="imputed")],
      cov.in.bins[, .(dataStart, dataEnd, count, set, type="raw")])
    ggplot()+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      geom_point(aes(
        dataEnd, count, color=type),
        data=show.dt)+
      geom_segment(aes(
        dataStart+0.5, count,
        xend=dataEnd+0.5, yend=count, color=type),
        data=show.dt)+
      facet_grid(set ~ .)+
      coord_cartesian(xlim=c(205000, 210000))
  }
  ## just repeat train data in valid regions
  n.folds <- 2
  cov.in.bins[, fold := ((bin.i-1) %% n.folds) + 1]
  cov.in.bins[, dataStart := ifelse(chromStart < binStart, binStart, chromStart)]
  cov.in.bins[, dataEnd := ifelse(chromEnd < binEnd, chromEnd, binEnd)]
  stopifnot(cov.in.bins[, dataStart[-1] == dataEnd[-.N]])
  gg <- ggplot()+
    theme_bw()+
    theme(panel.spacing=grid::unit(0, "lines"))+
    geom_point(aes(
      dataEnd, count),
      data=cov.in.bins)+
    facet_grid(fold ~ .)
  print(gg)
  gg+coord_cartesian(xlim=c(2, 3)*1e5)
  gg+coord_cartesian(xlim=c(200000, 225000))
  gg+coord_cartesian(xlim=c(205000, 210000))
  ## 1 2 3 4 1 2 3 4
  ## _ 2 3 _ _ 2 3 _
  ## 2 2 3 3 2 2 3 3
  ## 1 2 3 4 1 2 3 4
  ## 1 _ _ 4 1 _ _ 4
  ## 1 1 4 4 1 1 4 4
  ## 1 2 3 4 1 2 3 4
  ## 1 2 _ _ 1 2 _ _
  ## 1 2 2 1 1 2 2 1
  ## 1 2 3 4 1 2 3 4
  ## _ _ 3 4 _ _ 3 4
  ## 3 3 3 4 4 3 3 4
  ## _ 2 _ 2
  ## 2 2 2 2
  ## 1 _ 1 _
  ## 1 1 1 1
  for(validation.fold in 1:n.folds){
    cov.in.bins[, set := ifelse(fold==validation.fold, "validation", "train")]
    train.dt <- cov.in.bins[set=="train"]
    s <- if(validation.fold==2)1 else -1
    train.imp <- rbind(train.dt[, .(
      chrom,
      dataStart=dataStart+bin.size*s,
      dataEnd=dataEnd+bin.size*s,
      count
    )],
    train.dt[, .(
      chrom,
      dataStart,
      dataEnd,
      count
    )])[order(dataStart)]
    train.imp[, table(dataStart[-1] != dataEnd[-.N])]
    stopifnot(train.imp[, dataStart[-1] == dataEnd[-.N]])
    show.dt <- rbind(
      train.imp[, .(
        dataStart, dataEnd, count,
        set="train", type="imputed")],
      cov.in.bins[, .(dataStart, dataEnd, count, set, type="raw")])
    ggplot()+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      geom_point(aes(
        dataEnd, count, color=type),
        data=show.dt)+
      geom_segment(aes(
        dataStart+0.5, count,
        xend=dataEnd+0.5, yend=count, color=type),
        data=show.dt)+
      facet_grid(set ~ .)+
      coord_cartesian(xlim=c(205000, 210000))
    vdir <- file.path("figure-proof-of-concept", bin.size, validation.fold)
    dir.create(vdir, showWarnings=FALSE, recursive=TRUE)
    coverage.bedGraph <- file.path(vdir, "coverage.bedGraph")
    train.imp[, line := with(rle(count), rep(seq_along(lengths), lengths))]
    out.dt <- train.imp[, .(
      chromStart=as.integer(min(dataStart)),
      chromEnd=as.integer(max(dataEnd))
    ), by=.(chrom, count, line)]
    fwrite(
      out.dt[, .(chrom, chromStart, chromEnd, count)],
      coverage.bedGraph,
      sep="\t", col.names=FALSE, quote=FALSE)
    for(n.peaks in c(5, 10, 20, 40)){
      PeakSegDisk::sequentialSearch_dir(vdir, as.integer(n.peaks), verbose=TRUE)
    }
    segments.bed.vec <- Sys.glob(file.path(vdir, "*_segments.bed"))
    for(segments.i in seq_along(segments.bed.vec)){
      segments.bed <- segments.bed.vec[[segments.i]]
      seg.dt <- fread(
        segments.bed,
        col.names=c(
          "chrom",
          "segStart",
          "segEnd",
          "status",
          "mean"))
      ggplot()+
        theme_bw()+
        theme(panel.spacing=grid::unit(0, "lines"))+
        geom_point(aes(
          dataEnd, count),
          data=cov.in.bins)+
        geom_segment(aes(
          dataStart+0.5, count,
          xend=dataEnd+0.5, yend=count),
          data=cov.in.bins)+
        facet_grid(set ~ .)+
        geom_segment(aes(
          segStart, mean,
          xend=segEnd, yend=mean),
          color="green",
          data=seg.dt)+
        coord_cartesian(xlim=c(205000, 210000))
      seg.dt[, segStart1 := segStart+1L]
      setkey(seg.dt, segStart1, segEnd)
      cov.in.bins[, dataStart1 := dataStart +1L]
      setkey(cov.in.bins, dataStart1, dataEnd)
      segs.over.cov <- foverlaps(seg.dt, cov.in.bins, nomatch=0L)
      segs.over.cov[
      , overStart := ifelse(dataStart < segStart, segStart, dataStart)]
      segs.over.cov[
      , overEnd := ifelse(dataEnd < segEnd, dataEnd, segEnd)]
      loss.tsv <- sub("segments.bed", "loss.tsv", segments.bed)
      loss.dt <- fread(loss.tsv, col.names=PeakSegDisk::col.name.list$loss)
      segs.over.cov[, l1err := abs(mean - count)]
      segs.over.cov[, l2err := l1err^2]
      pen.dt.list[[paste(bin.size, validation.fold, segments.i)]] <- segs.over.cov[, {
        bases <- overEnd-overStart
        data.table(
          bin.size,
          validation.fold,
          penalty=loss.dt$penalty,
          peaks=loss.dt$peaks,
          bases=sum(bases),
          mae=sum(l1err*bases)/sum(bases),
          mse=sum(l2err*bases)/sum(bases),
          PoissonLoss=PeakSegOptimal::PoissonLoss(count, mean, bases)
        )
      }, by=.(set)]
    }
  }
}
(pen.dt <- do.call(rbind, pen.dt.list)[order(validation.fold, set, peaks)])


tall.dt <- melt(pen.dt, measure.vars=c("PoissonLoss", "mse", "mae"))
u.pen.vec <- unique(sort(tall.dt$penalty))
tall.approx <- tall.dt[, {
  L <- approx(log(penalty), value, log(u.pen.vec))
  with(L, data.table(penalty=exp(x), value=y))
}, by=.(bin.size, set, variable, validation.fold)][!is.nan(value)]
tall.stats <- tall.approx[, .(
  mean=mean(value),
  sd=sd(value)
), by=.(bin.size, set, variable, penalty)]
best.pen <- tall.stats[
  set=="validation",
  .SD[which.min(mean)],
  by=.(bin.size, variable)]

all.dir <- file.path("figure-proof-of-concept", "all")
dir.create(all.dir)
all.bg <- file.path(all.dir, "coverage.bedGraph")
fwrite(
  Mono27ac$coverage[, .(chrom, chromStart, chromEnd, count)],
  all.bg,
  sep="\t", col.names=FALSE, quote=FALSE)
err.dt.list <- list()
for(penalty in c(100, 500, 1000, 5000, 10000, 20000, 30000)){
  fit <- PeakSegDisk::PeakSegFPOP_dir(all.dir, penalty)
  peak.dt <- fit$segments[
    status=="peak",
    .(chrom,
      chromStart,
      chromEnd)]
  err.df <- PeakError::PeakErrorChrom(peak.dt, Mono27ac$labels)
  err.dt.list[[paste(penalty)]] <- with(
    err.df,
    data.table(
      penalty,
      set="test",
      variable="label errors",
      errors=sum(fp+fn)))
}
(err.dt <- do.call(rbind, err.dt.list))

set.colors <- c(
  test="black",
  validation="blue",
  train="red")
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(variable ~ bin.size, scales="free")+
  scale_color_manual(values=set.colors, breaks=names(set.colors))+
  scale_fill_manual(values=set.colors)+
  guides(fill="none")+
  geom_ribbon(aes(
    penalty, ymin=mean-sd, ymax=mean+sd, fill=set),
    data=tall.stats,
    alpha=0.5)+
  geom_line(aes(
    penalty, mean,
    color=set),
    size=1,
    data=tall.stats)+
  geom_point(aes(
    penalty, value,
    color=set),
    data=tall.dt)+
  geom_point(aes(
    penalty, errors, color=set),
    data=err.dt)+
  geom_vline(aes(
    xintercept=penalty, color=set),
    data=best.pen)+
  scale_x_log10()+
  ylab("")
png("figure-proof-of-concept.png", 10, 6, units="in", res=100)
print(gg)
dev.off()
