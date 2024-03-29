library(data.table)
data(neuroblastoma, package="neuroblastoma")
nb.list <- lapply(neuroblastoma, function(DF){
  setkey(data.table(DF), profile.id, chromosome)
})

seq.vec <- 1:nrow(nb.list$annotation)
seq.vec <- 1:10
cv.dt.list <- list()
rect.dt.list <- list()
err.dt.list <- list()
for(seq.i in seq.vec){
  ann <- nb.list$annotations[seq.i]
  problem.vars <- c("profile.id", "chromosome")
  meta <- ann[,problem.vars,with=FALSE]
  pro <- nb.list$profiles[ann]
  cv.fit <- pro[, binsegRcpp::binseg_normal_cv(
    logratio,
    position.vec=position)]
  sbs.fit <- wbs::sbs(pro$logratio)
  sbs.res <- data.table(sbs.fit$res)[order(-min.th)]
  greedy.fit <- binsegRcpp::binseg_normal(pro$logratio)
  ## sigma <- mad(diff(x)/sqrt(2))
  mad.var.est <- mad(diff(pro$logratio)/sqrt(2))
  ## th <- th.const * sigma * sqrt(2 * log(n))
  ## th.const = th/(sigma*sqrt(2*log(n)))
  sbs.res[, min.th.const := min.th/(mad.var.est*sqrt(2*log(nrow(pro))))]
  sbs.res[, ch := c(.N, rep(0, .N-1)), by=min.th]
  sbs.res[, Kmax := cumsum(ch)]
  greedy.selection <- penaltyLearning::modelSelection(
    greedy.fit$splits, "loss", "segments")
  rbind(
    greedy.fit$splits[-1, .(algo="greedy", change=end, param=segments)],
    sbs.res[, .(algo="th", change=cpt, param=min.th)])

  
  th.vec <- unique(sort(sbs.fit$res[,"min.th"]))
  sbs.models <- data.table(meta, th=th.vec)
  cpt.list <- wbs::changepoints(sbs.fit, th=th.vec)
  sbs.changes <- data.table(th.i=seq_along(th.vec))[, {
    th <- th.vec[[th.i]]
    maybe.na.vec <- cpt.list$cpt.th[[th.i]]
    change <- maybe.na.vec[!is.na(maybe.na.vec)]
    change.pos <- cv.fit$subtrain.borders[change+1]
    data.table(meta, th, change, change.pos)
  }, by=th.i]

  err.list <- penaltyLearning::labelError(
    sbs.models,
    ann,
    sbs.changes,
    problem.vars=problem.vars,
    model.vars="th",
    change.var="change.pos")
  err.sort <- err.list$model.errors[order(-th)]
  for(fx in c("fp","fn")){
    d <- diff(err.sort[[fx]])
    set(
      err.sort,
      j=paste0(fx,".diff"),
      value=c(0, d))
  }
  err.diff <- err.sort[fp.diff != 0 | fn.diff != 0]
  err.dt.list[[seq.i]] <- data.table(
    seq.i, err.diff)
  cv.dt.list[[seq.i]] <- data.table(
    seq.i, binseg.fit$cv)
  rect.dt.list[[seq.i]] <- binseg.fit$cv[, data.table(
    seq.i, max.times=times[1], second.times=times[2])]
}
err.dt <- do.call(rbind, err.dt.list)
cv.dt <- do.call(rbind, cv.dt.list)
rect.dt <- do.call(rbind, rect.dt.list)

library(ggplot2)
ggplot()+
  theme_bw()+
  geom_rect(aes(
    xmin=-Inf, xmax=Inf,
    ymin=second.times, ymax=max.times),
    data=rect.dt,
    fill="grey")+
  geom_point(aes(
    segments, times),
    data=cv.dt)+
  facet_grid(seq.i ~ .)
