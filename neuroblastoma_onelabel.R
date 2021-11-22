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
diff.dt.list <- list()
for(seq.i in seq.vec){
  ann <- nb.list$annotations[seq.i]
  problem.vars <- c("profile.id", "chromosome")
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
  greedy.selection <- data.table(penaltyLearning::modelSelection(
    greedy.fit$splits[, .(loss, segments)], "loss", "segments")
    )[, .(max.log.lambda, segments)]
  greedy.cpts <- greedy.selection[
    greedy.fit$splits[-1, .(cpt=end, segments)],
    on="segments"]
  penalty.vec <- c(AIC=2, BIC=log(nrow(pro)))
  for(penalty.name in names(penalty.vec)){
    penalty.value <- penalty.vec[[penalty.name]]
    set(
      greedy.cpts,
      j=penalty.name,
      value=greedy.cpts$max.log.lambda-log(penalty.value))
  }
  algo.list <- list(
    "greedy.cpts"=c("AIC","BIC"),
    "sbs.res"=c("min.th","min.th.const"))
  for(algo.name in names(algo.list)){
    cpt.dt <- get(algo.name)    
    cpt.dt[, cpt.pos := cv.fit$subtrain.borders[cpt+1] ]
    in.ann <- cpt.dt[ann$min < cpt.pos & cpt.pos < ann$max]
    param.name.vec <- algo.list[[algo.name]]
    for(param.name in param.name.vec){
      first.i <- which(is.finite(in.ann[[param.name]]))[1]
      first <- in.ann[first.i]
      diff.dt.list[[paste(seq.i, algo.name, param.name)]] <- data.table(
        ann[,.(profile.id,chromosome)],
        param.name, 
        param.value=first[[param.name]],
        diff.fp=ifelse(ann$ann=="normal", 1, 0),
        diff.tp=ifelse(ann$ann=="breakpoint", 1, 0))
    }
  }
  cv.dt.list[[seq.i]] <- data.table(
    seq.i, binseg.fit$cv)
  rect.dt.list[[seq.i]] <- binseg.fit$cv[, data.table(
    seq.i, max.times=times[1], second.times=times[2])]
}
diff.dt <- do.call(rbind, diff.dt.list)
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
