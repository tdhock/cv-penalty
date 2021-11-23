library(data.table)
data(neuroblastoma, package="neuroblastoma")
nb.list <- lapply(neuroblastoma, function(DF){
  setkey(data.table(DF), profile.id, chromosome)
})

seq.vec <- 1:10
seq.vec <- 1:nrow(nb.list$annotation)
err.dt.list <- list()
diff.dt.list <- list()
for(seq.i in seq.vec){
  cat(sprintf(
    "%4d / %4d sequences\n",
    seq.i,
    length(seq.vec)))
  ann <- nb.list$annotations[seq.i]
  meta <- ann[,.(profile.id,chromosome)]
  pro <- nb.list$profiles[ann]
  set.seed(1)
  cv.fit <- pro[, binsegRcpp::binseg_normal_cv(
    logratio,
    position.vec=position)]
  sbs.fit <- wbs::sbs(pro$logratio)
  cpt.fit <- changepoint::cpt.mean(
    pro$logratio, method="BinSeg",
    Q=nrow(pro)/2)
  not.end <- function(e)e[e<nrow(pro)]
  change.list <- list(
    binsegRcpp=coef(cv.fit)[, end],
    changepoint=cpt.fit@cpts,
    wbs=wbs::changepoints(sbs.fit)$cpt.th[[1]])
  model.dt <- data.table(meta, pkg=names(change.list))
  err.dt.list[[seq.i]] <- model.dt[, {
    change.i <- not.end(change.list[[pkg]])
    change.dt <- data.table(
      meta, pkg, change.pos=cv.fit$subtrain.borders[change.i+1])
    model.dt <- data.table(meta, pkg)
    err.list <- penaltyLearning::labelError(
      model.dt, ann, change.dt,
      problem.vars=names(meta),
      change.var="change.pos",
      model.vars="pkg")
    err.list$label.errors
  }, by=pkg]
  sbs.res <- data.table(sbs.fit$res)[order(-min.th)]
  greedy.fit <- pro[, binsegRcpp::binseg_normal(
    logratio,
    position.vec = position)]
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
    cpt.dt[, cpt.pos := greedy.fit$subtrain.borders[cpt+1] ]
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
}
out.list <- list(
  diff=do.call(rbind, diff.dt.list),
  err=do.call(rbind, err.dt.list))
saveRDS(out.list, "figure-neuroblastoma-one-label-data.rds")

