library(data.table)
library(ggplot2)
data.list <- readRDS(
  "figure-neuroblastoma-one-label-data.rds")
data.list$err[, lapply(.SD, sum), by=pkg, .SDcols=c("fp","fn")]
err.diffs <- data.list$diff
err.diffs[, .SD[c(1,.N)], by=param.name]
thresh.diffs <- err.diffs[, lapply(
  .SD, sum),
  by=.(param.name, param.value),
  .SDcols=paste0("diff.", c("fp","tp"))]
for(suffix in c("fp","tp")){
  s <- function(prefix)paste0(prefix,".",suffix)
  thresh.diffs[, s("cum") := cumsum(.SD[[s("diff")]]), by=param.name]
  thresh.diffs[, s("rate") := {
    cum.vec <- .SD[[s("cum")]]
    cum.vec/max(cum.vec)
  }]
}

roc.dt <- thresh.diffs[, .(
  FPR=c(0, rate.fp),
  TPR=c(0, rate.tp),
  param.value=c(Inf, param.value)
), by=param.name]

auc.dt <- roc.dt[, {
  dt <- .SD[order(param.value)]
  .(auc=WeightedROC::WeightedAUC(dt))
}, by=param.name]
ggplot()+
  theme_bw()+
  geom_path(aes(
    FPR, TPR, color=param.name),
    data=roc.dt)
