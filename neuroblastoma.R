library(data.table)
data(neuroblastoma, package="neuroblastoma")
nb.list <- lapply(neuroblastoma, function(DF){
  setkey(data.table(DF), profile.id, chromosome)
})

seq.vec <- 1:nrow(nb.list$annotation)
seq.vec <- 1:10
cv.dt.list <- list()
rect.dt.list <- list()
for(seq.i in seq.vec){
  ann <- nb.list$annotations[seq.i]
  pro <- nb.list$profiles[ann]
  binseg.fit <- pro[, binsegRcpp::binseg_normal_cv(
    logratio,
    position.vec=position)]
  sbs.fit <- wbs::sbs(pro$logratio)
  th.vec <- unique(sort(wbs.fit$res[,"min.th"]))
  sbs.models <- data.table(th=th.vec)
  cpt.list <- wbs::changepoints(sbs.fit, th=th.vec)
  sbs.changes <- data.table(th.i=seq_along(th.vec))[, {
    th <- th.vec[[th.i]]
    maybe.na.vec <- cpt.list$cpt.th[[th.i]]
    change <- maybe.na.vec[!is.na(maybe.na.vec)]
    change.pos <- binseg.fit$subtrain.borders[change+1]
    data.table(th, change, change.pos)
  }, by=th.i]
  cv.dt.list[[seq.i]] <- data.table(
    seq.i, binseg.fit$cv)
  rect.dt.list[[seq.i]] <- binseg.fit$cv[, data.table(
    seq.i, max.times=times[1], second.times=times[2])]
}
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
