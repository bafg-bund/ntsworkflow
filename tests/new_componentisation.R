

pl <- peaklist[[1]]
dl <- datenList[[1]]
pl_componentisation <- function(pl, sl, dl) {
  # make matrix of all EICs of each peak
  eicsl <- Map(function(m, l, r) xcms::rawEIC(dl, mzrange = c(m-0.01, m+0.01), scanrange = c(l-60, r+60)), pl$mz, pl$Leftendscan, pl$Rightendscan)
  eicsl <- lapply(eicsl, as.data.frame)
  eicsdf <- Reduce(function(x, y) merge(x, y, by = "scan", all = TRUE), eicsl)
  eicsdf$scan <- NULL
  stopifnot(ncol(eicsdf) == nrow(pl))
  colnames(eicsdf) <- 1:nrow(pl)
  eicsdf[is.na(eicsdf)] <- 0
  
  znorm <- function(ts) {
    ts.mean <- mean(ts)
    ts.dev <- sd(ts)
    (ts - ts.mean)/ts.dev
  }
  eicsm <- t(apply(eicsdf, 2, znorm))
  
  eicsd <- parallelDist::parDist(eicsm, method = "dtw", threads = 6, window.size = 2)
  #dbscan::kNNdistplot(eicsd, 2)
  dbscan_res <- dbscan::dbscan(eicsd, 35, 2)
  pl$Gruppe2 <- dbscan_res$cluster
}
pl$newGruppe <- newGruppe
View(subset(pl,,c(mz, RT,Gruppe,Gruppe2,newGruppe)))
changeTo0 <- as.numeric(names(table(pl$Gruppe)[table(pl$Gruppe) == 1]))
newGruppe <- pl$Gruppe
newGruppe[newGruppe %in% changeTo0] <- 0
View(data.frame(a = newGruppe, b = dbscan_res$cluster))
table(newGruppe)
ows <- function(ws) {
  eicsd <- parallelDist::parDist(eicsm, method = "dtw", threads = 6, window.size = ws)
  #dbscan::kNNdistplot(eicsd, 3)
  dbscan_res <- dbscan::dbscan(eicsd, 6, 2)
  
  tdf <- data.frame(a = newGruppe, b = dbscan_res$cluster)
  sum(
    tapply(tdf$b, tdf$a, function(x) length(unique(x))),
    tapply(tdf$a, tdf$b, function(x) length(unique(x)))
  )
}
optimise(ows, interval = c(0, 10), tol = 0.01)

oeps <- function(eps) {
  eicsd <- parallelDist::parDist(eicsm, method = "dtw", threads = 6, window.size = 2)
  #dbscan::kNNdistplot(eicsd, 3)
  dbscan_res <- dbscan::dbscan(eicsd, eps, 2)
  
  tdf <- data.frame(a = newGruppe, b = dbscan_res$cluster)
  sum(
    tapply(tdf$b, tdf$a, function(x) length(unique(x))),
    tapply(tdf$a, tdf$b, function(x) length(unique(x)))
  )
}

optimise(oeps, interval = c(1, 20), tol = 0.1)


oboth <- function(x) {
  ws <- x[1]
  ep <- x[2]
  eicsd <- parallelDist::parDist(eicsm, method = "dtw", threads = 6, window.size = ws)
  #dbscan::kNNdistplot(eicsd, 3)
  dbscan_res <- dbscan::dbscan(eicsd, ep, 2)
  
  tdf <- data.frame(a = pl$Gruppe, b = dbscan_res$cluster)
  sum(
    tapply(tdf$b, tdf$a, function(x) length(unique(x))),
    tapply(tdf$a, tdf$b, function(x) length(unique(x)))
  )
}
optim(c(0, 1), oboth, method = "SANN")




dbscan_res
pl$Gruppe2 <- dbscan_res$cluster
View(pl)
View(data.frame(eicsm[66,],eicsdf[,66]))
