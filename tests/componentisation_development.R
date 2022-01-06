

alig_componentisation(altable = grouped)

saveRDS(list(combiDist = combiDist, al = al), "~/projects/ntsworkflow/tests/test1.RDS")

View(altable[, c(1,2,4,grep("gruppe_", colnames(altable)))])

# calculate intergroup RT variation


which.max(unlist(by(al, al$Gruppe, function(x) sd(x$mean_RT), simplify = F))[-1])
unlist(by(al, al$Gruppe, function(x) sd(x$mean_RT), simplify = F))[34]
unlist(by(al, al$Gruppe, function(x) sd(x$mean_RT), simplify = F))


dbscan_res2 <- dbscan::dbscan(rtDistT, 0.01, 2)
altable[, "Gruppe"] <- dbscan_res2$cluster
al <- as.data.frame(altable)
View(al[c(90,390),])
View(as.matrix(rtDistT)[c(90,390), c(90,390)])

combiDist <- readRDS("tests/test1.RDS")$combiDist
al <- readRDS("tests/test1.RDS")$al
m <- as.matrix(combiDist)
m[lower.tri(m, diag = T)] <- 1
m2 <- m <= 0
gr <- apply(m2, 1, which)
counter <- 1
clusters <- numeric(length(gr))
allAssigments <- unname(unlist(gr))
# allAssigments[duplicated(allAssigments)][1]
gr[sapply(gr, is.element, el = 70)]
gr[["28"]]
# gr[["3"]]

#639 %in% allAssigments
for (grId in names(gr)) { # grId <- "122"
  
  grNum <- as.numeric(grId)
  if (grNum %in% allAssigments)
    next
  others <- gr[[grId]]
  if (length(others) > 0) {
    rows <- unname(c(grNum, gr[[grId]]))
    if (!all(clusters[rows] == 0))
      
    clusters[rows] <- counter
    counter <- counter + 1
  }
}
al$Gruppe <- clusters
View(al[al$Gruppe ==117, ])
View(al[al$Gruppe ==117, grep("gruppe_", colnames(al))])
View(al[c(28,  579), ])
View(m[c(28,  579),c(28,  579)])

View(al[c(70,  245), ])
View(al[c(70, 245), grep("gruppe_", colnames(al))])



dbscan_res <- dbscan::dbscan(combiDist, 0.1, 2)
al$Gruppe <- dbscan_res$cluster
