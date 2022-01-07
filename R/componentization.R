



#' @export
annotate_grouped_compononents <- function(sampleListLocal,
                                          peakListList,
                                          alignmentTable, numcores) {
  # browser()
  # prepare data
  alig <- as.data.frame(alignmentTable)
  stopifnot(length(peakListList) == nrow(sampleListLocal))
  pl <- do.call("rbind", peaklist)
  peakListList2 <- Map(function(x, s) {
    x <- x[,c(1,2,19:32,37)]
    x$row <- seq_len(nrow(x))
    x$sample <- s
    x
  }, peakListList, sampleListLocal$ID)
  pl <- do.call("rbind", peakListList2)
  
  ap <- alig[, grep("^PeakID_", colnames(alig))]
  apl <- split(ap, seq_len(nrow(ap)))
  apl <- lapply(apl, as.numeric)
  # loop through ap, in each row, collect all available annotations from pl
  rowAnnot <- parallel::mclapply(apl, function(x) { # x <- as.numeric(ap[1, ])
    subpl <- pl[pl$peak_id_all %in% x, ]
    # see if any of these peakIDs have Cl or Br links to them
    # the monoisotopic peak does not have an annotation, only the
    # 37Cl peak so you have to look for the row number of the 
    # current peakID in the whole table (Peaks are referenced with row
    # numbers not with fixed peakIDs in the peaklist (historic))
    anyCl <- Map(function(r, s) {  # r <- 1, s <- 1
      allCls <- peakListList[[s]][, c("Cl1","Cl2","Cl3","Cl4")]
      any(unlist(allCls) == r)
    }, subpl$row, subpl$sample)
    anyCl <- any(unlist(anyCl))
    anyBr <- Map(function(r, s) {  # r <- 1, s <- 1
      allCls <- peakListList[[s]][, c("Br1","Br2","Br3","Br4")]
      any(unlist(allCls) == r)
    }, subpl$row, subpl$sample)
    anyBr <- any(unlist(anyBr))
    if (anyCl && anyBr)
      return("Possible Cl and Br")
    if (anyCl)
      return("Possible Cl")
    if (anyBr)
      return("Possible Br")
    character(1)
  }, mc.cores = numcores)

  alig$clbr <- unlist(rowAnnot)
  res <- alig[alig$clbr != "", c("alignmentID", "mean_mz", "mean_RT", "clbr")]
  colnames(res) <- c("alignmentID", "mzData", "rtData", "name")
  res$rtData <- round(res$rtData / 60, 2)
  res$CAS <- NA
  res$mzDB <- NA
  res$rtDB <- NA
  res$expID <- NA
  res$formula <- NA
  res$SMILES <- NA
  res$adduct <- NA
  res$isotope <- NA
  res$score <- NA
  res$db_available <- NA
  res$CE <- NA
  res$CES <- NA
  # get sample with the highest intensity for each alignmentID
  res$sample <- vapply(res$alignmentID, function(aid) {
    intCols <- alig[alig$alignmentID == aid, grep("^Int_", colnames(alig)), drop = FALSE]
    te <- t(intCols)
    sampHigh <- row.names(te)[which.max(te[1,])]
    as.numeric(stringr::str_match(sampHigh, "_(\\d+)$")[,2])
  }, FUN.VALUE = numeric(1))
  res
}


#' @export
componentization_BfG <- function(Liste,
                                 daten,
                                 ppm = 10,
                                 Grenzwert_RT = 10,     
                                 Grenzwert_FWHM_left = 10,
                                 Grenzwert_FWHM_right = 10,
                                 Summe_all = 10,
                                 adjust_tolerance = TRUE) {
  
  Liste$C13 <- 0
  Liste$NaAddukt <- 0
  Liste$KAddukt <- 0
  Liste$NH4Addukt <- 0
  Liste$Cl1 <- 0
  Liste$Cl2 <- 0
  Liste$Cl3 <- 0
  Liste$Cl4 <- 0
  Liste$Br1 <- 0
  Liste$Br2 <- 0
  Liste$Br3 <- 0
  Liste$Br4 <- 0
  Liste$S1 <- 0
  Liste$S2 <- 0
  
  
  Liste$Gruppe <- 0
  
  #Grenzwert_RT <- 10
  #Grenzwert_FWHM_left <- 10
  #Grenzwert_FWHM_right <- 10
  #Summe_all <- 10
  
  Liste <- Liste[order(Liste$Intensity,decreasing = T), ] 
  i <- 1
  g <- 0
  groupleader <- NULL
  
  Grenzwert_RT[1:nrow(Liste)] <- Grenzwert_RT
  Grenzwert_FWHM_right[1:nrow(Liste)] <- Grenzwert_FWHM_right
  Grenzwert_FWHM_left[1:nrow(Liste)] <- Grenzwert_FWHM_left
  Summe_all[1:nrow(Liste)] <- Summe_all
  
  zaehler <- NULL
  
  int_min <- min(Liste$Intensity) #das beim PeakPicking gewählte Intensitäts-Minimum
 
  mzDiff13C <- 1.00335
  mzNa <- 21.98194
  
  while(any(Liste$Gruppe == 0)) {
    g <- g+1  #nächste Gruppe
    
    groupleader <- c(groupleader,which(Liste$Gruppe == 0)[1])
    
    
    Liste$Gruppe[groupleader[g]] <- g
    n <- 1
    abbruch <- FALSE
    while (!abbruch) {
      
      C13 <- which((abs(Liste$mz-Liste$mz[groupleader[g]]-mzDiff13C*n) < Liste$mz[groupleader[g]]/1000000*ppm) & (Liste$Intensity/Liste$Intensity[groupleader[g]]*100 <= Liste$mz[groupleader[g]]/12^n) &
                     abs(Liste$RT[groupleader[g]]-Liste$RT)<Grenzwert_RT[groupleader[g]] &                          # RT im Limit?
                     abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left)<Grenzwert_FWHM_left[groupleader[g]] &     # FWHM_left im Limit?
                     abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right)<Grenzwert_FWHM_right[groupleader[g]] &  # FWHM_right im Limit?
                     (abs(Liste$RT[groupleader[g]]-Liste$RT) +
                        abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left) +
                        abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right))<Summe_all[groupleader[g]])
      if (length(C13) > 1) { #es gibt mehrere Kandidaten
        C13 <- C13[which(Liste$Intensity[C13]/Liste$Intensity[groupleader[g]] > 0.01)] #ein Kohlenstoff sollte schon mindestens drin sein :)
        C13 <- C13[which.min(abs(Liste$RT[C13]-Liste$RT[groupleader[g]]))] #nimm den, der gemessen an der RT am besten passt
      }
      if (length(C13) < 1) abbruch <- TRUE
      Liste$C13[C13] <- groupleader[g]
      Liste$Gruppe[C13] <- g
      n <- n+1
    }
    suppressWarnings({
    if (adjust_tolerance) {
      Grenzwert_RT[groupleader[g]] <- max(abs(Liste$RT[(Liste$Gruppe == g) & (Liste$C13 == groupleader[g])]-Liste$RT[groupleader[g]]))*2
      
      if (Grenzwert_RT[groupleader[g]] == 0) 
        Grenzwert_RT[groupleader[g]] <- abs(daten@scantime[Liste$Scan[groupleader[g]]]-daten@scantime[Liste$Scan[groupleader[g]]+1])*2
      
      if (is.infinite(Grenzwert_RT[groupleader[g]])) 
        Grenzwert_RT[groupleader[g]] <- abs(daten@scantime[Liste$Scan[groupleader[g]]]-daten@scantime[Liste$Scan[groupleader[g]]+1])*2
      
      Grenzwert_FWHM_left[groupleader[g]] <- max(abs(Liste$FWHM_left[(Liste$Gruppe == g) & (Liste$C13 == groupleader[g])]-Liste$FWHM_left[groupleader[g]]))*2
      
      if (Grenzwert_FWHM_left[groupleader[g]] == 0) 
        Grenzwert_FWHM_left[groupleader[g]] <- abs(daten@scantime[Liste$Scan[groupleader[g]]]-daten@scantime[Liste$Scan[groupleader[g]]+1])*2
      
      if (is.infinite(Grenzwert_FWHM_left[groupleader[g]])) 
        Grenzwert_FWHM_left[groupleader[g]] <- abs(daten@scantime[Liste$Scan[groupleader[g]]]-daten@scantime[Liste$Scan[groupleader[g]]+1])*2
      
      Grenzwert_FWHM_right[groupleader[g]] <- max(abs(Liste$FWHM_right[(Liste$Gruppe == g) & (Liste$C13 == groupleader[g])]-Liste$FWHM_right[groupleader[g]]))*2
      
      if (Grenzwert_FWHM_right[groupleader[g]] == 0) 
        Grenzwert_FWHM_right[groupleader[g]] <- abs(daten@scantime[Liste$Scan[groupleader[g]]]-daten@scantime[Liste$Scan[groupleader[g]]+1])*2
      if (is.infinite(Grenzwert_FWHM_right[groupleader[g]])) 
        Grenzwert_FWHM_right[groupleader[g]] <- abs(daten@scantime[Liste$Scan[groupleader[g]]]-daten@scantime[Liste$Scan[groupleader[g]]+1])*2
      Summe_all[groupleader[g]] <- Grenzwert_RT[groupleader[g]]/2+Grenzwert_FWHM_left[groupleader[g]]/2+Grenzwert_FWHM_right[groupleader[g]]/2
    }
    })
  }
  
  
  
  groupleader <- groupleader[Liste$C13[groupleader] == 0]
  
  mzDiffCl <- 1.99705
  isotopePattern <- c(0.64,0.1)
  patternTol <- 0.05
  Liste$Cl2 <- 0
  
  
  for (g in 1:length(groupleader)) {
    #if (nrow(Liste[Liste$Gruppe == g,]) > 0) {  
    Cl2_Gruppe <- NULL
    for (n in 1:2) {
      if (Liste$Intensity[groupleader[g]]*isotopePattern[n] >= int_min & length(Cl2_Gruppe) == n-1) { 
        #die zu erwartendene Intensität muss über Threshold liegen und
        #alle vorherigen Isotope müssen gefunden worden sein
        # This line and others like it take a very long time for long peak lists
        Cl2_isotop <- which((abs(Liste$mz[groupleader]-Liste$mz[groupleader[g]]-mzDiffCl*n) < Liste$mz[groupleader[g]]/1000000*ppm) & (abs(Liste$Intensity[groupleader]/Liste$Intensity[groupleader[g]]-isotopePattern[n]) <= patternTol) &
                              abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader])<Grenzwert_RT[groupleader[g]] &                          # RT im Limit?
                              abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader])<Grenzwert_FWHM_left[groupleader[g]] &     # FWHM_left im Limit?
                              abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader])<Grenzwert_FWHM_right[groupleader[g]] &  # FWHM_right im Limit?
                              (abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader]) +
                                 abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader]) +
                                 abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader]))<Summe_all[groupleader[g]])
        if (length(Cl2_isotop) > 1) Cl2_isotop <- Cl2_isotop[which.min(abs(Liste$RT[groupleader[Cl2_isotop]]-Liste$RT[groupleader[g]]))]
        Cl2_Gruppe <- c(Cl2_Gruppe, Cl2_isotop)
      }
    }
    if (length(Cl2_Gruppe) > 0) {
      if (Liste$Cl2[groupleader[Cl2_Gruppe[1]]] == 0) { #wenn noch nicht als Cl2-Isotop zugeordnet
        Liste$Cl2[groupleader[Cl2_Gruppe]] <- groupleader[g]
        Liste$Gruppe[Liste$Gruppe == Liste$Gruppe[groupleader[Cl2_Gruppe[1]]]] <- Liste$Gruppe[groupleader[g]]
        if (length(Cl2_Gruppe) > 1) Liste$Gruppe[Liste$Gruppe == Liste$Gruppe[groupleader[Cl2_Gruppe[2]]]] <- Liste$Gruppe[groupleader[g]]
      } else { #wenn schon vorher als Cl2-Isotop zugeordnet, prüfen, welches besser passt
        if (abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader[Cl2_Gruppe]]) <= abs(Liste$RT[groupleader[Liste$Gruppe[groupleader[Cl2_Gruppe]]]]-Liste$RT[groupleader[Cl2_Gruppe]]) |
            abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader[Cl2_Gruppe]]) <= abs(Liste$FWHM_left[groupleader[Liste$Gruppe[groupleader[Cl2_Gruppe]]]]-Liste$FWHM_left[groupleader[Cl2_Gruppe]]) |
            abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader[Cl2_Gruppe]]) <= abs(Liste$FWHM_right[groupleader[Liste$Gruppe[groupleader[Cl2_Gruppe]]]]-Liste$FWHM_right[groupleader[Cl2_Gruppe]])) {
          #wenn RT, FWHM_left oder FWHM-right in der vorherigen Gruppe besser passen, bleibt es auch dort. Nur wenn alles besser passt kommt es in die neue Gruppe mit der niedrigeren Intensität.
          Liste$Cl2[groupleader[Cl2_Gruppe]] <- groupleader[g]
          Liste$Gruppe[Liste$Gruppe == Liste$Gruppe[groupleader[Cl2_Gruppe[1]]]] <- Liste$Gruppe[groupleader[g]]
          if (length(Cl2_Gruppe) > 1) Liste$Gruppe[Liste$Gruppe == Liste$Gruppe[groupleader[Cl2_Gruppe[2]]]] <- Liste$Gruppe[groupleader[g]]
        }
      }
    }
    #}
  }
  
  groupleader <- groupleader[Liste$Cl2[groupleader] == 0]
  
  mzDiffCl <- 1.99705
  isotopePattern <- 0.32
  patternTol <- 0.05
  Liste$Cl1 <- 0
  
  
  for (g in 1:length(groupleader)) {
    #if (nrow(Liste[Liste$Gruppe == g,]) > 0) {  
    Cl1_Gruppe  <- which((abs(Liste$mz[groupleader]-Liste$mz[groupleader[g]]-mzDiffCl) < Liste$mz[groupleader[g]]/1000000*ppm) & (abs(Liste$Intensity[groupleader]/Liste$Intensity[groupleader[g]]-isotopePattern) <= patternTol) &
                           abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader])<Grenzwert_RT[groupleader[g]] &                          # RT im Limit?
                           abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader])<Grenzwert_FWHM_left[groupleader[g]] &     # FWHM_left im Limit?
                           abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader])<Grenzwert_FWHM_right[groupleader[g]] &  # FWHM_right im Limit?
                           (abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader]) +
                              abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader]) +
                              abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader]))<Summe_all[groupleader[g]])
    if (length(Cl1_Gruppe) > 1) Cl1_Gruppe <- Cl1_Gruppe[which.min(abs(Liste$RT[groupleader[Cl1_Gruppe]]-Liste$RT[groupleader[g]]))]
    
    if (length(Cl1_Gruppe) > 0) {
      if (Liste$Cl1[groupleader[Cl1_Gruppe[1]]] == 0) { #wenn noch nicht als Cl1-Isotop zugeordnet
        Liste$Cl1[groupleader[Cl1_Gruppe]] <- groupleader[g]
        Liste$Gruppe[Liste$Gruppe == Liste$Gruppe[groupleader[Cl1_Gruppe[1]]]] <- Liste$Gruppe[groupleader[g]]
      } else { #wenn schon vorher als Cl2-Isotop zugeordnet, prüfen, welches besser passt
        if (abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader[Cl1_Gruppe]]) <= abs(Liste$RT[groupleader[Liste$Gruppe[groupleader[Cl1_Gruppe]]]]-Liste$RT[groupleader[Cl1_Gruppe]]) |
            abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader[Cl1_Gruppe]]) <= abs(Liste$FWHM_left[groupleader[Liste$Gruppe[groupleader[Cl1_Gruppe]]]]-Liste$FWHM_left[groupleader[Cl1_Gruppe]]) |
            abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader[Cl1_Gruppe]]) <= abs(Liste$FWHM_right[groupleader[Liste$Gruppe[groupleader[Cl1_Gruppe]]]]-Liste$FWHM_right[groupleader[Cl1_Gruppe]])) {
          #wenn RT, FWHM_left oder FWHM-right in der vorherigen Gruppe besser passen, bleibt es auch dort. Nur wenn alles besser passt kommt es in die neue Gruppe mit der niedrigeren Intensität.
          Liste$Cl1[groupleader[Cl1_Gruppe]] <- groupleader[g]
          Liste$Gruppe[Liste$Gruppe == Liste$Gruppe[groupleader[Cl1_Gruppe]]] <- Liste$Gruppe[groupleader[g]]
        }
      }
    }
    #}
  }
  
  
  groupleader <- groupleader[Liste$Cl1[groupleader] == 0]
  
  
  
  # Na-Addukte
  Liste$NaAddukt <- 0
  
  
  for (g in 1:length(groupleader)) {
    
    #if (nrow(Liste[Liste$Gruppe == g,]) > 0) {  
    NaAdduktGruppe <- which((abs(Liste$mz[groupleader]-Liste$mz[groupleader[g]]-mzNa) < Liste$mz[groupleader[g]]/1000000*ppm) &
                              abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader])<Grenzwert_RT[groupleader[g]] &                          # RT im Limit?
                              abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader])<Grenzwert_FWHM_left[groupleader[g]] &     # FWHM_left im Limit?
                              abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader])<Grenzwert_FWHM_right[groupleader[g]] &  # FWHM_right im Limit?
                              (abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader]) +
                                 abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader]) +
                                 abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader]))<Summe_all[groupleader[g]])
    if (length(NaAdduktGruppe) > 0)  {
      
      #wenn mehrere in Frage kommen, nimm den, der gemessen an der RT am besten passt:
      if (length(NaAdduktGruppe) > 1) NaAdduktGruppe <- NaAdduktGruppe[which.min(abs(Liste$RT[groupleader[NaAdduktGruppe]]-Liste$RT[groupleader[g]]))]
      
      if (Liste$NaAddukt[groupleader[NaAdduktGruppe]] == 0) { #wenn noch nicht als Na-Addukt zugeordnet
        Liste$NaAddukt[groupleader[NaAdduktGruppe]] <- groupleader[g]
        Liste$Gruppe[Liste$Gruppe == Liste$Gruppe[groupleader[NaAdduktGruppe]]] <- Liste$Gruppe[groupleader[g]]
      } else { #wenn schon vorher als 13C zugeordnet, prüfen, welches besser passt
        if (abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader[NaAdduktGruppe]]) <= abs(Liste$RT[groupleader[Liste$Gruppe[groupleader[NaAdduktGruppe]]]]-Liste$RT[groupleader[NaAdduktGruppe]]) |
            abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader[NaAdduktGruppe]]) <= abs(Liste$FWHM_left[groupleader[Liste$Gruppe[groupleader[NaAdduktGruppe]]]]-Liste$FWHM_left[groupleader[NaAdduktGruppe]]) |
            abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader[NaAdduktGruppe]]) <= abs(Liste$FWHM_right[groupleader[Liste$Gruppe[groupleader[NaAdduktGruppe]]]]-Liste$FWHM_right[groupleader[NaAdduktGruppe]])) {
          #wenn RT, FWHM_left oder FWHM-right in der vorherigen Gruppe besser passen, bleibt es auch dort. Nur wenn alles besser passt kommt es in die neue Gruppe mit der niedrigeren Intensität.
          Liste$NaAddukt[groupleader[NaAdduktGruppe]] <- groupleader[g]
          Liste$Gruppe[Liste$Gruppe == Liste$Gruppe[groupleader[NaAdduktGruppe]]] <- Liste$Gruppe[groupleader[g]]
        }
      }
      #}
    }
    
    #}
  }
  # Ende Na-Addukte
  groupleader <- groupleader[Liste$NaAddukt[groupleader] == 0]
  
  Liste$groupleader <- FALSE
  
  for (g in 1:length(groupleader)) {
    
    if (!any(Liste$groupleader[Liste$Gruppe == Liste$Gruppe[groupleader[g]]])) {  #nur wenn nicht schon einem anderen groupleader zugeordnet
      KandidatenGruppen <- which(Liste$Gruppe[groupleader] > Liste$Gruppe[groupleader[g]] &
                                   abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader])<Grenzwert_RT[groupleader[g]] &                # RT im Limit?
                                   abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader])<Grenzwert_FWHM_left[groupleader[g]] &     # FWHM_left im Limit?
                                   abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader])<Grenzwert_FWHM_right[groupleader[g]] &  # FWHM_right im Limit?
                                   (abs(Liste$RT[groupleader[g]]-Liste$RT[groupleader]) +
                                      abs(Liste$FWHM_left[groupleader[g]]-Liste$FWHM_left[groupleader]) +
                                      abs(Liste$FWHM_right[groupleader[g]]-Liste$FWHM_right[groupleader]))<Summe_all[groupleader[g]])
      Liste$Gruppe[Liste$Gruppe %in% Liste$Gruppe[groupleader[KandidatenGruppen]]] <- Liste$Gruppe[groupleader[g]]
      Liste$groupleader[groupleader[g]] <- TRUE 
      
    }
    
  }
  
  reihenfolge <- order(Liste$mz)
  reihenfolge <- order(reihenfolge)
  for (i in 1:nrow(Liste)) {
    if (Liste$C13[i] > 0) {
      Liste$C13[i] <- reihenfolge[Liste$C13[i]]
    } else {
      Liste$C13[i] <- 0
    }
    if (Liste$NaAddukt[i] > 0) {
      Liste$NaAddukt[i] <- reihenfolge[Liste$NaAddukt[i]]
    } else {
      Liste$NaAddukt[i] <- 0
    }
    if (Liste$Cl1[i] > 0) {
      Liste$Cl1[i] <- reihenfolge[Liste$Cl1[i]]
    } else {
      Liste$Cl1[i] <- 0
    }
    if (Liste$Cl2[i] > 0) {
      Liste$Cl2[i] <- reihenfolge[Liste$Cl2[i]]
    } else {
      Liste$Cl2[i] <- 0
    }
  }
  
  Liste <- Liste[order(Liste$mz), ] 
  
  return(Liste)
}

#' Second stage componentisation using the alignment table
#' 
#' Based on four criteria: rt, peak shape, intensity correlation and common mz differences.
#' Will replace group column in the alignment table. Distance matrices for all 
#' four criteria are combined and then DBSCAN clusters the features.
#' 
#' @param altable alignment table matrix ("grouped" in the app)
#' @param rttols rt tolerance in s
#' @param fracComponMatch fraction of samples with matching peak shape
#' @param mztol mz tolerance
#' @param pearsonCorr minimum pearson's r for comparing intensity trends
#' @param pol polarity, must be either "pos" or "neg" 
#'
#' @return Alignment table with the group column "Gruppe" replaced with the new values
#' @export
alig_componentisation <- function(altable, rttols = 3, fracComponMatch = 0.5, mztol = 0.005, pearsonCorr = 0.5, pol = "pos") {
  stopifnot(pol %in% c("pos", "neg"))
  # RT difference
  
  message("Computing RT distance matrix")
  rt_comp <- RcppXPtrUtils::cppXPtr(
    sprintf("double customDist(const arma::mat &A, const arma::mat &B) {
    double rt1 = as_scalar(A);
    double rt2 = as_scalar(B);
    double answer;
    if (std::abs(rt1 - rt2) <= %i) {
      answer = 0;
    } else {
      answer = 3;  // rt must match, extra penalty if not
    }
    return(answer);
  }", rttols), 
    depends = c("RcppArmadillo")
  )
  
  rtDistT <- parallelDist::parDist(altable[, "mean_RT", drop = F], 
                                   method = "custom", func = rt_comp)
  rtDist <- as.integer(rtDistT)
  attributes(rtDist) <- attributes(rtDistT)
  rm(rtDistT)
  
  # componentisation grouping from peaklist
  # the componentisation from the peaklist needs to be improved using dtw correlation.
  message("Computing peak shape componentisation distance matrix")
  group_comp <- RcppXPtrUtils::cppXPtr(
    sprintf("double customDist(const arma::mat &A, const arma::mat &B) {
      arma::rowvec x = A.row(0);
      arma::rowvec y = B.row(0);
      double answer = 3;
      // remove any entries with 0
      arma::uvec keepx = arma::find(x);
      arma::uvec keepy = arma::find(y);
      arma::uvec keep = arma::intersect(keepx, keepy);
      
      if (keep.n_elem < 3) {  // need to have at least 3 values to work with
        return(answer);
      }
      arma::vec z(keep.n_elem, arma::fill::zeros);
      z = x.elem(keep) - y.elem(keep);
  
      arma::uvec isSameGroup = arma::find(z == 0);
      double fractionSame = isSameGroup.n_elem / (double)z.n_elem;
      if (fractionSame > %.2f) {  // more than half of samples show same group
        answer = 0;
      } else if (fractionSame > %.2f) {  // can still be compensated for
        answer = 1;
      } else {
        answer = 3;  // higher penalty if never in similar component
      }
      // std::cout << fractionSame << std::endl;
      return(answer);
    }", fracComponMatch, fracComponMatch / 2), depends = c("RcppArmadillo")
  )
  all_gruppe <- altable[, grep("^gruppe_\\d+$", colnames(altable)), drop = F]
  groupDistT <- parallelDist::parDist(all_gruppe, method = "custom", func = group_comp)
  
  shapeDist <- as.integer(groupDistT)
  attributes(shapeDist) <- attributes(groupDistT)
  
  rm(groupDistT)
  
  # known mass differences
  message("Computing known adduct and isotopologue mz distance matrix")
  mass_comp <- switch(pol, 
                      pos = RcppXPtrUtils::cppXPtr(
                        sprintf("double customDist(const arma::mat &A, const arma::mat &B) {
                    double mz1 = as_scalar(A);
                    double mz2 = as_scalar(B);
                    double mzDiff = std::abs(mz1 - mz2);
                    double answer = 0;
                    arma::vec knownMasses = arma::zeros<arma::vec>(10);
                    knownMasses(0) = 1.003;  // 13C
                    knownMasses(1) = 21.982;
                    knownMasses(2) = 18.011;
                    knownMasses(3) = 23.075;
                    knownMasses(4) = 39.993;
                    knownMasses(5) = 22.985;
                    knownMasses(6) = 20.979;
                    knownMasses(7) = 45.057;
                    knownMasses(8) = 1.997;
                    knownMasses(9) = 2.005;
                    
                    arma::vec mzDiffDiff = arma::ones<arma::vec>(10);
                    for (int i = 0; i < knownMasses.n_elem; ++i) {
                      mzDiffDiff(i) = std::abs(knownMasses[i] - mzDiff); 
                    }
                    
                    arma::uvec isAdduct = arma::find(mzDiffDiff <= %.4f);
                    if (isAdduct.n_elem == 1) {
                      answer = 0;
                    } else {
                      answer = 1;
                    }

                    return(answer);
                  }", mztol), depends = c("RcppArmadillo")
                      ),
                      neg = RcppXPtrUtils::cppXPtr(
                        sprintf("double customDist(const arma::mat &A, const arma::mat &B) {
                        double mz1 = as_scalar(A);
                        double mz2 = as_scalar(B);
                        double mzDiff = std::abs(mz1 - mz2);
                        double answer = 0;
                        arma::vec knownMasses = arma::zeros<arma::vec>(8);
                        knownMasses(0) = 1.003;
                        knownMasses(1) = 1.997;
                        knownMasses(2) = 16.990;
                        knownMasses(3) = 43.989;
                        knownMasses(4) = 21.982;
                        knownMasses(5) = 46.005;
                        knownMasses(6) = 0.998;
                        knownMasses(7) = 15.987;
                        
                        arma::vec mzDiffDiff = arma::ones<arma::vec>(8);
                        for (int i = 0; i < knownMasses.n_elem; ++i) {
                          mzDiffDiff(i) = std::abs(knownMasses[i] - mzDiff); 
                        }
                        
                        arma::uvec isAdduct = arma::find(mzDiffDiff <= %.4f);
                        if (isAdduct.n_elem == 1) {
                          answer = 0;
                        } else {
                          answer = 1;
                        }
                        
                        return(answer);
                    }", mztol), depends = c("RcppArmadillo"))
  )
  mzDistT <- parallelDist::parDist(altable[, "mean_mz", drop = F],
                                   method = "custom", func = mass_comp)
  
  mzDist <- as.integer(mzDistT)
  attributes(mzDist) <- attributes(mzDistT)
  rm(mzDistT)
  
  # intenstiy correlation
  message("Computing intensity correlation distance matrix")
  pearsonCompCustom <- RcppXPtrUtils::cppXPtr(
    "double customDist(const arma::mat &A, const arma::mat &B) {
      double summe = 0;
      arma::rowvec X = A.row(0);
      arma::rowvec Y = B.row(0);
      // must have at least 5 non-zero intensities in at least one vector
      arma::uvec keepx = arma::find(X);
      arma::uvec keepy = arma::find(Y);
      if (keepx.n_elem < 5 && keepy.n_elem < 5) {
        return(summe);
      }
      for(int i = 0; i < X.n_elem; ++i) {
        summe += (X[i]  - arma::mean(X))*(Y[i] - arma::mean(Y));
      }
      return summe / (X.n_elem * arma::stddev(X) * arma::stddev(Y));
    }", depends = c("RcppArmadillo")
  )
  
  intens <- altable[, grep("^Int_", colnames(altable))]
  sumStat <- apply(intens, 1, max)
  intens <- sweep(intens, 1, sumStat, "/")
  intens <- as.matrix(intens)
  distM <- parallelDist::parDist(intens, method = "custom", func = pearsonCompCustom)
  
  corrDist <- ifelse(distM >= pearsonCorr, 0L, 1L)
  attributes(corrDist) <- attributes(distM)
  rm(distM)
  
  # combination
  message("Combining distance matrices and clustering")
  # if mz fits, you can afford to have one of the other criteria (corr or shape) fail.
  combiDist <- ifelse(mzDist == 1L, corrDist + shapeDist + rtDist, 
                      corrDist + shapeDist + rtDist - 1L)
  attributes(combiDist) <- attributes(corrDist)
  rm(mzDist, corrDist, shapeDist, rtDist)
  
  #dbscan::kNNdistplot(combiDist, 1)
  #browser()
  dbscan_res <- dbscan::dbscan(combiDist, 0.1, 2)
  # message(sprintf("DBSCAN clustering for %i objects.
  #                 The clustering contains %i cluster(s) and %i noise points.",
  #                 nrow(altable), 
  #                 length(table(dbscan_res$cluster))-1, 
  #                 table(dbscan_res$cluster)["0"]))
  altable[, "Gruppe"] <- dbscan_res$cluster
  #clusterSize <- unname(table(dbscan_res$cluster)[-1])
  #hist(clusterSize, breaks = seq(2, max(clusterSizes), 2), 
  #     main = "Histogram of cluster sizes", xlab = "Number of features per cluster")
  message("Completed 2nd stage componentisation")
  return(altable)
}
