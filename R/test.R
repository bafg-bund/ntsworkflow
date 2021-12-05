# test Toni

#' Replikate
#'
#'  Diese Funktion vereinigt Replikate in der Messung unter der Vorgabe, dass ein Peak in allen Replikaten vorkommen muss
#'  Dabei sollen Blanks nicht wereinigt werden und müssen deswegen mit "Anzahl_Blank" angegeben und am Ende der
#'  input-peaklist stehen
#'  
#'  Nur für Rep=2 bis Rep=4 programmiert!
#'
#' @param peaklist_r  aus der Non-target App
#' @param ppm erlaubte Massendifferenz zwischen 2 Peaks
#' @param DeltaRT erlaubt RT-differenz zwischen 2 Peaks
#' @param Rep Anzahl der zu einer Probe gehörenden Replikate (z.B. Triplikate, Rep=3)
#'
#' @return Peaklist mit vereinigten Replikaten und Blanks
#' @export
#' 
#'
#' @examples Messung von 5 Proben als Triplikat mit 3 Blanks
#'           Rep = 3, Anzahl_Proben = 15, Anzahl_Blank = 3
Replikate <- function(peaklist_r=peaklist, ppm=5, DeltaRT=3,Rep){
  
  # Welche und wieviele Messungen sind Probe und Blank?
  blanks <- NULL
  proben <- NULL
  for (i in (1:nrow(sampleList))) {
    if (sampleList[i, "sampleType"] == "Blank") blanks <- c(blanks,i)}
  
  for (i in (1:nrow(sampleList))) {
    if (sampleList[i, "sampleType"] == "Unknown") proben <- c(proben,i)}
  
  # Sicherungen 
  #if (length(proben) + length(blanks) != length(peaklist_r)){stop("Die Anzahl der Proben und Blanks stimmen nicht mit der peaklist überein!!!")}
  if ((length(proben)%%Rep)>0) {stop("Die Anzahl der Proben ist nicht ganzzalig durch die Zahl der Replikate teilbar!!!")}
  
  # Replikatefindung (nur Proben, keine Blanks)
  # passen mehrere Peaks aus Probe 2 zu einem Peak aus Probe 1, gewinnt der Peaks aus Probe 2 mit der ähnlicheren Intensität zu dem aus Probe 1
  ergebnis <- vector("list", length(peaklist_r))
  
  Rep1 <- 1
  
  while(Rep1<length(proben)){
    #for (i in i:(length(peaklist_r)/Rep)){
    
    switch(as.character(Rep),
           "2"={j <- 1
           for (j in j:nrow(peaklist_r[[proben[Rep1]]])){
             
             p2 <- which((abs(peaklist_r[[proben[Rep1]]][j,"mz"]-peaklist_r[[proben[Rep1]+1]][,"mz"]) < (peaklist_r[[proben[Rep1]]][j,"mz"]*ppm/1000000) & (abs(peaklist_r[[proben[Rep1]]][j,"RT"]-peaklist_r[[proben[Rep1]+1]][,"RT"])<DeltaRT)))
             
             if (length(p2)>1){
               int_p1 <- peaklist_r[[proben[Rep1]]][j,"Intensity"]
               p2 <-which.min(peaklist_r[[proben[Rep1]]][j,"Intensity"]-peaklist_r[[proben[Rep1]+1]][p2,"Intensity"])}
             
             if (length(p2)>0){
               if (p2>0){  
                 ergebnis[[proben[Rep1]]]   <- rbind(ergebnis[[proben[Rep1]]],peaklist_r[[proben[Rep1]]][j,])
                 ergebnis[[proben[Rep1]+1]] <- rbind(ergebnis[[proben[Rep1]+1]],peaklist_r[[proben[Rep1]+1]][p2,]) }}
           }},
           
           
           "3"={j <- 1
           for (j in j:nrow(peaklist_r[[proben[Rep1]]])){
             
             p2 <- which((abs(peaklist_r[[proben[Rep1]]][j,"mz"]-peaklist_r[[proben[Rep1]+1]][,"mz"]) < (peaklist_r[[proben[Rep1]]][j,"mz"]*ppm/1000000) & (abs(peaklist_r[[proben[Rep1]]][j,"RT"]-peaklist_r[[proben[Rep1]+1]][,"RT"])<DeltaRT)))
             p3 <- which((abs(peaklist_r[[proben[Rep1]]][j,"mz"]-peaklist_r[[proben[Rep1]+2]][,"mz"]) < (peaklist_r[[proben[Rep1]]][j,"mz"]*ppm/1000000) & (abs(peaklist_r[[proben[Rep1]]][j,"RT"]-peaklist_r[[proben[Rep1]+2]][,"RT"])<DeltaRT)))
             
             if (length(p2)>1){
               int_p1 <- peaklist_r[[proben[Rep1]]][j,"Intensity"]
               p2 <-which.min(peaklist_r[[proben[Rep1]]][j,"Intensity"]-peaklist_r[[proben[Rep1]+1]][p2,"Intensity"])}
             
             if (length(p3)>1){
               int_p1 <- peaklist_r[[proben[Rep1]]][j,"Intensity"]
               p3 <-which.min(peaklist_r[[proben[Rep1]]][j,"Intensity"]-peaklist_r[[proben[Rep1]+1]][p3,"Intensity"])}
             
             if (length(p2)>0 && length(p3)>0){
               if (p2>0 && p3>0){  
                 ergebnis[[proben[Rep1]]]   <- rbind(ergebnis[[proben[Rep1]]],peaklist_r[[proben[Rep1]]][j,])
                 ergebnis[[proben[Rep1]+1]] <- rbind(ergebnis[[proben[Rep1]+1]],peaklist_r[[proben[Rep1]+1]][p2,])
                 ergebnis[[proben[Rep1]+2]] <- rbind(ergebnis[[proben[Rep1]+2]],peaklist_r[[proben[Rep1]+2]][p3,]) }}
           }}, 
           
           
           "4"={j <- 1
           for (j in j:nrow(peaklist_r[[proben[Rep1]]])){
             
             p2 <- which((abs(peaklist_r[[proben[Rep1]]][j,"mz"]-peaklist_r[[proben[Rep1]+1]][,"mz"]) < (peaklist_r[[proben[Rep1]]][j,"mz"]*ppm/1000000) & (abs(peaklist_r[[proben[Rep1]]][j,"RT"]-peaklist_r[[proben[Rep1]+1]][,"RT"])<DeltaRT)))
             p3 <- which((abs(peaklist_r[[proben[Rep1]]][j,"mz"]-peaklist_r[[proben[Rep1]+2]][,"mz"]) < (peaklist_r[[proben[Rep1]]][j,"mz"]*ppm/1000000) & (abs(peaklist_r[[proben[Rep1]]][j,"RT"]-peaklist_r[[proben[Rep1]+2]][,"RT"])<DeltaRT)))
             p4 <- which((abs(peaklist_r[[proben[Rep1]]][j,"mz"]-peaklist_r[[proben[Rep1]+3]][,"mz"]) < (peaklist_r[[proben[Rep1]]][j,"mz"]*ppm/1000000) & (abs(peaklist_r[[proben[Rep1]]][j,"RT"]-peaklist_r[[proben[Rep1]+3]][,"RT"])<DeltaRT)))
             
             if (length(p2)>1){
               int_p1 <- peaklist_r[[proben[Rep1]]][j,"Intensity"]
               p2 <-which.min(peaklist_r[[proben[Rep1]]][j,"Intensity"]-peaklist_r[[proben[Rep1]+1]][p2,"Intensity"])}
             
             if (length(p3)>1){
               int_p1 <- peaklist_r[[proben[Rep1]]][j,"Intensity"]
               p3 <-which.min(peaklist_r[[proben[Rep1]]][j,"Intensity"]-peaklist_r[[proben[Rep1]+1]][p3,"Intensity"])}
             
             if (length(p4)>1){
               int_p1 <- peaklist_r[[proben[Rep1]]][j,"Intensity"]
               p4 <-which.min(peaklist_r[[proben[Rep1]]][j,"Intensity"]-peaklist_r[[proben[Rep1]+1]][p4,"Intensity"])}
             
             if (length(p2)>0 && length(p3)>0 && length(p4)>0){
               if (p2>0 && p3>0 && p4>0){  
                 ergebnis[[proben[Rep1]]]   <- rbind(ergebnis[[proben[Rep1]]],peaklist_r[[proben[Rep1]]][j,])
                 ergebnis[[proben[Rep1]+1]] <- rbind(ergebnis[[proben[Rep1]+1]],peaklist_r[[proben[Rep1]+1]][p2,])
                 ergebnis[[proben[Rep1]+2]] <- rbind(ergebnis[[proben[Rep1]+2]],peaklist_r[[proben[Rep1]+2]][p3,])
                 ergebnis[[proben[Rep1]+3]] <- rbind(ergebnis[[proben[Rep1]+3]],peaklist_r[[proben[Rep1]+3]][p4,]) }}
           }},
           
           {stop("mehr als 4 Replikate wurden leider nicht programmiert")})
    
    Rep1 <- Rep1+Rep
  }
  
  # Anhängen der Blanks an ihren vorherigen Platz in der Peaklist 
  if (length(blanks)>0){
    
    for (i in 1:length(blanks)){
      ergebnis[blanks[i]] <- peaklist_r[blanks[i]]}}
  
  
  return(ergebnis) 
  
}  




#' Gruppierung_Summe_Grenzwert
#' 
#' ordnet die Peaks in Gruppen (beginnend bei den intensivsten)
#' die Gruppen sollten neben dem [M+H]+ Peak die Addukte, Isotope und Fragmente enthalten
#' Die Spalte "Gruppe" wird an die Peaklist angehangen
#' 
#' Für das Alignment wird die Funktion grouping_BfG2 benötigt
#'
#' @param peaklist_df Peaklist aus der Non-target App, oder aus der Funktion Replikate, egal, hauptsache Peaklist
#' @param mit_RT soll die RT für die Gruppierung mit benutzt weden (TRUE), oder nur die FWHM_left und FWHM_right (FALSE)
#' @param Grenzwert_RT erlaubte RT-differenz zwischen 2 Peaks
#' @param Grenzwert_FWHM_left erlaubte FWHM_left-differenz zwischen 2 Peaks
#' @param Grenzwert_FWHM_right erlaubte FWHM_right-differenz zwischen 2 Peaks
#' @param Summe_all (bei mit_RT=TRUE) erlaubte maximale Summe aus RT-diff, FWHM_left-diff und FWHM_right-diff zwischen 2 Peaks 
#' @param Summe_FWHM (bei mit_RT=FALSE) erlaubte maximale Summe aus FWHM_left-diff und FWHM_right-diff zwischen 2 Peaks 
#'
#' @return Peaklist mit der Information der Gruppenzuordnung
#' @export
#'
#' @examples
Gruppierung_SG <- function(peaklist_df, 
                           mit_RT=TRUE, Grenzwert_RT=1.50, 
                           Grenzwert_FWHM_left=0.5, 
                           Grenzwert_FWHM_right=1, Summe_all=2.5, Summe_FWHM=1.5) {
  
  
  Liste <- peaklist_df
  Liste$Gruppe <- 1000000
  Liste <- Liste[order(Liste$Intensity,decreasing = T), ] 
  i <- 1
  g <- 0
  
  if (mit_RT == TRUE) {
    
    while(any(Liste$Gruppe > 10000)) {
      
      g <- g+1  #nächste Gruppe
      
      i <- which(Liste$Gruppe > 10000)[1]   #nächstes Ungruppiertes
      
      passende <- which((Liste$Gruppe > 10000) &                                            # nur die, die noch nicht gruppiert wurden
                          abs(Liste$RT[i]-Liste$RT)<Grenzwert_RT &                          # RT im Limit?
                          abs(Liste$FWHM_left[i]-Liste$FWHM_left)<Grenzwert_FWHM_left &     # FWHM_left im Limit?
                          abs(Liste$FWHM_right[i]-Liste$FWHM_right)<Grenzwert_FWHM_right &  # FWHM_right im Limit?
                          (abs(Liste$RT[i]-Liste$RT) +
                             abs(Liste$FWHM_left[i]-Liste$FWHM_left) +
                             abs(Liste$FWHM_right[i]-Liste$FWHM_right))<Summe_all)
      Liste$Gruppe[passende] <- g
      
    }
    
    
  } else {
    
    while(any(Liste$Gruppe > 10000)) {
      
      g <- g+1  #nächste Gruppe
      
      i <- which(Liste$Gruppe > 10000)[1]   #nächstes Ungruppiertes
      
      passende <- which((Liste$Gruppe > 10000) &                                            # nur die, die noch nicht gruppiert wurden
                          abs(Liste$FWHM_left[i]-Liste$FWHM_left)<Grenzwert_FWHM_left &     # FWHM_left im Limit?
                          abs(Liste$FWHM_right[i]-Liste$FWHM_right)<Grenzwert_FWHM_right &  # FWHM_right im Limit?
                          (abs(Liste$RT[i]-Liste$RT) +
                             abs(Liste$FWHM_left[i]-Liste$FWHM_left) +
                             abs(Liste$FWHM_right[i]-Liste$FWHM_right))<Summe_all)
      Liste$Gruppe[passende] <- g
      
    }

  }
  
  Liste <- Liste[order(Liste$mz, decreasing = F), ]  
  Liste
  
  
}



#' grouping_BfG_ohne_datenlist
#' 
#' Funktion für das Alignment der Ergebnispeaklist der Gruppierungsfunktion Gruppierung_SG
#' Die Spalte Gruppe wird für jede Probe in die Alignmenttabelle integriert
#' der MS2 Vergleich ist deaktiviert, da die Peaklist in der Funktion Replikate gekürzt wirde und ein Verweis in
#' die Datenlist der Non-target App fehlschlägt
#'
#' @param peaklists Peaklisten aus der Funktion Gruppierung_SG
#' @param ppm erlaubte Massendifferenz zwischen 2 Peaks
#' @param DeltaRT erlaubte RT-differenz zwischen 2 Peaks
#' @param ms2ppm (DEAKTIVIERT) erlaubte Massendifferenz im MS2-Spektrum zwischen 2 Peaks
#'
#' @return Alignmenttabelle mit zusätzlicher Gruppeninformation
#' @export
#'
#' @examples
grouping_BfG_ohne_datenlist <- function(peaklists, ppm=10, DeltaRT=10, ms2ppm=10){
  #browser()
  ende <- FALSE
  spaltenlaenge <- 1 #nrow der l?ngsten Peaklist
  SuchStart <- 1
  
  ergebnis <- matrix(, nrow = 0, ncol = length(peaklists)*6)
  masse <- 0
  
  ii <- 0
  ms2vergleich <- vector()
  
  for (i in 1:length(peaklists)) {
    if (spaltenlaenge < nrow(peaklists[[i]])) 
      spaltenlaenge <- nrow(peaklists[[i]])
    peaklists[[i]][, "Scan"] <- 0
  }
  
  while (ende == FALSE) {
    ende <- TRUE
    iii <- SuchStart - 1
    jjj <<- iii
    
    ii <- ii + 1
    ergebnis <- rbind(ergebnis,0)
    abbruch <- FALSE
    ms2spectra <- list()
    
    # KJ: gets the first peak from the first peaklist as mz and rt
    while (abbruch == FALSE) {
      iii <- iii + 1
      for (i in 1:length(peaklists)) {
        if (iii <= nrow(peaklists[[i]])) {
          if ((peaklists[[i]][iii, "Scan"] == 0) & (abbruch == FALSE)) {
            masse <- peaklists[[i]][iii, "mz"]
            RT <- peaklists[[i]][iii, "RT"]
            abbruch <- TRUE
            ende <- FALSE
          }  
        }
        if (iii > spaltenlaenge) 
          abbruch <- TRUE
      }
      SuchStart <- iii
      
    }#end while abbruch == false
    
    # get that mz and rt of the "highest" peak in all files
    kandidaten <- list()
    highest_i <- 1
    for (i in 1:length(peaklists)) {
      kandidaten[[i]] <- which((abs(peaklists[[i]][, "mz"] - masse) < (masse*ppm/1000000)) & 
                                 (abs(peaklists[[i]][, "RT"] - RT) < DeltaRT) & 
                                 (peaklists[[i]][, "Scan"] == 0))
      #order by intensity:
      kandidaten[[i]] <- kandidaten[[i]][order(peaklists[[i]][kandidaten[[i]], "Intensity"], decreasing=TRUE)]
      if (length(kandidaten[[highest_i]]) > 0) {
        if (length(kandidaten[[i]]) > 0) {
          if (peaklists[[i]][kandidaten[[i]][1], "Intensity"] > peaklists[[highest_i]][kandidaten[[highest_i]][1], "Intensity"]) 
            highest_i <- i
        } 
      } else {
        highest_i <- i
      }
    }  
    
    masse <- peaklists[[highest_i]][kandidaten[[highest_i]][1],1]
    RT <- peaklists[[highest_i]][kandidaten[[highest_i]][1],2]
    
    #browser(expr = !is.na(masse) && isTRUE(all.equal(masse, 237.1016, tolerance = 0.005, scale = 1)))
    # KJ only check files in which this peak was found
    zuPruefende_i <- which(lengths(kandidaten) > 0)
    
    for (i in zuPruefende_i) {
      # if this is not the sample with the highest peak, check candidates again according to new mz and rt
      if (i != highest_i) {
        kandidaten[[i]] <- which((abs(peaklists[[i]][, "mz"]-masse) < (masse*ppm/1000000)) & 
                                   (abs(peaklists[[i]][, "RT"]-RT)<DeltaRT) & (peaklists[[i]][, "Scan"] == 0))
      }  
      
      if (length(kandidaten[[i]]) > 0) {
        passendeZeile <- kandidaten[[i]][which.min(abs(peaklists[[i]][kandidaten[[i]], "RT"] - 
                                                         peaklists[[highest_i]][kandidaten[[highest_i]][1], "RT"]))]
        ergebnis[ii,i*6-5] <- peaklists[[i]][passendeZeile, "peak_id_all"]  # peak ID from peaklist
        ergebnis[ii,i*6-4] <- peaklists[[i]][passendeZeile, "mz"]
        ergebnis[ii,i*6-3] <- peaklists[[i]][passendeZeile, "RT"]
        ergebnis[ii,i*6-2] <- peaklists[[i]][passendeZeile, "Intensity"]
        ergebnis[ii,i*6-1] <- peaklists[[i]][passendeZeile, "MS2scan"]
        ergebnis[ii,i*6] <- peaklists[[i]][passendeZeile, "Gruppe"]
        peaklists[[i]][passendeZeile, "Scan"] <- 1
      }
      
      #      if ((ergebnis[ii,i*6-1] > 0)) {
      #        RTMS2 <- datenList[[i]]@msnRt[peaklists[[i]][passendeZeile, "MS2scan"]]
      #        if ((RTMS2 > peaklists[[i]][passendeZeile, "LeftendRT"]) & 
      #            (RTMS2 < peaklists[[i]][passendeZeile, "RightendRT"])) {
      #          
      #          ms2spectra[[i]] <- xcms::getMsnScan(datenList[[i]], ergebnis[ii,i*6-1])
      #          if (any(is.na(ms2spectra[[i]]))) 
      #            ms2spectra[[i]] <- NULL
      #        } else {
      #          ms2spectra[[i]] <- NULL
      #          
      #          ergebnis[ii,i*6-1] <- 0
      #        }
      #      }  
    }  
    
    #remove empty ms2 spectra:
    ms2spectra <- Filter(Negate(function(x) is.null(unlist(x))), ms2spectra)
    
    #if there are at least two samples with ms2 spectra, compare them:
    if (length(ms2spectra) > 1) {
      ms2vergleich[ii] <- spektraVergleichen(ms2spectra,ms2ppm)
    } else {
      ms2vergleich[ii] <- 0
    }
    
    
  }#end while ende == false
  
  
  mean_mz <- vector()
  mean_RT <- vector()
  Gruppe <- numeric(nrow(ergebnis))
  for (i in 1:nrow(ergebnis)) {
    mean_mz[i] <- mean(ergebnis[i,which(ergebnis[i,seq(2,ncol(ergebnis),by=6)] > 0)*6-4])
    mean_RT[i] <- mean(ergebnis[i,which(ergebnis[i,seq(3,ncol(ergebnis),by=6)] > 0)*6-3])
  }
  ergebnis <- cbind(mean_mz,mean_RT,ms2vergleich,Gruppe,ergebnis)
  
  spaltennamen <- c("mean_mz","mean_RT","MS2Fit", "Gruppe")
  for (i in 1:(ncol((ergebnis)-4)/6)) {
    spaltennamen <- c(spaltennamen,
                      paste0("peak_id_all_",as.character(i)),
                      paste0("mz_",as.character(i)),
                      paste0("RT_",as.character(i)),
                      paste0("Int_",as.character(i)),
                      paste0("ms2scan_",as.character(i)),
                      paste0("gruppe_",as.character(i)))  # componentization information
  }
  colnames(ergebnis) <- spaltennamen
  
  # Vereinheitlichung der Gruppen
  gruppenzaehler <- 1
  
  for (i in 1:nrow(ergebnis)){
    if (ergebnis[i,"Gruppe"]==0){
      # max_Int = col with highest Int in row i
      max_Int <- grep("Int_", colnames(ergebnis)) [which.max(ergebnis[i,grep("Int_", colnames(ergebnis))])]
      max_gruppe <- paste0("gruppe_",stringr::str_match(colnames(ergebnis)[max_Int], "Int_(\\d+)$")[,2])
      # b = rows with Gruppe == 0 AND same groups found in first samples of a
      b <- which(ergebnis[,"Gruppe"]==0 & ergebnis[,max_gruppe]==ergebnis[i,max_gruppe])
      ergebnis[b,"Gruppe"] <- gruppenzaehler
      gruppenzaehler <- gruppenzaehler+1}}
  
  # Eliminieren von ein-Peak-Gruppen  (Gruppe wird 0)
  for (i in 1:max(ergebnis[,"Gruppe"])){
    if (length(which(ergebnis[,"Gruppe"]==i))==1) ergebnis[which(ergebnis[,"Gruppe"]==i),"Gruppe"] <-0}
  
  # Schließen der entstandenen Zwischenräume
  V1 <- sort(unique(ergebnis[,"Gruppe"])[-(which(unique(ergebnis[,"Gruppe"])==0))])
  V2 <- seq_along(V1) 
  for (i in 1:length(V1)){
    ergebnis[which(ergebnis[,"Gruppe"]==V1[i]),"Gruppe"] <- V2[i]}
  
  
  ergebnis <- ergebnis[!is.na(ergebnis[,1]), ]
  return(ergebnis)
}



#' blankCorrection_variable_headerList
#' 
#' Blankabzug nach Vorbild der orginalen Non-target App Funktion mit Berücksichtigung der Gruppenzugehörigkeit
#'  
#' Ist ein Peak im Blank wird IN JEDER Probe je nach seiner Intensität in der Gruppe entschieden:
#'   Int_max der Gruppe: -> die ganze Gruppe wird in dieser Probe entfernt
#'   nicht Int_max der Gruppe: -> nur dieser Peak wird in der Probe entfernt
#'
#' @param grouptable Alignmenttabelle der grouping_BfG2 Funktion
#' @param head headerlist aus der Non-target App
#' @param intensityFactor Faktor für die Intensitätsdifferenz zwischen 2 Peaks (Vgl. grouping_BfG2)
#'
#' @return Blankkorrigierte Alignmenttabelle
#' @export
#'
#' @examples
blankCorrection_variable_headerList <- function(grouptable, intensityFactor = 10,head=headerList){
  stop("Headerlist no longer exists, needs to be rewritten using sampleList")
  blanks <- NULL
  intensityrows <- grep("Int_", colnames(grouptable)) #Einbindung der Gruppenspalte erledigt
  grouprows <- grep("gruppe_", colnames(grouptable))
  ergebnis <- matrix(,ncol=ncol(grouptable),nrow=0)
  
  for (i in (1:length(head))) {
    if (head[[i]]$sampleType == "Blank") blanks <- c(blanks,i)
  }
  
  
  #Erstellt Matrix mit den intensivsten Peaks jeder Probe, jeder Gruppe in den "Unknown" Proben
  grouptable_Int_max <- matrix(1000000,ncol=(2*length(grouprows[-blanks])),nrow=max(grouptable[,grouprows[1:(length(grouprows[-blanks]))]]))
  
  for (i in 1: (length(grouprows[-blanks]))){
    
    for (j in 1:max(grouptable[,grouprows[i]])){
      
      if (0!=length(which(grouptable[,grouprows[i]]==j))){
        grouptable_Int_max[j,((2*i)-1)] <- max(c(grouptable[which(grouptable[,grouprows[i]]==j),intensityrows[i]]))}
      grouptable_Int_max[j,(2*i)] <- j
      
    } #for j
  } #for i
  
  
  #Absuchen der grouptable nach Blanks
  for (i in (1:nrow(grouptable))) {
    #print(paste0(round((i/nrow(grouptable)*100),2),"% Blankabzug"))
    
    if (any(grouptable[i,intensityrows] > intensityFactor*max(grouptable[i,intensityrows[blanks]]))) {
      #grouptable[i,intensityrows] <- grouptable[i,intensityrows]-max(grouptable[i,intensityrows[blanks]])
      ergebnis <- rbind(ergebnis,grouptable[i,])
      
    } else {
      
      
      # Reihe ist ein Vector, welche anzeigt welcher Peak in den Proben 
      #      1.) im Blank gefunden wurde und
      #      2.) Gruppenintensivster Peak ist
      #      mit grouptable[i,[intensityrows[Reihe]]] bekommt man sie angezeigt
      
      Probe <- 1
      for (Probe in Probe:length(intensityrows[-blanks])){    #Zähler, der durch die Anzahl an Proben(Col in der Grouptable)(ohne blank geht)
        
        Reihe <- which(grouptable[i,intensityrows[Probe]] == grouptable_Int_max[,(2*Probe)-1])    #Ist der Peak aus Probe auch ein Gruppenchef
        # Reihe ist die Gruppe des Peaks in Probe Probe
        
        if (length(Reihe)>0){
          
          t<-i   #Laufvariable t bekommt den Wert von i, da nicht die ganze grouptable durchsucht werden muss, sondern nur der noch kommende Teil
          for (t in t:nrow(grouptable)){
            if (grouptable[t,grouprows[Probe]]==Reihe){grouptable[t,c((intensityrows[Probe]-3):(intensityrows[Probe]+2))] <- 0}}
          
          f<-1
          if (length(ergebnis)>0){
            for (f in f:nrow(ergebnis)){
              if (ergebnis[f,grouprows[Probe]]==Reihe){ergebnis[f,c((intensityrows[Probe]-3):(intensityrows[Probe]+2))] <- 0}}}
          
          
        } #if length(Reihe)
        
      } #for Probe
      
      #Sind die betroffenen Peaks die intensivsten in den betroffenen Gruppen?
      # der logische Vektor gr_max gibt an, ob der untersuchte Peak der intensivste der Gruppe ist
      
      #Was passiert wenn der betroffene Peak das Int-max der Gruppe ist?
      # Wenn gr_max TRUE ist wird die ganze Gruppe entfernt (peakID,mz,RT,Int,ms2scan, Gruppe <-0)
      # Wenn gr_max FALSE passiert an dieser Stelle nichts und der Peak wird wie in Christians Orginal einfach nicht in
      # die Ergebnisliste übertragen
      
    } #else
    
    
  } # for
  
  if (0<length(which(rowSums(ergebnis[,-c(1,2,3)])==0))){
    ergebnis <- ergebnis[-c(which(rowSums(ergebnis[,-c(1,2,3)])==0)),]}
  
  # Eliminieren von ein-Peak-Gruppen  (Gruppe wird 0)
  for (i in 1:max(ergebnis[,"Gruppe"])){
    if (length(which(ergebnis[,"Gruppe"]==i))==1) ergebnis[which(ergebnis[,"Gruppe"]==i),"Gruppe"] <-0}
  
  # Schließen der entstandenen Zwischenräume
  V1 <- sort(unique(ergebnis[,"Gruppe"])[-(which(unique(ergebnis[,"Gruppe"])==0))])
  V2 <- seq_along(V1) 
  for (i in 1:length(V1)){
    ergebnis[which(ergebnis[,"Gruppe"]==V1[i]),"Gruppe"] <- V2[i]}
  
  
  return(ergebnis)   
}



#' Einleiter
#' 
#' Funktion zur Ermittlung "neuer" Peaks einer Probe im Vergleich zur vorherigen Probe
#' 
#' Hinweis: wenn man das Ergebnis noch alignen will sollte man die Werte für ppm und DeltaRT in beiden Funktionen angleichen,
#'           da sonst u.U. Werte alignt werden, welche laut Einleiter eigentlich entfernt sein sollten!
#' 2. Hinweis: Es werden NUR aufeinander folgende Proben verglichen!
#'             Wenn ein Peak in der vorherigen Probe (im Code: i-1) nicht vorkommt, aber in 2 Proben davor (sinngemäß i-2) schon, 
#'             wird er als "neu" in die Ergebnisliste aufgenommen
#'             Alternative -> Funktion Einleiter_alles_vorher
#'
#' @param peaklist_e Peaklist aus der Non-target App oder auch aus Gruppierung_SG
#' @param Proben_Auswahl_Vector Vektor welche Proben untersucht werden sollen
#' @param ppm erlaubte Massendifferenz zwischen 2 Peaks
#' @param DeltaRT erlaubte RT-differenz zwischen 2 Peaks
#' @param Int_Lim Intensitätsthreshold ab wann Peaks untersucht werden sollen
#'                Peaks unterhalb des Thresholds kommen nicht in die Ergebnisliste
#' @param intensityFactor Faktor für die Intensitätsdifferenz zwischen 2 Peaks (Vgl. grouping_BfG2)
#'
#' @return Peaklist mit der Information welche Peaks zwischen 2 Proben unterschiedlich sind, Einleiter halt
#' @export
#'
#' @examples
Einleiter <- function(peaklist_e = new_peaklist, Proben_Auswahl_Vector = c(1:length(new_peaklist)),ppm=10, DeltaRT=10, Int_Lim=10, intensityFactor =10){
  
  ProbenList <- peaklist_e[Proben_Auswahl_Vector]
  ergebnis <- vector("list", length(Proben_Auswahl_Vector))
  ergebnis[1] <- ProbenList[1]
  
  i<-2
  for (i in i:length(ProbenList)){
    #print(100*i/length(Proben_Auswahl_Vector))
    j<-1
    
    for (j in j:nrow(ProbenList[[i]])){
      
      
      if (ProbenList[[i]][j,3] > Int_Lim){
        if(is.na(ProbenList[[i-1]][j,3])){ProbenList[[i-1]][j,3] <-0}
        Reihe <- which((abs(ProbenList[[i]][j,1]-ProbenList[[i-1]][,1]) < (ProbenList[[i]][j,1]*ppm/1000000) & (abs(ProbenList[[i]][j,2]-ProbenList[[i-1]][,2])<DeltaRT)))
        if (0 == length(Reihe)
            || ProbenList[[i]][j,3] > ProbenList[[i-1]][Reihe,3]*intensityFactor
        ){
          
          ergebnis[[i]] <- rbind(ergebnis[[i]],ProbenList[[i]][j,])
          
        }
      }
    }
  }
  
  
  return(ergebnis)
  
  #dann die Probe von der vorherigen abziehen um zu ermiteln wieviele Peaks dazugekommen sind (nrow)
  
} 



#' Einleiter_alles_vorher
#' 
#' Funktion zur Ermittlung "neuer" Peaks einer Probe im Vergleich zu ALLEN vorherigen Proben in der Peaklist
#' Hinweis: wenn man das Ergebnis noch alignen will sollte man die Werte für ppm und DeltaRT in beiden Funktionen angleichen,
#'          da sonst u.U. Werte alignt werden, welche laut Einleiter eigentlich entfernt sein sollten!
#'          
#' kann bei vielen Proben zu leeren Listen führen, wenn alle Peaks einer Probe bereits in vorherigen Proben gemessen wurden
#'  -> führt auf die Weise zum Absturz der Alignmentfunktion
#'
#' @param peaklist_e Peaklist aus der Non-target App oder auch aus Gruppierung_SG
#' @param Proben_Auswahl_Vector Vektor welche Proben untersucht werden sollen
#' @param ppm erlaubte Massendifferenz zwischen 2 Peaks
#' @param DeltaRT erlaubte RT-differenz zwischen 2 Peaks
#' @param Int_Lim Intensitätsthreshold ab wann Peaks untersucht werden sollen
#'                Peaks unterhalb des Thresholds kommen nicht in die Ergebnisliste
#' @param intensityFactor Faktor für die Intensitätsdifferenz zwischen 2 Peaks (Vgl. grouping_BfG2)
#'
#' @return Peaklist mit der Information welche Peaks zwischen einer Proben und allen vorherigen unterschiedlich sind
#' @export
#'
#' @examples
Einleiter_alles_vorher <- function(peaklist_e = new_peaklist, Proben_Auswahl_Vector = c(1:length(new_peaklist)),ppm=10, DeltaRT=10, Int_Lim=10, intensityFactor =10){
  
  ProbenList <- peaklist_e[Proben_Auswahl_Vector]
  ergebnis <- vector("list", length(Proben_Auswahl_Vector))
  ergebnis[1] <- ProbenList[1]
  
  u<-0
  
  for (i in 2:length(ProbenList)){

    for (j in 1:nrow(ProbenList[[i]])){
      
      
      if (ProbenList[[i]][j,3] > Int_Lim){
        x<-1  # x sind alle Proben vor der Probe i
        u <-0 # wenn u 0 bleibt ist der Peak in i in den anderen Proben nicht, oder intensityFactor niedriger enthalten und geht in das Ergebnis über. Wenn u 1 wird, nicht.
        for (x in x:(i-1)){
          
          if(is.na(ProbenList[[x]][j,3])){ProbenList[[x]][j,3] <-0}
          Reihe <- which((abs(ProbenList[[i]][j,1]-ProbenList[[x]][,1]) < (ProbenList[[i]][j,1]*ppm/1000000) & (abs(ProbenList[[i]][j,2]-ProbenList[[x]][,2])<DeltaRT)))
          
          if (0 == length(Reihe)
              || ProbenList[[i]][j,3] > ProbenList[[x]][Reihe,3]*intensityFactor)
          {u <- u} else {u <- u+1}
        }
        if (u==0){
          ergebnis[[i]] <- rbind(ergebnis[[i]],ProbenList[[i]][j,])}
        
      }
    }
  }
  
  
  return(ergebnis)
  
  #dann die Probe von der vorherigen abziehen um zu ermiteln wieviele Peaks dazugekommen sind (nrow)
  
} 



#' get Trends
#' 
#' Trenderkennung mit Angabe der Start und Endpunkte (Proben) des Intensitätsanstiegs
#'
#' @param grouptable      von Alignment oder Blankabzug
#' @param Int_Columns     Int Spalten der grouptable
#' @param sn              Signal/ Noise minimum
#' @param Int_Threshold   min Int der Maxima
#'
#' @return In welchen Proben haben Peaks ihre Intensitätsmaxima und in welchen
#'         dazugehörigen Proben beginnt der Anstieg der Intensität
#' @export
#'
#' @examples
getTrends <- function(grouptable, Int_Columns=grep("Int_", colnames(grouptable)),sn=5,Int_Threshold=20) {
  #groupedColums sind die Int der Alignmenttabelle
  
  #trends <- grouped[,seq(12,ncol(grouped),by = 5)] 
  trends <- grouptable[,Int_Columns] 
  #sn <- 5 
  
  resultTable <- NULL
  
  for (i in 1:nrow(trends)) {
    singleTrend <- cbind(seq(1:ncol(trends)),trends[i,])
    singleTrend[singleTrend[,2]==0,2] <- 1
    derivative <- cbind(singleTrend[2:nrow(singleTrend),1], diff(singleTrend[,2]))
    maxima <- derivative[(derivative[(2:(nrow(derivative)-1)-1),2] > 0) & (derivative [(2:(nrow(derivative)-1)),2] <= 0),1]
    maxima_start <- maxima
    intensity <- singleTrend[maxima,2]
    
    if (length(maxima) > 0) {
      for (n in 1:length(maxima)) {
        
        maxima_start[n] <- max(singleTrend[which(derivative[1:(maxima[n]-1),2] <= 0),1])+1
        if (is.infinite(maxima_start[n])) {maxima_start[n] <- singleTrend[1,1]}
        
        if (singleTrend[maxima[n],2] < sn*singleTrend[maxima_start[n],2]) maxima[n] <- 0
        
      }
    }
    #Christians Version
    maxima_start <- maxima_start[maxima[]>0]
    intensity <- intensity[maxima[]>0]
    maxima <- maxima[maxima[]>0]
    
    #if (length(maxima) > 0) resultTable <- rbind(resultTable,cbind(i,maxima_start,maxima,intensity))
    
    #meine erste Version mit Gruppenangabe
    #    Row <- i                                                  # Reihe in der Alignmenttabelle
    #    Sample_maxima_start <- maxima_start[maxima[]>0]           # Anfang des Int_Anstiege
    #    Int_maxima_start <- grouptable[i,(1+(6*maxima_start))]    # Int des Anfangs des Int_Anstieges, in den meisten Fällen 0
    #    Group_maxima_start <- grouptable[i,(3+(6*maxima_start))]  # Gruppe des Anfangs des Int_Anstieges, in den meisten Fällen 0
    #    Sample_maxima <- maxima[maxima[]>0]                       # Probe des Int_max
    #    Int_maxima <- intensity[maxima[]>0]                       # Int_max
    #    Group_maxima <- grouptable[i,(3+(6*maxima))]              # Gruppe des Int_max
    #    
    #    if (length(maxima) > 0) resultTable <- rbind(resultTable,cbind(Row,Sample_maxima_start,Int_maxima_start,Group_maxima_start,Sample_maxima,Int_maxima, Group_maxima))
    
    #meine zweite Version mit Gruppenangabe (ohne Werte die eh 0 sind)
    Row <- i                                                  # Reihe in der Alignmenttabelle
    mean_mz <- grouptable[[i,"mean_mz"]]   
    mean_RT <- grouptable[[i,"mean_RT"]]  
    Sample_maxima_start <- maxima_start[maxima[]>0]           # Anfang des Int_Anstiege
    Sample_maxima <- maxima[maxima[]>0]                       # Probe des Int_max
    Int_maxima <- intensity[maxima[]>0]                       # Int_max
    #Group_maxima <- grouptable[i,(4+(6*maxima))]              # Gruppe des Int_max
    Group_maxima <- grouptable[i,sprintf("gruppe_%i",Sample_maxima)]

    if (length(maxima) > 0) resultTable <- rbind(resultTable,cbind(Row,mean_mz,mean_RT,Sample_maxima_start,Sample_maxima,Int_maxima, Group_maxima))
    
  }
  
  for (i in nrow(resultTable):1){
    if(resultTable[i,"Int_maxima"]<Int_Threshold){resultTable <- resultTable[-i,]}
  }
  
  
  return(resultTable)
}

#' create_log_file
#'
#' @param no parameters needed
#'
#' @return log_file with used settings for peakpicking, alignment, blankcorrectio,normalization, filtering, annotation
#' @export
#'
#' @examples
create_log_file <- function () {
  
  # browser()
  
  if(!exists("log_file")){ 
    log_file <- list()
    log_file[6][[1]] <- list()}
  
  # Peakpicking Settings="Unknown"
  # Funktionsaufruf in app.R unter PP single sample und PP multiple sample 
  if(any(sampleList$sampleType=="Unknown") & 
     any(ls(envir=parent.frame())=="adjust_tolerance")) {
      sample <- which(sampleList$sampleType=="Unknown")
      Settings_PP_sample <- list()
      Settings_PP_sample$Peakpicking_Settings <- "Unknown"
      Settings_PP_sample$Mass_Range <- peakPickSettings[[sample[1]]]$massrange
      Settings_PP_sample$mz_Step <- peakPickSettings[[sample[1]]]$mz_step
      Settings_PP_sample$RT_Range_min <- peakPickSettings[[sample[1]]]$rtrange/60
      Settings_PP_sample$Peakwidth_s <- peakPickSettings[[sample[1]]]$peakwidth
      Settings_PP_sample$Noise_Scans <- peakPickSettings[[sample[1]]]$NoiseScans
      Settings_PP_sample$SN <- peakPickSettings[[sample[1]]]$sn
      Settings_PP_sample$Int_Threshold <- peakPickSettings[[sample[1]]]$int_threshold
      Settings_PP_sample$MaxNumPeaks <- peakPickSettings[[sample[1]]]$maxNumPeaks
      
      adjust_tolerance <- get ("adjust_tolerance", parent.frame())
      if (adjust_tolerance) {Settings_PP_sample$Dynamic_Tolerance_componentization <- TRUE
      } else {Settings_PP_sample$Dynamic_Tolerance_componentization <- FALSE
              Settings_PP_sample$mz_margin_Componentization_ppm <- peakPickSettings[[sample[1]]]$ppm
              Settings_PP_sample$RT_Tolerance_Componentization <- peakPickSettings[[sample[1]]]$RT_Tol}
      
      log_file[[1]] <- Settings_PP_sample
  }
  
  # Peakpicking Settings="Blank"
  # Funktionsaufruf in app.R unter PP single sample und PP multiple sample
  if(any(sampleList$sampleType=="Blank") &
     any(ls(envir=parent.frame())=="adjust_tolerance")){
    blanksample <- which(sampleList$sampleType=="Blank")
    Settings_PP_blank <- list()
    Settings_PP_blank$Peakpicking_Settings <- "Blank"
    Settings_PP_blank$Mass_Range <- peakPickSettings[[blanksample[1]]]$massrange
    Settings_PP_blank$mz_Step <- peakPickSettings[[blanksample[1]]]$mz_step
    Settings_PP_blank$RT_Range_min <- peakPickSettings[[blanksample[1]]]$rtrange
    Settings_PP_blank$Peakwidth_s <- peakPickSettings[[blanksample[1]]]$peakwidth
    Settings_PP_blank$NoiseScans <- peakPickSettings[[blanksample[1]]]$NoiseScans
    Settings_PP_blank$SN <- peakPickSettings[[blanksample[1]]]$sn
    Settings_PP_blank$Int_Threshold <- peakPickSettings[[blanksample[1]]]$int_threshold
    Settings_PP_blank$MaxNumPeaks <- peakPickSettings[[blanksample[1]]]$maxNumPeaks
    
      adjust_tolerance <- get ("adjust_tolerance", parent.frame())
      if (adjust_tolerance) {Settings_PP_blank$Dynamic_Tolerance_componentization <- TRUE
       } else {Settings_PP_blank$Dynamic_Tolerance_componentization <- FALSE
               Settings_PP_blank$mz_margin_Componentization_ppm <- peakPickSettings[[blanksample[1]]]$ppm
               Settings_PP_blank$RT_Tolerance_Componentization <- peakPickSettings[[blanksample[1]]]$RT_Tol}
    
    log_file[[2]] <- Settings_PP_blank
  }
  
  # Alignment Settings
  # Funktionsaufruf in app.R unter Alignment
  if(any(ls(envir=parent.frame())=="ppm_dev") &
     any(ls(envir=parent.frame())=="DeltaRT")){
    ppm_dev <- get ("ppm_dev", parent.frame())
    DeltaRT <- get ("DeltaRT", parent.frame())
    log_file[[3]] <- list(Alignment_Settings="",
                          Mass_Deviation_ppm=ppm_dev,
                          RT_Deviation_s=DeltaRT)
    # bei erneutem Alignment wird eine neue Alignmenttabelle erstellt,daher...
    if(length(log_file)>3){ log_file[[4]] <- list()
                            log_file[[5]] <- list()
                            log_file[[6]] <- list()}
    if(length(log_file)>6) log_file[[7]][[1]] <- "Durch neue Datenbearnbeitung ist die Annotation evtl. nicht mehr aktuell und sollte wiederholt werden"
    
    }
  
  # Blankcorrection Settings
  # Funktionsaufruf in app.R unter Blank Correction
  if(any(ls(envir=parent.frame())=="intensityFactor") &
     any(ls(envir=parent.frame())=="deleteGrouped")){
    intensityFactor <- get ("intensityFactor", parent.frame())
    deleteGrouped <- get ("deleteGrouped", parent.frame())
    grouped_before <- get ("nrow_grouped_before_blankcorrection", parent.frame())
    number_of_removed_features <- grouped_before - nrow(grouped)
    log_file[[4]] <- list(Blankcorrection_Settings="",
                          IntensityFactor=intensityFactor,
                          DeleteGrouped=deleteGrouped,
                          NumRemoveFeat=number_of_removed_features)
    # bei erneuter Blankcorrection sind evtl. vorhandene Filterschritte überspielt, daher...
    if(length(log_file)>4){ log_file[[5]] <- list()
                            log_file[[6]] <- list()}
    if(length(log_file)>6)  log_file[[7]][[1]] <- "Durch neue Datenbearnbeitung ist die Annotation evtl. nicht mehr aktuell und sollte wiederholt werden"
    }
  
  # Normalize Settings 
  # Funktionsaufruf in app.R unter Normalize
  if(any(ls(envir=parent.frame())=="intCols_normalize")){
    mean_mz_chosen_feature <- get ("mean_mz_chosen_feature", parent.frame())
    mean_RT_chosen_feature <- get ("mean_RT_chosen_feature", parent.frame())/60
    row_number <- get ("row_number", parent.frame())
      log_file[[5]] <- list(normalization_settings="",
                           Mean_mz_chosen_Feature=mean_mz_chosen_feature,
                           Mean_RT_chosen_Feature_min=mean_RT_chosen_feature,
                           RowNumber_chosen_Feature=row_number)
  } 
  
  # Filter Settings
  # Remove rows with fewer than X rows
      if(any(ls(envir=parent.frame())=="minDetections")){
        minDetections <- get ("minDetections", parent.frame())
        grouped_before <- get ("nrow_grouped_before_filter", parent.frame())
        number_of_removed_features <- grouped_before - nrow(grouped)
        log_file[[6]][[1]] <- list(Filter_Settings="Remove rows with fewer than X detections",
                                   MinSampleDetect=minDetections,
                                   NumRemoveFeat=number_of_removed_features)
        if(length(log_file)>6) log_file[[7]][[1]] <- "Durch neue Datenbearnbeitung ist die Annotation evtl. nicht mehr aktuell und sollte wiederholt werden"
        }    
  
  # Remove features which are not the highest feature in their component in files x:y
      if(any(ls(envir=parent.frame())=="leader_first_sample") &
        any(ls(envir=parent.frame())=="leader_last_sample")){
        leader_first_sample <- get ("leader_first_sample", parent.frame())
        leader_last_sample <- get ("leader_last_sample", parent.frame())
        grouped_before <- get ("nrow_grouped_before_filter", parent.frame())
          number_of_removed_features <- grouped_before - nrow(grouped)
          log_file[[6]][[2]] <- list(Filter_Settings="Remove features which are not the highest feature in their component",
                                     First_sample=leader_first_sample,
                                     Last_sample=leader_last_sample,
                                     NumRemoveFeat=number_of_removed_features)
      if(length(log_file)>6) log_file[[7]][[1]] <- "Durch neue Datenbearnbeitung ist die Annotation evtl. nicht mehr aktuell und sollte wiederholt werden"
      }
  
  # Remove features which are not found in replicate injections
  if(any(ls(envir=parent.frame())=="rep_first_sample") &
     any(ls(envir=parent.frame())=="rep_last_sample") &
     any(ls(envir=parent.frame())=="rep_number_replicates") &
     any(ls(envir=parent.frame())=="in_at_least_samples")){
    rep_first_sample <- get ("rep_first_sample", parent.frame())
    rep_last_sample <- get ("rep_last_sample", parent.frame())
    number_replicates <- get ("rep_number_replicates", parent.frame())
    in_at_least_samples <- get ("in_at_least_samples", parent.frame())
    grouped_before <- get ("nrow_grouped_before_filter", parent.frame())
    number_of_removed_features <- grouped_before - nrow(grouped)
    log_file[[6]][[3]] <- list(Filter_Settings="Remove features which are not found in replicate injections",
                               Replicate_first_sample=rep_first_sample,
                               Replicate_last_sample=rep_last_sample,
                               Number_Replicates=number_replicates,
                               In_at_least_samples=in_at_least_samples,
                               NumRemoveFeat=number_of_removed_features)
    if(length(log_file)>6) log_file[[7]][[1]] <- "Durch neue Datenbearnbeitung ist die Annotation evtl. nicht mehr aktuell und sollte wiederholt werden"
      }
  
  # Get average intensity of features in replicates
  if(any(ls(envir=parent.frame())=="repAve_first_sample") &
     any(ls(envir=parent.frame())=="repAve_last_sample") &
     any(ls(envir=parent.frame())=="repAve_number_replicates")){
    repAve_first_sample <- get ("repAve_first_sample", parent.frame())
    repAve_last_sample <- get ("repAve_last_sample", parent.frame())
    number_replicates <- get ("repAve_number_replicates", parent.frame())
    log_file[[6]][[4]] <- list(Filter_Settings="Get average intensity of features in replicates",
                               RepAve_first_sample=repAve_first_sample,
                               RepAve_last_sample=repAve_last_sample,
                               Number_Replicates=number_replicates)
    if(length(log_file)>6) log_file[[7]][[1]] <- "Durch neue Datenbearnbeitung ist die Annotation evtl. nicht mehr aktuell und sollte wiederholt werden"
      }
  
  # Remove rows where features are found in fewer than X consecutive files
  #browser()
  if(any(ls(envir=parent.frame())=="consec_first_sample") &
     any(ls(envir=parent.frame())=="consec_last_sample") &
     any(ls(envir=parent.frame())=="number_consecutives")){
    consec_first_sample <- get ("consec_first_sample", parent.frame())
    consec_last_sample <- get ("consec_last_sample", parent.frame())
    number_consecutives <- get ("number_consecutives", parent.frame())
    grouped_before <- get ("nrow_grouped_before_filter", parent.frame())
    number_of_removed_features <- grouped_before - nrow(grouped)
    log_file[[6]][[5]] <- list(Filter_Settings="Remove rows where features are found in fewer than X consecutive files",
                               Consec_first_sample=consec_first_sample,
                               Consec_last_sample=consec_last_sample,
                               Number_consecutives=number_consecutives,
                               NumRemoveFeat=number_of_removed_features)
    if(length(log_file)>6) log_file[[7]][[1]] <- "Durch neue Datenbearnbeitung ist die Annotation evtl. nicht mehr aktuell und sollte wiederholt werden"
      }
  
  # Annotation Settings
  # Funktionsaufruf in app.R unter Annotation
  if(any(ls(envir=parent.frame())=="threshold_score") &
     any(ls(envir=parent.frame())=="mztolu") &
     any(ls(envir=parent.frame())=="rttol") &
     any(ls(envir=parent.frame())=="polarity") &
     any(ls(envir=parent.frame())=="CE") &
     any(ls(envir=parent.frame())=="CES") &
     any(ls(envir=parent.frame())=="mztolu_ms2") &
     any(ls(envir=parent.frame())=="rtoffset")){
    threshold_score <- get ("threshold_score", parent.frame())
    mztolu <- get ("mztolu", parent.frame())
    rttol <- get ("rttol", parent.frame())
    polarity <- get ("polarity", parent.frame())
    CE <- get ("CE", parent.frame())
    CES <- get ("CES", parent.frame())
    mztolu_ms2 <- get ("mztolu_ms2", parent.frame())
    rtoffset <- get ("rtoffset", parent.frame())
    groupedWithAnnotation <- get ("groupedWithAnnotation", parent.frame())
    log_file[[7]] <- list(Annotation_Settings="",
                          Threshold_Score=threshold_score,
                          mz_Tolerance_mDa=mztolu*1000,
                          RT_Tolerance_min=rttol,
                          Polarity=polarity,
                          CE=CE,
                          CES=CES,
                          mz_Tolerance_MS2_mDa=mztolu_ms2*1000,
                          RT_Offset_min=rtoffset,
                          NumAnnoCompounds=length(unique(groupedWithAnnotation$name))-1)
  }
  
  return(log_file)
}


