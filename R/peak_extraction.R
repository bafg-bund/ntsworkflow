
# Peak picking for database screening ####



#' Christian's FindPeaks algorithm
#'
#' @param i
#' @param daten
#' @param mz_step
#' @param rt_min_scan
#' @param rt_max_scan
#' @param sn
#' @param int_threshold
#' @param peak_NoiseScans
#' @param precursormzTol
#'
#' @return
#' @export
#'
#' @examples
FindPeaks_SingleXIC2 <- function(i,
                                 daten,
                                 mz_step,
                                 rt_min_scan,
                                 rt_max_scan,
                                 sn,
                                 int_threshold,
                                 peak_NoiseScans,
                                 precursormzTol) {

  peaklist_singleXIC <- matrix( , nrow = 0, ncol = 15)
  XIC <- xcms::rawEIC(daten, mzrange = c(i,i+mz_step))
  #browser()


  #maxPeaksPerSignal <- 60/(max(daten@scantime)/length(daten@scantime))/5
  maxPeaksPerSignal <- 20
  deleteZeroIntensityNoise <- FALSE

  #XIC mit Scan und Intensity in eine besser lesbare Tabelle umschreiben:
  non_derived <- cbind(XIC$scan,XIC$intensity)

  #if (median(non_derived[,2]) > 0) {
  #1st derivative of XIC:
  deriv1 <- cbind(XIC$scan[2:length(XIC$scan)], diff(XIC$intensity))

  #2nd derivative of XIC:
  #deriv2 <- cbind(XIC$scan[3:length(XIC$scan)], diff(diff(XIC$intensity)))
  #Scan-Nummern der Maxima aufschreiben: (Bedingung f?r Maxima: 2.Ableitung < 0, 1.Ableitung geht durch 0 von negativ (scan vorher < 0) zu positiv (dieser scan < 0); und Intensity > Threshold) - Es wird der vorhergehende Scan aufgeschrieben, da durch Ableitungsbildung (immer zwischen zwei Punkten) das Maximum "verrutscht"
  #maxima <- deriv2[(deriv2[1:(nrow(deriv2)-1),2]<0) & (deriv1[(2:(nrow(deriv1)-1)-1),2] > 0) & (deriv1[(2:(nrow(deriv1)-1)),2] < 0) & (non_derived[3:(nrow(non_derived)-1),2]>int_threshold),1]-1

  #Vereinfachung: 2.Ableitung wird nicht gebraucht, da Bedingung "1.Ableitung geht durch 0 von positiv (scan vorher > 0) zu negativ (dieser scan < 0)" dies schon implementiert
  #maxima <- deriv1[(deriv1[(2:(nrow(deriv1)-1)-1),2] > 0) & (deriv1[(2:(nrow(deriv1)-1)),2] <= 0) & (non_derived[3:(nrow(non_derived)-1),2]>int_threshold),1]
  maxima <- deriv1[(deriv1[(2:(nrow(deriv1)-1)-1),2] > 0) & (deriv1[(2:(nrow(deriv1)-1)),2] <= 0),1]
  #int_threshold2 <- quantile(non_derived[,2])[4]
  int_threshold2 <- quantile(non_derived[,2],probs=0.99)
  if ((int_threshold2 == 0) | (int_threshold2 > int_threshold)) {
    maxima <- maxima[non_derived[maxima,2] >= int_threshold]
  } else {
    maxima <- maxima[non_derived[maxima,2] >= int_threshold2]
  }

  #Nur die Maxima im gew?nschten RT-Bereich behalten. (Es wurde vorher absichtlich der gesamte RT-Bereich abgesucht, um Peaks nicht abzuschneiden - macht von der Rechenzeit nicht viel Unterschied)
  maxima <- maxima[(maxima > rt_min_scan) & (maxima < rt_max_scan)]
  left_end <- maxima
  right_end <- maxima
  noiselevel <- maxima
  #noiselevelafter <- maxima
  #noisedeviation <- maxima
  #noisedeviationafter <- maxima
  amountofpeaks <- rep(1,length(maxima))
  #einen Bereich des XIC-Spektrums um jeweiliges Maximum (+/- peak_NoiseScans), um noiselevel zu bestimmen
  #diejenigen Maxima unterhalb des sn (signal/noise)-Levels werden auf 0 gesetzt
  if (length(maxima) > 0) { #wird nur gemacht, wenn auch Maxima eingetragen wurden
    for (n in 1:length(maxima)) {
      if (maxima[n]<=peak_NoiseScans) {
        noise <- non_derived[1:(maxima[n]+peak_NoiseScans),]
        left_end[n] <- max(noise[which(diff(noise[1:maxima[n],2]) <= 0),1])+1
        right_end[n] <- min(noise[(which(diff(noise[(maxima[n]+1):(maxima[n]+peak_NoiseScans),2]) >= 0)+maxima[n]),1])
        if (is.infinite(left_end[n])) {left_end[n] <- noise[1,1]}
        if (is.infinite(right_end[n])) {right_end[n] <- noise[peak_NoiseScans*2+1,1]}
        #noise_afterpeak <- non_derived[(right_end[n]+1):(right_end[n]+peak_NoiseScans),]
        if (deleteZeroIntensityNoise) noise <- noise[noise[,2] != 0,,drop=FALSE]
        noise <- non_derived[1:(left_end[n]-1),,drop=FALSE]

      } else {
        if (maxima[n]>(length(non_derived[,1])-peak_NoiseScans)) {
          noise <- non_derived[(maxima[n]-peak_NoiseScans):length(non_derived[,1]),]
          left_end[n] <- max(noise[which(diff(noise[1:peak_NoiseScans,2]) <= 0),1])+1
          right_end[n] <- min(noise[(which(diff(noise[(peak_NoiseScans+1):nrow(noise),2]) >= 0)+peak_NoiseScans),1])
          if (is.infinite(left_end[n])) {left_end[n] <- noise[1,1]}
          if (is.infinite(right_end[n])) {right_end[n] <- noise[nrow(noise),1]}
          noise <- non_derived[(left_end[n]-peak_NoiseScans):(left_end[n]-1),,drop=FALSE]
          if (deleteZeroIntensityNoise) noise <- noise[noise[,2] != 0,,drop=FALSE]
        }
        else {
          noise <- non_derived[(maxima[n]-peak_NoiseScans):(maxima[n]+peak_NoiseScans),]


          #Anfang und Ende des Peaks ?ber 1.Ableitung bestimmen:
          #absolute Scan-Nummer:
          left_end[n] <- max(noise[which(diff(noise[1:peak_NoiseScans,2]) <= 0),1])+1
          right_end[n] <- min(noise[(which(diff(noise[(peak_NoiseScans+2):(peak_NoiseScans*2+1),2]) >= 0)+peak_NoiseScans+1),1])

          if (is.infinite(left_end[n])) {left_end[n] <- noise[1,1]}
          if (is.infinite(right_end[n])) {right_end[n] <- noise[peak_NoiseScans*2+1,1]}


          if (left_end[n] > peak_NoiseScans) {
            noise <- non_derived[(left_end[n]-peak_NoiseScans):(left_end[n]-1),,drop=FALSE]
          } else {
            noise <- non_derived[1:(left_end[n]-1),,drop=FALSE]
          }

          #if ((right_end[n]+peak_NoiseScans) > length(non_derived[,1])) {
          #  noise_afterpeak <- non_derived[(right_end[n]+1):length(non_derived[,1]),]
          #} else {
          #  noise_afterpeak <- non_derived[(right_end[n]+1):(right_end[n]+peak_NoiseScans),]
          #}

          #Werte mit Intensity=0 rauswerfen:
          if (deleteZeroIntensityNoise) noise <- noise[noise[,2] != 0,,drop=FALSE]
          #noise_afterpeak <- noise_afterpeak[noise_afterpeak[,2] != 0,,drop=FALSE]
        }
      }


      noiselevel[n] <- mean(noise[,2])
      noiselevel[is.na(noiselevel)] <- mean(non_derived[,2])
      #USP EP definition: !!!
      #noisedeviation[n] <- (max(noise[,2])-min(noise[,2]))/2

      #noiselevelafter[n] <- mean(noise_afterpeak[,2])
      #noiselevelafter[is.na(noiselevelafter)] <- mean(non_derived[,2])
      #noisedeviationafter[n] <- (max(noise_afterpeak[,2])-min(noise_afterpeak[,2]))/2

      #?berpr?fen, ob zwei Peaks zusammen geh?ren (right_end des vorherigen und left_end dieses Maximums gleich oder überlappend)
      if ((n > 1) & (maxima[n] > 0)){
        if ((right_end[n-1] >= (left_end[n])) & (maxima[n-1] > 0))  {
          #if (non_derived[left_end[n],2] > noiselevel[n-1]) {
          #if ((non_derived[left_end[n],2] > non_derived[maxima[n],2]/2) | (non_derived[left_end[n],2] > non_derived[maxima[n-1],2]/2)) {
          if ((non_derived[left_end[n],2]-noiselevel[n-1]) > (min(c(non_derived[maxima[n],2],non_derived[maxima[n-1],2]))-noiselevel[n-1])/2) {
            #if (((non_derived[left_end[n],2] > non_derived[maxima[n],2]/2) | (non_derived[left_end[n],2] > non_derived[maxima[n-1],2]/2)) & ((min(c(non_derived[maxima[n],2],non_derived[maxima[n-1],2]))-noiselevel[n-1])>(max(c(non_derived[maxima[n],2],non_derived[maxima[n-1],2]))-noiselevel[n-1])/10)) {
            if (min(c(non_derived[maxima[n],2],non_derived[maxima[n-1],2])) < max(c(non_derived[maxima[n-1],2],non_derived[maxima[n],2]))/2) amountofpeaks[which.min(c(non_derived[maxima[n-1],2],non_derived[maxima[n],2]))+n-2] <- 0
            #in diesem Fall wird das Maximum dem intensiveren Peak zugeordnet:
            maxima[n] <- maxima[which.max(non_derived[c(maxima[n-1],maxima[n]),2])+n-2]
            amountofpeaks[n] <- amountofpeaks[n]+amountofpeaks[n-1]
            #das Noiselevel (links des Peaks) und die Position des linken Rands wird vom vorherigen Maximum ?bernommen:
            noiselevel[n] <- noiselevel[n-1]
            left_end[n] <- left_end[n-1]
            #Das vorherige Maximum auf 0 setzen, damit es sp?ter gel?scht wird:
            maxima[n-1] <- 0
            left_end[n-1] <- 0
            right_end[n-1] <- 0
            noiselevel[n-1] <- 0
            #noiselevelafter[n-1] <- 0
            #noisedeviation[n-1] <- 0
            #noisedeviationafter[n-1] <- 0
            amountofpeaks[n-1] <- 0
          }
          else {
            #noisedeviation[n] <- noisedeviation[n-1]
            #noisedeviationafter[n-1] <- noisedeviation[n-1]
            noiselevel[n] <- noiselevel[n-1]
            #noiselevelafter[n-1] <- noiselevel[n-1]
          }
        }
      }
    }
    #nur diejenigen Maxima behalten, die nicht auf 0 gesetzt wurden:
    left_end <- left_end[maxima[]>0]
    right_end <- right_end[maxima[]>0]
    noiselevel <- noiselevel[maxima[]>0]
    #noiselevelafter <- noiselevelafter[maxima[]>0]
    #noisedeviation <- noisedeviation[maxima[]>0]
    #noisedeviationafter <- noisedeviationafter[maxima[]>0]
    amountofpeaks <- amountofpeaks[maxima[]>0]
    maxima <- maxima[maxima[]>0]

    #und nun diejenigen rauswerfen, die unter int_threshold sind: (geht leider irgendwie nicht in einem Schritt)
    maxima[(non_derived[maxima,2] < int_threshold)] <- 0
    left_end <- left_end[maxima[]>0]
    right_end <- right_end[maxima[]>0]
    noiselevel <- noiselevel[maxima[]>0]
    #noiselevelafter <- noiselevelafter[maxima[]>0]
    #noisedeviation <- noisedeviation[maxima[]>0]
    #noisedeviationafter <- noisedeviationafter[maxima[]>0]
    amountofpeaks <- amountofpeaks[maxima[]>0]
    maxima <- maxima[maxima[]>0]
  }



  #check if there are peak clusters (most probably noise):
  if (length(maxima) > 1) {
    peaksectionstartid <- 1
    for (n in 2:length(maxima)) {
      #go through all maxima until 1) you find one, which is not directly neighboring the one before or 2) the larger peak is higher than 2 times the smaller one or 3) we reached the last maximum:
      if ((right_end[n-1] < left_end[n]) | (min(c(non_derived[maxima[n],2],non_derived[maxima[n-1],2])) < max(c(non_derived[maxima[n-1],2],non_derived[maxima[n],2]))/2) | (n==length(maxima))) {
        #check if the sum of all peaks from the selection before (n-1) is too high (one maxima might already consist of several peaks). Exception: one single signal with a lot of peaks is fine
        if ((sum(amountofpeaks[peaksectionstartid:(n-1)]) > maxPeaksPerSignal) & (peaksectionstartid != n)) {
          maxima[peaksectionstartid:(n-1)] <- 0
        }
        #next round start from here:
        peaksectionstartid <- n
      }
    }
    left_end <- left_end[maxima[]>0]
    right_end <- right_end[maxima[]>0]
    noiselevel <- noiselevel[maxima[]>0]
    #noiselevelafter <- noiselevelafter[maxima[]>0]
    amountofpeaks <- amountofpeaks[maxima[]>0]
    maxima <- maxima[maxima[]>0]
  }

  #now that we know, which peaks belong together and which ones are baseline, we cut out the ones which do not belong to the selecte RT range:
  #left_end <- left_end[(maxima > rt_min_scan) & (maxima < rt_max_scan)]
  #right_end <- right_end[(maxima > rt_min_scan) & (maxima < rt_max_scan)]
  #noiselevel <- noiselevel[(maxima > rt_min_scan) & (maxima < rt_max_scan)]
  #noiselevelafter <- noiselevelafter[(maxima > rt_min_scan) & (maxima < rt_max_scan)]
  #amountofpeaks <- amountofpeaks[(maxima > rt_min_scan) & (maxima < rt_max_scan)]
  #maxima <- maxima[(maxima > rt_min_scan) & (maxima < rt_max_scan)]

  noisedeviation <- maxima

  if (length(maxima) > 0) {
    for (n in 1:length(maxima)) {
      if ((right_end[n]+peak_NoiseScans) > nrow(non_derived)) {
        noise <- cbind(non_derived[(left_end[n]-peak_NoiseScans):(nrow(non_derived)),],0)
      } else {
        if ((left_end[n]-peak_NoiseScans) < 1) {
          noise <- cbind(non_derived[1:(right_end[n]+peak_NoiseScans),],0)
        } else {
          noise <- cbind(non_derived[(left_end[n]-peak_NoiseScans):(right_end[n]+peak_NoiseScans),],0)
        }
      }

      # TODO: no need to always go through all m maxima. Only check for those in the peak_NoiseScans - range
      noise <- noise[-((left_end[n]:right_end[n])-left_end[n]+peak_NoiseScans+1),]

      for (m in 1:length(maxima)) {
        #for (m in which((maxima < max(noise[,1])) & (maxima > min(noise[,1])))) {
        #for (m in which(((maxima < max(noise[,1])) & (maxima > min(noise[,1]))) | ((right_end < max(noise[,1])) & (right_end > min(noise[,1]))) | ((left_end < max(noise[,1])) & (left_end > min(noise[,1]))))) {
        noise[noise[,1] %in% non_derived[(left_end[m]:right_end[m]),1],3] <- 1
      }
      noise <- noise[noise[,3]==0,,drop=FALSE]
      if (deleteZeroIntensityNoise) noise <- noise[noise[,2] != 0,,drop=FALSE]
      if (nrow(noise) == 0) {noise <- non_derived[c(left_end[n],right_end[n]),]}

      noiselevel[n] <- mean(noise[,2])

      #USP EP definition:
      #noisedeviation[n] <- (max(noise[,2])-min(noise[,2]))/2
      noisedeviation[n] <- IQR(noise[,2])

      if (((non_derived[maxima[n],2]-noiselevel[n]) < noisedeviation[n]*sn) | ((non_derived[maxima[n],2]-noiselevel[n]) < int_threshold)) {
        maxima[n] <- 0
        left_end[n] <- 0
        right_end[n] <- 0
        noiselevel[n] <- 0
        #noiselevelafter[n] <- 0
        amountofpeaks[n] <- 0
        noisedeviation[n] <- 0
        #noisedeviationafter[n] <- 0
      }

    }

    left_end <- left_end[maxima[]>0]
    right_end <- right_end[maxima[]>0]
    noiselevel <- noiselevel[maxima[]>0]
    #noiselevelafter <- noiselevelafter[maxima[]>0]
    noisedeviation <- noisedeviation[maxima[]>0]
    #noisedeviationafter <- noisedeviationafter[maxima[]>0]
    amountofpeaks <- amountofpeaks[maxima[]>0]
    maxima <- maxima[maxima[]>0]
  }




  if (length(maxima) > 0) {

    #exactmasses <- NULL
    #peak_intens <- NULL
    exactmass <- maxima
    peak_intens <- maxima

    for (n in 1:length(maxima)) {
      #das Spektrum ?ffnen, um die exakte Masse und Intensit?t abzulesen:
      spektrum <- xcms::getScan(daten,maxima[n])
      spektrum <- spektrum[(spektrum[,1] >= i) & (spektrum[,1] <= i+mz_step),,drop=FALSE]
      if (length(spektrum)>2) { #bei mehreren m?glichen exakten Massen wird die intensivste aus der Liste genommen
        #exactmasses <- c(exactmasses,spektrum[spektrum[,2] == max(spektrum[,2]),1])
        #peak_intens <- c(peak_intens,spektrum[spektrum[,2] == max(spektrum[,2]),2]-noiselevel[n])
        exactmass[n] <- spektrum[spektrum[,2] == max(spektrum[,2]),1]
        peak_intens[n] <- spektrum[spektrum[,2] == max(spektrum[,2]),2]-noiselevel[n]
      }
      if (length(spektrum)==2) { #gibt es nur eine passende exakte Masse im Bereich, wird die Liste zum Vector:
        #exactmasses <- c(exactmasses,spektrum[1])
        #peak_intens <- c(peak_intens,spektrum[2]-noiselevel[n])
        exactmass[n] <- spektrum[1]
        peak_intens[n] <- spektrum[2]-noiselevel[n]
      }
      if (length(spektrum)==0) { #gibt es keine passende exakte Masse im Bereich, wird 0 vermerkt:
        #exactmasses <- c(exactmasses,0)
        #peak_intens <- c(peak_intens,0)
        exactmass[n] <- 0
        peak_intens[n] <- 0
      }



      #calculate FWHM:
      peakspektrum <- non_derived[right_end[n]:left_end[n],,drop=FALSE]
      #FWHM_left_scan <- min(peakspektrum[peakspektrum[,2] > (peak_intens[n]/2),1])
      FWHM_left_scan <- min(peakspektrum[(peakspektrum[,2]-noiselevel[n]) > (peak_intens[n]/2),1])
      if (FWHM_left_scan == 1) FWHM_left_scan <- 2
      #FWHM_right_scan <- max(peakspektrum[peakspektrum[,2] > (peak_intens[n]/2),1])
      FWHM_right_scan <- max(peakspektrum[(peakspektrum[,2]-noiselevel[n]) > (peak_intens[n]/2),1])
      if (FWHM_right_scan == nrow(non_derived)) FWHM_right_scan <- FWHM_right_scan-1
      Steigung <- (non_derived[FWHM_left_scan,2] - non_derived[FWHM_left_scan-1,2])/(daten@scantime[FWHM_left_scan]-daten@scantime[FWHM_left_scan-1])
      FWHM_left <- daten@scantime[FWHM_left_scan-1]+(peak_intens[n]/2-non_derived[FWHM_left_scan-1,2]+noiselevel[n])/Steigung
      Steigung <- (non_derived[FWHM_right_scan+1,2] - non_derived[FWHM_right_scan,2])/(daten@scantime[FWHM_right_scan+1]-daten@scantime[FWHM_right_scan])
      FWHM_right <- daten@scantime[FWHM_right_scan]+(peak_intens[n]/2-non_derived[FWHM_right_scan,2]+noiselevel[n])/Steigung

      #calculate area:
      peakArea <- flux::auc(x = daten@scantime[peakspektrum[,1]], y = peakspektrum[,2],thresh = noiselevel[n])

      #Get MS2 scan:
      ms2scan <- 0 #in case there is no MS2 we make sure that there will be an entry
      #check which MS2 scans were performed with a suitable parent mass:
      ms2candidates <- which(abs(daten@msnPrecursorMz-exactmass[n])<exactmass[n]/1000000*precursormzTol)
      #select the candidate, whose RT is closest to the one of the MS1 maximum:
      #if (length(ms2candidates) > 0) ms2scan <- ms2candidates[which(abs(daten@msnRt[ms2candidates]-daten@scantime[maxima[n]]) == min(abs(daten@msnRt[ms2candidates]-daten@scantime[maxima[n]])))]
      if (length(ms2candidates) > 0) ms2scan <- ms2candidates[which.min(abs(daten@msnRt[ms2candidates]-daten@scantime[maxima[n]]))]



      if ((peak_intens[n] >= int_threshold) & (peak_intens[n] >= noisedeviation[n]*sn)) peaklist_singleXIC <- rbind(peaklist_singleXIC,c(exactmass[n],daten@scantime[maxima[n]],peak_intens[n],maxima[n],daten@scantime[left_end[n]],daten@scantime[right_end[n]],left_end[n],right_end[n],noisedeviation[n],peakArea,FWHM_left,FWHM_right,noiselevel[n],i,ms2scan))
    }
  }
  colnames(peaklist_singleXIC) <- c("exactmass","scantime","peak_intens","maxima","scantimeleft_end","scantimeright_end","left_end","right_end","noisedeviation","peakArea","FWHM_left","FWHM_right","noiselevel","i","ms2scan")

  #} #end if median
  return(peaklist_singleXIC)
}

