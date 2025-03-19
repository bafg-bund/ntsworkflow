

#include <RcppArmadillo.h>
#include "peakPicking.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Find peaks in an ion chromatogram
//'
//' @description Peak finding algorithm using apex detection by 1st derivative, 
//' an iterative search method and no chromatogram smoothing. The method is 
//' published in Dietrich, C., Wick, A., & Ternes, T. A. (2021). 
//' Open source feature detection for non‐target LC‐MS analytics. 
//' Rapid Communications in Mass Spectrometry, e9206. https://doi.org/10.1002/rcm.9206 
//'
//' @param mz m/z of current ion chromatogram (Da)
//' @param mz_step binning width used to extract chromatogram (da)
//' @param eic extracted ion chromatogram (just intensities)
//' @param scantime Scan time of each index in eic (in s) 
//' @param minIntensity Minimum intensity for peak-picking
//' @param sn Minimum signal-to-noise ratio
//' @param noisescans Number of scans before and after peak to determine noise
//' @param peakwidth_min Minimum width of a peak
//' @param peakwidth_max Maximum width of a peak
//' @param maxPeaksPerSignal Maximum number of sub-peaks (direction changes) within a peak
//' 
//' @return A numeric matrix of class Rcpp::numericMatrix.
//'  rows: peaks found. cols: 16 peak descriptors.
//'  col 1: 0 (placeholder for m/z)
//'  col 2: Retention time of peak (s)
//'  col 3: 0 (Placeholder for peak intensity)
//'  col 4: Intensity found in chromatogram 
//'  col 5: Scan number of peak apex
//'  col 6: Scantime of peak start
//'  col 7: Scantime of peak end
//'  col 8: Scan number of peak start
//'  col 9: Scan number of peak end
//'  col 10: UNKNOWN noisedeviation
//'  col 11: Peak area
//'  col 12: Left RT of peak at half height (s)
//'  col 13: Right RT of peak at half height (s)
//'  col 14: Baseline of peak (intensity)
//'  col 15: m/z of this chromatogram
//'  col 16: 0 (placeholder for ms2 scan number)
//' @export
// [[Rcpp::export]]
NumericMatrix pickPeaksOneEicCpp(
    double mz, 
    double mz_step, 
    std::vector<double> eic,  
    std::vector<double> scantime, 
    double minIntensity, 
    int sn, 
    int noisescans, 
    double peakwidth_min, 
    double peakwidth_max,
    int maxPeaksPerSignal
  ) {
  
  
  
  double minIntensityChecked = checkMinIntensity(eic, minIntensity);
  std::vector<double> eicDerivative = computeDerivative(eic);
  std::vector<int> apexLocs = getApexIndices(eic, eicDerivative, minIntensityChecked);
  
  int anzahlmaxima = apexLocs.size();
  
  // Determine peak boundaries
  // --Loop through each local maximum to merge adjacent maxima--
  // 
  // If the intensity of the valley between the current maximum
  // and the previous one is lower than 1/2 of the lower of the two
  // maxima, then the current maximum is merged with the previous one
  // and the new maximum gets the intensity of the higher maximum.
  // The merging is done by setting the previous peak's location to
  // zero (this way it is marked for deletion later) 
  
  std::vector<int> startLoc(anzahlmaxima, 0);	 
  std::vector<int> endLoc(anzahlmaxima, 0);
  std::vector<double> noiselevel(anzahlmaxima, 0);
  std::vector<int> peaksInPeak(anzahlmaxima, 0);
  double intensity = 0;
  int noisecounter = 0;
  int indLook = 0;
  
  for (int apexNum = 0; apexNum < anzahlmaxima; ++apexNum) { 
    
    startLoc[apexNum] = getPeakStartLoc(apexLocs, eicDerivative, apexNum);
    endLoc[apexNum] = getPeakEndLoc(apexLocs, eicDerivative, apexNum);
    
    // 1st noiselevel calculation based on mean intensity in front of peak
    // Compute average intensity in the area before the peak by summing all
    // intensities in a loop and then dividing by the number of iterations
    intensity = 0;
    noisecounter = 0;  
    indLook = startLoc[apexNum];
    while ((indLook > (startLoc[apexNum]-noisescans)) && (indLook > 0)) {
      indLook--;
      noisecounter++;
      intensity += eic.at(indLook);
    }
    
    noiselevel[apexNum] = intensity/noisecounter;
    
    peaksInPeak[apexNum] = 1;
    // Merging of consecutive local maxima (Check if two peaks belong together) 
    // do not do this for the first peak or any peak at the start of the eic
    if ((apexNum > 0) && (apexLocs[apexNum] > 0)) {
      // determine if the current peak and the previous peak are next to each 
      // other (and we are not at the beginning)
      if ((endLoc[apexNum-1] >= startLoc[apexNum]) && (apexLocs[apexNum-1] > 0)) {
        // Determine if the valley between the current and previous 
        // peak is higher than half the height of the smaller peak
        // if this is the case the peaks belong together and will be joined
        double noisePeak1 = noiselevel[apexNum-1];
        double valleyBetween = eic[startLoc[apexNum]] - noisePeak1;
        double heightSmaller = std::min(eic[apexLocs[apexNum]], eic[apexLocs[apexNum-1]]) - noisePeak1;
        if (valleyBetween > heightSmaller / 2) {
          // peaks will be joined
          
          
          // Determine peaksInPeak
          // Determine if the smaller of the two peaks is higher than half the 
          // higher of the two peaks
          double heightHigher = std::max(eic[apexLocs[apexNum-1]],eic[apexLocs[apexNum]]) - noisePeak1;
          if (heightSmaller > heightHigher / 2) {
            // it is clear that it's part of same peak - so peaksInPeak is the sum from both peaks
            peaksInPeak[apexNum] = peaksInPeak[apexNum] + peaksInPeak[apexNum-1];
          } else {
            // if the previous peak is deemed to be noisy (too many peaksInPeak), then it is deleted
            // by setting the left end to the current left end
            if (peaksInPeak[apexNum-1] > maxPeaksPerSignal) {
              noiselevel[apexNum-1] = eic[apexLocs[apexNum-1]];
              startLoc[apexNum-1] = startLoc[apexNum];
            }
            
            // Here we look to see if the intensity has 'levelled off' by 
            // comparing to the next 2 peaks
            // this basically ends a tailing after it is no longer decreasing
            if ((eic[apexLocs[apexNum-1]] > eic[apexLocs[apexNum]]) && (anzahlmaxima > apexNum+2)) {
              std::vector<double> heights = {eic[apexLocs[apexNum]], eic[apexLocs[apexNum+1]], eic[apexLocs[apexNum+2]]};
              auto smallestElem = std::min_element(heights.begin(), heights.end());
              auto highestElem = std::max(heights.begin(), heights.end());
              if ( *smallestElem > *highestElem / 2 ) {
                // The intensity has levelled off, so the current peak marks the end of the current peak
                // delete the current peak by setting endLoc to previous peak's endLoc
                endLoc[apexNum] = endLoc[apexNum-1];
                peaksInPeak[apexNum] = peaksInPeak[apexNum-1];
              }
            }
          }
          
          // If the previous peak is the higher peak, that becomes the new apex
          if (eic[apexLocs[apexNum-1]] > eic[apexLocs[apexNum]]) {
            apexLocs[apexNum] = apexLocs[apexNum-1];
          }
          
          // The noise level and peak start become that of first peak
          noiselevel[apexNum] = noiselevel[apexNum-1];
          startLoc[apexNum] = startLoc[apexNum-1];
          // delete the first of the two peaks by setting values to 0
          apexLocs[apexNum-1] = 0;
          startLoc[apexNum-1] = 0;
          endLoc[apexNum-1] = 0;
          noiselevel[apexNum-1] = 0;
          peaksInPeak[apexNum-1] = 0;
        } else {
          // peaks are not joined
          // if the previous peak is clean, set the current peak's noise level to that of the previous peak
          // this way you get a better determination of noise (since there is no peak in the way)
          if (peaksInPeak[apexNum-1] < maxPeaksPerSignal) {
            noiselevel[apexNum] = noiselevel[apexNum-1];
          }
        }
      }
    }
  }

	 
	int j = 0;

	
	// delete maxima that have been set to 0, those below intensity threshold and too narrow ones
	j = 0;
	for (int i = 0; i < anzahlmaxima; ++i) { 
	  if ((apexLocs[i] > 0) && (eic[apexLocs[i]] > minIntensity) && 
         (scantime.at(endLoc[i])-scantime.at(startLoc[i]) > peakwidth_min)) {
	  	apexLocs[j] = apexLocs[i];
      startLoc[j] = startLoc[i];
      endLoc[j] = endLoc[i];
      noiselevel[j] = noiselevel[i];
      peaksInPeak[j] = peaksInPeak[i];
      j++;
	  }
	}
	anzahlmaxima = j;
	apexLocs.resize(anzahlmaxima); 
	
	// 190722 KJ
	// in some cases peaksInPeak can be larger then 10...
	// To correct for this set all back to 10
	for (int k = 0; k < anzahlmaxima; ++k) {
	  if (peaksInPeak.at(k) > maxPeaksPerSignal) {
	    peaksInPeak.at(k) = maxPeaksPerSignal;
	  }
	}
	
	int peaksectionstartid = 0;
	int peaksInPeak_sum = 0;
	
	// check if there are peak clusters
	for (int i = 1; i < anzahlmaxima; ++i) { 
		peaksInPeak_sum += peaksInPeak[i-1];
		if ((endLoc[i-1] < startLoc[i]) || 
          (std::min(eic[apexLocs[i]],eic[apexLocs[i-1]]) < 
            std::max(eic[apexLocs[i-1]],eic[apexLocs[i]])/2) || 
              ((apexLocs[i]-apexLocs[i-1]) > noisescans) || 
                (i == anzahlmaxima-1)) {
			if ((peaksInPeak_sum > maxPeaksPerSignal) && (peaksectionstartid != i)) {
				for (int n = peaksectionstartid; n < i; n++) {
					apexLocs.at(n) = 0;
				}
			}
			peaksectionstartid = i;
			peaksInPeak_sum = 0;	
		}
	}
	
	/* delete apexLocs that have been set to 0 */
	j = 0;
	for (int i = 0; i < (anzahlmaxima); ++i) { 
	  if (apexLocs[i] > 0) {
	  	apexLocs[j] = apexLocs[i];
      startLoc[j] = startLoc[i];
      endLoc[j] = endLoc[i];
      noiselevel[j] = noiselevel[i];
      peaksInPeak[j] = peaksInPeak[i];
      j++;
	  }
	}
	anzahlmaxima = j;
	apexLocs.resize(anzahlmaxima); 
	
			
	/* accurate noise calculation */
	std::vector<double> noisedeviation(anzahlmaxima,0);
	std::vector<double> noiseRegion(noisescans*2,0);
	int otherMaximum;
	
	for (int i = 0; i < anzahlmaxima; ++i) { 
		intensity = 0;
		noisecounter = 0;  
		std::fill(noiseRegion.begin(), noiseRegion.end(), 0);
		    	
  	/* check if there are peaks in front of this one, otherwise add the 
		respective scan to the noiseRegion */
		j = startLoc[i];
		otherMaximum = i-1;
		while ((j > (startLoc[i]-noisescans)) && (j > 0)) {
  		if ((otherMaximum >= 0) && (endLoc.at(otherMaximum) >= j)) {
  		/* if there is another peak, jump to the beginning of that peak*/
			j = startLoc.at(otherMaximum);
			otherMaximum -= 1;
  		} else {
  			j--;
				noiseRegion.at(noisecounter) = eic.at(j);
			  noisecounter++;
				intensity += eic.at(j);
			}
		}
  				
  	/* check if there are peaks after this one, otherwise add the respective scan to the noiseRegion*/
		j = endLoc[i];
		otherMaximum = i+1;
		/*while ((j < (endLoc[i]+noisescans)) && (j < eic.nrow()-1)) {*/
		while ((j < (endLoc[i]+noisescans)) && ((unsigned)j < eic.size()-1)) {
  		if ((otherMaximum < anzahlmaxima) && (startLoc.at(otherMaximum) <= j)) {
  			/* if there is another peak, jump to the beginning of that peak*/
  			j = endLoc.at(otherMaximum);
  			otherMaximum += 1;
  		} else {
  			j++;
  			noiseRegion.at(noisecounter) = eic.at(j);
				noisecounter++;
				intensity += eic.at(j);
			}
		}
		
		std::sort(noiseRegion.begin(), noiseRegion.begin()+noisecounter);
		noisedeviation[i] = noiseRegion.at(noisecounter*0.9)-noiseRegion.at(noisecounter*0.1);
    noiselevel[i] = intensity/noisecounter;
	} 
	
	/* delete apexLocs that are below S/N treshold */
	j = 0;
	for (int i = 0; i < anzahlmaxima; ++i) { 
	  if (eic[apexLocs[i]]-noiselevel[i] >= noisedeviation[i]*sn) {	
	  	apexLocs[j] = apexLocs[i];
      startLoc[j] = startLoc[i];
      endLoc[j] = endLoc[i];
      noiselevel[j] = noiselevel[i];
      noisedeviation[j] = noisedeviation[i];
      peaksInPeak[j] = peaksInPeak[i];
      j++;
	  }
	}
	anzahlmaxima = j;
	apexLocs.resize(anzahlmaxima); 
	
	
	// -- Determine FWHM -- 
	// We determine the RT of the left end of half peak heigh and the right end 
	// of half peak height. This is done by starting at the left end of the peak
	// and working backwards and vice versa for the right end. The actual position
	// is likely to be between two scans, so the slope computed between the two
	// scans and the corresponding scan time is computed
	std::vector<double> FWHM_left(anzahlmaxima,0);
	std::vector<double> FWHM_right(anzahlmaxima,0);
	double slope = 0;
	// Loop through all apexLocs
	for (int i = 0; i < anzahlmaxima; ++i) { 
	  
	  // Find the position of half height on the left side by starting at the left
	  // end of the peak and moving right until the intensity is no longer less
	  // than half of the peak apex intensity. j is then the scan number of half height
	  // We must make sure that j never goes beyond the end of the spectrum 
	  j = startLoc[i];
	  while ((unsigned)j < eic.size()-1 && 
          (eic[j]-noiselevel[i] < (eic[apexLocs[i]]-noiselevel[i])/2)) {
		  j++;
	  }
	  
	  // Compute the slope at this position and use this to compute the scantime
	  // between scans
	  slope = (eic[j]-eic[j-1])/(scantime[j]-scantime[j-1]);
	  if (slope == 0) {
	  	FWHM_left[i] = scantime[j];
	  } else {
	  	FWHM_left[i] = scantime[j-1]+(eic[apexLocs[i]]/2 - eic[j-1] + noiselevel[i])/slope;
	  }
	  
	  // Find the position of half height on the right side analogously to the 
	  // left side. Make sure j never goes beyond the beginning of the eic
	  j = endLoc[i];
	  while (j > 0 && (eic[j]-noiselevel[i] < (eic[apexLocs[i]]-noiselevel[i])/2)) {
		  j--;
	  }
    
    // j+1 might be longer than eic, in which case no slope can be computed
	  if ((unsigned)j+1 >= eic.size()) {
	    slope = 0;
	  } else {
	    slope = (eic[j+1]-eic[j])/(scantime[j+1]-scantime[j]);
	  }
	  
	  if (slope == 0) {
	  	FWHM_right[i] = scantime[j];
	  } else {
	 	  FWHM_right[i] = scantime[j]+(eic[apexLocs[i]]/2-eic[j]+noiselevel[i])/slope;
	  }
	}
	
	// -- Delete peaks that are too broad (2 x FWHM > peakwidth_max) and calculate area --
	j = 0;
	std::vector<double> area(anzahlmaxima, 0);
	
	for (int i = 0; i < anzahlmaxima; ++i) { 
	  if ((FWHM_right[i]-FWHM_left[i])*2 <= peakwidth_max) {
	  	apexLocs[j] = apexLocs[i];
      startLoc[j] = startLoc[i];
      endLoc[j] = endLoc[i];
      noiselevel[j] = noiselevel[i];
      noisedeviation[j] = noisedeviation[i];
      FWHM_left[j]= FWHM_left[i];
      FWHM_right[j] = FWHM_right[i];
      peaksInPeak[j] = peaksInPeak[i];
      for (int ii = startLoc[i]; ii < endLoc[i]; ++ii) {
      	if ((eic[ii] >= noiselevel[i]) && (eic.at(ii+1) >= noiselevel[i])) {
      	  area[j] += ((eic[ii]-noiselevel[i]) + (eic.at(ii+1)-noiselevel[i])) * 
      	    (scantime.at(ii+1) - scantime[ii])/2;
      	}
      }
      j++;
	  }
	}
	anzahlmaxima = j;
	apexLocs.resize(anzahlmaxima); 
   
  NumericMatrix ergebnis(anzahlmaxima, 16);


  for(int i = 0; i < anzahlmaxima; ++i) { 
    ergebnis(i,0)  = 0;                                    // placeholder m/z
    ergebnis(i,1)  = scantime.at(apexLocs.at(i));            // retention time
    ergebnis(i,2)  = 0;                                    // placeholder intensity
    ergebnis(i,3)  = eic.at(apexLocs.at(i))-noiselevel.at(i);// Intensity found in chromatogram
    ergebnis(i,4)  = apexLocs.at(i)+1;                       // scan number of peak
    ergebnis(i,5)  = scantime.at(startLoc.at(i));
    ergebnis(i,6)  = scantime.at(endLoc.at(i));
    ergebnis(i,7)  = startLoc.at(i)+1;
    ergebnis(i,8)  = endLoc.at(i)+1;
    ergebnis(i,9)  = noisedeviation.at(i);
    ergebnis(i,10) = area.at(i);                           // peak area
    ergebnis(i,11) = FWHM_left.at(i);                      // left RT of peak at half height (s)
    ergebnis(i,12) = FWHM_right.at(i);                     // right RT of peak at half height (s)
    ergebnis(i,13) = noiselevel.at(i);                     // Baseline
    ergebnis(i,14) = mz;                                   // m/z of this chromatogram
    ergebnis(i,15) = 0;                                    // placeholder ms2 scan number
	}
     
  return(ergebnis);
 
}

// The min intensity for the apexLocs calculation must be at most the 90%-tile intensity
double checkMinIntensity(const std::vector<double>& eic, double oldMinIntensity) {
  std::vector<double> eicIntensitiesSorted(eic.size(), 0);
  eicIntensitiesSorted.assign(eic.begin(), eic.end());
  std::sort(eicIntensitiesSorted.begin(), eicIntensitiesSorted.end());
  
  double newMinIntensity = eicIntensitiesSorted[eicIntensitiesSorted.size() * 0.9];
  if ((newMinIntensity == 0) || (newMinIntensity > oldMinIntensity)) {
    return oldMinIntensity;
  } else {
    return newMinIntensity;
  }
}


std::vector<double> computeDerivative(const std::vector<double>& eic) {
  std::vector<double> derivative((eic.size()-1));
  derivative[0] = eic[1] - eic[0];
  for(int i = 1; (unsigned)i < (eic.size()-1); ++i) {
    derivative[i] = eic[i+1] - eic[i];
  }
  return derivative;
}

std::vector<int> getApexIndices(const std::vector<double>& eic, const std::vector<double>& eicDerivative, double minIntensityChecked) {
  std::vector<int> crossingPoints(eic.size(), 0);
  int count = 0;
  for(int i = 1; (unsigned)i < (eic.size()-1); ++i) {
    if ((eicDerivative[(i-1)] > 0) && (eicDerivative[i] <= 0) && (eic[i] >= minIntensityChecked)) {
      crossingPoints.at(count) = i;
      count++;
    }
  }
  crossingPoints.resize(count);
  return crossingPoints;
}

int getPeakStartLoc(const std::vector<int>& apexLocs, const std::vector<double>& eicDerivative, const int apexNum) {
  int indLook = apexLocs[apexNum]-1;
  while ((eicDerivative[indLook] > 0) && (indLook > 0)) {
    indLook--;
  }
  return indLook+1; 
}

int getPeakEndLoc(const std::vector<int>& apexLocs, const std::vector<double>& eicDerivative, const int apexNum) {
  int indLook = apexLocs[apexNum];
  while (((unsigned)indLook < eicDerivative.size()) && (eicDerivative[indLook] < 0)) {
    indLook++;
  }
  return indLook;
}

// Copyright 2016-2025 Bundesanstalt für Gewässerkunde
// This file is part of ntsworkflow


