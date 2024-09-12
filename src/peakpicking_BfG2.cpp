


// Copyright 2016-2024 Bundesanstalt für Gewässerkunde
// This file is part of ntsworkflow
// ntsworkflow is free software: you can redistribute it and/or modify it under the 
// terms of the GNU General Public License as published by the Free Software 
// Foundation, either version 3 of the License, or (at your option) any 
// later version.
// 
// ntsworkflow is distributed in the hope that it will be useful, but WITHOUT ANY 
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along 
// with ntsworkflow. If not, see <https://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>
using namespace Rcpp;

//' Find peaks in an ion chromatogram
//'
//' @description Peak finding algorithm using maxima detection by 1st derivative, 
//' an iterative search method and no chromatogram smoothing. The method is 
//' published in Dietrich, C., Wick, A., & Ternes, T. A. (2021). 
//' Open source feature detection for non‐target LC‐MS analytics. 
//' Rapid Communications in Mass Spectrometry, e9206. https://doi.org/10.1002/rcm.9206 
//'
//' @param mz m/z of current ion chromatogram (Da)
//' @param mz_step binning width used to extract chromatogram (da)
//' @param XIC ion chromatogram (intensities)
//' @param scantime Scan time of each index in XIC (in s) 
//' @param min_intensity Minimum intensity for peak-picking
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
//' 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix peakPickingBfGC(
    double mz, 
    double mz_step, 
    std::vector<double> XIC,  
    std::vector<double> scantime, 
    double min_intensity, 
    int sn, 
    int noisescans, 
    double peakwidth_min, 
    double peakwidth_max,
    int maxPeaksPerSignal
  ) {
  
   
  std::vector<double> derivative((XIC.size()-1));
  std::vector<int> maxima(XIC.size(), 0);
  int anzahlmaxima = 0;
  int j = 0;
  
  double min_intensity2 = 0.1;
  
  // determine if the min intensity is set too high, at least for the first maxima calculation
  // Set the min_intensity2 to 90%-tile intensity
  std::vector<double> intensities(XIC.size(), 0);
  intensities.assign(XIC.begin(), XIC.end());
  std::sort(intensities.begin(), intensities.end());
  min_intensity2 = intensities[intensities.size() * 0.9];
  if ((min_intensity2 == 0) || (min_intensity2 > min_intensity)) {
  	min_intensity2 = min_intensity;
  }
  
  // Detection of local maxima within the chromatogramm
  // calculate the derivative and where the derivative crosses 0 from positive to negative
  derivative[0] = XIC[1] - XIC[0];
  for(int i = 1; (unsigned)i < (XIC.size()-1); ++i) {
	  derivative[i] = XIC[i+1] - XIC[i];
	  if ((derivative[(i-1)] > 0) && (derivative[i] <= 0) && (XIC[i] >= min_intensity2)) {
	 	  maxima.at(anzahlmaxima) = i;
	 	  anzahlmaxima++;
	  }
  }
  
  // Determine peak boundaries
  maxima.resize(anzahlmaxima);
  std::vector<int> left_end(anzahlmaxima, 0);	 
  std::vector<int> right_end(anzahlmaxima, 0);
  std::vector<double> noiselevel(anzahlmaxima, 0);
  std::vector<int> amountofpeaks(anzahlmaxima, 0);
  double intensity = 0;
  int noisecounter = 0;

  // --Loop through each local maximum to merge adjacent maxima--
  // 
  // If the intensity of the valley between the current maximum
  // and the previous one is lower than 1/2 of the lower of the two
  // maxima, then the current maximum is merged with the previous one
  // and the new maximum gets the intensity of the higher maximum.
  // The merging is done by setting the previous peak's intensity to
  // zero (this way it is marked for deletion later)
  for (int i = 0; i < anzahlmaxima; ++i) { 
    
    // Find start of peak (left_end) 
    j = maxima[i]-1;
    while ((derivative[j] > 0) && (j > 0)) {
      j--;
    }
    left_end[i] = j+1; 
    
    // find end of peak (right_end)
    j = maxima[i];
    while (((unsigned)j < derivative.size()) && (derivative[j] < 0)) {
      j++;
    }
    
    right_end[i] = j; 
    
    // 1st noiselevel calculation based on mean intensity in front of peak
    intensity = 0;
    noisecounter = 0;  
    
    // Compute average intensity in the area before the peak by summing all
    // intensities in a loop and then dividing by the number of iterations
    j = left_end[i];
    while ((j > (left_end[i]-noisescans)) && (j > 0)) {
      j--;
      noisecounter++;
      intensity += XIC.at(j);
    }
    
    noiselevel[i] = intensity/noisecounter;
    
    amountofpeaks[i] = 1;
    // Merging of consecutive local maxima
    // Check if two peaks belong together 
    if ((i > 0) && (maxima[i] > 0)) {
      // determine if the current peak and the previous peak are next to each 
      // other (and we are not at the beginning)
      if ((right_end[i-1] >= left_end[i]) && (maxima[i-1] > 0)) {
        // Determine if the valley between the current and previous 
        // peak is higher than half the height of the smallest peak
        // if this is the case the peaks belong together and will be joined
        if ((XIC[left_end[i]]-noiselevel[i-1]) > 
              (std::min(XIC[maxima[i]],XIC[maxima[i-1]])-noiselevel[i-1])/2) {
          // Determine if the smallest of the two peaks is higher than half the 
          // higher of the two peaks
          if (std::min(XIC[maxima[i-1]],XIC[maxima[i]])-noiselevel[i-1] > 
                (std::max(XIC[maxima[i-1]],XIC[maxima[i]])-noiselevel[i-1])/2) {
            amountofpeaks[i] = amountofpeaks[i]+amountofpeaks[i-1];
          } else {
            if (amountofpeaks[i-1] > maxPeaksPerSignal) {
              noiselevel[i-1] = XIC[maxima[i-1]];
              // if the previous peak is deemed to be noisy, then it is deleted
              // by setting the left end to the current left end
              left_end[i-1] = left_end[i];
            }
            
            // Here we look to see if the intensity has 'levelled off' by 
            // looking at the next 2 peaks
            // this basically ends a tailing after it is no longer decreasing
            if ((XIC[maxima[i-1]] > XIC[maxima[i]]) && (anzahlmaxima > i+2)) {
              double lowInt = std::min({XIC[maxima[i]], XIC[maxima[i+1]], XIC[maxima[i+2]]});
              double halfHighInt = (std::max({XIC[maxima[i]], XIC[maxima.at(i+1)], XIC[maxima.at(i+2)]}))/2;
              if ( lowInt > halfHighInt ) {
                right_end[i] = right_end[i-1];
                amountofpeaks[i] = amountofpeaks[i-1];
              }
            }
          }
          
          if (XIC[maxima[i-1]] > XIC[maxima[i]]) {
            maxima[i] = maxima[i-1];
          }
          
          noiselevel[i] = noiselevel[i-1];
          left_end[i] = left_end[i-1];
          maxima[i-1] = 0;
          left_end[i-1] = 0;
          right_end[i-1] = 0;
          noiselevel[i-1] = 0;
          amountofpeaks[i-1] = 0;
        } else { 
          if (amountofpeaks[i-1] < maxPeaksPerSignal) {
            noiselevel[i] = noiselevel[i-1];
          }
        }
      }
    }
  }

	 


	// min_intensity = 1;
	
	/* delete maxima that have been set to 0, those below intensity threshold and
  too narrow ones*/
	j = 0;
	for (int i = 0; i < anzahlmaxima; ++i) { 
	  if ((maxima[i] > 0) && (XIC[maxima[i]] > min_intensity) && 
         (scantime.at(right_end[i])-scantime.at(left_end[i]) > peakwidth_min)) {
	  	maxima[j] = maxima[i];
      left_end[j] = left_end[i];
      right_end[j] = right_end[i];
      noiselevel[j] = noiselevel[i];
      amountofpeaks[j] = amountofpeaks[i];
      j++;
	  }
	}
	anzahlmaxima = j;
	maxima.resize(anzahlmaxima); 
	
	// 190722 KJ
	// in some cases amountofpeaks can be larger then 10...
	// To correct for this set all back to 10
	for (int k = 0; k < anzahlmaxima; ++k) {
	  if (amountofpeaks.at(k) > maxPeaksPerSignal) {
	    amountofpeaks.at(k) = maxPeaksPerSignal;
	  }
	}
	
	int peaksectionstartid = 0;
	int amountofpeaks_sum = 0;
	
	/* check if there are peak clusters */
	for (int i = 1; i < anzahlmaxima; ++i) { 
		amountofpeaks_sum += amountofpeaks[i-1];
		if ((right_end[i-1] < left_end[i]) || 
          (std::min(XIC[maxima[i]],XIC[maxima[i-1]]) < 
            std::max(XIC[maxima[i-1]],XIC[maxima[i]])/2) || 
              ((maxima[i]-maxima[i-1]) > noisescans) || 
                (i == anzahlmaxima-1)) {
			if ((amountofpeaks_sum > maxPeaksPerSignal) && (peaksectionstartid != i)) {
				for (int n = peaksectionstartid; n < i; n++) {
					maxima.at(n) = 0;
				}
			}
			peaksectionstartid = i;
			amountofpeaks_sum = 0;	
		}
	}
	
	/* delete maxima that have been set to 0 */
	j = 0;
	for (int i = 0; i < (anzahlmaxima); ++i) { 
	  if (maxima[i] > 0) {
	  	maxima[j] = maxima[i];
      left_end[j] = left_end[i];
      right_end[j] = right_end[i];
      noiselevel[j] = noiselevel[i];
      amountofpeaks[j] = amountofpeaks[i];
      j++;
	  }
	}
	anzahlmaxima = j;
	maxima.resize(anzahlmaxima); 
	
			
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
		j = left_end[i];
		otherMaximum = i-1;
		while ((j > (left_end[i]-noisescans)) && (j > 0)) {
  		if ((otherMaximum >= 0) && (right_end.at(otherMaximum) >= j)) {
  		/* if there is another peak, jump to the beginning of that peak*/
			j = left_end.at(otherMaximum);
			otherMaximum -= 1;
  		} else {
  			j--;
				noiseRegion.at(noisecounter) = XIC.at(j);
			  noisecounter++;
				intensity += XIC.at(j);
			}
		}
  				
  	/* check if there are peaks after this one, otherwise add the respective scan to the noiseRegion*/
		j = right_end[i];
		otherMaximum = i+1;
		/*while ((j < (right_end[i]+noisescans)) && (j < XIC.nrow()-1)) {*/
		while ((j < (right_end[i]+noisescans)) && ((unsigned)j < XIC.size()-1)) {
  		if ((otherMaximum < anzahlmaxima) && (left_end.at(otherMaximum) <= j)) {
  			/* if there is another peak, jump to the beginning of that peak*/
  			j = right_end.at(otherMaximum);
  			otherMaximum += 1;
  		} else {
  			j++;
  			noiseRegion.at(noisecounter) = XIC.at(j);
				noisecounter++;
				intensity += XIC.at(j);
			}
		}
		
		std::sort(noiseRegion.begin(), noiseRegion.begin()+noisecounter);
		noisedeviation[i] = noiseRegion.at(noisecounter*0.9)-noiseRegion.at(noisecounter*0.1);
    noiselevel[i] = intensity/noisecounter;
	} 
	
	/* delete maxima that are below S/N treshold */
	j = 0;
	for (int i = 0; i < anzahlmaxima; ++i) { 
	  if (XIC[maxima[i]]-noiselevel[i] >= noisedeviation[i]*sn) {	
	  	maxima[j] = maxima[i];
      left_end[j] = left_end[i];
      right_end[j] = right_end[i];
      noiselevel[j] = noiselevel[i];
      noisedeviation[j] = noisedeviation[i];
      amountofpeaks[j] = amountofpeaks[i];
      j++;
	  }
	}
	anzahlmaxima = j;
	maxima.resize(anzahlmaxima); 
	
	
	/*calculate FWHM*/
	std::vector<double> FWHM_left(anzahlmaxima,0);
	std::vector<double> FWHM_right(anzahlmaxima,0);
	double slope = 0;
	for (int i = 0; i < anzahlmaxima; ++i) { 
	  j = left_end[i];
	  while (XIC[j]-noiselevel[i] < (XIC[maxima[i]]-noiselevel[i])/2) {
		  j++;
	  }
	  slope = (XIC[j]-XIC[j-1])/(scantime[j]-scantime[j-1]);
	  if (slope == 0) {
	  	FWHM_left[i] = scantime[j];
	  } else {
	  	FWHM_left[i] = scantime[j-1]+(XIC[maxima[i]]/2 - XIC[j-1] + noiselevel[i])/slope;
	  }
	  
	  j = right_end[i];
	  while (XIC[j]-noiselevel[i] < (XIC[maxima[i]]-noiselevel[i])/2) {
		  j--;
	  }
    
    // j+1 might be longer than XIC
	  if ((unsigned)j+1 >= XIC.size()) {
	    slope = 0;
	  } else {
	    slope = (XIC[j+1]-XIC[j])/(scantime[j+1]-scantime[j]);
	  }
	  
	  if (slope == 0) {
	  	FWHM_right[i] = scantime[j];
	  } else {
	 	  FWHM_right[i] = scantime[j]+(XIC[maxima[i]]/2-XIC[j]+noiselevel[i])/slope;
	  }
	}
	
	/* delete too broad peaks (2 x FWHM > peakwidth_max) and calculate area*/
	j = 0;
	std::vector<double> area(anzahlmaxima, 0);
	
	for (int i = 0; i < anzahlmaxima; ++i) { 
	  if ((FWHM_right[i]-FWHM_left[i])*2 <= peakwidth_max) {
	  	maxima[j] = maxima[i];
      left_end[j] = left_end[i];
      right_end[j] = right_end[i];
      noiselevel[j] = noiselevel[i];
      noisedeviation[j] = noisedeviation[i];
      FWHM_left[j]= FWHM_left[i];
      FWHM_right[j] = FWHM_right[i];
      amountofpeaks[j] = amountofpeaks[i];
      for (int ii = left_end[i]; ii < right_end[i]; ++ii) {
      	if ((XIC[ii] >= noiselevel[i]) && (XIC.at(ii+1) >= noiselevel[i])) {
      	  area[j] += ((XIC[ii]-noiselevel[i]) + (XIC.at(ii+1)-noiselevel[i])) * 
      	    (scantime.at(ii+1) - scantime[ii])/2;
      	}
      }
      j++;
	  }
	}
	anzahlmaxima = j;
	maxima.resize(anzahlmaxima); 
   
  NumericMatrix ergebnis(anzahlmaxima, 16);


  for(int i = 0; i < anzahlmaxima; ++i) { 
    ergebnis(i,0)  = 0;                                    // placeholder m/z
    ergebnis(i,1)  = scantime.at(maxima.at(i));            // retention time
    ergebnis(i,2)  = 0;                                    // placeholder intensity
    ergebnis(i,3)  = XIC.at(maxima.at(i))-noiselevel.at(i);// Intensity found in chromatogram
    ergebnis(i,4)  = maxima.at(i)+1;                       // scan number of peak
    ergebnis(i,5)  = scantime.at(left_end.at(i));
    ergebnis(i,6)  = scantime.at(right_end.at(i));
    ergebnis(i,7)  = left_end.at(i)+1;
    ergebnis(i,8)  = right_end.at(i)+1;
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

