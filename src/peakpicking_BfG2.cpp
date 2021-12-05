#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
NumericMatrix peakPickingBfGC(double mz, double mz_step, std::vector<double> XIC,  
                              std::vector<double> scantime, double min_intensity, int sn, 
                              int noisescans, double peakwidth_min, double peakwidth_max,
                              int maxPeaksPerSignal) {
  
  std::vector<double> derivative((XIC.size()-1));
  std::vector<int> maxima(XIC.size(),0);
  int anzahlmaxima=0;
  int j = 0;
  /*int noisescans = 100;*/
  double min_intensity2 = 0.1;
  // maxPeaksPerSignal, this depends heavily on how noisy the chromatogram is and should be adjustable
  //int maxPeaksPerSignal = 100;

  // determine if the min intensity is set too high, at least for the first maxima calculation
  std::vector<double> intensities(XIC.size(),0);
  intensities.assign(XIC.begin(),XIC.end());
  std::sort(intensities.begin(),intensities.end());
  min_intensity2 = intensities[intensities.size()*0.9];
  if ((min_intensity2 == 0) || (min_intensity2 > min_intensity)) {
  	min_intensity2 = min_intensity;
  }
  
  // calculate the derivative and where the derivative crosses 0 from positive to negative
  derivative[0] = XIC[1]-XIC[0];
  for(int i = 1; i < (XIC.size()-1); ++i) {
	 derivative[i] = XIC[i+1]-XIC[i];
	 if ((derivative[(i-1)] > 0) && (derivative[i] <= 0) && (XIC[i] >= min_intensity2)) {
	 	maxima[anzahlmaxima] = i;
	 	anzahlmaxima++;
	 }
  }
  
  maxima.resize(anzahlmaxima);
  std::vector<int> left_end(anzahlmaxima,0);	 
  std::vector<int> right_end(anzahlmaxima,0);
  std::vector<double> noiselevel(anzahlmaxima,0);
  std::vector<int> amountofpeaks(anzahlmaxima,0);
  double intensity=0;
  int noisecounter=0;
  	 
  for(int i = 0; i < (anzahlmaxima); ++i) { 
  	/* find start of peak (left_end) */
	j = maxima[i]-1;
  	while ((derivative[j] > 0) && (j > 0)) {
  		j--;
   	}
  	left_end[i] = j+1; 
  	/* find end of peak (right_end) */
	j = maxima[i];
   	while ((derivative[j] < 0) && (j < XIC.size())) {
  		j++;
   	}
  	right_end[i] = j; 
  	
	/* 1st noiselevel calculation based on mean intensity around the peak */
	intensity = 0;
	noisecounter = 0;  
		    	
  	j = left_end[i];
  	while ((j > (left_end[i]-noisescans)) && (j > 0)) {
  		j--;
  		noisecounter++;
   		intensity += XIC[j];
  	}
  	
  
  	noiselevel[i] = intensity/noisecounter;
  	
	amountofpeaks[i] = 1;
  	/* check if two peaks belong together */
  	if ((i > 0) && (maxima[i] > 0)) {
  	  // determine if the current peak and the previous peak are next to each other (and we are not at the beginning)
  		if ((right_end[i-1] >= left_end[i]) && (maxima[i-1] > 0)) {
  		  // Determine if the valley between the current and previous peak is higher than half the height of the smallest peak
  		  // if this is the case the peaks belong together and will be joined
  			if ((XIC[left_end[i]]-noiselevel[i-1]) > (std::min(XIC[maxima[i]],XIC[maxima[i-1]])-noiselevel[i-1])/2) {
  			  // Determine if the smallest of the two peaks is higher than half the higher of the two peaks
  				if (std::min(XIC[maxima[i-1]],XIC[maxima[i]])-noiselevel[i-1] > (std::max(XIC[maxima[i-1]],XIC[maxima[i]])-noiselevel[i-1])/2) {
  					amountofpeaks[i] = amountofpeaks[i]+amountofpeaks[i-1];
  				} else {
  					if (amountofpeaks[i-1] > maxPeaksPerSignal) {
  						noiselevel[i-1] = XIC[maxima[i-1]];
  					  // if the previous peak is deemed to be noisy, then it is deleted by setting the left end to the current left end
              left_end[i-1] = left_end[i];
  					}

  					// here we look to see if the intensity has 'levelled off' by looking at the next 2 peaks
  					// this basically ends a tailing after it is no longer decreasing
  					if ((XIC[maxima[i-1]] > XIC[maxima[i]]) && (anzahlmaxima > i+1)) {
  						if (std::min({XIC[maxima[i]],XIC[maxima[i+1]],XIC[maxima[i+2]]}) > (std::max({XIC[maxima[i]],XIC[maxima[i+1]],XIC[maxima[i+2]]}))/2) {
  							right_end[i] = right_end[i-1];
  							amountofpeaks[i] = amountofpeaks[i-1];
  						}
  					}
  				}
  				if (XIC[maxima[i-1]] > XIC[maxima[i]]) maxima[i] = maxima[i-1];
  				noiselevel[i] = noiselevel[i-1];
            	left_end[i] = left_end[i-1];
            	maxima[i-1] = 0;
            	left_end[i-1] = 0;
            	right_end[i-1] = 0;
            	noiselevel[i-1] = 0;
            	amountofpeaks[i-1] = 0;
  			} else { 
           		 if (amountofpeaks[i-1] < maxPeaksPerSignal) noiselevel[i] = noiselevel[i-1];
          		}
  		}
  	}
  	
  	}
	
	// min_intensity = 1;
	
	/* delete maxima that have been set to 0, those below intensity threshold and too narrow ones*/
	j = 0;
	for(int i = 0; i < (anzahlmaxima); ++i) { 
	  if ((maxima[i] > 0) && (XIC[maxima[i]] > min_intensity) && (scantime[right_end[i]]-scantime[left_end[i]] > peakwidth_min)) {
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
	// in some cases amountofpeaks can be larger then 10... dunno why...
	// To correct for this set all back to 10
	for (int k = 0; k < anzahlmaxima; ++k) {
	  if (amountofpeaks[k] > maxPeaksPerSignal) {
	    amountofpeaks[k] = maxPeaksPerSignal;
	  }
	}
	
	
	int peaksectionstartid = 0;
	int amountofpeaks_sum = 0;
	
	/* check if there are peak clusters */
	for(int i = 1; i < (anzahlmaxima); ++i) { 
		amountofpeaks_sum += amountofpeaks[i-1];
		if ((right_end[i-1] < left_end[i]) || (std::min(XIC[maxima[i]],XIC[maxima[i-1]]) < std::max(XIC[maxima[i-1]],XIC[maxima[i]])/2) || ((maxima[i]-maxima[i-1]) > noisescans) || (i==anzahlmaxima-1)) {
			if ((amountofpeaks_sum > maxPeaksPerSignal) && (peaksectionstartid != i)) {
				
				for(int n = peaksectionstartid; n < i; n++) {
					maxima[n] = 0;
				}
			}
			peaksectionstartid = i;
			amountofpeaks_sum = 0;	
		}
		
	}
	
	/* delete maxima that have been set to 0 */
	j = 0;
	for(int i = 0; i < (anzahlmaxima); ++i) { 
	  if (maxima[i] > 0)  {
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
	
	for(int i = 0; i < (anzahlmaxima); ++i) { 
		intensity = 0;
		noisecounter = 0;  
		std::fill(noiseRegion.begin(),noiseRegion.end(),0);
		    	
  		/* check if there are peaks in front of this one, otherwise add the respective scan to the noiseRegion*/
		j = left_end[i];
  		otherMaximum = i-1;
  		while ((j > (left_end[i]-noisescans)) && (j > 0)) {
	  		if ((otherMaximum >= 0) && (right_end[otherMaximum] >= j)) {
	  			/* if there is another peak, jump to the beginning of that peak*/
				j = left_end[otherMaximum];
				otherMaximum -= 1;
	  		} else {
	  			j--;
  				noiseRegion[noisecounter] = XIC[j];
				noisecounter++;
  				intensity += XIC[j];
  			}
  		}
  				
  		/* check if there are peaks after this one, otherwise add the respective scan to the noiseRegion*/
		j = right_end[i];
  		otherMaximum = i+1;
  		/*while ((j < (right_end[i]+noisescans)) && (j < XIC.nrow()-1)) {*/
  		while ((j < (right_end[i]+noisescans)) && (j < XIC.size()-1)) {
	  		if ((otherMaximum < anzahlmaxima) && (left_end[otherMaximum] <= j)) {
	  			/* if there is another peak, jump to the beginning of that peak*/
				j = right_end[otherMaximum];
				otherMaximum += 1;
	  		} else {
	  			j++;
	  			noiseRegion[noisecounter] = XIC[j];
  				noisecounter++;
  				intensity += XIC[j];
  			}
  		}
  		std::sort(noiseRegion.begin(),noiseRegion.begin()+noisecounter);
  		noisedeviation[i] = noiseRegion[noisecounter*0.9]-noiseRegion[noisecounter*0.1];
	    noiselevel[i] = intensity/noisecounter;
	} 
	
	/* delete maxima that are below S/N treshold */
	j = 0;
	for(int i = 0; i < (anzahlmaxima); ++i) { 
	  if (XIC[maxima[i]]-noiselevel[i] >= noisedeviation[i]*sn)  {	
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
	for(int i = 0; i < (anzahlmaxima); ++i) { 
	  j = left_end[i];
	  while (XIC[j]-noiselevel[i] < (XIC[maxima[i]]-noiselevel[i])/2) {
		j++;
	  }
	  slope = (XIC[j]-XIC[j-1])/(scantime[j]-scantime[j-1]);
	  if (slope == 0) {
	  	FWHM_left[i] = scantime[j];
	  } else {
	  	FWHM_left[i] = scantime[j-1]+(XIC[maxima[i]]/2-XIC[j-1]+noiselevel[i])/slope;
	  }
	  
	  
	  j = right_end[i];
	  while (XIC[j]-noiselevel[i] < (XIC[maxima[i]]-noiselevel[i])/2) {
		j--;
	  }
	  slope = (XIC[j+1]-XIC[j])/(scantime[j+1]-scantime[j]);
	  if (slope == 0) {
	  	FWHM_right[i] = scantime[j];
	  } else {
	 	FWHM_right[i] = scantime[j]+(XIC[maxima[i]]/2-XIC[j]+noiselevel[i])/slope;
	 }
	}
	
	
	
	/* delete too broad peaks (2 x FWHM > peakwidth_max) and calculate area*/
	j = 0;
	std::vector<double> area(anzahlmaxima,0);
	
	for(int i = 0; i < (anzahlmaxima); ++i) { 
	  if ((FWHM_right[i]-FWHM_left[i])*2 <= peakwidth_max)  {
	  	maxima[j] = maxima[i];
        left_end[j] = left_end[i];
        right_end[j] = right_end[i];
        noiselevel[j] = noiselevel[i];
        noisedeviation[j] = noisedeviation[i];
        FWHM_left[j]= FWHM_left[i];
        FWHM_right[j] = FWHM_right[i];
        amountofpeaks[j] = amountofpeaks[i];
        for (int ii = left_end[i]; ii < right_end[i]; ++ii) {
        	if ((XIC[ii] >= noiselevel[i]) && (XIC[ii+1] >= noiselevel[i])) area[j] += ((XIC[ii]-noiselevel[i])+(XIC[ii+1]-noiselevel[i]))*(scantime[ii+1]-scantime[ii])/2;
        }
        j++;
	  }
	}
	anzahlmaxima = j;
	maxima.resize(anzahlmaxima); 
	   
  
   
   NumericMatrix ergebnis(anzahlmaxima,16);
   
   for(int i = 0; i < (anzahlmaxima); ++i) { 
	ergebnis(i,0) = 0;
	ergebnis(i,1) = scantime[maxima[i]];
	ergebnis(i,2) = 0;
	ergebnis(i,3) = XIC[maxima[i]]-noiselevel[i];
	ergebnis(i,4) = maxima[i]+1;
	ergebnis(i,5) = scantime[left_end[i]];
	ergebnis(i,6) = scantime[right_end[i]];
	ergebnis(i,7) = left_end[i]+1;
	ergebnis(i,8) = right_end[i]+1;
	ergebnis(i,9) = noisedeviation[i];
	ergebnis(i,10) = area[i];
	ergebnis(i,11) = FWHM_left[i];
	ergebnis(i,12) = FWHM_right[i];
	ergebnis(i,13) = noiselevel[i];
	ergebnis(i,14) = mz;
	ergebnis(i,15) = 0; /*ms2scan*/
	}
   
   
  //Rcout << mz << std::endl; 
     
  return(ergebnis);
  
  
}

