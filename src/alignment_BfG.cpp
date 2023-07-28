#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Copyright 2016-2023 Bundesanstalt für Gewässerkunde
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


double sum(NumericVector a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += a[i];
  }
  return s;
}

double mean(NumericVector a)
{
  return sum(a) / a.size();
}

double sqsum(NumericVector a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += pow(a[i], 2);
  }
  return s;
}

double stdev(NumericVector nums)
{
  double N = nums.size();
  return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

// [[Rcpp::export]]
double pearsoncoeff(NumericVector X, NumericVector Y)
{
  double summe = 0;
  for(int i = 0; i < X.size(); ++i) {
    summe += (X[i] - mean(X))*(Y[i] - mean(Y));
  }
  return summe / (X.size()*stdev(X)* stdev(Y));
}

// [[Rcpp::export]]
LogicalVector correlates_with(NumericMatrix aligned_intensities, int row, double coefficient)
{
  LogicalVector ergebnis(aligned_intensities.nrow());
  for(int i = 0; i < aligned_intensities.nrow(); ++i) {
    if (pearsoncoeff(aligned_intensities(row-1,_),aligned_intensities(i,_)) > coefficient) {
      ergebnis[i] = true;
    } else {
      ergebnis[i] = false;
    }
  }	
  return(ergebnis);
}

// [[Rcpp::export]]
NumericMatrix correlation(NumericMatrix aligned_intensities)
{
  NumericMatrix ergebnis(aligned_intensities.nrow(),aligned_intensities.nrow());
  for(int j = 0; j < aligned_intensities.nrow(); ++j) {
    for(int i = 0; i < aligned_intensities.nrow(); ++i) {
      ergebnis(j,i) = pearsoncoeff(aligned_intensities(j,_),aligned_intensities(i,_));
    }
    ergebnis(j,j) = 0;
  }
  return(ergebnis);
}

// [[Rcpp::export]]
IntegerMatrix alignmentBfGC(List peaklistR, int mz_dev, int DeltaRT, int mz_dev_unit) {  // change ppm_dev to mz_dev
  NumericMatrix peaklistn = as<NumericMatrix>(peaklistR[0]);


  
  int peaklistelements = peaklistR.size();
  int suchstart = 0;
  int spaltenlaenge = 0;
  int maximalezeilen = 0;

  for(int i = 0; i < peaklistelements; ++i) {
    peaklistn = as<NumericMatrix>(peaklistR[i]);
  	if (spaltenlaenge < peaklistn.nrow()) spaltenlaenge = peaklistn.nrow();

  	maximalezeilen += peaklistn.nrow();
  }
  
  IntegerMatrix ergebnis(maximalezeilen,peaklistelements);
  LogicalMatrix erledigt(spaltenlaenge,peaklistelements);
  LogicalMatrix kandidaten(spaltenlaenge,peaklistelements);
  std::fill(ergebnis.begin(),ergebnis.end(),0);
  std::fill(erledigt.begin(),erledigt.end(),false);
  std::fill(kandidaten.begin(),kandidaten.end(),false);
   
  int iii = 0;
  bool abbruch = false;
  bool ende = false;
  double masse = 0;
  double masse_mz_dev = 0;
  int RT=0;
  int zeile=0;
  int ergebniszeile = 0;
  NumericVector highest_i(2);
  int anzahlkandidaten=0;
  int kandidat_zeile=0; 
  
  highest_i[0] = 0;
  highest_i[1] = 0;
  
  while (ende == false) { 
  /* for (int j = 0; j < maximalezeilen; ++j) {  */

  iii = suchstart;
  abbruch = false;
  ende = true;
  // Schleife sucht Elmente die 
   while (abbruch == false) { 
   	for(int i = 0; i < peaklistelements; ++i) {
  		peaklistn = as<NumericMatrix>(peaklistR[i]);
	  	if (iii < peaklistn.nrow()) {
	  		if ((erledigt(iii,i) == false) and (abbruch == false)) {
	  			masse = peaklistn(iii,0);
	  			RT = peaklistn(iii,1);
	  			highest_i[0] = i;
	  			highest_i[1] = peaklistn(iii,2);
	  			abbruch = true;
	  			ende = false;
			}	
	  	}
	  	if (iii > spaltenlaenge) abbruch = true;
  	}
  	suchstart = iii;
  	iii += 1;
  } 
  
  // calculate the mz tolerance
  if (mz_dev_unit == 1) {
    masse_mz_dev = masse*mz_dev*2/1000000;
  } else {
    masse_mz_dev = (double)mz_dev*2/1000;  // must be cast as double
  }
  

  
  for(int i = 0; i < peaklistelements; ++i) {
  peaklistn = as<NumericMatrix>(peaklistR[i]);
  zeile = suchstart-1;
    // Schleife zählt die Zeilen mit der Bedingung der Massentoleranz
    while ((peaklistn(zeile,0) < masse-masse_mz_dev) and (zeile < peaklistn.nrow())) {
    	zeile += 1;
    }
    while ((peaklistn(zeile,0) <= masse+masse_mz_dev) and (zeile < peaklistn.nrow())) {
      
    	if ((peaklistn(zeile,1) >= RT-DeltaRT*2) and (peaklistn(zeile,1) <= RT+DeltaRT*2)) {
    	  
    	  if (erledigt(zeile,i) == false) {
    	  	kandidaten(zeile,i) = true;
      		if (peaklistn(zeile,2) > highest_i[1]) {
      			highest_i[0] = i;
      			highest_i[1] = peaklistn(zeile,2);
      			masse = peaklistn(zeile,0);
      			RT = peaklistn(zeile,1);
      		}
      		} 
    	
    	}
  	zeile += 1;
  }
  }
   
   // calculate the mz tolerance according to this new mass
   if (mz_dev_unit == 1) {
     masse_mz_dev = masse*mz_dev/1000000;
   } else {
     masse_mz_dev = (double)mz_dev/1000;
   }
  
  //  std::cout << masse_mz_dev << std::endl;
  
  for(int i = 0; i < peaklistelements; ++i) {
  peaklistn = as<NumericMatrix>(peaklistR[i]);
  anzahlkandidaten = std::count(kandidaten(_,i).begin(),kandidaten(_,i).end(),true);


  
    if (anzahlkandidaten > 0) {
  	zeile = 0;
    	for(int kandidatennr = 0; kandidatennr < anzahlkandidaten; ++kandidatennr) {
    		while (kandidaten(zeile,i) == false) {
    			zeile += 1;
    		}
    		if (kandidatennr == 0) {
    			kandidat_zeile = zeile;
    		} else {
    			if ((abs(peaklistn(zeile,1)-RT) <= abs(peaklistn(kandidat_zeile,1)-RT)) and (abs(peaklistn(zeile,1)-RT) <= DeltaRT) and (abs(peaklistn(zeile,0)-masse) <= masse_mz_dev)) kandidat_zeile = zeile;
    		}
    		kandidaten(zeile,i) = false;
    		zeile += 1;
    	}
  	
    	if ((abs(peaklistn(kandidat_zeile,1)-RT) <= DeltaRT) and (abs(peaklistn(kandidat_zeile,0)-masse) <= masse_mz_dev)) {
    		ergebnis(ergebniszeile,i) = kandidat_zeile+1;
    		erledigt(kandidat_zeile,i) = true;
    	}
    }
  }
  ergebniszeile += 1;
   } 
  
  
  return ergebnis(Range(0,ergebniszeile-2),_);
  
}

