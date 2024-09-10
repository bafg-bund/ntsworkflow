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
  
  // take the first peaklist
  NumericMatrix peaklistn = as<NumericMatrix>(peaklistR[0]);
  
  int peaklistelements = peaklistR.size();
  int suchstart = 0;
  // The longest peaklist
  int spaltenlaenge = 0;
  int maximalezeilen = 0;
  
  // Sum of all rows in all peaklists
  for (int i = 0; i < peaklistelements; ++i) {
    peaklistn = as<NumericMatrix>(peaklistR[i]);
    if (spaltenlaenge < peaklistn.nrow()) {
      spaltenlaenge = peaklistn.nrow();
    }
    
    maximalezeilen += peaklistn.nrow();
  }
  
  // ergebnis holds results by reference (pl-row, pl-index)
  // numbers on the same row mean they have been aligned.
  // the maximum number of rows would be reached if no feature could be
  // aligned to any to any other feature, so this would be the total number
  // of features
  IntegerMatrix ergebnis(maximalezeilen, peaklistelements);
  // records if features are still "open" or if they have been grabbed
  // already, by reference (pl row, pl index) 
  // can be at the most the number of features in the largest
  // peaklist by total number of peaklists
  LogicalMatrix erledigt(spaltenlaenge, peaklistelements);
  // holds the current candidates by reference (pl-row, pl-index)
  LogicalMatrix kandidaten(spaltenlaenge, peaklistelements);
  std::fill(ergebnis.begin(), ergebnis.end(),0);
  std::fill(erledigt.begin(), erledigt.end(),false);
  std::fill(kandidaten.begin(), kandidaten.end(),false);
  
  // This counter 
  int iii = 0;
  bool abbruch = false;
  bool ende = false;
  double masse = 0;
  double masse_mz_dev = 0;
  int RT = 0;
  int zeile = 0;
  int ergebniszeile = 0;
  // For intensity, the first value holds the index of the pl in which the
  // intensity is found, the second the intensity itself
  NumericVector highest_i(2);
  int anzahlkandidaten = 0;
  int kandidat_zeile = 0; 
  
  highest_i[0] = 0;
  highest_i[1] = 0;
  
  while (ende == false) { 
    
    
    iii = suchstart;
    abbruch = false;
    ende = true;
    // Loop to grab the next available feature which has not yet been processed
    while (abbruch == false) { 
      // Loop though each peaklist
      for (int i = 0; i < peaklistelements; ++i) {
        peaklistn = as<NumericMatrix>(peaklistR[i]);
        if (iii < peaklistn.nrow()) {
          if ((erledigt(iii, i) == false) and (abbruch == false)) {
            masse = peaklistn(iii, 0);
            RT = peaklistn(iii, 1);
            highest_i[0] = i;
            highest_i[1] = peaklistn(iii, 2);
            abbruch = true;
            ende = false;
          }	
        }
        // If we have passed the largest peaklist nrow, can also stop
        if (iii > spaltenlaenge) {
          abbruch = true;
        }
      }
      suchstart = iii;
      iii += 1;
    } 
    
    // calculate the mz tolerance for this mass
    if (mz_dev_unit == 1) {
      masse_mz_dev = masse*mz_dev*2/1000000;
    } else {
      masse_mz_dev = (double)mz_dev*2/1000;  // must be cast as double
    }
    
    // Loop through all peaklists and get candidates for this mass
    // Find the candidate with the highest intensity, this becomes
    // the reference feature for further comparisons
    for (int i = 0; i < peaklistelements; ++i) {
      peaklistn = as<NumericMatrix>(peaklistR[i]);
      zeile = suchstart-1;
      // Schleife zählt die Zeilen mit der Bedingung der Massentoleranz
      // goes through all masses until it reaches the highest mass less than
      // mass minus tolerance
      while ((peaklistn(zeile,0) < masse-masse_mz_dev) and (zeile < peaklistn.nrow())) {
        zeile += 1;
      }
      // The loop continutes until it finds the highest mass <= mass plus tolerance
      while ((peaklistn(zeile,0) <= masse+masse_mz_dev) and (zeile < peaklistn.nrow())) {
        
        // If the feature is within the 2*RT tolerance
        if ((peaklistn(zeile,1) >= RT-DeltaRT*2) and (peaklistn(zeile,1) <= RT+DeltaRT*2)) {
          // and still "open"
          if (erledigt(zeile,i) == false) {
            // it is recorded as one of the candidates
            // if it is the candidate with the highest intensity, 
            // this features is taken as the "true" mass and RT for any
            // further comparisons
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
    
    // At this stage we have a number of possible candidates
    
    // Calculate the mz tolerance according to this new reference mass
    // (the candidate with the highest intensity)
    if (mz_dev_unit == 1) {
      masse_mz_dev = masse*mz_dev/1000000;
    } else {
      masse_mz_dev = (double)mz_dev/1000;
    }
    
    // Loop through each peaklist
    for (int i = 0; i < peaklistelements; ++i) {
      peaklistn = as<NumericMatrix>(peaklistR[i]);
      // Get the total number of candidates in this peaklist
      anzahlkandidaten = std::count(kandidaten(_,i).begin(),kandidaten(_,i).end(),true);
      
      if (anzahlkandidaten > 0) {
        zeile = 0;
        for (int kandidatennr = 0; kandidatennr < anzahlkandidaten; ++kandidatennr) {
          // Get the first row where a candidate is found (in this peaklist)
          while (kandidaten(zeile, i) == false) {
            zeile += 1;
          }
          // first iteration we just set the kandidat_zeile to the current row
          // so this is a start point for the next comparison
          if (kandidatennr == 0) {
            kandidat_zeile = zeile;
          } else {
            // In proceeding candidates we check if this candidate is closer to
            // the current reference RT than the preceeding candidate (and the tolerances are met), 
            // if so, this becomes the new candidate for comparison
            if ((abs(peaklistn(zeile,1)-RT) <= abs(peaklistn(kandidat_zeile,1)-RT)) and 
                  (abs(peaklistn(zeile,1)-RT) <= DeltaRT) and 
                  (abs(peaklistn(zeile,0)-masse) <= masse_mz_dev)) {
              kandidat_zeile = zeile;
            }
          }
          kandidaten(zeile,i) = false;
          zeile += 1;
        }
        // At the end of this loop we have kandidat_zeile, which is the best option
        // in the current peaklist (closest to reference)
        
        if ((abs(peaklistn(kandidat_zeile,1)-RT) <= DeltaRT) and 
              (abs(peaklistn(kandidat_zeile,0)-masse) <= masse_mz_dev)) {
          // +1 to convert to R type index counting
          ergebnis(ergebniszeile,i) = kandidat_zeile+1; 
          // This feature has no been "taken" so it is no longer checked
          // when going to the next feature or when looking for candidates
          erledigt(kandidat_zeile,i) = true;
        }
      }
    }
    
    ergebniszeile += 1;
  } 
  
  // Return all columns but only rows that were filled
  return ergebnis(Range(0,ergebniszeile-2),_);
  
}

