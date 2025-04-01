

#ifndef PEAKPICKING_H
#define PEAKPICKING_H

#include <Rcpp.h>

double checkMinIntensity(const std::vector<double>& eic, double oldMinIntensity);
std::vector<double> computeDerivative(const std::vector<double>& eic);
std::vector<int> getApexIndices(const std::vector<double>& eic, const std::vector<double>& eicDerivative, double minIntensityChecked);
int getPeakStartLoc(const std::vector<int>& apexLocs, const std::vector<double>& eicDerivative, const int apexNum);
int getPeakEndLoc(const std::vector<int>& apexLocs, const std::vector<double>& eicDerivative, const int apexNum);
double getNoiseLevel(const std::vector<int>& startLoc, const std::vector<double>& eic, const int noisescans, int apexNum);
#endif