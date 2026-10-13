# Package index

## All functions

- [`Report-class`](http://r.bafg.de/reference/Report-class.md)
  [`Report`](http://r.bafg.de/reference/Report-class.md) : A reference
  class for database screening.
- [`RunGenForm()`](http://r.bafg.de/reference/RunGenForm.md) : Run
  GenForm with various Commands
- [`alig_componentisation()`](http://r.bafg.de/reference/alig_componentisation.md)
  : Second stage componentisation using the alignment table
- [`alignmentBfGC()`](http://r.bafg.de/reference/alignmentBfGC.md) :
  Align peaks in separate peak lists by m/z and retention time tolerance
- [`alignment_BfG_cpp()`](http://r.bafg.de/reference/alignment_BfG_cpp.md)
  : Alignment of peaks from difference samples
- [`alignment_BfG_cpp_par()`](http://r.bafg.de/reference/alignment_BfG_cpp_par.md)
  : Align peaks with parallel processing
- [`annotateCleanMSMS()`](http://r.bafg.de/reference/annotateCleanMSMS.md)
  : Annotate Cleaned MSMS
- [`annotate_grouped()`](http://r.bafg.de/reference/annotate_grouped.md)
  : Annotate alignment table with compound names by spectral library
  screening
- [`annotate_grouped_components()`](http://r.bafg.de/reference/annotate_grouped_components.md)
  : Annotate alignment table with adduct and isotopologue information
- [`annotate_grouped_mz_rt()`](http://r.bafg.de/reference/annotate_grouped_mz_rt.md)
  : Annotate alignment table based only on m/z and RT comparison
- [`averageReplicates()`](http://r.bafg.de/reference/averageReplicates.md)
  : Get mean intensities for replicate injections
- [`blankCorrection()`](http://r.bafg.de/reference/blankCorrection.md) :
  Remove features also found in the blanks
- [`calc_ndp()`](http://r.bafg.de/reference/calc_ndp.md) : Calculate
  Normalized Dot Product Similarity Score Between MS2 Spectra
- [`clean_and_normalize_ms2()`](http://r.bafg.de/reference/clean_and_normalize_ms2.md)
  : Filter and normalize MS2 spectra by removing precursor peaks and
  noise.
- [`componentization_BfG()`](http://r.bafg.de/reference/componentization_BfG.md)
  : Group features into components
- [`correct_formula()`](http://r.bafg.de/reference/correct_formula.md) :
  Combine a formula with the adduct
- [`correlates_with_r()`](http://r.bafg.de/reference/correlates_with_r.md)
  : Find Pearson's correlation between intensity trends
- [`create_log_file()`](http://r.bafg.de/reference/create_log_file.md) :
  Create a list for logging when processing with the Non-Target App
- [`dbGetSpectrum()`](http://r.bafg.de/reference/dbGetSpectrum.md) :
  Extraction of a MS2-Spectrum for a given ExperimentID from spectral
  library
- [`filter_matching_ms2()`](http://r.bafg.de/reference/filter_matching_ms2.md)
  : Compare MS2 spectra to a reference database and record the best
  match and similarity score.
- [`filter_matching_rt_mz()`](http://r.bafg.de/reference/filter_matching_rt_mz.md)
  : Match sample data by precursor mass and RT of csl compound within
  tolerances
- [`fixNoLossSpace()`](http://r.bafg.de/reference/fixNoLossSpace.md) :
  Fix Whitespace Error in GenForm Loss Annotation
- [`forIdentExport()`](http://r.bafg.de/reference/forIdentExport.md) :
  Export data to For-Ident.
- [`get_mass()`](http://r.bafg.de/reference/get_mass.md) : Calculate the
  mass/charge (m/z) from molecular formula and adduct form
- [`hello`](http://r.bafg.de/reference/hello.md) : A simple function
  doing little
- [`keepConsecutive()`](http://r.bafg.de/reference/keepConsecutive.md) :
  Keep only rows found in at least X consecutive samples
- [`keepReps()`](http://r.bafg.de/reference/keepReps.md) : Only retain
  peaks found in all replicates
- [`keep_reps_regex()`](http://r.bafg.de/reference/keep_reps_regex.md) :
  Replicate feature based on regular expression
- [`label_separate_peaks_by_rt()`](http://r.bafg.de/reference/label_separate_peaks_by_rt.md)
  : Group and label separate chromatographic peaks (A-Z) based on
  retention time, keeping only the most intense from each group.
- [`loadReport()`](http://r.bafg.de/reference/loadReport.md) : Load a
  previously saved Report
- [`mergeReport()`](http://r.bafg.de/reference/mergeReport.md) : Merge
  two Report objects
- [`ms2_search()`](http://r.bafg.de/reference/ms2_search.md) : Library
  screening without prior peak-picking
- [`new_peaklist()`](http://r.bafg.de/reference/new_peaklist.md) : Add
  additional peak and sample ID information to peaklists
- [`ntsworkflow`](http://r.bafg.de/reference/ntsworkflow.md) :
  ntsworkflow: A package for non-target-screening data processing.
- [`onlyGroupLeaders()`](http://r.bafg.de/reference/onlyGroupLeaders.md)
  : Keep only features which are the group leaders
- [`optimumMzStep()`](http://r.bafg.de/reference/optimumMzStep.md) : Get
  a recommendation for the m/z-step parameter (binning width)
- [`pickPeaksMzRange()`](http://r.bafg.de/reference/pickPeaksMzRange.md)
  : Peak-picking for an m/z range
- [`pickPeaksOneEic()`](http://r.bafg.de/reference/pickPeaksOneEic.md) :
  Peakpeaking algorithm using C++
- [`pickPeaksOneEicCpp()`](http://r.bafg.de/reference/pickPeaksOneEicCpp.md)
  : Find peaks in an ion chromatogram
- [`rcpparma_hello_world()`](http://r.bafg.de/reference/rcpparma_hello_world.md)
  [`rcpparma_outerproduct()`](http://r.bafg.de/reference/rcpparma_hello_world.md)
  [`rcpparma_innerproduct()`](http://r.bafg.de/reference/rcpparma_hello_world.md)
  [`rcpparma_bothproducts()`](http://r.bafg.de/reference/rcpparma_hello_world.md)
  : Set of functions in example RcppArmadillo package
- [`removeRare()`](http://r.bafg.de/reference/removeRare.md) : Remove
  features found in less than x samples
- [`runPeakPicking()`](http://r.bafg.de/reference/runPeakPicking.md) :
  Start the shiny peak picking app from the R console
- [`search_ms1_for_precursor()`](http://r.bafg.de/reference/search_ms1_for_precursor.md)
  : Search MS1 scans near the MS2 retention time for a precursor ion
  matching the expected m/z.
- [`summarize_groups()`](http://r.bafg.de/reference/summarize_groups.md)
  : Summarize componentization groups into one column
- [`wind()`](http://r.bafg.de/reference/wind.md) : Simple function to
  create window tolerance
