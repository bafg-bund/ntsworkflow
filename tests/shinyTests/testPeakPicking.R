
# 1. Test current code ####
# 1a. Set cores to 1 and "Add sample" the files in tests/testthat/fixtures/mzXML-test-files
# 1b. Set m/z range to 450-460 and run peak picking for all files
# 1c. Align all files with standard settings
# 1d. Annotate files using CSL in tests/testthat/fixture/CSL_olmesartan-d6. Look for annotation of olmesartan-d6

devtools::load_all()
runPeakPicking()

rm(list = ls())

# 2. Test installed code ####
# 2a. Install ntsworkflow 
# 2b. Run tests 1a-1d again.
library(ntsworkflow)
runPeakPicking()
