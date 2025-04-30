
# Test that peakpicking works by loading the file tests/shinyTests/fixtures/loadExampleDataRhine2016.RDS
# and testing run peak picking for all samples.

devtools::load_all()
runPeakPicking()
