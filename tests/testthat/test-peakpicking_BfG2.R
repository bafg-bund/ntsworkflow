test_that("peakPickingBfGC returns maxima without crashing for a given XIC at mz: 252.33", {
  # Load the fixture data
  XIC_252 <- readRDS("~/projects/ntsworkflow/tests/testthat/fixture/XIC_252.rds")
  rawData_scantime <- readRDS("~/projects/ntsworkflow/tests/testthat/fixture/rawData_scantime.rds")
  
  # Call the function
  maxima <- peakPickingBfGC(
    mz = 352.31, 
    mz_step = 0.02, 
    XIC = XIC_252, 
    scantime = rawData_scantime, 
    min_intensity = 1, 
    sn = 1, 
    noisescans = 30, 
    peakwidth_min = 5, 
    peakwidth_max = 60, 
    maxPeaksPerSignal = 10
  )
  
  # Expectations 
  expect_type(maxima.is.matrix()) 

})