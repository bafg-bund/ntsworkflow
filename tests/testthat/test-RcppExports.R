test_that("Pick-picking cpp functions returns maxima without crashing for problematic XIC", {
  # Load the fixture data
  testXic <- readRDS(test_path("fixtures", "xic_352.33_2019mq1_RcppExports1.RDS"))
    
  testScantime <- readRDS(test_path("fixtures", "scantime_2019mq1_RcppExports1.RDS.RDS"))
  
  # Call the function
  maxima <- pickPeaksOneEicCpp(
    mz = 352.31, 
    mz_step = 0.02, 
    eic = testXic, 
    scantime = testScantime, 
    minIntensity = 1, 
    sn = 1, 
    noisescans = 30, 
    peakwidth_min = 5, 
    peakwidth_max = 60, 
    maxPeaksPerSignal = 10
  )
  
  # Expectations 
  expect_true(is.matrix(maxima)) 
  expect_equal(nrow(maxima),1)
  expect_equal(ncol(maxima),16)
  expect_equal(maxima[1,5], 1370)

})

