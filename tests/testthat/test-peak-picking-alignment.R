

test_that("peakpicking works for Olmesartan-d6 test file", {
  rd <- xcms::xcmsRaw(test_path("fixtures", "RH_pos_20230812.mzXML"), includeMSn = T)
  mz_min <- 450
  mz_max <- 460
  mz_step <- 0.02
  s <- seq(mz_min, mz_max, by = mz_step*0.5)
  x <- lapply(s, function(ii) {
    peakpicking_BfG_cpp(
      i = ii,
      rawData = rd, 
      mz_step = 0.02, 
      rt_min_scan = 1,
      rt_max_scan = 235, 
      sn = 3, 
      int_threshold = 10, 
      NoiseScans = 30, 
      peakwidth_min = 5, 
      peakwidth_max = 60,
      precursormzTol = 20, 
      maxPeaksPerSignal = 10
    )
  })
  y <- purrr::discard(x, is.null)
  expect_equal(length(y), 3, info = "There should be 3 peaks found between
               m/z 450 and 460, all other XICs should return NULL")
  z <- purrr::keep(x, is.null)
  expect_equal(length(z), 998, info = "There should be 3 peaks found between
               m/z 450 and 460, all other XICs should return NULL")
  f <- function(mi, ma) {
    FindPeaks_BfG(
      daten = rd, 
      mz_min = mi, 
      mz_max = ma, 
      mz_step = 0.02,
      rt_min = 120,
      rt_max = 1200,
      sn = 3,
      int_threshold = 10,
      peak_NoiseScans = 30,
      precursormzTol = 20,
      peakwidth_min = 5,
      peakwidth_max = 60,
      maxPeaksPerSignal = 10
    )
  }
  dfPeakList <- f(450, 460)
  expect_true(is.data.frame(dfPeakList))
  expect_equal(nrow(dfPeakList), 3)
  dfPeakList <- f(450, 451)
  expect_true(is.matrix(dfPeakList), info = "If no peaks found, returns logi 
              matrix with 0 rows")
  expect_equal(nrow(dfPeakList), 0, info = "If no peaks found, returns logi 
               matrix with 0 rows")
  
  
})