
test_that("Peak-picking for individual EICs of olmesartanD6TestFile returns 3 peaks", {
  mz_min <- 450
  mz_max <- 460
  mz_step <- 0.02
  s <- seq(mz_min, mz_max, by = mz_step*0.5)
  x <- lapply(s, function(ii) {
    pickPeaksOneEic(
      i = ii,
      rawData = olmesartanD6TestFile, 
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
})


test_that("Peak-picking for an mz range olmesartanD6TestFile runs", {
  
  f <- function(mi, ma) {
    pickPeaksMzRange(
      daten = olmesartanD6TestFile, 
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
  expect_true(is.data.frame(dfPeakList), info = "If no peaks found, returns df with 0 rows")
  expect_equal(nrow(dfPeakList), 0, info = "If no peaks found, returns df with 0 rows")
  
})

test_that("File with no peaks returns an emptly peaklist", {
  rd <- xcms::xcmsRaw(test_path("fixtures", "mzXML-test-files", "RH_pos_20220602_no_peaks.mzXML"), includeMSn = T)
  mi <- 100
  ma <- 1000
  mz_step <- 0.02
  
  x <- pickPeaksMzRange(
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
  
  expect_true(is.data.frame(x), info = "If no peaks found, returns df with 0 rows")
  expect_equal(nrow(x), 0, info = "If no peaks found, returns df with 0 rows")
  
})







