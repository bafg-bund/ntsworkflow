test_that("ms2_search finds Olmesartan-d6 in sample files", {
  # Prepare paths for CSL file and sample files
  db_path <- testthat::test_path("fixtures", "CSL_olmesartan-d6.db")
  data_dir <- testthat::test_path("fixtures", "mzXML-test-files")
  mzxml_files <- list.files(data_dir, pattern = "\\.mzXML$", full.names = TRUE)
  
  # Run m2_search
  results <- ms2_search(
    data_path = mzxml_files,
    db_path = db_path
  )
  
  # Test return class
  expect_s3_class(results, "data.frame")

  # Test if the expected compound was found (twice)
  expect_true("Olmesartan-d6" %in% results$comp_name)
  expect_equal(length(results$comp_name), 2)
  
  # Test if no compounds were found for negative control
  path_expected_no_comp <- "RH_pos_20220602_no_peaks.mzXML"
  expect_false(path_expected_no_comp %in% results$samp)
})
