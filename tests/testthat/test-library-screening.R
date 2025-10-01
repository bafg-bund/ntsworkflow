# Integration tests ####

test_that("ms2_search finds Olmesartan-d6 in sample files", {
  # Prepare paths for CSL file and sample files
  csl_path <- testthat::test_path("fixtures", "CSL_olmesartan-d6.db")
  data_dir <- testthat::test_path("fixtures", "mzXML-test-files")
  mzxml_files <- list.files(data_dir, pattern = "\\.mzXML$", full.names = TRUE)
  
  results <- ms2_search(
    data_path = mzxml_files,
    csl_path = csl_path
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

test_that("Multiple peaks are labelled by letters according to RT", {
  csl_path <- test_path("fixtures", "csl-test-files", "CSL_v25.6.db")
  data_dir <- test_path("fixtures", "mzXML-test-files", "lib-screening-multi-peaks")
  mzxml_files <- list.files(data_dir, pattern = "\\.mzXML$", full.names = TRUE)
  
  results <- ms2_search(
    data_path = mzxml_files,
    csl_path = csl_path,
    rttolm = 2,
    compounds = "Tebuconazole"
  )
  expect_equal(nrow(results), 2)
  expect_equal(results$peak, c("A", "B"))
  expect_lt(results[results$peak == "A", "rt"], results[results$peak == "B", "rt"])
})

test_that("ms2_search finds two out of three available compounds in sample files", {
  # Prepare paths for CSL file and sample files
  csl_path <- testthat::test_path("fixtures", "CSL_olm-d6_beza-d4_2Ctp.db")
  data_dir <- testthat::test_path("fixtures", "mzXML-test-files")
  mzxml_files <- list.files(data_dir, pattern = "\\.mzXML$", full.names = TRUE)
  
  # Run m2_search
  results <- ms2_search(
    data_path = mzxml_files,
    csl_path = csl_path
  )
  
  # Test return class
  expect_s3_class(results, "data.frame")
  
  # Test if two out of three compounds were found
  expect_true("Olmesartan-d6" %in% results$comp_name)
  expect_true("(2-Carboxyethyl)triphenylphosphonium" %in% results$comp_name)
  expect_false("Bezafibrate-d4" %in% results$comp_name)
  
  # Test if the expected compounds were found twice each
  expect_equal(length(results$comp_name), 4)  
  
  # Test if no compounds were found for the negative control
  path_expected_no_comp <- "RH_pos_20220602_no_peaks.mzXML"
  expect_false(path_expected_no_comp %in% results$samp)
})


# Unit tests ####


test_that("get_csl_data returns csl data correctly", {
  # Path to mock test database
  test_csl_path <- testthat::test_path("fixtures", "CSL_olm-d6_beza-d4_2Ctp.db")
  
  # Prepare test inputs
  chromatography <- "bfg_nts_rp1"
  pol <- "pos"
  CE_s <- 30:40
  CES_s <- 0:15
  instr <- "LC-ESI-QTOF TripleTOF 5600 SCIEX"
  compounds <- NULL  # Don't filter by compound
  
  # Run function
  csl_data <- get_csl_data(
    csl_path = test_csl_path,
    chromatography = chromatography,
    pol = pol,
    CE_s = CE_s,
    CES_s = CES_s,
    instr = instr,
    compounds = compounds
  )
  
  # Assert that output type is correct and expected tables exist
  expect_type(csl_data, "list")
  expect_named(csl_data, c("suspects", "exptbl", "comptbl", "fragtbl"))
  
  # Assert that suspects are as expected
  expect_true("Olmesartan-d6" %in% csl_data$suspects$name)
  expect_true("(2-Carboxyethyl)triphenylphosphonium" %in% csl_data$suspects$name)
  expect_false("Bezafibrate-d4" %in% csl_data$suspects$name)
})


test_that("get_csl_data filters suspects correctly", {
  # Path to mock test database
  test_csl_path <- testthat::test_path("fixtures", "CSL_olm-d6_beza-d4_2Ctp.db")
  
  # Prepare test inputs
  chromatography <- "bfg_nts_rp1"
  pol <- "pos"
  CE_s <- 30:40
  CES_s <- 0:15
  instr <- "LC-ESI-QTOF TripleTOF 5600 SCIEX"
  compounds <- "(2-Carboxyethyl)triphenylphosphonium"  # Compound filter
  
  # Run function
  csl_data <- get_csl_data(
    csl_path = test_csl_path,
    chromatography = chromatography,
    pol = pol,
    CE_s = CE_s,
    CES_s = CES_s,
    instr = instr,
    compounds = compounds
  )
  
  # Assert that suspects are as expected
  expect_true("(2-Carboxyethyl)triphenylphosphonium" %in% csl_data$suspects$name)
  expect_false("Olmesartan-d6" %in% csl_data$suspects$name)
  expect_false("Bezafibrate-d4" %in% csl_data$suspects$name)
})

test_that("filter_matching_rt_mz returns correct number of entries for 
            different mz/RT combinations", {
  mock_lib_entries <- tribble(
    ~csl_comp_rt, ~csl_comp_mz, ~expected_nrow,
    6.7,          335,          1,
    1.5,          335,          0,
    8.3,          200,          0,
    3.3,          100,          2,
    NA,           100,          2
  )
  mock_peak_data <- data.frame(
    mz = c(100.0, 100.1, 450.5, 335.4),
    rt = c(3.5, 3.4, 1.5, 6.6) * 60,
    index = 1:4
  )
  test_mock_entries <- function(lib_entry) {
    result <- filter_matching_rt_mz(mock_peak_data, lib_entry["csl_comp_rt"], lib_entry["csl_comp_mz"], mztolu = 0.5, rtoffset = 0, rttolm = 1)
    expect_equal(nrow(result), unname(lib_entry["expected_nrow"]))  
  }
  apply(mock_lib_entries, 1, test_mock_entries)
})


