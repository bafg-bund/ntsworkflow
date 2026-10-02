

test_that("Merging a report with no IS results does not result in a warning", {
  tempDir <- withr::local_tempdir()

  # Prepare IS table with incorrect RT
  newIsTable <- file.path(tempDir, "temp_is_table_pos.csv")
  file.copy(test_path("fixtures", "IS_table_pos.csv"), newIsTable)
  readLines(newIsTable)
  system(glue::glue("sed -i '3s/7/1/' {newIsTable}"))
  
  test <- Report$new()
  test$changeSettings("pol", "pos")
  test$addDB(F, test_csl_path())
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFiles(F, c(
    test_path("fixtures", "mzXML-test-files", "RH_pos_20230812.mzXML")
  ))
  test$process_all(comp_names = "Olmesartan-d6")
  test$clearAndSave(F, file.path(tempDir, "test1.report"))
  
  test2 <- Report$new()
  test2$changeSettings("pol", "pos")
  test2$addDB(F, test_csl_path())
  test2$addIS(F, newIsTable)
  test2$addRawFiles(F, c(
    test_path("fixtures", "mzXML-test-files", "RH_pos_20230813.mzXML")
  ))
  test2$process_all(comp_names = "Olmesartan-d6")
  # replace IS table to match test1
  test2$IS <- data.frame()
  test2$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test2$clearAndSave(F, file.path(tempDir, "test2.report"))
  
  expect_no_warning(mergedReport <- mergeReport(test, test2))
  numIsTotal <- nrow(mergedReport$ISresults)
  numPeaksTotal <- nrow(mergedReport$peakList)
  expect_equal(numIsTotal, 1)
  expect_equal(numPeaksTotal, 2)
  
  expect_no_warning(mergedReportReversed <- mergeReport(test2, test))
  numIsTotalReversed <- nrow(mergedReportReversed$ISresults)
  numPeaksTotalReversed <- nrow(mergedReportReversed$peakList)
  expect_equal(numIsTotalReversed, 1)
  expect_equal(numPeaksTotalReversed, 2)
  expect_true(all.equal(mergedReport$ISresults, mergedReportReversed$ISresults))
  unlink(tempDir, recursive = TRUE)
})
