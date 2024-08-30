test_that("dbas processing can find IS in test files, these can be saved without
          removing files from memory", {
  test <- Report$new()
  test$changeSettings("pol", "pos")
  test$addDB(F, test_path("fixtures", "CSL_olmesartan-d6.db"))
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFiles(F, c(
    test_path("fixtures", "RH_pos_20230812.mzXML"),
    test_path("fixtures", "RH_pos_20230813.mzXML")
  ))
  test$process_all()
  expect_equal(nrow(test$peakList), 2)
  expect_contains(test$peakList$comp_name, "Olmesartan-d6")
  pth <- withr::local_tempfile(fileext = ".report")
  test$clearAndSave(F, pth, clearData = F)
  expect_true(file.exists(pth))
  expect_s4_class(test$rawData[[1]], "xcmsRaw") 
  test2 <- loadReport(F, pth)
  expect_contains(test2$peakList$comp_name, "Olmesartan-d6")
  file.remove(pth)
})

test_that("Processing an empty file returns an empty report", {
  test <- Report$new()
  test$changeSettings("pol", "pos")
  test$addDB(F, test_path("fixtures", "CSL_olmesartan-d6.db"))
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFiles(F, test_path("fixtures", "RH_pos_20220602_no_peaks.mzXML"))
  test$process_all()
  expect_equal(nrow(test$peakList), 0)
  expect_equal(nrow(test$ISresults), 0)
})

