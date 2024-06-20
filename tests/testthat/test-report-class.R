test_that("dbas processing can find IS in test files", {
  test <- Report$new()
  test$changeSettings("pol", "pos")
  test$addDB(F, "/scratch/nts/MS2_db_v11.db")
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFiles(F, c(
    test_path("fixtures", "RH_pos_20230812.mzXML"),
    test_path("fixtures", "RH_pos_20230813.mzXML")
  ))
  test$process_all()
  expect_equal(nrow(test$peakList), 6)
  expect_contains(test$peakList$comp_name, "Olmesartan-d6")
  pth <- withr::local_tempfile(fileext = ".report")
  test$clearAndSave(F, pth)
  file.exists(pth)
  test2 <- loadReport(F, pth)
  expect_contains(test2$peakList$comp_name, "Olmesartan-d6")
})
