test_that("dbas processing can find IS in test files, these can be saved without
          removing files from memory", {
  test <- Report$new()
  test$changeSettings("pol", "pos")
  test$addDB(F, test_csl_path)
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFiles(F, c(
    test_path("fixtures", "mzXML-test-files", "RH_pos_20230812.mzXML"),
    test_path("fixtures", "mzXML-test-files", "RH_pos_20230813.mzXML")
  ))
  test$process_all(comp_names = "Olmesartan-d6")
  expect_no_warning(test$plotMS2(1))
  expect_equal(nrow(test$peakList), 2)
  expect_contains(test$peakList$comp_name, "Olmesartan-d6")
  test$reIntegrate()
  expect_equal(test$integRes$samp[1], "RH_pos_20230812.mzXML")
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
  test$addDB(F, test_csl_path)
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFiles(F, test_path("fixtures", "mzXML-test-files", "RH_pos_20220602_no_peaks.mzXML"))
  test$process_all()
  expect_equal(nrow(test$peakList), 0)
  expect_equal(nrow(test$ISresults), 0)
})

test_that("False positive in multiple indicies can be deleted", {
  test <- loadReport(F, test_path("fixtures", "testReport.report"))
  test$deleteFP("Olmesartan-d6")
  expect_equal(nrow(test$peakList), 0)
  expect_equal(nrow(test$integRes), 0)
  expect_equal(nrow(test$MS2), 0)
  expect_equal(nrow(test$MS1), 0)
  expect_equal(nrow(test$EIC), 0)
  
  test <- loadReport(F, test_path("fixtures", "testReport.report"))
  test$deleteFP("Olmesartan-d6", 1)
  expect_equal(test$peakList$samp, "RH_pos_20230813.mzXML")
  expect_equal(nrow(test$peakList), 1)
  expect_equal(nrow(test$integRes), 1)
  expect_gt(nrow(test$MS2), 1)
  expect_gt(nrow(test$MS1), 1)
  expect_gt(nrow(test$EIC), 1)
  
  test <- loadReport(F, test_path("fixtures", "testReport.report"))
  test$deleteFP("Olmesartan-d6", 2)
  expect_equal(test$peakList$samp, "RH_pos_20230812.mzXML")
  expect_equal(nrow(test$peakList), 1)
  expect_equal(nrow(test$integRes), 1)
  expect_gt(nrow(test$MS2), 1)
  expect_gt(nrow(test$MS1), 1)
  expect_gt(nrow(test$EIC), 1)
  
  test <- loadReport(F, test_path("fixtures", "testReport.report"))
  test$deleteFP("Olmesartan-d6", c(1,2))
  expect_equal(nrow(test$peakList), 0)
  expect_equal(nrow(test$integRes), 0)
  expect_equal(nrow(test$MS2), 0)
  expect_equal(nrow(test$MS1), 0)
  expect_equal(nrow(test$EIC), 0)
})

test_that("Reintegration runs without error", {
  # Prepare test files for reintegration
  test <- Report$new()
  test$changeSettings("pol", "pos")
  test$addDB(F, test_csl_path)
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFilesDir(F, test_path("fixtures", "mzXML-test-files"))
  test$process_all(comp_names = "Olmesartan-d6")
  test$reIntegrate()
  expect_equal(nrow(test$integRes), 2)
})


test_that("Reintegration runs even with strange file names", {
  # Prepare test files with strange names
  tDir <- tempSaveDir <- withr::local_tempdir() 
  oFiles <- list.files(test_path("fixtures", "mzXML-test-files"), pattern = "*.mzXML", full.names = TRUE)
  nFiles <- file.path(tDir, paste0("foobar.!, +^#_", basename(oFiles)))
  file.copy(oFiles, nFiles)
  
  # Prepare test files for reintegration
  test <- Report$new()
  test$changeSettings("pol", "pos")
  test$addDB(F, test_csl_path)
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFiles(F, nFiles)
  suppressMessages(test$process_all(comp_names = "Olmesartan-d6"))
  expect_no_error(test$reIntegrate())
})

test_that("Blank correction also removes peaks from the reintegration table", {
  testFiles <- c(
    test_path("fixtures", "mzXML-test-files", "RH_pos_20230812.mzXML"),
    test_path("fixtures", "mzXML-test-files", "RH_pos_20230813.mzXML")
  )
  test <- Report$new()
  test$changeSettings("pol", "pos")
  test$addDB(F, test_csl_path)
  test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
  test$addRawFiles(F, testFiles)
  suppressMessages(test$process_all(comp_names = "Olmesartan-d6"))
  suppressMessages(test$reIntegrate())
  expect_equal(nrow(test$peakList), 2)
  expect_equal(nrow(test$integRes), 2)
  test$deleteBackground(2, 1)
  expect_equal(nrow(test$peakList), 0)
  expect_equal(nrow(test$integRes), 0)
  
  # Detection only in integRes (gap-filled) and not in peaklist is removed
  suppressMessages(test$reprocess(comp_names = "Olmesartan-d6"))
  suppressMessages(test$reIntegrate())
  newIntegRes <- test$integRes[c(1, 2, 2), ]
  newIntegRes[3, "samp"] <- "foobar"
  test$rawFiles <- c(test$rawFiles, "foobar")
  test$integRes <- newIntegRes
  test$deleteBackground(2:3, 1)
  expect_equal(nrow(test$integRes), 0)
  expect_equal(nrow(test$peakList), 0)
  test$rawFiles <- test$rawFiles[1:2]
  
  # Comp found only in blank in peaklist and found in blank and sample in integRes is deleted from both
  suppressMessages(test$reprocess(comp_names = "Olmesartan-d6"))
  suppressMessages(test$reIntegrate())
  newPeakList <- test$peakList
  newIntegRes <- test$integRes
  newPeakList <- newPeakList[c(1, 1, 2), ]
  newPeakList[1, "peakID"] <- 4
  newPeakList[1, "comp_name"] <- "foobar"
  newIntegRes <- newIntegRes[c(1,1,2,2), ]
  newIntegRes[1, "comp_name"] <- "foobar"
  newIntegRes[4, "comp_name"] <- "foobar"
  test$integRes <- newIntegRes
  test$peakList <- newPeakList
  test$deleteBackground(2, 1)
  expect_equal(nrow(test$integRes), 0)
  expect_equal(nrow(test$peakList), 0)
  
  # Detection with area much higher in sample than in blank is kept
  suppressMessages(test$reprocess(comp_names = "Olmesartan-d6"))
  suppressMessages(test$reIntegrate())
  newPeakList <- test$peakList
  newIntegRes <- test$integRes
  newPeakList[2, "int_a"] <- 1e5
  newPeakList[2, "int_h"] <- 1e5
  newIntegRes[2, "int_a"] <- 1e5
  newIntegRes[2, "int_a"] <- 1e5
  test$integRes <- newIntegRes
  test$peakList <- newPeakList
  test$deleteBackground(2, 1)
  expect_equal(nrow(test$peakList), 1)
  expect_equal(nrow(test$integRes), 1)
  
  # no detections in blank, no peak removed.
  suppressMessages(test$reprocess(comp_names = "Olmesartan-d6"))
  suppressMessages(test$reIntegrate())
  newPeakList <- test$peakList[2, ]
  newIntegRes <- test$integRes[2, ]
  test$integRes <- newIntegRes
  test$peakList <- newPeakList
  suppressMessages(test$deleteBackground(2, 1))
  expect_equal(nrow(test$peakList), 1)
  expect_equal(nrow(test$integRes), 1)
})

test_that("mzML as well as mzXML files can be added", {
  test <- Report$new()
  test$addRawFiles(F, c(
    test_path("fixtures", "mzXML-test-files", "RH_pos_20230812.mzXML"),
    test_path("fixtures", "mzXML-test-files", "RH_pos_20230813.mzXML")
  ))
  expect_length(test$rawFiles, 2)
  
  test <- Report$new()
  test$addRawFilesDir(F, test_path("fixtures", "mzXML-test-files"))
  expect_length(test$rawFiles, 3)
  
  test <- Report$new()
  test$addRawFiles(F, list.files(test_path("fixtures", "mzML-test-files"), full.names = TRUE))
  expect_length(test$rawFiles, 2)
  
  test <- Report$new()
  test$addRawFilesDir(F, test_path("fixtures", "mzML-test-files"))
  expect_length(test$rawFiles, 2)
})

test_that("mzML files can be processed", {
  test <- Report$new()
  test$addRawFilesDir(F, test_path("fixtures", "mzML-test-files"))
  test$addDB(F, test_csl_path)
  test$addIS(F, test_path("fixtures", "istd-tables-dbas", "bfg_nts_rp1_neochema_is_table_pos.csv"))
  suppressMessages(
    test$process_all(comp_names = "Bezafibrate-d6")
  )
  expect_contains(test$peakList$comp_name, "Bezafibrate-d6")
  expect_contains(test$ISresults$IS, "Bezafibrate-d6")
  expect_no_error(suppressMessages(
    test$reIntegrate()
  ))
})
  