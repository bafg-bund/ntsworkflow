test_that("annotate_grouped returns annotation table with expected compound", {
  # Prepare inputs
  csl_path <- testthat::test_path("fixtures", "csl-test-files", "CSL_v25.7.db")
  sampleList <- readRDS(testthat::test_path("fixtures", "annotation", "sampleListLocal.RDS"))
  peaklist <- readRDS(testthat::test_path("fixtures", "annotation", "peakListList.RDS"))
  grouped <- readRDS(testthat::test_path("fixtures", "annotation", "alignmentTable.RDS"))
  datenList <- readRDS(testthat::test_path("fixtures", "annotation", "datenListLocal.RDS"))
    
  # Run the function
  annotationTable <- annotate_grouped(
    sampleListLocal = sampleList,
    peakListList = peaklist,
    alignmentTable = grouped,
    db_path = csl_path,
    datenListLocal = datenList
  )
  
  # Test return class
  expect_s3_class(annotationTable, "data.frame")
  
  # Test if the expected compound was found
  expect_true("Olmesartan-d6" %in% annotationTable$name)
})