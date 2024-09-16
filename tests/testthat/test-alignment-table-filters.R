

test_that("filter replicates from alignement table", {
  g <- readRDS(test_path("fixtures", "alignment-table-filters", "keepReps-grouped.RDS"))
  expect_equal(nrow(g), 66)
  # Consider all samples expect blanks (see sampleList below)
  g2 <- keepReps(g, 3:8, 3, 3)
  expect_equal(nrow(g2), 35)
  
  # Consider only the first set of replicates
  g3 <- keepReps(g, 3:5, 3, 3)
  expect_equal(nrow(g3), 66)
  
  # Consider only the second set of replicates
  g4 <- keepReps(g, 6:8, 3, 3)
  expect_equal(nrow(g4), 66)
  # To get the group table see "A batch with blanks has compound annotations"
  # test in ntsportal
  
  # View(g2[, grepl("^Int", colnames(g2))])
  # View(g3[, grepl("^Int", colnames(g2))])
  # View(g4[, grepl("^Int", colnames(g2))])
  # sl <- test_path("fixtures", "alignment-table-filters", "keepReps-sampleList.RDS")
  # View(readRDS(sl))
})

test_that("Filter replicates from alignment based on regexp", {
  
  alignment <- readRDS(test_path("fixtures", "alignment-table-filters", "keepReps-grouped.RDS"))
  sl <- test_path("fixtures", "alignment-table-filters", "keepReps-sampleList.RDS")
  sampleList <- readRDS(sl)
  al <- keep_reps_regex(
    alignment = alignment, samples = 3:8, sampleList = sampleList, 
    regexp = "0[123](_pos)", least = 3
  )
  expect_equal(nrow(al), 35)
  
  expect_message(
    al2 <- keep_reps_regex(
      alignment = alignment, samples = 3:8, sampleList = sampleList, 
      regexp = "0[123](_pos)", least = 4
    )
  )
  expect_equal(nrow(al2), 35)
  
  
  # View(al[, grepl("^Int", colnames(al))])
  
  # Get test data
  # alignment <- readRDS(test_path("fixtures", "alignment-table-filters", "keepReps-grouped.RDS"))
  # sl <- test_path("fixtures", "alignment-table-filters", "keepReps-sampleList.RDS")
  # sampleList <- readRDS(sl)
  # regexp <- "0[123](_pos)"
  # samples <- 3:8
})