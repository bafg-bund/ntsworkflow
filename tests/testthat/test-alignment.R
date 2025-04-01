

# Testing alignment ####
test_that("Alignment does not produce out-of-bounds warnings", {
  
  plPath <- test_path("fixtures", "peaklist_alignment_peak-picking-alignment3.RDS")
  pl <- readRDS(plPath)
  expect_no_warning(x <- alignment_BfG_cpp(pl, 5, 20, "mDa"))
  expect_equal(nrow(x), 289)
  expect_equal(ncol(x), 52)
})