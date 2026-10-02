

# Prepare testReport.report object ####
library(testthat)
devtools::load_all()
test <- Report$new()
test$changeSettings("pol", "pos")
test$addDB(F, test_csl_path)
test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
test$addRawFiles(F, c(
  test_path("fixtures", "mzXML-test-files", "RH_pos_20230812.mzXML"),
  test_path("fixtures", "mzXML-test-files", "RH_pos_20230813.mzXML")
))
test$process_all(comp_names = "Olmesartan-d6")
test$reIntegrate()
test$clearAndSave(F, test_path("fixtures", "testReport.report"))
