

# Prepare testReport.report object ####
library(testthat)
test <- Report$new()
test$changeSettings("pol", "pos")
test$addDB(F, test_path("fixtures", "CSL_olmesartan-d6.db"))
test$addIS(F, test_path("fixtures", "IS_table_pos.csv"))
test$addRawFiles(F, c(
  test_path("fixtures", "mzXML-test-files", "RH_pos_20230812.mzXML"),
  test_path("fixtures", "mzXML-test-files", "RH_pos_20230813.mzXML")
))
test$process_all()
test$reIntegrate()
test$clearAndSave(F, test_path("fixtures", "testReport.report"))
