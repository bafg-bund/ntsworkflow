


library(ntsworkflow)

test <- Report$new()

test$settings
test$changeSettings("pol", "pos")
test$addDB(F, "/scratch/nts/MS2_db_v9.db")
test$addIS(F, "/scratch/nts/IS_table_pos.csv")
test$addRawFiles(F, c(
  "/srv/cifs-mounts/g2/G/G2/HRMS/Messdaten/koblenz/wasser/2019/201904/pos/RH_pos_20190401.mzXML",
  "/srv/cifs-mounts/g2/G/G2/HRMS/Messdaten/koblenz/wasser/2019/201904/pos/RH_pos_20190402.mzXML"))
test$process_all()
test$reIntegrate()

View(test$peakList)
View(test$integRes)
test$view()
