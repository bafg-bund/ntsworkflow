


library(ntsworkflow)

test <- Report$new()

test$settings
test$changeSettings("pol", "pos")
test$addDB(F, "/scratch/nts/MS2_db_v9.db")
test$addIS(F, "/scratch/nts/IS_table_pos.csv")
test$addRawFiles(F, "/srv/cifs-mounts/g2/G/G2/HRMS/Messdaten/koblenz/wasser/2019/201904/pos/RH_pos_20190401.mzXML")
test$process_all()


View(test$peakList)
test$view()
