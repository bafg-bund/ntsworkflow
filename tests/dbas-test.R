


library(ntsworkflow)

test <- Report$new()

test$settings
test$changeSettings("pol", "pos")
test$addDB(F, "/scratch/nts/MS2_db_v9.db")
test$addIS(F, "~/projects/ntsautoeval/IS_table_pos.csv")
test$addRawFiles(F, "/home/Jewell/HRMS_Z/Messdaten/koblenz/wasser//2018/201802/pos/RH_pos_20180226.mzXML")
test$process_all()

test$peakList
test$view()
