


library(ntsworkflow)

test <- Report$new()

test$settings
test$changeSettings("pol", "neg")
test$addDB(F, "/scratch/nts/MS2_db_v11.db")
test$addIS(F, "/scratch/nts/IS_table_neg.csv")
test$addRawFiles(F, "/srv/cifs-mounts/g2/G/G2/HRMS/Messdaten/koblenz/wasser/2018/201803/neg/RH_neg_MQ201803_IS.mzXML")
test$process_all()
test$reIntegrate()

View(test$peakList)
View(test$integRes)
test$view()
