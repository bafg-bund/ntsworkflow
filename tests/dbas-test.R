


library(ntsworkflow)

test <- Report$new()

test$settings
test$changeSettings("pol", "pos")
test$addDB(F, "/scratch/nts/MS2_db_v11.db")
test$addIS(F, "/scratch/nts/IS_table_pos.csv")
test$addRawFiles(F, "/srv/cifs-mounts/g2/G/G2/4-Projekte_G2/4.3-Drittmittel-Projekte_G2/4.3.3-NTS_HessenIII/4.3.3.4-Daten/20230801_Hessen_III/mzXML/pos/040-31-05-2023_02-117-SP_Gersprenz_pos.mzXML")
test$process_all()
test$reIntegrate()

View(test$peakList)
View(test$integRes)
test$view()

runPeakPicking()
