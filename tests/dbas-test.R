


library(ntsworkflow)

test <- Report$new()

test$settings
test$changeSettings("pol", "pos")
test$addDB(F, "~/sqlite_local/MS2_db_v9.db")
test$addIS(F, "~/projects/ntsautoeval/IS_table_pos.csv")
test$addRawFiles(F, "~/messdaten/bimmen/schwebstoff/rhb_pos/B18_pos2.mzXML")
test$process_all(comp_names = "Avobenzone")

test$peakList
test$view()
