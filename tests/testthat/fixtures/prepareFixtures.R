

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


# Prepare mzML testfiles with current MSconvert version
# only works on windows
system('"C:\\Users\\Jewell\\AppData\\Local\\Apps\\ProteoWizard 3.0.26183.2740ae3 64-bit\\msconvert.exe" --zlib --numpressLinear --numpressPic --filter "peakPicking vendor msLevel=1-" --filter "scanTime [666,690]" --filter "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState> File:"""^<SourcePath^>""", NativeID:"""^<Id^>""""  --outdir "C:\\Users\\Jewell\\Documents\\rprojects\\ntsworkflow\\tests\\testthat\\fixtures\\mzML-test-files" "Z:\\Projekte\\ntsportal\\intern\\msrawfiles\\bfg_lab\\ringversuch\\fants_2026\\tof\\pos\\2026-08-30_sample10_rhein_spike_tof_injrep1_pos.wiff" "Z:\\Projekte\\ntsportal\\intern\\msrawfiles\\bfg_lab\\ringversuch\\fants_2026\\tof\\pos\\2026-08-30_sample10_rhein_spike_tof_injrep2_pos.wiff"')
