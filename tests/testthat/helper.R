

# Since opening the file takes time, only do it once
olmesartanD6TestFile <- xcms::xcmsRaw(test_path("fixtures", "mzXML-test-files", "RH_pos_20230812.mzXML"), includeMSn = T)

test_csl_path <- test_path("fixtures", "csl-test-file", "CSL_v26.2_monde.db")

