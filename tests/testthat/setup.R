
# Since opening the file takes time, only do it once
olmesartanD6TestFile <- xcms::xcmsRaw(test_path("fixtures", "mzXML-test-files", "RH_pos_20230812.mzXML"), includeMSn = T)
