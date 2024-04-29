

library(ntsworkflow)

test <- loadReport()
test$peakList


debug(ms2_search)
test$reprocess(1, "Bezafibrate-d4")



test2 <- Report(test)
test2$reprocess(1, "Bezafibrate-d4")
