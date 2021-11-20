## ----setup, include = FALSE---------------------------------------------------
library(ggplot2)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(ntsworkflow)

## -----------------------------------------------------------------------------
test <- Report$new()

## ---- eval = FALSE------------------------------------------------------------
#  test$addRawFiles()

## -----------------------------------------------------------------------------
listOfFiles <- c(
    "example1.mzXML",
    "example2.mzXML"
  )
test$addRawFiles(dialog = FALSE, normalizePath(listOfFiles))

## -----------------------------------------------------------------------------
test$rawFiles

## -----------------------------------------------------------------------------
# With prompt: test$addIS()
test$addIS(dialog = FALSE, "IS_table_pos.csv")

## -----------------------------------------------------------------------------
test$IS

## ---- eval=FALSE--------------------------------------------------------------
#  test$settings

## ---- results='asis', echo=FALSE----------------------------------------------
values <- sapply(test$settings, function(x) if (length(x) > 1) Reduce(paste, x) else paste(x))
values[1] <- "D:/Example/sqlite/MS2db.db"
knitr::kable(data.frame(Name = names(test$settings), Value = values), row.names = FALSE,
             caption = "Settings used by the suspect search algorithm")

## -----------------------------------------------------------------------------
test$addDB(dialog = FALSE, "example_database.db")

## ---- message=FALSE, warning=FALSE--------------------------------------------
test$process_all()

## ---- results="asis", echo=FALSE----------------------------------------------
knitr::kable(test$peakList[1:2, c(1:2,22,12:13)], caption="Excerpt of the peak list")

## ---- results="asis", echo=FALSE----------------------------------------------
knitr::kable(test$ISresults[, 1:5], caption = "IS integration results")

## ---- fig.cap="EIC plot from peak with ID 21, lamotrigine in this case.", fig.width=5, fig.height=3----
test$plotEIC(1)  # you must provide the peakID, which is found in the peak list

## ---- fig.cap="MS^1^ plot from peak with ID 21, lamotrigine in this case.", fig.width=5, fig.height=3----
test$plotMS1(1)

## ---- fig.cap="MS^2^ plot from peak with ID 21, lamotrigine in this case.", fig.width=5, fig.height=3----
test$plotMS2(1)

## ---- eval=FALSE--------------------------------------------------------------
#  test$deleteFP("Benzyl-dimethyl-hexadecylammonium")

## ---- eval=FALSE--------------------------------------------------------------
#  test$falsePos <- test$falsePos[test$falsePos != "Guanosine"]

## ---- eval = FALSE------------------------------------------------------------
#  test$deleteBackground(3, 2)  # will not delete IS peaks
#  test$deleteBackground(3, 2, includeIS = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  test$peakList <- test$peakList[!is.na(test$peakList$int_a), ]

## ---- eval=FALSE--------------------------------------------------------------
#  test$reIntegrate()

## ---- eval=FALSE--------------------------------------------------------------
#  test$saveSettings()

## ---- eval=FALSE--------------------------------------------------------------
#  test$clearAndSave()
#  
#  # Or, if you want to give it another name or location
#  test$clearAndSave("D:\\testFolder\\example_test")

## ---- eval=FALSE--------------------------------------------------------------
#  newName <- ntsworkflow::loadReport("example_test.report")

## ---- eval=FALSE--------------------------------------------------------------
#  write.csv(test$peakList, file = "new_peaklist.csv")

