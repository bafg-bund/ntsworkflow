

library(shiny)
library(shinyBS)
library(DT)
library(tcltk2)
library(shinyFiles)
batchFiles <<- NULL
batchFilesRAM <<- NULL
batchFilesSampleType <<- NULL
sampleTypes <<- c("Unknown","Blank","Standard")
#globalwd <- getwd()
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width = 2,
      actionButton("BatchProcess","Batch Process"),
      bsModal(
        "BatchProcessPopup", "Batch Process", "BatchProcess", size = "large", 
        DT::dataTableOutput("BatchProcessTable"),
        shinyFilesButton("BatchProcessAddFiles","Add Files", "Please select a data file", TRUE),
        actionButton("BatchProcessRemoveFiles","Remove Files"),
        actionButton("BatchProcessStart","Go"),
        fluidRow(column(6, checkboxInput("batchBlank", "Pick blanks at 1/x intensity threshold and 1/y sn")),
                 column(1, h5("x:", align = "right")),
                 column(2, numericInput("batchBlankIntFactor", NULL, 10, 1, NA, step = 1, "100%")),
                 column(1, h5("y:", align = "right")),
                 column(2, numericInput("batchBlankSnFactor", NULL, 2, 1, NA, step = 1, "100%"))),
        checkboxInput("batchAlign", "Align peaks after peak-picking"),
        fluidRow(column(6, checkboxInput("saveBatch", "Save after processing")),
                 column(6, textInput("saveBatchName", NULL, 
                                     format(lubridate::now(), "%Y%m%d"),"100%", 
                                     "Name (default: date)")))
        
      )
    ),
    mainPanel(
      fluidRow(
        column(
          6, selectInput("sampleType", "Sample type:", sampleTypes, sampleTypes[1])
        ),
        column(6, textInput("blankRegex", "Blank regex:", 
                            placeholder = "e.g., '_MQ_' or '_blank_'"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Batch add files ####
  home <- c(WD = globalwd, Home = fs::path_home(), "R Installation" = R.home(), getVolumes("(E:)")())
  home2 <- c(Home = fs::path_home())
  
  # For some reason on server the batch process add sample drop down does not work, 
  # the problem of launching shinyFiles from within a modal is known. The issue was evidently fixed
  # but not here
  # This is a problem with firefox. Chrome apparently does not have this issue
  # using home2 so that home appears first
  shinyFileChoose(input, 'BatchProcessAddFiles', roots=home2, session=session,
                  filetypes=c('mzML', 'mzXML'))
  
  observeEvent(input$BatchProcessAddFiles, {
    req(is.list(input$BatchProcessAddFiles))
    
    dateiInfo <- parseFilePaths(home, input$BatchProcessAddFiles)
    dateiInfo <- as.data.frame(lapply(dateiInfo, as.character), stringsAsFactors = FALSE)
    additionalFiles <- dateiInfo$datapath
    if (length(additionalFiles) > 0) 
      pfad <<- dirname(additionalFiles[1])
    batchFiles <<- c(batchFiles,additionalFiles)
    batchFilesRAM <<- c(batchFilesRAM,rep(FALSE,length(additionalFiles)))
    
    # Automatically detect blanks and mark these in the table as such
    types <- if (input$blankRegex == "") {
      rep(1, length(additionalFiles)) 
    } else {
      ifelse(grepl(input$blankRegex, basename(additionalFiles)), 2, 1)
    }
    # mark files with sample types
    batchFilesSampleType <<- c(batchFilesSampleType, types)
    
    batchTable <- cbind(
      dirname(batchFiles),
      basename(batchFiles),
      round(file.size(batchFiles)/1000000,1),
      sampleTypes[batchFilesSampleType],
      batchFilesRAM
    )
    colnames(batchTable) <- c("Dir","File","Size","SampleType","RAM")
    output$BatchProcessTable <- DT::renderDataTable(
      DT::datatable(
        batchTable,
        selection=list(target="cell"),
        options = list(
          columnDefs = list(
            list(
              targets = 0,
              render = JS("function(data, type, row, meta) {","return type === 'display' && data.length > 10 ?","'<span title=\"' + data + '\">' + data.substr(0, 7) + '...</span>' : data;","}"))
          )
        ),
        callback = JS('table.page("next").draw(false);')
      )
    )
  })
  
  observeEvent(input$BatchProcessTable_cell_clicked, {
    selectedBatchFiles <<- input$BatchProcessTable_cells_selected
    if (length(selectedBatchFiles) > 0) {
      if (any(selectedBatchFiles[,2] == 4)) {
        batchFilesRAM[selectedBatchFiles[which(selectedBatchFiles[,2] == 4),1]] <<- !batchFilesRAM[selectedBatchFiles[which(selectedBatchFiles[,2] == 4),1]]
        batchTable <- cbind(dirname(batchFiles),basename(batchFiles),round(file.size(batchFiles)/1000000,1),sampleTypes[batchFilesSampleType],batchFilesRAM)
        colnames(batchTable) <- c("Dir","File","Size","SampleType","RAM")
        output$BatchProcessTable <- DT::renderDataTable(
          DT::datatable(batchTable,
                        selection=list(target="cell"),
                        options = list(
                          columnDefs = list(
                            list(targets = 0, render = JS("function(data, type, row, meta) {","return type === 'display' && data.length > 10 ?","'<span title=\"' + data + '\">' + data.substr(0, 7) + '...</span>' : data;","}"))
                          )
                        ),
                        callback = JS('table.page("next").draw(false);')
          )
        )
      }
      
      if (any(selectedBatchFiles[,2] == 3)) {
        batchFilesSampleType[selectedBatchFiles[selectedBatchFiles[,2] == 3,1]] <<- batchFilesSampleType[selectedBatchFiles[selectedBatchFiles[,2] == 3,1]]+1
        
        if (batchFilesSampleType[selectedBatchFiles[selectedBatchFiles[,2] == 3,1]] > length(sampleTypes)) 
          batchFilesSampleType[selectedBatchFiles[selectedBatchFiles[,2] == 3,1]] <<- 1
        
        batchTable <-
          cbind(
            dirname(batchFiles),
            basename(batchFiles),
            round(file.size(batchFiles) / 1000000, 1),
            sampleTypes[batchFilesSampleType],
            batchFilesRAM
          )
        colnames(batchTable) <- c("Dir","File","Size","SampleType","RAM")
        output$BatchProcessTable <- DT::renderDataTable(
          DT::datatable(batchTable,
                        selection=list(target="cell"),
                        options = list(
                          columnDefs = list(
                            list(targets = 0, render = JS("function(data, type, row, meta) {","return type === 'display' && data.length > 10 ?","'<span title=\"' + data + '\">' + data.substr(0, 7) + '...</span>' : data;","}"))
                          )
                        ),
                        callback = JS('table.page("next").draw(false);')
          )
        )
      }
    }
  })
  
  # Batch remove files ####
  observeEvent(input$BatchProcessRemoveFiles, {
    if (length(selectedBatchFiles) > 0) {
      batchFiles <<- batchFiles[-selectedBatchFiles[,1]]
      batchFilesRAM <<- batchFilesRAM[-selectedBatchFiles[,1]]
      batchFilesSampleType <<- batchFilesSampleType[-selectedBatchFiles[,1]]
      batchTable <-
        cbind(
          dirname(batchFiles),
          basename(batchFiles),
          round(file.size(batchFiles) / 1000000, 1),
          batchFilesSampleType,
          batchFilesRAM
        )
      colnames(batchTable) <- c("Dir","File","Size","SampleType","RAM")
      output$BatchProcessTable <- DT::renderDataTable(
        DT::datatable(batchTable,
                      selection=list(target="cell"),
                      options = list(
                        columnDefs = list(
                          list(targets = 0, render = JS("function(data, type, row, meta) {","return type === 'display' && data.length > 10 ?","'<span title=\"' + data + '\">' + data.substr(0, 7) + '...</span>' : data;","}"))
                        )
                      ),
                      callback = JS('table.page("next").draw(false);')
        )
      )
    }  
  })
}

shinyApp(ui = ui, server = server)
