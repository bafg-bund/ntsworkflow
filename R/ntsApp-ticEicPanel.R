

ticPanelUi <- function(id) {
  tagList(
    plotlyOutput(NS(id, "ticPlot")),
    hr(),
    plotlyOutput(NS(id, "ms1Plot"))
  )
}

ticPanelServer <- function(id, msFileList, msFileSelected) {
  stopifnot(is.reactive(msFileSelected))
  
  moduleServer(id, function(input, output, session) {

    ticTbl <- reactive({
      req(!is.null(msFileSelected()), length(msFileList) > 0)
      tic <- tibble(
        time = round(msFileList[[msFileSelected()]]@scantime/60, 2),
        intensity = msFileList[[msFileSelected()]]@tic
      )
    })
    
    output$ticPlot <- renderPlotly({
      ticPlot <- ggplot(ticTbl(), aes(time, intensity)) + 
        geom_line() + 
        labs(x = "Time (min)", y = "Intentsity (cps)")
      ggplotly(ticPlot, source = "ticPlot", dynamicTicks = TRUE) %>% 
        event_register("plotly_click")
    })
    
    spectrum <- reactive({
      clickData <- event_data("plotly_click", source = "ticPlot")
      req(!is.null(msFileSelected()), length(msFileList) > 0, !is.null(clickData))
      position <- which(abs(msFileList[[msFileSelected()]]@scantime-clickData$x*60)==min(abs(msFileList[[msFileSelected()]]@scantime-clickData$x*60)))
      as_tibble(xcms::getScan(msFileList[[msFileSelected()]], position))
    })
    
    output$ms1Plot <- renderPlotly({
      ms1Plot <- ggplot(spectrum()) + 
        geom_segment(aes(mz, xend = mz, y = 0, yend = intensity)) +
        labs(x = "m/z", y = "Intensity (cps)")
      ggplotly(ms1Plot, tooltip = c("x", "intensity"), dynamicTicks = TRUE)
    })
  })
}

eicPanelUi <- function(id) {
  tagList(
    fluidRow(
      column(width = 3,numericInput(NS(id, "eicMass"),"XIC m/z:",value = 100,width = "120px",step = 0.001)),
      column(width = 9,numericInput(NS(id, "mzTolerance"),"+/- (ppm)",value = 10,width = "100px",step = 5))
    ),
    plotlyOutput(NS(id,"eicPlot"))
  )
}

eicPanelServer <- function(id, msFileList, msFileSelected) {
  stopifnot(is.reactive(msFileSelected))
  moduleServer(id, function(input, output, session) {
    eicTbl <- reactive({
      req(!is.null(msFileSelected()), length(msFileList) > 0)
      delta <- input$eicMass * input$mzTolerance / 1000000 
      eicList <- xcms::rawEIC(
        msFileList[[msFileSelected()]],
        mzrange = c(input$eicMass - delta / 2, input$eicMass + delta / 2)
      )
      rt <- msFileList[[msFileSelected()]]@scantime/60
      eic <- as_tibble(eicList)
      eic$rt <- rt[pluck(eic, "scan")]
      eic
    })
    
    output$eicPlot <- renderPlotly({
      p <- ggplot(eicTbl(), aes(rt, intensity)) +
        geom_line() + labs(x = "Retention time (min)", y = "Intensity (cps)")
      ggplotly(p)
    })
  })
}


# Copyright 2025 Bundesanstalt für Gewässerkunde (Federal Institute of Hydrology)
# This file is part of ntsworkflow
