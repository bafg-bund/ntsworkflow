# Copyright 2016-2023 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow
# ntsworkflow is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any 
# later version.
# 
# ntsworkflow is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along 
# with ntsworkflow. If not, see <https://www.gnu.org/licenses/>.


# Non-Target App
# written by: Christian Dietrich, Kevin Jewell, Toni Köppe
# Last update: 2023-07-11

library(shiny)
library(shinyBS)
library(DT)
library(tcltk2)
library(shinyFiles)
library(foreach)
library(ggplot2)
library(ntsworkflow)

#### Set-up ####


if (exists("datenList") && length(datenList) > 0) {
  message("Using available data from current environment, click load from 
          environment to use this data. To start a new session, close the app
          and clear the environment.")
} else {
  message("Starting new session.")
  sampleList <<- data.frame(ID = numeric(), File = character(), 
                            sampleType = character(), optMzStep = numeric(), 
                            RAM = logical(), deleted = logical(), 
                            stringsAsFactors = F)
  datenList <<- list()
  peaklist <<- list()
  peakPickSettings <<- list()
  grouped <<- NULL
  annotationTable <<- NULL
  groupedWithAnnotation <<- NULL
  multiHitsTable <<- NULL
  selected<<- NULL
  selectedcell <<- NULL
  # For the MS2 overlay in the alignment tab
  MS2 <<- list() 
  # MS2 selection to be shown in the overlay
  MS2selected <<- 1 
  MS2_max_intensity <<- vector()
  peaklist_NeutralLoss <<- list()
  peaklist_FragmentIons <<- list()
  peaklist_FragmentIonsSubset <<- NULL
  peaklist_NeutralLossSubset <<- NULL
  FragmentAnalysis_PeaklistRows <<- NULL
  selectedFragmentIon <<- 0
  selectedNeutralLoss <<- 0
  
  batchFiles <<- NULL
  batchFilesRAM <<- NULL
  batchFilesSampleType <<- NULL
  
  selectedBatchFiles <<- NULL
  
  sampleTypes <<- c("Unknown","Blank","Standard")
  
  neutralLossTable <<- NULL
  FragmentIonTable <<- NULL
  
  pfad <<- getwd()
  
  counter <<- 0
}


#### UI ####
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width = 2,
      #### Add Sample, Save, Load ####
      shinyFilesButton("addSample", "Add Sample", "Please select a data file", TRUE),
      actionButton("BatchProcess","Batch Process"),
      bsModal("BatchProcessPopup","Batch Process","BatchProcess", size = "large", 
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
              
      ),
      hr(),
      DT::dataTableOutput("sampleTable"),
      shinyFilesButton("load", "Load", "Please select an RDS file", FALSE),
      shinySaveButton("save","Save", "Save file as...", filetype = list(RDS = ".RDS")),
      actionButton("loadEnv", "Load env.")
    ),
    
    #### Sample info ####
    mainPanel(width = 10,
      tabsetPanel(
        type = "tabs",
        tabPanel(
          "Sample info",
          fluidRow(column(width = 5, strong("Sample ID:")),
                   column(width = 7, textOutput("sampleID"))),
          fluidRow(column(width = 5, strong("MS1-Scans:")),
                   column(width = 7, textOutput("ms1scans"))),
          fluidRow(column(width = 5, strong("Scan times:")),
                   column(width = 7, textOutput("ms1RT"))),
          hr(),
          fluidRow(column(width = 5, strong("MS2-Scans:")),
                   column(width = 7, textOutput("ms2scans"))),
          fluidRow(column(width = 5, strong("Scan times:")),
                   column(width = 7, textOutput("ms2RT"))),
          hr(),
          fluidRow(column(
            width = 7, textInput("samplePath", label = NULL, width = "200%")
          ),
          column(
            width = 5, actionButton("changeSamplePath", "Change path")
          )),
          fluidRow(
            column(width = 3, actionButton("removeFromRAM", "Remove selected file from memory")),
            column(width = 3, actionButton("reloadToRAM", "Reload selected file to memory")),
            column(width = 3, actionButton("deleteSample", "Delete selected file"))
          ),
          hr(),
          fluidRow(
            column(width = 3, actionButton("removeAllFromRAM", "Remove all files from memory")),
            column(width = 3, actionButton("reloadAllToRAM", "Reload multiple files to memory", width = "100%"),
                   textInput("reloadChoice", NULL, width = "100%", 
                             placeholder = "Choose files e.g.: x:y (default: all files)"))
          ),
          hr(),
          fluidRow(
            column(
              6, selectInput("sampleType", "Sample type:", sampleTypes, sampleTypes[1])
            ),
            column(6, textInput("blankRegex", "Blank regex:", 
                               placeholder = "e.g., '_MQ_' or '_blank_'"))
          ),
          
          hr(),
          fluidRow(
            column(4, numericInput("numcores", 
                                   "Number of parallel threads to use:", 2, 1, 1, 1)),
            column(3, textOutput("memUsage"), offset = 2),
            column(3, actionButton("updateMem", "Update"))
          )
        ),
        tabPanel(
          "TIC",
          plotOutput(
            "TICplot",
            brush = brushOpts(
              id = "TICplot_brush",
              resetOnNew = TRUE,
              delay = 6000
            ),
            click = "TICplot_click",
            dblclick = "TICplot_dblclick"
          ),
          hr(),
          plotOutput(
            "MSPlot",
            brush = brushOpts(
              id = "MSplot_brush",
              resetOnNew = TRUE,
              delay = 6000
            ),
            click = "MSplot_click",
            dblclick = "MSplot_dblclick"
          )
        ),
        tabPanel(
          "XIC",
          fluidRow(
            column(
              width = 3,
              numericInput(
                "EIC_mass",
                "XIC m/z:",
                value = 100,
                width = "120px",
                step = 0.001
              )
            ),
            column(
              width = 9,
              numericInput(
                "EIC_mass_delta",
                "+/- (ppm)",
                value = 10,
                width = "100px",
                step = 5
              )
            )
          ),
          plotOutput(
            "XICplot",
            brush = brushOpts(
              id = "XICplot_brush",
              resetOnNew = TRUE,
              delay = 6000
            ),
            click = "XICplot_click",
            dblclick = "XICplot_dblclick"
          )
        ),
        #### Peak Picking tab ####
        tabPanel("Peak Picking",
                 tabsetPanel(
                   type = "tabs",
                   tabPanel(
                     "Parameters",
                     wellPanel(fluidRow(
                       column(
                         width = 9,
                         sliderInput(
                           "PeakPick_massrange",
                           "Mass range m/z:",
                           width = "100%",
                           min = 0,
                           max = 1200,
                           value = c(100, 200),
                           step = 10
                         )
                       ),
                       column(
                         width = 2,
                         numericInput(
                           "PeakPick_mzstep",
                           "m/z step:",
                           value = 0.02,
                           min = 0.005,
                           max = 1,
                           step = 0.005
                         ),
                         textOutput("mzsteprecommendation")
                       )
                     )),
                     column(
                       width = 9,
                       sliderInput(
                         "PeakPick_RTrange",
                         "RT range (min):",
                         width = "100%",
                         min = 0,
                         max = 30,
                         value = c(2, 20),
                         step = 0.1
                       ),
                       fluidRow(
                         column(
                           width = 2,
                           numericInput(
                             "PeakPick_IntTresh",
                             "Min. Intensity:",
                             value = 10,
                             min = 0,
                             step = 1
                           )
                         ),
                         column(
                           width = 2,
                           numericInput(
                             "PeakPick_SN",
                             "S/N:",
                             value = 3,
                             min = 1,
                             step = 1
                           )
                         ),
                         column(
                           width = 2,
                           numericInput(
                             "PeakPick_NoiseScans",
                             "Noise (scans):",
                             value = 30,
                             min = 10,
                             step = 10
                           )
                         ),
                         column(2, numericInput("maxPeaks","Max. peaks per peak",10,5,NA,1)),
                         column(
                           width = 4,
                           sliderInput(
                             "PeakPick_PeakWidthRange",
                             "Peak width (sec):",
                             min = 0.1,
                             max = 120.0,
                             value = c(5.0, 60.0),
                             ticks = FALSE,
                             round = FALSE,
                             step = .1
                           )
                         )
                       ),
                       wellPanel(fluidRow(
                         column(width = 5, h4("Componentization")),
                         checkboxInput(
                           "Componentization_Dynamic_Tolerance",
                           "Automatically adjust RT tolerances based on 13C peak of each component",
                           value = FALSE
                         ),
                         column(2, numericInput("componentization_ppm",  # PeakPick_ppm
                             "m/z tol. (ppm)", 5, 0, step = 10)),
                         column(2, numericInput("componentization_rt_tol",  # PeakPick_RT_Tol
                             "RT tol. apex (s)", 1.5, .1, step = .1)),
                         column(2, numericInput("componentization_rt_tol_l",
                                                "RT tol. left (s)",
                                                .5, .1, step = .1)),
                         column(3, numericInput("componentization_rt_tol_r",
                                                "RT tol. right (s)",
                                                1, .1, step = .1)),
                         column(3, numericInput("componentization_rt_tol_sum",
                                                "RT tol. sum (s)",
                                                2.5, .1, step = .1))
                       ),
                       fluidRow(
                         "Köppe, T., K. S. Jewell, C. Dietrich, A. Wick and T. A. 
                         Ternes (2020). 'Application of a non-target workflow for 
                         the identification of specific contaminants using the 
                         example of the Nidda river basin.' Water Research 178: 115703."
                       )
                      )
                     ),
                     column(
                       width = 3,
                       wellPanel(
                         h4("MS2 Assignment"),
                         numericInput(
                           "PPMS2Fit_ppm",
                           "MS2 precursor m/z Tolerance (ppm)",
                           value = 20,
                           min = 0,
                           step = 10
                         )
                       ),
                       wellPanel(
                         radioButtons(
                           "PPRadioSingleAll",
                           "",
                           choices = c("Selected File" = "one", "All Files" = "all"),
                           selected = "one"
                         ),
                         actionButton("PickPeaks", "Pick Peaks", width = "100%")
                       )
                     ),
                     hr()
                   ),
                   #### Data tab ####
                   tabPanel("Data",
                            fluidRow(
                              column(
                                width = 6,
                                fluidRow(
                                  actionButton("PeakPick_ShowAll", "All", width = "10%"),
                                  actionButton("PeakPick_GetChlorinated", "Cl", width = "10%"),
                                  actionButton("PeakPick_GetBrominated", "Br", width = "10%"),
                                  actionButton("PeakPick_HideFalsePositives", 
                                               "Hide False Positives", width = "30%"),
                                  actionButton("peakpick_random", "Random Sort", width = "20%")
                                ),
                                DT::dataTableOutput("PeakPickingTable"),
                                tags$script(
                                  ' $(document).on("keydown", function (e) {Shiny.onInputChange("lastkeypresscode", [e.keyCode,Math.random()]);switch(e.keyCode){case 38: case 40: case 33: case 34: e.preventDefault()};});'
                                ),
                                shinySaveButton("PickPeaksTableExport","Export .csv", 
                                                "Save file as...", filetype = list(csv = "csv"))
                              ),
                              column(
                                width = 6,
                                plotOutput(
                                  "PeakPickXIC",
                                  brush = brushOpts(
                                    id = "PeakPickXIC_brush",
                                    resetOnNew = TRUE,
                                    delay = 6000
                                  ),
                                  dblclick = "PeakPickXIC_dblclick"
                                ),
                                plotOutput("PeakPickMS"),
                                plotOutput("PeakPickMS2",
                                           click = "PeakPickMS2_click"),
                                bsModal(
                                  "PeakPickMS2popup",
                                  "MS2 data",
                                  "PeakPickMS2_click",
                                  size = "small",
                                  DT::dataTableOutput("PeakPickMS2Table")
                                )
                              )
                            )),
                   #### Components tab ####
                   tabPanel("Components",  
                            fluidRow(
                              column(width=6, 
                                     DT::dataTableOutput("ComponentsTable"),
                                     tags$script(' $(document).on("keydown", function (e) {Shiny.onInputChange("Componentskeypresscode", [e.keyCode,Math.random()]);switch(e.keyCode){case 38: case 40: case 33: case 34: e.preventDefault()};});'),
                                     actionButton("ComponentsTableExport","Export to csv")),
                              column(width=6, plotOutput("ComponentsXIC",
                                                         brush = brushOpts(id="ComponentsXIC_brush",resetOnNew = TRUE, delay = 6000),
                                                         dblclick = "ComponentsXIC_dblclick"),
                                     DT::dataTableOutput("ComponentsTable2")
                              )
                            )
                   ),
                   #### GenForm tab ####
                   tabPanel(
                     "GenForm", 
                     fluidRow(  # settings
                       column(2, textInput("genformIon", "Ion type", "", placeholder = "-e, +e, -H, +H or +Na")),
                       column(2, numericInput("genformMS1tol", "MS1 tol ppm", 5, 1, 1000, 1)),
                       column(2, numericInput("genformMS2tol", "MS2 tol ppm", 15, 1, 1000, 1)),
                       column(2, textInput("genformElements", "Elements", "CHBrClFINOPSSi", placeholder = "CHBrClFINOPSSi")),
                       column(2, textInput("genformFf", "Fuzzy Formula", "", placeholder = "e.g.: C0-6H0-20O1-3")),
                       column(2, style = "margin-top: 25px;", actionButton("genformGo", "Compute"))
                     ),
                     fluidRow(  # results
                       column(12, DT::dataTableOutput("genFormResTab"))
                     ),
                     fluidRow(  # citation
                       "Meringer, M., Reinker, S., Zhang, J., & Muller, A. (2011). MS/MS data improves 
                      automated determination of molecular formulas by mass spectrometry. 
                      MATCH Commun. Math. Comput. Chem, 65(2), 259-290. & 'GenformR' (Schymanski, E.)"
                     ))
                 )),
        #### Alignment tab ####
        tabPanel(
          "Alignment",
            tabsetPanel(
              type = "tabs",  
                tabPanel(
                  "Table",            
                    fluidRow(
                      column(1, numericInput("Alignment_deltamz", "m/z tolerance:",value = 10,
                               min = 1, max = 100, step = 1)),
                      column(1, radioButtons("aligMzTolType", NULL, c("ppm", "mDa"), "ppm", width = "100%")),
                      column(2, numericInput("Alignment_deltaRT", "RT tolerance (s):", value = 30, 
                                              min = 1, max = 120, step = 1)),
                      
                      column(2, numericInput("blankCorrFac", "Blank factor:", value = 10, step = 0.1)),
                      column(2,numericInput("trendFac", "Similar trends r:", NA, .1, 1, .1)),
                      
                      column(2, actionButton("AlignPeaks", "Align Features", width = "100%"),
                             actionButton("Normalize", "Normalize", width = "100%"),
                             actionButton("BlankCorrection", "Blank Correction", width = "100%")),
                      column(2, checkboxInput("alignLeaders", "Align only group leaders", value = FALSE),
                             checkboxInput("parallelAlign", "Parallel alignment", value = FALSE),
                             checkboxInput("deleteGrouped", "Also remove grouped in blank correction", value = FALSE))
                      ),
                    fluidRow(column(6, DT::dataTableOutput("AlignmentTable"),
                        shinySaveButton("groupedExport", "Export .csv", "Save file as...", 
                                        filetype = list(csv = "csv")),
                        
                        h4("Further alignment table filters"),
                        # remove rows with less than x number of peaks
                        hr(),
                        fluidRow(
                          column(6, p("Remove rows with fewer than", align = "right")),
                          column(2, numericInput("minDetections", NULL, 1, 1, NA, 1, "100%")),
                          column(2, p("detections")),
                          column(2, actionButton("removeRare", "Go", width = "100%"))
                        ),
                        hr(),
                        # remove rows which are never groupleader
                        fluidRow(
                          column(
                            8, 
                            p("Remove features which are not the highest feature in their component in files", 
                              align = "right")
                          ),
                          column(2, textInput("leaderSamples", NULL, width = "100%", placeholder = "x:y")),
                          column(2, actionButton("onlyGroupLeaders", "Go", width = "100%"))),
                        # keep only peaks found in replicate injections
                        hr(),
                        fluidRow(
                          p("Remove features which are not found in replicate injections (intensity
                            will be set to zero)", align = "center")
                        ),
                        fluidRow(
                          column(4, textInput("repSamples", "Files to consider:", width = "100%", 
                                              placeholder = "x:y")),
                          column(3, numericInput("repNum", "No. of replicates:", 1,1,NA,1,"100%")),
                          column(3, numericInput("repLeast", "In at least:",1,1,NA,1,"100%")),
                          column(2, actionButton("keepReps", "Go", width = "100%"))
                        ),
                        hr(),
                        # average intensities of replicate injections and remove replicates
                        fluidRow(p("Get average intensity of features in replicates, 
                                   ignoring zero intensity features (replicate files will be removed)", 
                                   align = "center")),
                        fluidRow(
                          column(4, textInput("repAveSamples", "Files to consider:", width = "100%", 
                                              placeholder = "x:y")),
                          column(3, numericInput("repAveNum", "No. of replicates:", 1,1,NA,1,"100%")),
                          column(3, ""),
                          column(2, actionButton("aveReps", "Go", width = "100%"))
                        ),
                        hr(),
                        fluidRow(
                          p("Remove rows where features are found in fewer than X consecutive files.")
                        ),
                        fluidRow(
                          column(4, textInput("consecSamples", "Files to consider:", width = "100%",
                                              placeholder = "x:y")),
                          column(3, numericInput("consecNum", "No. of consecutive:", 1,1,NA,1,"100%")),
                          column(3, ""),
                          column(2, actionButton("consecGo", "Go", width = "100%"))
                        ),
                        hr()
                        ),
                        column(
                          6, 
                          plotOutput("AlignmentXICs", dblclick = "aligXIC_dblclick", 
                                     brush = brushOpts(id = "aligXIC_brush", resetOnNew = T)),
                          plotOutput("AlignmentTrend"),
                          plotOutput("AlignmentMS2", click = "AlignmentMS2_click"))
                    )
                ),
                tabPanel(
                  "Componentisation 2",
                  fluidRow(p("Second stage componentisation based on four criteria: 
                              rt, peak shape, intensity correlation and common mz differences.
                             Session must include five or more raw data files.
                             Will replace group column in the alignment table.")),
                  fluidRow(
                    column(2, numericInput("componen2mztol", "m/z tolerance (mDa)", 5, 1, 1000, 1)),
                    column(2, numericInput("componen2rttol", "RT tolerance (s)", 3, 1, 60, 1)),
                    column(2, numericInput("componen2fracShape", "Min. frac. peak shape matches", 0.5, 0.1, 1, 0.05)),
                    column(2, numericInput("componen2corr", "Min. Pearson's r for inten. trend", 0.5, 0.1, 1, 0.05)),
                    column(2, radioButtons("componen2pol", "Polarity", choices = c("pos", "neg"), selected = "pos", inline = T)),
                    column(2, actionButton("componen2go", "Compute", width = "100%"))
                  ),
                  fluidRow(
                    column(6, fluidRow(textOutput("componen2text")), fluidRow(DT::dataTableOutput("componen2table"))),
                    column(6, plotOutput("componen2hist"))
                  )
                ),
                tabPanel(
                  "Overview", 
                  fluidRow(
                    plotOutput("overview", click = "overview_click", dblclick = "overview_dblclick", 
                               hover = hoverOpts(id = "overviewHover"), 
                               brush = brushOpts(id = "overview_brush", resetOnNew = TRUE)),
                    plotOutput("intTrendOverview"),
                    DT::dataTableOutput("overviewSelected")
                  )
                ),
                tabPanel(
                  "Highest intensities",
                  fluidRow(
                    plotOutput("intRangePl", click = "intRangePl_click",dblclick = "intRangePl_dblclick", 
                               hover = hoverOpts(id = "highIntHover"),
                               brush = brushOpts(id = "intRangePl_brush", resetOnNew = TRUE)),
                    plotOutput("intTrendHighInt"),
                    DT::dataTableOutput("intRangeSelected")
                  )
                ),
                tabPanel(
                  "Cluster analysis",
                  fluidRow(plotOutput("intDendro", hover = hoverOpts(id ="dendHover"),
                                      brush = brushOpts(id = "dendroBrush", resetOnNew = TRUE),
                                      dblclick = "dendroDblClick", click = "dendroClick")),
                  
                  fluidRow(plotOutput("intTrendDendro", height = "300px")),
                  fluidRow(DT::dataTableOutput("dendroSelected"))         
                ),
                tabPanel(
                  "Similar trends",   
                  fluidRow(
                    column(6, DT::dataTableOutput("SimilarTrendsTable")),
                    column(6,
                           h4(textOutput("simTrenHead"), align = "center"),
                           plotOutput("SimilarTrends", hover = hoverOpts(id = "simTrenHover")),
                           h5(textOutput("whichSimTren"), align = "center")
                           ) 
                  )     
                )
              )
        ),
        #### Annotation tab ####
        tabPanel(
          "Annotation",
          br(),
          fluidRow(
            column(1, numericInput("annotMzTolmDa", "MS1 tol (mDa)", value = 5, min = 1, step = 1)),
            column(1, numericInput("annotRtTolM", "RT tol (min)", value = 1, min = 0.1, step = 0.1)),
            column(1, numericInput("annotDpWind", "MS2 tol (mDa)", value = 15, min = 1, step = 1)),
            column(1, numericInput("annotThresh", "DP threshold", value = 500, min = 10, 
                                   max = 1000, step = 10)),
            column(1, numericInput("annotRtOffset", "Ret time offset", value = 0, step = 0.1)),
            column(1, numericInput("annotIntCut", "Int. cut rel.", value = 0, min = 0, max = 0.99, step = 0.01)),
            column(1, radioButtons("annotPol", "Pol.", choices = c("pos", "neg"), selected = "pos")),
            column(2, sliderInput("annotCE", "CE", min = 10, max = 110, value = c(30,40), 
                                  ticks = FALSE, step = 10)),
            column(1, sliderInput("annotCES", "CES", min = 0, max = 15, value = c(0,15), 
                                  ticks = FALSE, step = 1)),
            column(1, shinyFilesButton("ms2Db", label="Choose DB", 
                                       title = "Please select a database", multiple = FALSE),
                   textOutput("dbLoc")),
            column(1, actionButton("dbHelp", "?"))
          ),
          fluidRow(
            column(3, hr()),
            column(2, checkboxInput("annotAppend", "Append annotation")),
            column(1, p(strong("Ret. times:"))),
            column(1, selectInput("annotChromMeth", NULL, "unknown", multiple = FALSE)),
            column(1, p(strong("Spec. Source:"))),
            column(1, selectInput("annotExpSource", NULL, "unknown", multiple = TRUE)),
            column(2, actionButton("annotGo", "Annotate", width = "100%")),
            column(1, shinySaveButton("annotationTableExport", "Export .csv", "Save file as...",
                                      filetype = list(csv = "csv")))
          ),
          fluidRow(
            DT::dataTableOutput("AlignmentAnnotTable")
          ),
          fluidRow(
            DT::dataTableOutput("multiHitsTable"),
            column(12, hr())
          ),
          fluidRow(
            column(6, plotOutput("annotationXIC")),
            column(6, plotOutput("annotationMS1"))
          ),
          fluidRow(
            column(6, plotOutput("annotationMS2")),
            column(6, plotOutput("annotationTrend"))
          ),
          fluidRow(
            "Jewell, K. S., U. Kunkel, B. Ehlig, F. Thron, M. Schlüsener, C. 
            Dietrich, A. Wick and T. A. Ternes (2019). 'Comparing mass, 
            retention time and MS2 spectra as criteria for the automated 
            screening of small molecules in aqueous environmental samples 
            analyzed by LC-QToF-MS/MS.' Rapid Communications in Mass Spectrometry 34: e8541."
          )
        ) 
       
      )
    )
  )
)

server <- function(input, output, session) {
  
  #### Set-up ####
  options(shiny.maxRequestSize = 3000 * 1024^2)
  
  `%notin%` <- function(x, y) !(x %in% y)
  
  # Genform is a command line program
  GENFORM_LOC <- "/usr/local/bin/"
  
  # reactive values
  TICplot_ranges <- reactiveValues(x = NULL, y = NULL)
  TICplot_RT <- reactiveValues(x = 1)
  XICplot_ranges <- reactiveValues(x = NULL, y = NULL)
  PPXIC_ranges <- reactiveValues(x = NULL, y = NULL)
  ComponentXIC_ranges <- reactiveValues(x = NULL, y = NULL)
  FAXIC_ranges <- reactiveValues(x = NULL, y = NULL)
  masse <- reactiveValues(mz = 100, y = 100)
  MSplot_ranges <- reactiveValues(x = NULL, y = NULL)
  
  # update ui
  
  # Number of threads needs to be kept low due to very high memory demands
  observeEvent(NA, updateNumericInput(
    session, "numcores",
    value = min(round(parallel::detectCores()/3), 6),
    max = parallel::detectCores() - 1
    ), 
    once = TRUE)
  
  # Plotting functions ####
  
  MSplot_plotten <- function() {
    output$MSPlot <- renderPlot({
      plot(
        spektrum,
        type = "h",
        xlim = MSplot_ranges$x,
        ylim = MSplot_ranges$y,
        main = "MS Spectrum",
        xlab = "m/z",
        ylab = "Intensity",
        xaxs = "i",
        yaxs = "i"
      )
      text(x=masse$mz,y=masse$y,label=round(masse$mz,4),pos=4,offset=0.2)
      
    })
  }
  
  XICplot_plotten <- function() {
    output$XICplot <- renderPlot({
      delta <- input$EIC_mass*input$EIC_mass_delta/1000000 
      if (length(datenList[[selected]]@env$intensity) == 1) {
        plot(x=0,y=0)
      } else {
        XIC <- xcms::rawEIC(
          datenList[[selected]],
          mzrange = c(input$EIC_mass - delta / 2, input$EIC_mass + delta / 2)
        )
        plot(
          x = datenList[[selected]]@scantime / 60,
          y = XIC$intensity,
          type = "l",
          xlim = XICplot_ranges$x,
          ylim = XICplot_ranges$y,
          main = "XIC",
          xlab = "RT",
          ylab = "Intensity",
          xaxs = "i",
          yaxs = "i"
        )
      }  
      counter <<- counter + 1
    })
  }
  
  PeakPickXIC_plotten <- function(zeile, dataSel){
     
    req(nrow(peaklist[[dataSel]]) >= 1)
    mz_step <- input$PeakPick_mzstep 
    noiseRTleft <- 0
    noiseRTright <- length(datenList[[dataSel]]@scantime)
    
    if ((peaklist[[dataSel]]$Leftendscan[zeile]-peakPickSettings[[dataSel]]$NoiseScans) > 0) 
      noiseRTleft <- datenList[[dataSel]]@scantime[peaklist[[dataSel]]$Leftendscan[zeile]-peakPickSettings[[dataSel]]$NoiseScans]/60
    
    if ((peaklist[[dataSel]]$Rightendscan[zeile]+peakPickSettings[[dataSel]]$NoiseScans) < length(datenList[[dataSel]]@scantime)) 
      noiseRTright <- datenList[[dataSel]]@scantime[peaklist[[dataSel]]$Rightendscan[zeile]+peakPickSettings[[dataSel]]$NoiseScans]/60
    
    infotext <- ""
    
    if (!peaklist[[dataSel]]$RealPeak[zeile]) 
      infotext <- "FALSE POSITIVE"
    
    XIC <- xcms::rawEIC(
      datenList[[dataSel]],
      mzrange = c(
        peaklist[[dataSel]]$mz[zeile] - mz_step / 2, 
        peaklist[[dataSel]]$mz[zeile] + mz_step / 2
      )
    )
    
    plot(
      x = datenList[[dataSel]]@scantime / 60,
      y = XIC$intensity,
      type = "l",
      xlim = PPXIC_ranges$x,
      ylim = PPXIC_ranges$y,
      main = paste0(
        "XIC (",
        round(peaklist[[dataSel]]$mz[zeile] - mz_step / 2, 4),
        "-",
        round(peaklist[[dataSel]]$mz[zeile] + mz_step / 2, 4),
        ")"
      ),
      xlab = "RT",
      ylab = "Intensity",
      xaxs = "i",
      yaxs = "i"
    )
    lines(
      x = c(peaklist[[dataSel]]$LeftendRT[zeile] / 60, peaklist[[dataSel]]$LeftendRT[zeile] / 60),
      y = c(
        peaklist[[dataSel]]$Baseline[zeile],
        (peaklist[[dataSel]]$Intensity[zeile]) + peaklist[[dataSel]]$Baseline[zeile]
      ),
      type = "l",
      col = "red"
    )
    lines(
      x = c(peaklist[[dataSel]]$RightendRT[zeile] / 60, peaklist[[dataSel]]$RightendRT[zeile] / 60),
      y = c(
        peaklist[[dataSel]]$Baseline[zeile],
        (peaklist[[dataSel]]$Intensity[zeile]) + peaklist[[dataSel]]$Baseline[zeile]
      ),
      type = "l",
      col = "red"
    )
    lines(
      x = c(noiseRTleft, noiseRTright),
      y = c(peaklist[[dataSel]]$Baseline[zeile], peaklist[[dataSel]]$Baseline[zeile]),
      type = "l",
      col = "red"
    )
    lines(
      x = c(noiseRTleft, noiseRTright),
      y = c(
        peaklist[[dataSel]]$Baseline[zeile] + peaklist[[dataSel]]$NoiseDeviation[zeile] / 2,
        peaklist[[dataSel]]$Baseline[zeile] + peaklist[[dataSel]]$NoiseDeviation[zeile] / 2
      ),
      type = "l",
      col = "blue"
    )
    lines(
      x = c(noiseRTleft, noiseRTright),
      y = c(
        peaklist[[dataSel]]$Baseline[zeile] - peaklist[[dataSel]]$NoiseDeviation[zeile] / 2,
        peaklist[[dataSel]]$Baseline[zeile] - peaklist[[dataSel]]$NoiseDeviation[zeile] / 2
      ),
      type = "l",
      col = "blue"
    )
    lines(
      x = c(peaklist[[dataSel]]$FWHM_left[zeile] / 60, peaklist[[dataSel]]$FWHM_left[zeile] / 60),
      y = c(
        peaklist[[dataSel]]$Baseline[zeile],
        (peaklist[[dataSel]]$Intensity[zeile] / 2) + peaklist[[dataSel]]$Baseline[zeile]
      ),
      type = "l",
      col = "green"
    )
    lines(
      x = c(peaklist[[dataSel]]$FWHM_right[zeile] / 60, peaklist[[dataSel]]$FWHM_right[zeile] / 60),
      y = c(
        peaklist[[dataSel]]$Baseline[zeile],
        (peaklist[[dataSel]]$Intensity[zeile] / 2) + peaklist[[dataSel]]$Baseline[zeile]
      ),
      type = "l",
      col = "green"
    )
    mtext(
      infotext,
      side = 3,
      line = 0,
      col = "red"
    )
  }
  
  ComponentXIC_plotten <- function(zeile) {
    req(selectedR())
    selected <- selectedR()
    mz_step <- input$PeakPick_mzstep 
    noiseRTleft <- 0
    noiseRTright <- length(datenList[[selected]]@scantime)
    
    if ((peaklist[[selected]]$Leftendscan[zeile]-peakPickSettings[[selected]]$NoiseScans) > 0) 
      noiseRTleft <- datenList[[selected]]@scantime[peaklist[[selected]]$Leftendscan[zeile]-peakPickSettings[[selected]]$NoiseScans]/60
    
    if ((peaklist[[selected]]$Rightendscan[zeile]+peakPickSettings[[selected]]$NoiseScans) < length(datenList[[selected]]@scantime)) 
      noiseRTright <- datenList[[selected]]@scantime[peaklist[[selected]]$Rightendscan[zeile]+peakPickSettings[[selected]]$NoiseScans]/60
    
    XIC <- list()
    Komponenten <- which(peaklist[[selected]]$Gruppe == peaklist[[selected]]$Gruppe[zeile])
    Komponenten <- Komponenten[order(peaklist[[selected]]$Intensity[Komponenten])]
    farben <- rainbow(length(Komponenten))
    
    for (i in 1:length(Komponenten)) {
      XIC[[i]] <- xcms::rawEIC(
        datenList[[selected]],
        mzrange = c(
          peaklist[[selected]]$mz[Komponenten[i]] - mz_step / 2,
          peaklist[[selected]]$mz[Komponenten[i]] + mz_step / 2
        )
      )
    }
    plot(
      x    = datenList[[selected]]@scantime / 60,
      y    = XIC[[1]]$intensity,
      type = "l",
      xlim = ComponentXIC_ranges$x,
      ylim = ComponentXIC_ranges$y,
      xlab = "RT",
      ylab = "Intensity",
      xaxs = "i",
      yaxs = "i",
      col  = farben[1],
      axes = FALSE
    )
      if (length(XIC) > 1) {
        for (i in 2:length(XIC)) {
          par(new=TRUE)
          plot(
            x    = datenList[[selected]]@scantime / 60,
            y    = XIC[[i]]$intensity,
            type = "l",
            xlim = ComponentXIC_ranges$x,
            ylim = ComponentXIC_ranges$y,
            xlab = "",
            ylab = "",
            xaxs = "i",
            yaxs = "i",
            col  = farben[i],
            axes = FALSE
          )
        }
      }
      axis(side=1)
      axis(side=2)
  }
  
  PeakPickMSPlotten <- function(zeile) {
    if (nrow(peaklist[[selected]]) < 1)
      return(NULL)
    spektrum_PP <- xcms::getScan(datenList[[selected]],peaklist[[selected]]$Scan[zeile])
    output$PeakPickMS <- renderPlot({
       
      spektrum_PP <- as.data.frame(spektrum_PP)
      masse <- peaklist[[selected]]$mz[zeile]
      RT <- peaklist[[selected]]$RT[zeile]
      inten <- peaklist[[selected]]$Intensity[zeile]
      spektrum_PP <- spektrum_PP[spektrum_PP$mz > masse - 15 & spektrum_PP$mz < masse + 15, ]
      plotme <- ggplot(spektrum_PP, aes(mz, intensity, label = round(mz,4))) +
        geom_segment(aes(x = mz, xend = mz, y = 0, yend = intensity),
                     stat = "identity", linewidth = .5, alpha = .5) +
        theme_bw(base_size = 14) +
        ggtitle(paste0("MS1 of ", round(masse, 4), " @ ", round(RT/60,1), " min")) +
        geom_text(data = spektrum_PP[spektrum_PP$intensity > inten / 10, ], 
                  check_overlap = TRUE, vjust = -0.5) +
        geom_vline(xintercept = masse, color = "red", alpha = 0.2) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, inten + inten/10), 
                        xlim = c(-10, 10) + masse) +
        ylab("Intensity") +
        xlab("m/z (u)")
      
      infotext <- ""
      
      if (peaklist[[selected]]$C13[zeile] > 0) 
        infotext <- paste0(infotext, "! 13C of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$C13[zeile]],4)," !")
      
      if (peaklist[[selected]]$NaAddukt[zeile] > 0) 
        infotext <- paste0(infotext,"! Na adduct of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$NaAddukt[zeile]],4)," !")
      
      if (peaklist[[selected]]$NH4Addukt[zeile] > 0) 
        infotext <- paste0(infotext,"! NH4 adduct of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$NH4Addukt[zeile]],4)," !")
      
      if ((peaklist[[selected]]$Cl1[zeile] > 0) & (peaklist[[selected]]$Cl1[zeile] != zeile)) 
        infotext <- paste0(infotext,"! Cl Isotope (1 Cl) of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$Cl1[zeile]],4)," !")
      
      if (peaklist[[selected]]$Cl1[zeile] == zeile) 
        infotext <- paste0(infotext,"! Shows Cl Isotope Pattern (1 Cl) !")
      
      if ((peaklist[[selected]]$Cl2[zeile] > 0) & (peaklist[[selected]]$Cl2[zeile] != zeile)) 
        infotext <- paste0(infotext,"! Cl Isotope (2 Cl) of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$Cl2[zeile]],4)," !")
      
      if (peaklist[[selected]]$Cl2[zeile] == zeile) 
        infotext <- paste0(infotext,"! Shows Cl Isotope Pattern (2 Cl) !")
      
      if ((peaklist[[selected]]$Cl3[zeile] > 0) & (peaklist[[selected]]$Cl3[zeile] != zeile)) 
        infotext <- paste0(infotext,"! Cl Isotope (3 Cl) of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$Cl3[zeile]],4)," !")
      
      if (peaklist[[selected]]$Cl3[zeile] == zeile) 
        infotext <- paste0(infotext,"! Shows Cl Isotope Pattern (3 Cl) !")
      
      if ((peaklist[[selected]]$Cl4[zeile] > 0) & (peaklist[[selected]]$Cl4[zeile] != zeile))  
        infotext <- paste0(infotext,"! Cl Isotope (4 Cl) of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$Cl4[zeile]],4)," !")
      
      if (peaklist[[selected]]$Cl4[zeile] == zeile) 
        infotext <- paste0(infotext,"! Shows Cl Isotope Pattern (4 Cl) !")
      
      if (peaklist[[selected]]$Br1[zeile] > 0) 
        infotext <- paste0(infotext,"! Br Isotope (1 Br) of ",
                           round(peaklist[[selected]]$mz[peaklist[[selected]]$Br1[zeile]],4)," !")
      
      if (peaklist[[selected]]$Br2[zeile] > 0) 
        infotext <- paste0(infotext,"! Br Isotope (2 Br) of ",
                           round(peaklist[[selected]]$mz[peaklist[[selected]]$Br2[zeile]],4)," !")
      
      if (peaklist[[selected]]$Br3[zeile] > 0) 
        infotext <- paste0(infotext,"! Br Isotope (3 Br) of ",
                           round(peaklist[[selected]]$mz[peaklist[[selected]]$Br3[zeile]],4)," !")
      
      if (peaklist[[selected]]$Br4[zeile] > 0) 
        infotext <- paste0(infotext,"! Br Isotope (4 Br) of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$Br4[zeile]],4)," !")
      
      if (peaklist[[selected]]$KAddukt[zeile] > 0) 
        infotext <- paste0(infotext,"! K adduct of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$KAddukt[zeile]],4)," !")
      
      if ((peaklist[[selected]]$S1[zeile] > 0) & (peaklist[[selected]]$S1[zeile] != zeile)) 
        infotext <- paste0(infotext,"! S Isotope (1 S) of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$S1[zeile]],4)," !")
      
      if ((peaklist[[selected]]$S2[zeile] > 0) & (peaklist[[selected]]$S2[zeile] != zeile)) 
        infotext <- paste0(infotext,"! S Isotope (2 S) of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$S2[zeile]],4)," !")
      
      if (peaklist[[selected]]$S1[zeile] == zeile) 
        infotext <- paste0(infotext,"! Shows S Isotope Pattern (1 S) !")
      
      if (peaklist[[selected]]$S2[zeile] == zeile) 
        infotext <- paste0(infotext,"! Shows S Isotope Pattern (2 S) !")
      
      if (peaklist[[selected]]$InSourceFragmentOf[zeile] > 0) 
        infotext <- paste0(infotext,"! Source Fragment of ",round(peaklist[[selected]]$mz[peaklist[[selected]]$InSourceFragmentOf[zeile]],4)," !")
      
      plotme <- plotme + annotate("text", x = masse - 10, y = Inf, label = infotext, vjust = 1.1, hjust = -.1, 
                                  colour = "red", alpha = .7)
      plotme
    })
    
  }
  
  ms2spektrumR <- reactiveVal() 
  observe({
    input$PeakPickingTable_cell_clicked
    input$sampleTable_cell_clicked
    zeile <- input$PeakPickingTable_rows_selected 
    req(zeile)
    req(!is.null(selectedR()))
    req(!is.null(peaklist[[selectedR()]]))
    req(nrow(peaklist[[selectedR()]]) >= 1)
    
    if (zeile > nrow(peaklist[[selectedR()]]))
      zeile <- nrow(peaklist[[selectedR()]])
    
    if (peaklist[[selectedR()]]$MS2scan[zeile] > 0) {
      mz <- peaklist[[selectedR()]]$mz[zeile]
      specTable <- xcms::getMsnScan(datenList[[selectedR()]],peaklist[[selectedR()]]$MS2scan[zeile])
      specTable <- specTable[specTable[, "mz"] < mz + 1, ]
      ms2spektrumR(specTable)
    } else {
      ms2spektrumR(NULL)
    }
  })  
 
  PeakPickMS2Plotten <- function(zeile) {
    
      req(zeile)
      req(ms2spektrumR())
      req(selectedR())
      selected <- selectedR()
      masse <- datenList[[selected]]@msnPrecursorMz[peaklist[[selected]]$MS2scan[zeile]]
      RT <- datenList[[selected]]@msnRt[peaklist[[selected]]$MS2scan[zeile]]
      
      ms2spektrum_ <- as.data.frame(ms2spektrumR())
      comp_mz <- peaklist[[selected]][zeile, "mz"]
      ms2spektrum_ <- ms2spektrum_[ms2spektrum_$mz < comp_mz + 10, ]
      plotme <- ggplot(ms2spektrum_, aes(mz, intensity, label = round(mz,4))) +
        geom_segment(aes(x = mz, xend = mz, y = 0, yend = intensity),
                     stat = "identity", linewidth = .5, alpha = .5) +
        theme_bw(base_size = 14) +
        ggtitle(paste0("MS2 of ",round(masse,4)," @ ",round(RT/60,1)," min")) +
        geom_text(data = ms2spektrum_[ms2spektrum_$intensity > 0.01, ], check_overlap = TRUE, vjust = -0.5) +
        geom_vline(xintercept = comp_mz, color = "red", alpha = 0.2) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, max(ms2spektrum_$intensity)*1.1), xlim = c(0, comp_mz + 5)) +
        ylab("Intensity") +
        xlab("m/z (u)")
      if ((RT < peaklist[[selected]][zeile, "LeftendRT"]) || (RT > peaklist[[selected]][zeile, "RightendRT"])) {
        plotme +  
          annotate("text", x = -Inf, y = 1.05, label = "MS2 scan outside chromatographic peak", 
                   vjust = 0, hjust = -.1, fontface = "bold.italic", alpha = .5)
      }
      plotme

  }
  
  alignmentTrendPlotten <- function(zeile, titel = "", tabelle = grouped) {
    # get sample numbers
    intCols <- colnames(tabelle)[grep("^Int_", colnames(tabelle))]
    sampNumbers <- as.numeric(stringr::str_match(intCols, "_(\\d+)$")[,2])
    
    inten <- data.frame(samp = as.factor(sampNumbers), 
                        int = as.numeric(tabelle[zeile, intCols])) 
    
    ggplot(inten, aes(samp, int)) + geom_point(shape = 1) + 
      scale_y_continuous(limits = c(0, NA)) + xlab("Sample") + ylab("Inten. (cps)") +
      ggtitle(titel) + theme_bw(14)
  }
  
  aligXICranges <- reactiveValues(x = NULL, y = NULL)
  
  alignmentXICsPlotten <- function(ID) {
    XIC <- list()
    sl <- sampleListR()
    validate(need(any(sl$RAM), "Load sample to RAM to display EICs"))
    # get samples to plot from grouped table
    gn <- colnames(grouped)
    intColNames <- gn[grep("^Int_", gn)]
    files <- as.numeric(stringr::str_match(intColNames, "^Int_(.*)$")[,2])
    files <- files[files %in% sl[sl$RAM & !sl$deleted, "ID"]]
    # which row has this ID?
    zeile <- which(grouped[, "alignmentID"] == ID)
    for (i in files) {
      mzstep <- peakPickSettings[[i]]$mz_step
      if (length(datenList[[i]]@env$intensity) > 1) {
        XIC[[i]] <- xcms::rawEIC(datenList[[i]],
                           mzrange = c(grouped[zeile, paste0("mz_", i)] - mzstep / 2,
                                       grouped[zeile, paste0("mz_", i)] + mzstep / 2),
                           rtrange = peakPickSettings[[i]]$rtrange)
        XIC[[i]]$scanTime <- round(datenList[[i]]@scantime[XIC[[i]]$scan]/60, 2)
      } else {
        XIC[[i]] <- list()
        XIC[[i]]$intensity <- rep(0, length(datenList[[1]]@scantime))
        XIC[[i]]$scanTime <- round(datenList[[i]]@scantime/60, 2)
      }  
    }
     
    XIC <- compact(XIC)
    # make one table from XIC
    XIC <- lapply(XIC, as.data.frame)
    XIC <- Map(function(df, num) transform(df, samp = num), XIC, files)
    XIC <- do.call("rbind", XIC)
    XIC$samp <- as.factor(XIC$samp)
    XIC <- XIC[XIC$intensity != 0,]  # remove row with intensity 0 to reduce length and increase speed
    # overlay plot
    ggplot(XIC, aes(scanTime, intensity, color = samp)) + geom_line() + guides(color = "none") +
      xlab("Time (min.)") + ylab("Inten. (cps)") + 
      geom_vline(xintercept = grouped[zeile, "mean_RT"]/60, color = "red", alpha = .2) + 
      theme_bw(14) + coord_cartesian(xlim = aligXICranges$x, ylim = aligXICranges$y, expand = F)
  }
  
  observeEvent(input$aligXIC_dblclick, {
    brush <- input$aligXIC_brush
    if (!is.null(brush)) {
      aligXICranges$x <- c(brush$xmin, brush$xmax)
      aligXICranges$y <- c(brush$ymin, brush$ymax)
    } else {
      aligXICranges$x <- NULL
      aligXICranges$y <- NULL
    }
    output$AlignmentXICs <- renderPlot(alignmentXICsPlotten(alignmentIdSelected()))
  })
  
  aligMs2R <- reactiveValues(options = 0, current = 1)
  
  alignmentMS2Plotten <- function(ID, selection = 0) {
    # which row has this ID?
    zeile <- which(grouped[, "alignmentID"] == ID)
    sl <- sampleListR()
    validate(need(any(sl$RAM), "Load sample to RAM to display MS2s"))
    MS2 <- list()
    gn <- colnames(grouped)
    intColNames <- gn[grep("^Int_", gn)]
    files <- as.numeric(stringr::str_match(intColNames, "^Int_(.*)$")[,2])
    files <- files[files %in% sl[sl$RAM & !sl$deleted, "ID"]]
    # browser()
    for (i in files) {
      if ((grouped[zeile, paste0("ms2scan_", i)] > 0) &&
          (length(datenList[[i]]@env$intensity) > 1)) {
        MS2[[i]] <- xcms::getMsnScan(datenList[[i]], grouped[zeile, paste0("ms2scan_", i)])
        MS2_max_intensity[i] <- max(MS2[[i]][, 2])
      } else {
        MS2[i] <- list(NULL)
      }
    }
    # get samples which had MS2 available
    samplesWithMS2 <- which(Negate(sapply)(MS2, is.null))
    if (length(samplesWithMS2) == 0)
      aligMs2R$options <- 0 else aligMs2R$options <- c(0, samplesWithMS2)
    req(samplesWithMS2)
    MS2 <- compact(MS2)
    MS2 <- lapply(MS2, as.data.frame)
    
    MS2 <- Map(function(df, num) transform(df, samp = num), MS2, samplesWithMS2)
    MS2 <- do.call("rbind", MS2)
    massRange <- c(min(MS2$mz), grouped[zeile, "mean_mz"] + 1)
    intMax <- max(MS2$intensity)
    # browser()
    if (selection != 0) 
      MS2 <- MS2[MS2$samp == selection, ]
     
    MS2$samp <- as.factor(MS2$samp)
    # 10 most intense fragments for text labels
    ggplot(MS2, aes(mz, intensity, color = samp)) +
      geom_segment(aes(x = mz, xend = mz, y = 0, yend = intensity),
                   stat = "identity", size = .5, alpha = .5) +
      geom_text(aes(label = round(mz,4)), alpha = .6, nudge_y = intMax*.02, check_overlap = T) +
      guides(color = "none") + xlab("m/z") + ylab("Inten. (cps)") +
      scale_x_continuous(limits = massRange) + scale_y_continuous(limits = c(NA, intMax*1.05)) +
      theme_bw(base_size = 14) + annotate("text", -Inf, Inf, vjust = 2, hjust = -.3, label = selection, color = "grey20")
  }
  
  observeEvent(input$AlignmentMS2_click, {
     
    aligMs2R$current <- aligMs2R$current + 1
    if (aligMs2R$current > length(aligMs2R$options))
      aligMs2R$current <- 1
    output$AlignmentMS2 <- renderPlot(
      alignmentMS2Plotten(alignmentIdSelected(), aligMs2R$options[aligMs2R$current]))
  })
  
  # Store table used to make plot in order to use this for hover
  simTrenInten <- reactiveVal()  
  
  SimilarTrends_plotten <- function(ID, highlight = NULL) {
    zeile <- which(grouped[, "alignmentID"] == ID)
    sl <- stringr::str_match(colnames(grouped), "^Int_(\\d+)$")[,2]
    sl <- as.numeric(sl[!is.na(sl)])
    ids <- c(ID, similarTrends())
    inten <- foreach(x = ids, .combine = "rbind") %do% {
      a <- grouped[grouped[,"alignmentID"] == x, grep("^Int_", colnames(grouped))]
      b <- as.factor(sl)
      stopifnot(length(a) == length(b))
      data.frame(samp = b,
                 alignmentID = x,
                 int = a)
    }
    
    # Create highlighted trend
    if (is.null(highlight) || Negate(is.element)(highlight, similarTrends()))
      highlight <- ID
    hlrows <- inten[inten$alignmentID == highlight, ]
    
    simTrenInten(inten)
    inten$alignmentID <- as.factor(inten$alignmentID)
    hlrows$alignmentID <- as.factor(hlrows$alignmentID)
    output$SimilarTrends <- renderPlot({
      ggplot(inten, aes(samp, int, color = alignmentID, group = alignmentID)) + geom_point(shape = 1) + 
        geom_line() + scale_y_continuous(limits = c(0, NA)) + xlab("Sample") + ylab("Inten. (cps)") +
        theme_bw(14) + guides(colour = "none") + geom_point(data = hlrows, shape = 19, size = 3)
    })
  }
  
  annotationMS1Plotten <- function(zeile) {
    thisID <- groupedWithAnnotation[zeile, "alignmentID"]
    dataId <- annotationTable[annotationTable$alignmentID == thisID, "sample", drop = T]
    formula <- annotationTable[annotationTable$alignmentID == thisID, "formula", drop = T]
    adduct <- annotationTable[annotationTable$alignmentID == thisID, "adduct", drop = T]
    if (groupedWithAnnotation[zeile, "multHits"]) {
      zeileMult <- input$multiHitsTable_rows_selected
      dataId <- dataId[zeileMult]
      formula <- formula[zeileMult]
      adduct <- adduct[zeileMult]
    }
    
    req(length(dataId) == 1, length(formula) == 1, length(adduct) == 1)
    
    peakId <- grouped[zeile, paste0("PeakID_", dataId)]
    thisPl <- peaklist[[dataId]]
    scanNumber <- thisPl[thisPl$peak_id_all == peakId, "Scan"]
    spec <- xcms::getScan(datenList[[dataId]], scanNumber)
    comp_mz <- thisPl$mz[thisPl$peak_id_all == peakId]
    RT <- thisPl$RT[thisPl$peak_id_all == peakId]
    comp_int <- thisPl$Intensity[thisPl$peak_id_all == peakId]
    spec <- as.data.frame(spec)
    spec <- spec[spec$mz > comp_mz - 15 & spec$mz < comp_mz + 15, ]
    colnames(spec) <- c("mz", "int")

    # Compute theoretical isotopic distribution
    if (!is.na(formula)) {
      form_charge <- ntsworkflow::correct_formula(formula, adduct)
      mf <- rcdk::get.formula(form_charge$form, charge = form_charge$charge)
      formLab <- mf@string
      isotope_spec <- rcdk::get.isotopes.pattern(mf, minAbund = 0.01)
      isotope_spec <- as.data.frame(isotope_spec)
      colnames(isotope_spec) <- c("mz","int")
      # mirror intensity of MS1 spec
      isotope_spec$int <- isotope_spec$int * -max(comp_int)
    } else {
      isotope_spec <- data.frame(mz = 0, int = 0)
      formLab <- "Unknown"
    }
    
    ggplot(spec, aes(mz, int, label = round(mz,4))) +
      geom_segment(aes(x = mz, xend = mz, y = 0, yend = int),
                   stat = "identity", size = .5, alpha = 0.5) +
      geom_segment(data=isotope_spec, aes(x = mz, xend = mz, y = 0, yend = int),
                   stat = "identity", size = .5, alpha = 0.5)+
      theme_bw() +
      geom_text(data = spec[spec$int > comp_int*0.2, ], check_overlap = TRUE, vjust = -0.5) +
      geom_text(data = isotope_spec[isotope_spec$int < -comp_int*0.2, ],
                check_overlap = TRUE, vjust = 1.5) +
      geom_vline(xintercept = comp_mz, color = "red", alpha = 0.2) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(xlim = c(comp_mz - 2, comp_mz + 8),
                      ylim = c(-1.2*max(comp_int), 1.2*max(comp_int))) +
      ylab("Intensity (abs.)") +
      xlab("m/z (u)") +
      geom_hline(yintercept = 0, color = "blue") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = .1, fill = "orange") +
      annotate("text", x = -Inf, y = -Inf, label = "Calculated-Spec", vjust = -1, hjust = -.1,
               fontface = "bold.italic", alpha = .5) +
      annotate("text", x = -Inf, y = Inf, label = "Data-Spec", vjust = 1.5, hjust = -.1,
               fontface = "bold.italic", alpha = .5) +
      annotate("text", x = Inf, y = Inf, label = formLab, vjust = 1.5, hjust = 1.1)
  }
  annotationMS2Plotten <- function(zeile) {
    # Get annotationTable row as a list to access the values later
    thisID <- groupedWithAnnotation[zeile, "alignmentID"]
    atr <- as.list(annotationTable[annotationTable$alignmentID == thisID, ])
    # if we have multiple hits, these values are the sample anyway, so just take 1st
    dataId <- atr[["sample"]][1]
    comp_mz <- atr[["mzData"]][1]
    
    # If we have multiple hits, need to get the correct values from multiHitsTable
    if (groupedWithAnnotation[zeile, "multHits"]) {
      zeileMult <- input$multiHitsTable_rows_selected
      dbExpId <- atr[["expID"]][zeileMult] 
      score <- atr[["score"]][zeileMult]
      name <- atr[["name"]][zeileMult]
    } else {
      dbExpId <- atr[["expID"]] 
      score <- atr[["score"]]
      name <- atr[["name"]]
    }
    
    req(length(dataId) == 1)
    req(length(dbExpId) == 1)
    # get ms2 spec from data
    dataSpec <- xcms::getMsnScan(datenList[[dataId]], grouped[zeile, paste0("ms2scan_", dataId)])
    dataSpec <- as.data.frame(dataSpec)
    colnames(dataSpec) <- c("mz", "int")
    # get ms2 spec from database
    dbPath <- as.character(parseFilePaths(home, input$ms2Db)$datapath)
    if (length(dbPath) == 0) {
      shiny::showNotification("No database selected", type = "warning")
      dbSpec <- data.frame(mz = 0, int = 0) 
    } else if (grepl("\\.db$", dbPath) && !is.na(dbExpId)) {
      db <- DBI::dbConnect(RSQLite::SQLite(), dbPath)
      dbSpec <- ntsworkflow::dbGetSpectrum(db, dbExpId)
      DBI::dbDisconnect(db)
    } else if (grepl("\\.yaml$", dbPath)) {
      yamldb <- yaml::read_yaml(dbPath)
      frags <- yamldb[[name]]$fragments
      if (is.list(frags)) {
        dbSpec <- as.data.frame(do.call("rbind", frags))
        colnames(dbSpec) <- c("mz", "int")
      } else if (is.vector(frags)) {
        dbSpec <- data.frame(mz = frags, int = 1)
      } else {
        dbSpec <- data.frame(mz = 0, int = 0) 
      }
    } else {
      dbSpec <- data.frame(mz = 0, int = 0) 
      shiny::showNotification("No spectra available", type = "warning")
    }
    shiny::req(is.data.frame(dbSpec))
    
    # get relative intensities
    dataSpec <- dataSpec[dataSpec$mz < comp_mz + 0.015, ]
    dataSpec$int <- dataSpec$int / max(dataSpec$int)
    dbSpec$int <- dbSpec$int / max(dbSpec$int)
    dbSpec$int <- dbSpec$int * -1
    
    ggplot(dataSpec, aes(mz, int, label = round(mz,4))) +
      geom_segment(aes(x = mz, xend = mz, y = 0, yend = int),
                   stat = "identity", size = .5, alpha = .5) +
      geom_segment(data = dbSpec, aes(x = mz, xend = mz, y = 0, yend = int),
                   stat = "identity", size = .5, alpha = .5) +
      theme_bw() +
      geom_text(data = dataSpec[dataSpec$int > 0.01, ], check_overlap = TRUE, vjust = -0.5) +
      geom_text(data = dbSpec[dbSpec$int < -0.01, ], check_overlap = TRUE, vjust = 1.5) +
      geom_vline(xintercept = comp_mz, color = "red", alpha = 0.2) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = c(-1.2, 1.2), xlim = c(0, comp_mz + 5)) +
      ylab("Intensity (relative)") +
      xlab("m/z (u)") +
      geom_hline(yintercept = 0, color = "blue") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = .1,
               fill = "orange") +
      annotate("text", x = -Inf, y = -1.1, label = "DB-Spec", vjust = 0.5, hjust = -.1,
               fontface = "bold.italic", alpha = .5) +
      annotate("text", x = -Inf, y = 1.05, label = "Data-Spec", vjust = 0, hjust = -.1,
               fontface = "bold.italic", alpha = .5) +
      annotate("text", x = -Inf, y = .95, label = paste("Score:", round(score)), vjust = 0.5,
               hjust = -.1, fontface = "bold.italic", alpha = .5, colour = "darkorange4")
  }
  
  
  
  
  # Table fill functions: ####
  
  FillPeakPickTable <- function(filterMaske = list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL), 
                                HighlightType = 0, orderType = list(list(0,"asc")), 
                                displayStart = 0, TableLength = 10, selectedRow = 1) {
    
    outputtable <- data.frame(
      mz=double(),RT=double(),Intensity=double(),Area=double(),SN=double(),Cl=numeric(),
      Br=numeric(),HighlightType=numeric(),RealPeak = numeric())
    
    if (length(peaklist) >= selected) {
      if (!is.null(peaklist[[selected]])) {
        if (nrow(peaklist[[selected]]) > 0) {
          if (length(HighlightType) == 1) HighlightType <- !peaklist[[selected]]$RealPeak
          
          outputtable <- as.data.frame(
            cbind(peaklist[[selected]]$mz,
                  round(peaklist[[selected]]$RT/60,1),
                  round(peaklist[[selected]]$Intensity,0),
                  round(peaklist[[selected]]$Area),
                  peaklist[[selected]]$SN,(peaklist[[selected]]$Cl1>0) | 
                    (peaklist[[selected]]$Cl2>0) | (peaklist[[selected]]$Cl3>0) | 
                    (peaklist[[selected]]$Cl4>0),(peaklist[[selected]]$Br1>0) | 
                    (peaklist[[selected]]$Br2>0),HighlightType,peaklist[[selected]]$RealPeak))
          colnames(outputtable) <- c("mz","RT","Intensity","Area","SN","Cl","Br","HighlightType","RealPeak")
        }
      }
    }
    
    
    outputtable$random <- runif(nrow(outputtable))
    
    columnDefs_list <- list(list(type = "num", width = "20%",className="dt-right", targets = 0),
                              list(type = "num", width = "10%",className="dt-right", targets = 1),
                              list(type = "num", width = "20%",className="dt-right", targets = 2),
                              list(type = "num", width = "20%",className="dt-right", targets = 3),
                              list(type = "num", width = "20%",className="dt-right", targets = 4),
                              list(visible=FALSE, targets = 5),
                              list(visible=FALSE, targets = 6),
                              list(visible=FALSE, targets = 7),
                              list(visible=FALSE, targets = 8),
                              list(visible=FALSE, targets = 9)
    )
    output$PeakPickingTable <- DT::renderDataTable(
      DT::datatable(outputtable,
                    colnames = c("m/z","RT","Intensity","Area","S/N","Cl","Br","HighlightType","RealPeak", "random"),
                    rownames = FALSE,
                    selection = list(mode = 'single', selected = selectedRow),
                    filter = list(position = "bottom", plain = TRUE),
                    width = "200px",
                    options = list(
                      dom = 'ltip',
                      autowidth = TRUE,
                      stateSave = TRUE,
                      stateDuration = -1,
                      pageLength = TableLength,
                      columnDefs = columnDefs_list,
                      searchCols = filterMaske,
                      order = orderType,
                      displayStart = displayStart
                    )
      )
      %>% formatRound(columns = "mz",digits = 4)
      %>% formatRound(columns = "RT",digits = 1)
      %>% formatRound(columns = "Intensity",digits = 0)
      %>% formatRound(columns = "SN",digits = 0)
      %>% formatStyle(columns = "HighlightType", target = "row", 
                      backgroundColor = styleEqual(c(1,2),c("yellow","red")))
      
    )
  }
  
  
  FillAlignmentTable <- function(gr = grouped) {
     
    # Get general intensity, look at first intensity column
    firstIntCol <- grep("^Int_", colnames(gr))[1]
    upperQ <- quantile(gr[, firstIntCol], type = 1, na.rm = TRUE, probs = 0.75) 
    intDigits <- ifelse(upperQ < 10, 2, 0)
    grouped2 <- cbind(gr[, c("mean_mz", "mean_RT", "alignmentID", "Gruppe"), drop=F],
                      round(gr[, grep("^Int_", colnames(gr)), drop = FALSE], intDigits))
    
    grouped2[,"mean_RT"] <- grouped2[,"mean_RT"] / 60

    spaltennamen <- c("m/z","RT","ID", "Group", colnames(grouped2)[-(1:4)])

    # Merge with annotation table if there is one, ignore multi hits entries, just take first one
    if (!is.null(annotationTable)) {
       
      grouped2 <- as.data.frame(grouped2)
      
      if("samplename" %in% colnames(annotationTable)){
        
        Label <- paste(annotationTable$samplename,annotationTable$sample_type,sep=", ")
        annot <- annotationTable[, c("alignmentID", "name" )]
        annot <- cbind(annot,Label)
        spaltennamen <- c(spaltennamen,"name","label")
        
      } else {
        annot <- annotationTable[, c("alignmentID", "name")]
        spaltennamen <- c(spaltennamen,"name")
      }
      
      notDup <- !duplicated(annotationTable$alignmentID)
      annot <- annot[notDup, ]
      grouped2 <- merge(grouped2, annot, all.x = TRUE, by = "alignmentID")
      
    }
    # Remove ID column
    grouped2 <- subset(grouped2, , -alignmentID)
    spaltennamen <- spaltennamen[spaltennamen != "ID"]
    filterMask <- vector("list", length = ncol(grouped2))
    output$AlignmentTable <- DT::renderDataTable(
      DT::datatable(grouped2,
                    colnames = spaltennamen,
                    rownames = FALSE,
                    selection = 'single',
                    width = "100%",
                    filter = list(position = "bottom", plain = TRUE),
                    options = list(
                      scrollX = TRUE,
                      dom = 'ltip',
                      searchCols = filterMask
                      )
      )
      %>% formatRound(columns = "mean_mz",digits = 4)
      %>% formatRound(columns = "mean_RT",digits = 1)
    )
  }

  # function to get intensity from grouped table
  getIntensity <- function(row, ofTable) {
    fileNum <- ofTable[row, "sample"]
    if (!is.na(fileNum)) 
      grouped[row, paste0("Int_", fileNum)] else NA
  }
  
  similarTrends <- reactiveVal()
  FillSimilarTrendsTable <- function(ID) {
    req(grouped)
    req(is.numeric(grouped[,seq(8,ncol(grouped), by = 6)]))
    req(!is.na(input$trendFac))
    req(is.matrix(grouped[,seq(8,ncol(grouped), by = 6)]))
    # browser()
    zeile <- which(grouped[, "alignmentID"] == ID)
    req(zeile >= 1)
    similarTrendRows <- which(correlates_with_r(grouped[,grep("^Int_", colnames(grouped))], zeile, input$trendFac))
    similarTrends(grouped[similarTrendRows, "alignmentID"])
    grouped2 <- cbind(
      grouped[similarTrendRows, c("mean_mz", "mean_RT", "Gruppe"), drop = F],
      round(grouped[similarTrendRows,grep("^Int_", colnames(grouped))],0)
    )
    grouped2[, "mean_RT"] <- grouped2[, "mean_RT"]/60
    spaltennamen <- c("m/z","RT", "Compon.", colnames(grouped2)[-1:-3])
    ti <- sprintf("Similar to m/z: %.4f Da, RT %.1f min.", grouped[zeile, "mean_mz"], 
                  grouped[zeile, "mean_RT"]/60)
    if ("name" %in% colnames(aligR()))
      ti <- paste(ti, aligR()[zeile, "name"], sep = ", ")
    output$simTrenHead <- renderText(ti)
    output$SimilarTrendsTable <- DT::renderDataTable(
      DT::datatable(grouped2,
                    colnames = spaltennamen,
                    rownames = FALSE,
                    selection = 'single',
                    width = "100%",
                    options = list(scrollX = TRUE)
      )
      %>% formatRound(columns = "mean_mz",digits = 4)
      %>% formatRound(columns = "mean_RT",digits = 1)
    )
  }
  
  FillAnnotationTable <- function() {
     
    req(annotationTable)
    groupedPlus <- as.data.frame(grouped)
    
    # include column indicating duplicates in the annotationTable
    dupAZ <- annotationTable[duplicated(annotationTable$alignmentID), 
                             "alignmentID", drop = TRUE]
    # remove duplicates from the annotationTable
    annotationTable2 <- annotationTable[!duplicated(annotationTable$alignmentID), ]
    groupedPlus <- merge(groupedPlus, annotationTable2, by = "alignmentID", all.x = TRUE)
    groupedPlus$multHits <- groupedPlus$alignmentID %in% dupAZ
    
    if("samplename" %in% colnames(groupedPlus)){
      
      groupedPlusRed <- subset(
        groupedPlus, , c(mean_mz, mean_RT, name, samplename, sample_type, score,
                         sample, multHits, alignmentID, mzDB, rtDB)
      )
      
    } else {
    
      groupedPlusRed <- subset(
        groupedPlus, , c(mean_mz, mean_RT, name, CAS, score, sample, multHits, 
                         alignmentID, mzDB, rtDB)
      )
    }
    
    groupedPlusRed$mean_RT <- groupedPlusRed$mean_RT / 60
    groupedPlusRed$mzDev <- groupedPlusRed$mean_mz - groupedPlusRed$mzDB
    groupedPlusRed$rtDev <- groupedPlusRed$mean_RT - groupedPlusRed$rtDB
    groupedPlusRed$mzDB <- NULL
    groupedPlusRed$rtDB <- NULL
    
    groupedPlusRed$highest_int <- vapply(
      seq_len(nrow(groupedPlusRed)), 
      getIntensity, 
      numeric(1), 
      ofTable = groupedPlusRed
    )
    
    stopifnot(nrow(groupedPlusRed) == nrow(grouped))
     
    groupedWithAnnotation <<- groupedPlusRed
    output$AlignmentAnnotTable <- DT::renderDataTable(
      DT::datatable(groupedPlusRed,
                    rownames = FALSE,
                    selection = 'single',
                    width = "100%",
                    options = list(scrollX = TRUE)
      )
      %>% formatRound(columns = "mean_mz", digits = 4)
      %>% formatRound(columns = "mean_RT", digits = 2)
      %>% formatRound(columns = "mzDev", digits = 4)
      %>% formatRound(columns = "rtDev", digits = 2)
      %>% formatRound(columns = "highest_int", digits = 2)
      %>% formatStyle("multHits", target = "row", backgroundColor = styleEqual(1,"tan"))
    )
  }
  clearAnnotationTable <- function() {
    output$AlignmentAnnotTable <- DT::renderDataTable(data.frame(empty = "empty"))
  }
  
  fillMultiHitsTable <- function() {
    az <- groupedWithAnnotation[input$AlignmentAnnotTable_rows_selected, "alignmentID"]
    multiHits <- annotationTable[annotationTable$alignmentID == az, ]
    groupedPlus <- as.data.frame(grouped)
    multiHits <- merge(groupedPlus, multiHits, by = "alignmentID")
    multiHits$multHits <- TRUE
    
    if("samplename" %in% colnames(multiHits)){
      
      multiHits <- multiHits[, c("mean_mz", "mean_RT", "name", "samplename", "sample_type", "score", "sample", "multHits",
                                 paste0("Int_", multiHits[1, "sample"]))]
      colnames(multiHits)[9] <- "highest_int"
      
    } else {
    
    
    multiHits <- multiHits[, c("mean_mz", "mean_RT", "name", "CAS", "score", "sample", "multHits",
                               paste0("Int_", multiHits[1, "sample"]))]
    colnames(multiHits)[8] <- "highest_int"
    }
    
    multiHits$mean_RT <- multiHits$mean_RT / 60
    
    multiHitsTable <<- multiHits
    
    output$multiHitsTable <- DT::renderDataTable(
      DT::datatable(multiHits,
                    rownames = FALSE,
                    selection = list(mode="single", selected=1),
                    width = "100%",
                    options = list(dom = 't')
                    
      )  
      %>% formatRound(columns = "mean_mz", digits = 4)
      %>% formatRound(columns = "mean_RT", digits = 2)
      %>% formatRound(columns = "highest_int", digits = 2)
    )
    
  }
  
  # Other functions: ####
  
  compact <- function(x) {
    Filter(Negate(is.null), x)
  }
  
  deleteSample <- function(zeile) {
    datenList[[zeile]] <<- NULL
    peaklist[[zeile]] <<- NULL
    peakPickSettings[[zeile]] <<- NULL
    sampleList[zeile, "RAM"] <<- FALSE
    sampleList[zeile, "deleted"] <<- TRUE
    selected <<- NULL
    selectedR(NULL)
  }
  
  ppProgress <- function() {
    ppProgressFile <- file.path(globalwd, paste0("Peak-picking-progress-", Sys.Date(), ".txt"))
    cat("Peak-picking progress indicator:\n", file = ppProgressFile)
    ppProgressFile  # return file name
  }
  
  ppProgressInc <- function(completed, filen) {
    progressText <- paste(completed, date(), "completed\n")
    cat(progressText, file = filen, append = T)
  }
  
  # Load environment ####
  observeEvent(input$loadEnv, {
    if (any(c("sampleList", "datenList") %notin% ls(globalenv())))
      showNotification("No data found in environment")
    req(all(c("sampleList", "datenList") %in% ls(globalenv())))
    loadp <- shiny::Progress$new(session, 0, 1)
    loadp$set(1, message = "Loading data")
    selected <<- NULL
    selectedR(NULL)
    if (!is.null(grouped))
      FillAlignmentTable()
    if (!is.null(annotationTable))
      FillAnnotationTable()
    sampleListR(sampleList)
    loadp$close()
  })
  
  # Memory usage ####
  output$memUsage <- renderText({
    input$updateMem
    x <- system2('free', args = "-m",stdout=TRUE)
    used <- as.numeric(strsplit(x[2], " +")[[1]][3])
    tot <- as.numeric(strsplit(x[2], " +")[[1]][2])
    percentUsed <- round(100 * used / tot)
    paste0("Current memory usage: ", percentUsed, "%")
  })
  
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
  
  # Batch peak picking and alignment ####
  observeEvent(input$BatchProcessStart, {
    # Get peak picking settings
    mz_min <- isolate(input$PeakPick_massrange[1])
    mz_max <- isolate(input$PeakPick_massrange[2])
    mz_step <- isolate(input$PeakPick_mzstep)
    rt_min <- isolate(input$PeakPick_RTrange[1])*60
    rt_max <- isolate(input$PeakPick_RTrange[2])*60
    PeakWidth_min <- isolate(input$PeakPick_PeakWidthRange[1])
    PeakWidth_max <- isolate(input$PeakPick_PeakWidthRange[2])
    peak_NoiseScans <- isolate(input$PeakPick_NoiseScans)
    sn <- isolate(input$PeakPick_SN)
    int_threshold <- isolate(input$PeakPick_IntTresh)
    precursormzTol <- isolate(input$PPMS2Fit_ppm)
    componentization_ppm <- isolate(input$componentization_ppm)
    componentization_rt_tol <- isolate(input$componentization_rt_tol)
    componentization_rt_tol_l <- isolate(input$componentization_rt_tol_l)
    componentization_rt_tol_r <- isolate(input$componentization_rt_tol_r)
    componentization_rt_tol_sum <- isolate(input$componentization_rt_tol_sum)
    maxNumPeaks <- isolate(input$maxPeaks)
    
    # Save settings
    settings <- list()
    settings$massrange <- c(mz_min,mz_max)
    settings$mz_step <- mz_step
    settings$rtrange <- c(rt_min,rt_max)
    settings$peakwidth <- c(PeakWidth_min,PeakWidth_max)
    settings$NoiseScans <- peak_NoiseScans
    settings$sn <- sn
    settings$int_threshold <- int_threshold
    settings$precursormzTol <- precursormzTol
    settings$ppm <- componentization_ppm
    settings$RT_Tol <- componentization_rt_tol
    settings$PPTableRow <- 1  
    settings$orderType <- list(list(0,"asc"))
    settings$TableLength <- 10
    settings$displayStart <- 0
    settings$maxNumPeaks <- maxNumPeaks
    
    
    sampleAmountOffset <- length(datenList)
    
    req(length(batchFiles) > 0) 

    batchp <- Progress$new(session, 0, 1)
    batchp$set(0.1, message = "Overall progress", detail = "Set-up")
    
    # Append data to sampleList
    # optMzStep will be added after peakpicking
    newIds <- if(length(sampleList$ID) == 0) 
      seq_along(batchFiles) else max(sampleList$ID) + seq_along(batchFiles)
    newRows <- Map(function(newID, filep, inRam, sampTypeNum) {
      data.frame(ID = newID, File = filep, RAM = inRam, 
                 deleted = F, sampleType = sampleTypes[sampTypeNum], 
                 optMzStep = NA, stringsAsFactors = F )
    }, newIds, batchFiles, batchFilesRAM, batchFilesSampleType)
    sampleList <<- do.call("rbind", append(list(sampleList), newRows))
    sampleListR(sampleList)  # update reactiveVal
    
    # foreach goes through batch files and returns a list with all the results
    # as a list (headerlist, datenList, peaklist, settings)
    batchp$set(0.3, detail = "Peak-picking")
    
    cl <- parallel::makeForkCluster(input$numcores)
    doParallel::registerDoParallel(cl)
    mcoptions <- list(preschedule = FALSE)
    
    combiRes <- foreach(filep = batchFiles, 
                        batchSampleType = batchFilesSampleType, 
                        inRam = batchFilesRAM, .export = "input",
                        .options.multicore = mcoptions) %dopar% {
      datenx <- xcms::xcmsRaw(filep, includeMSn = TRUE)
      headerEntry <- list()
      headerEntry$mz_step_recommendation <- ntsworkflow::optimumMzStep(datenx, 0.90)
      headerEntry$sampleType <- sampleTypes[batchSampleType]
      headerEntry$samplingDate <- as.Date(file.mtime(filep))
      headerEntry$samplingSite <- ""
      headerEntry$samplingPosition <- 0
      
      ppsettingsx <- settings
      
      # If sample is a blank, it should be picked at 1/10th intensity and 1/2 sn ratio,
      # to ensure that all peaks near threshold are captured in blank.
      int_threshold_i <- int_threshold
      sn_i <- sn
      if (headerEntry$sampleType == "Blank" && isolate(input$batchBlank)) {
        int_threshold_i <- int_threshold / isolate(input$batchBlankIntFactor)
        sn_i <- sn / isolate(input$batchBlankSnFactor)
        ppsettingsx$int_threshold <- int_threshold_i
        ppsettingsx$sn <- sn_i
      }
      pl <- ntsworkflow::FindPeaks_BfG(
        daten = datenx, 
        mz_min, 
        mz_max, 
        mz_step,
        rt_min,
        rt_max,
        sn_i,
        int_threshold_i,
        peak_NoiseScans,
        precursormzTol,
        PeakWidth_min,
        PeakWidth_max,
        maxNumPeaks
      )
      
      # Remove samples from RAM
      if (!inRam) {
        datenx@env$intensity <- 0
        datenx@env$mz <- 0
        datenx@env$profile <- 0
        datenx@env$msnIntensity <- 0
        datenx@env$msnMz <- 0
      }  
      gc() 
      
      # Collect the results in a list
      resultsList <- list(dat = datenx, headE = headerEntry, 
                          ppset = ppsettingsx, peakl = pl)
      resultsList
    }
    
    # extract data from workers and save to global environment
    newdat <- lapply(combiRes, function(x) x[["dat"]])
    newheadE <- lapply(combiRes, function(x) x[["headE"]])
    newppset <- lapply(combiRes, function(x) x[["ppset"]])
    newpeakl <- lapply(combiRes, function(x) x[["peakl"]])
    datenList <<- append(datenList, newdat)
    peakPickSettings <<- append(peakPickSettings, newppset)
    peaklist <<- append(peaklist, newpeakl)
    
    sampleList[is.na(sampleList$optMzStep), "optMzStep"] <<- vapply(
      newheadE, function(x) x[["mz_step_recommendation"]], numeric(1))
    sampleListR(sampleList)

    # Componentization for all peaklists at the same time (parallel)
    batchp$set(0.6, detail = "Componentization")
    
    peaklist <<- foreach(peaklistx = peaklist, datenx = datenList,
                         .export = "input") %dopar% {
      
      if (is.null(peaklist) || nrow(peaklistx) == 0)
        return(peaklistx)
      
      newpl <- ntsworkflow::componentization_BfG(
        peaklistx, 
        daten = datenx, 
        ppm = componentization_ppm, 
        Grenzwert_RT = componentization_rt_tol, 
        Grenzwert_FWHM_left = componentization_rt_tol_l, 
        Grenzwert_FWHM_right = componentization_rt_tol_r, 
        Summe_all = componentization_rt_tol_sum,
        adjust_tolerance = isolate(input$Componentization_Dynamic_Tolerance))
      
      newpl$RealPeak <- TRUE
      newpl
    }
    
    # For logging
    adjust_tolerance <- input$Componentization_Dynamic_Tolerance
    log_file <<- create_log_file()
    
    # Start alignment if requested
    if (input$batchAlign) {
      batchp$set(0.8, detail = "Alignment")
      if (is.null(alignPeaksEvent()))
        alignPeaksEvent(1) else alignPeaksEvent(alignPeaksEvent() + 1)
    }
    
    # Save results if requested
    if (input$saveBatch) {
      batchp$set(0.9, detail = "Saving")
      if (is.null(saveLater()))
        saveLater(1) else saveLater(saveLater() + 1)
    }
    
    batchp$close()
    parallel::stopCluster(cl)
    rm(cl)
  })
  
  # Add Sample ####
  home3 <- c(WD = globalwd, Home = fs::path_home())
  shinyFileChoose(input, 'addSample', roots=home3, session=session,
                  filetypes=c('mzML', 'mzXML'))
  
  sampleListR <- reactiveVal(NULL)
  
  observeEvent(input$addSample, {
    req(is.list(input$addSample))
    dateiInfo <- parseFilePaths(home, input$addSample)
    dateiInfo <- as.data.frame(lapply(dateiInfo, as.character), stringsAsFactors = FALSE)
    dateiAlle <- dateiInfo$datapath
    req(dateiAlle)
    req(all(grepl("\\.mzX?ML$", dateiAlle)))  # all files must be mzML or mzXML
    pfad <<- dirname(dateiAlle[1])
    datSeq <- nrow(sampleList) + 1:length(dateiAlle)
    loadp <- shiny::Progress$new(session, 0, 1)
    loadp$set(0.4, message = "Loading...")
    
    # parallel data loading
    datenList[datSeq] <<- parallel::mclapply(dateiAlle, xcms::xcmsRaw, includeMSn = TRUE,
                                   mc.preschedule = TRUE, mc.cores = input$numcores)
    
    br <- input$blankRegex
    createRow <- function(newId, dateiPfad, datRaw) {
      bl <- if (br == "") FALSE else grepl(br, basename(dateiPfad))
      data.frame(ID = newId, File = dateiPfad, 
                           sampleType = sampleTypes[ifelse(bl, 2, 1)],
                           optMzStep = ntsworkflow::optimumMzStep(datRaw, 0.90),
                           RAM = T, 
                           deleted = F, 
                           stringsAsFactors=F)
    
    }
    newRowsList <- Map(createRow, datSeq, dateiAlle, datenList[datSeq])
    newRows <- do.call("rbind", newRowsList)
      
    sampleList <<- rbind(sampleList, newRows)
    loadp$set(0.5)
    
    sampleListR(sampleList)
    
    loadp$close()
  })
  
  # SampleTable Output ####
  output$sampleTable <- DT::renderDataTable({
     
    req(nrow(sampleListR()) >= 1)
    sampleList_output <- sampleListR()
    sampleList_output[, "File"] <- basename(sampleList_output[, "File"])
    sampleList_output <- sampleList_output[!sampleList_output$deleted, ]
    sampleList_output <- sampleList_output[, c("ID", "File", "RAM")]
    req(nrow(sampleList_output) >= 1)
    DT::datatable(sampleList_output, 
                  selection = 'single',
                  colnames = c("ID","File", "RAM"),
                  rownames = FALSE, 
                  options = list(
                    dom = 't',
                    paging = FALSE, 
                    searching = FALSE, 
                    autowidth = TRUE,
                    columnDefs = list(
                      list(targets = 1, render = JS("function(data, type, row, meta) {",
                                                    "return type === 'display' && data.length > 20 ?",
                                                    "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                    "}")),
                    list(targets = 2, visible = F)
                    )
                  ),
                  callback = JS('table.page("next").draw(false);')
    ) %>% DT::formatStyle(
      "RAM", 
      target = "row", 
      backgroundColor = styleEqual(c(1, 0), c("#66FF66", "#FF6666"))
    )
  })
  
  # SampleTable observer ####
  selectedR <- reactiveVal()
  observeEvent(input$sampleTable_cell_clicked, {
    
    if (!is.null(selected)) {
      if (length(peakPickSettings) >= selected) {
        TableState <- isolate(input$PeakPickingTable_state)
        peakPickSettings[[selected]]$orderType <<- TableState$order
        peakPickSettings[[selected]]$TableLength <<- TableState$length
        zeilen <<- isolate(input$PeakPickingTable_rows_all)
        angezeigteZeilen <<- isolate(input$PeakPickingTable_rows_current)
        peakPickSettings[[selected]]$displayStart <<- which(zeilen == angezeigteZeilen[1])-1
      }  
    }  
    currentRow <- input$sampleTable_rows_selected
    if (!is.null(currentRow)) {
      # Changed selection so that you select the right sample id, not the row of the table
      sl <- sampleListR()[!sampleListR()$deleted, ]
      selected <<- sl[input$sampleTable_rows_selected, "ID"]
      selectedR(sl[input$sampleTable_rows_selected, "ID"])
    }
    
    if (length(selected) == 1) {
      output$sampleID <- renderText(selected)
    
      updateSelectInput(session, "sampleType", choices = sampleTypes,
                        selected = sampleListR()[sampleListR()$ID == selectedR(), "sampleType"])

      
      updateTextInput(session, "samplePath", label = NULL, value = unname(sampleList[selected, "File"]))
      
      output$ms1scans <- renderText(length(datenList[[selected]]@tic))
      output$ms2scans <- renderText(length(datenList[[selected]]@msnRt))
      output$ms1RT <-
        renderText(c(round(min(
          datenList[[selected]]@scantime / 60
        ), 1), "-", round(max(
          datenList[[selected]]@scantime / 60
        ), 1), "min"))
      output$ms2RT <-
        renderText(c(round(min(
          datenList[[selected]]@msnRt / 60
        ), 1), "-", round(max(
          datenList[[selected]]@msnRt / 60
        ), 1), "min"))
      
      output$mzsteprecommendation <- renderText(paste("Recomm.: ", sampleListR()[sampleListR()$ID == selectedR(), "optMzStep"]))
      output$TICplot <- renderPlot({
        plot(x = datenList[[selected]]@scantime/60, y = datenList[[selected]]@tic, type = "l", xlim = TICplot_ranges$x, ylim = TICplot_ranges$y,main = "TIC",xlab="RT",ylab="Intensity",xaxs = "i",yaxs="i")
        lines(x = TICplot_RT$x, y = max(datenList[[selected]]@tic), type = "h", col = "red")
      })
      XICplot_plotten()   
      
      output$PeakPickXIC <- renderPlot(NULL)
      output$PeakPickMS <- renderPlot(NULL)
      
      # PeakPicking Table fill
      updateSliderInput(session,
                        "PeakPick_massrange",
                        min = datenList[[selected]]@mzrange[1],
                        max = datenList[[selected]]@mzrange[2])
      updateSliderInput(session, "PeakPick_RTrange", 
                        max = round(max(datenList[[selected]]@scantime) / 60))
      if (length(peakPickSettings) >= selected) { 
        settings <- peakPickSettings[[selected]]
        if (!is.null(settings$massrange)) {
          updateSliderInput(session, "PeakPick_massrange",value = settings$massrange)
          updateNumericInput(session, "PeakPick_mzstep", value = settings$mz_step)
          updateSliderInput(session, "PeakPick_RTrange",value = settings$rtrange/60)    
          updateSliderInput(session, "PeakPick_PeakWidthRange" ,value = settings$peakwidth)
          updateNumericInput(session, "PeakPick_NoiseScans", value = settings$NoiseScans)
          updateNumericInput(session, "PeakPick_SN", value = settings$sn)
          updateNumericInput(session, "PeakPick_IntTresh", value = settings$int_threshold)
          updateNumericInput(session, "PPMS2Fit_ppm", value = settings$precursormzTol)
          updateNumericInput(session, "componentization_ppm", value = settings$ppm)
          updateNumericInput(session, "componentization_rt_tol", value = settings$RT_Tol)
          updateNumericInput(session, "maxPeaks", value = settings$maxPeaks)
          updateCheckboxInput(
            session, "Componentization_Dynamic_Tolerance", 
            value = settings$componentization_dynamic_tolerance
          )
          
          if (settings$PPTableRow > 0) {  
            FillPeakPickTable(orderType = settings$orderType, TableLength = settings$TableLength, selectedRow = settings$PPTableRow, displayStart = settings$displayStart)
            
            if (sampleList[selected, "RAM"]) {
              output$PeakPickXIC <- renderPlot(PeakPickXIC_plotten(settings$PPTableRow, selectedR()))
              PeakPickMSPlotten(settings$PPTableRow)
              output$PeakPickMS2 <- renderPlot(PeakPickMS2Plotten(settings$PPTableRow))
            }
          } else {
            FillPeakPickTable(orderType = settings$orderType, TableLength = settings$TableLength)
            
            output$ComponentsXIC <- renderPlot(ComponentXIC_plotten(1))
          }
        } else {
          FillPeakPickTable()
          output$ComponentsXIC <- renderPlot(ComponentXIC_plotten(1))
        }
      }

    } 
  })
  
  # Load ####
  shinyFileChoose(input, 'load', roots=home3, session=session, filetypes = "RDS")
  observeEvent(input$load, {
    req(is.list(input$load))
    
    loadp <- shiny::Progress$new(session, 0, 1)
    loadp$set(value = 1, message = "Loading...", detail = "Reading file")
    
    sl <- readRDS(as.character(parseFilePaths(home, input$load)$datapath)) 
    if (!("sampleList" %in% names(sl)) ||
        !is.data.frame(sl[["sampleList"]]) ||
        nrow(sl[["sampleList"]]) == 0)
      showNotification("The file is of the wrong format or is corrupt.", 
                       type = "warning")
    req("sampleList" %in% names(sl))
  
    
    if ("sampleList" %in% names(sl)) {
      sampleList <<- sl[["sampleList"]]
      sampleListR(sampleList)
    }
    
    if ("datenList" %in% names(sl))
      datenList <<- sl[["datenList"]]
    
    if ("peaklist" %in% names(sl))
      peaklist <<- sl[["peaklist"]]
    
    if ("peakPickSettings" %in% names(sl))
      peakPickSettings <<- sl[["peakPickSettings"]]
    
    if ("log_file" %in% names(sl)) 
      log_file <<- sl[["log_file"]]
    
    if ("grouped" %in% names(sl))
      grouped <<- sl[["grouped"]]
    if (!is.null(grouped))
      FillAlignmentTable()
    
    if ("clusteringData" %in% names(sl))
      clusteringData <<- sl[["clusteringData"]]
    
    if ("annotationTable" %in% names(sl))
      annotationTable <<- sl[["annotationTable"]]
    if (!is.null(annotationTable))
      FillAnnotationTable()
    
    loadp$set(detail = "Reading data files")
    
    ids <- sampleList[sampleList$RAM & !sampleList$deleted, "ID"]
    paths <- sampleList[sampleList$RAM & !sampleList$deleted, "File"]

    datenList[ids] <<- parallel::mclapply(paths, xcms::xcmsRaw, includeMSn=TRUE,
                                          mc.preschedule = FALSE, 
                                          mc.cores = input$numcores)
   
    # Add some other necessary variables
    sampleTypes <<- c("Unknown","Blank","Standard")
    loadp$close()   
  })
  
  # Save ####
  save_data <- function(speicherort) {
    
    savep <- shiny::Progress$new(session, 0, 1)
    savep$set(1, message = "Saving...")
    save_var <- list()
    if (exists("log_file")) {
      log_file_save <- jsonlite::toJSON(log_file, pretty = TRUE)
      write(log_file_save,file=paste0(gsub(".RDS","",speicherort,fixed=TRUE),"_log_file.json"))
      save_var[["log_file"]] <- log_file
    }
    
    save_var[["sampleList"]] <- sampleList
    save_var[["peaklist"]] <- peaklist
    save_var[["peakPickSettings"]] <- peakPickSettings
    save_var[["grouped"]] <- grouped
    save_var[["annotationTable"]] <- annotationTable
    if (exists("clusteringData"))
      save_var[["clusteringData"]] <- clusteringData
    
    datenList_save <- list()
    for (i in sampleList[!sampleList$deleted, "ID"]) {
      # Need to copy all stuff individually. copying the whole list and 
      # deleting the env also delete the env of the original list
      text_xcmsRaw <- "xcmsRaw"
      attr(text_xcmsRaw, "package") <- "xcms"  # need to define class's package as attribute
      datenList_save[[i]] <- new(text_xcmsRaw)
      datenList_save[[i]]@tic <- datenList[[i]]@tic
      datenList_save[[i]]@scantime <- datenList[[i]]@scantime
      datenList_save[[i]]@scanindex <- datenList[[i]]@scanindex
      datenList_save[[i]]@polarity <- datenList[[i]]@polarity
      datenList_save[[i]]@acquisitionNum <- datenList[[i]]@acquisitionNum
      datenList_save[[i]]@profmethod <- datenList[[i]]@profmethod
      datenList_save[[i]]@profparam <- datenList[[i]]@profparam
      datenList_save[[i]]@mzrange <- datenList[[i]]@mzrange
      datenList_save[[i]]@gradient <- datenList[[i]]@gradient
      datenList_save[[i]]@msnScanindex <- datenList[[i]]@msnScanindex
      datenList_save[[i]]@msnAcquisitionNum <- datenList[[i]]@msnAcquisitionNum
      datenList_save[[i]]@msnPrecursorScan <- datenList[[i]]@msnPrecursorScan
      datenList_save[[i]]@msnLevel <- datenList[[i]]@msnLevel
      datenList_save[[i]]@msnRt <- datenList[[i]]@msnRt
      datenList_save[[i]]@msnPrecursorMz <- datenList[[i]]@msnPrecursorMz
      datenList_save[[i]]@msnPrecursorIntensity <- datenList[[i]]@msnPrecursorIntensity
      datenList_save[[i]]@msnPrecursorCharge <- datenList[[i]]@msnPrecursorCharge
      datenList_save[[i]]@msnCollisionEnergy <- datenList[[i]]@msnCollisionEnergy
      datenList_save[[i]]@filepath <- datenList[[i]]@filepath
      datenList_save[[i]]@scanrange <- datenList[[i]]@scanrange
      datenList_save[[i]]@mslevel <- datenList[[i]]@mslevel
      datenList_save[[i]]@env$intensity <- 0
      datenList_save[[i]]@env$mz <- 0
      datenList_save[[i]]@env$profile <- 0
      datenList_save[[i]]@env$msnIntensity <- 0
      datenList_save[[i]]@env$msnMz <- 0
    }
    save_var[["datenList"]] <- datenList_save
    saveRDS(save_var, speicherort)
    savep$close()
  }
  
  shinyFileSave(input, 'save', roots=home, session=session)
  observeEvent(input$save, {
    req(is.list(input$save))
    speicherort <- as.character(parseSavePath(home, input$save)$datapath)
    req(length(speicherort) > 0)
    save_data(speicherort)
  })
  
  # ignoreInit is needed otherwise will set to unknown when loading
  observeEvent(input$sampleType, {
    if (length(selectedR()) == 1) {
      #headerList[[selectedR()]]$sampleType <<- input$sampleType  # deprecated
      sampleList[sampleList$ID == selectedR(), "sampleType"] <<- input$sampleType
      sampleListR(sampleList)
    }
      
  }, ignoreInit = TRUE)
  
  # Change path ####
  observeEvent(input$changeSamplePath, {
    oldPath <- sampleListR()[sampleListR()$ID == selectedR(), "File"]
    if (length(selectedR()) == 1)  
      newPath <- input$samplePath  # the user must manually change the sample path in the text box
    req(newPath)
    if (!file.exists(newPath))
      shiny::showNotification("File not found", type = "error")
    req(file.exists(newPath))
    if (newPath == oldPath)
      shiny::showNotification("Path has not been changed", type = "error")
    req(newPath != oldPath)
    
    withProgress(message = 'Loading', value = 1, {
      newData <- xcms::xcmsRaw(newPath,includeMSn = TRUE)
    })
    if (length(datenList[[selected]]@tic) == length(newData@tic) &&
          all(datenList[[selected]]@tic == newData@tic)) {
      sampleList[selected, "File"] <<- newPath
      sampleList[selected, "RAM"] <<- TRUE
      datenList[[selected]] <<- newData
      sampleListR(sampleList)
    } else {
      shiny::showNotification("The replacement is not the same as the original", type = "error")
    }
  })
  
  # Remove and reload sample from RAM ####
  removeFromRamFun <- function(fileID) {
    datenList[[fileID]]@env$intensity <<- 0
      datenList[[fileID]]@env$mz <<- 0
      datenList[[fileID]]@env$profile <<- 0
      datenList[[fileID]]@env$msnIntensity <<- 0
      datenList[[fileID]]@env$msnMz <<- 0
      sampleList[fileID, "RAM"] <<- FALSE
      sampleListR(sampleList)
  }
  
  observeEvent(input$removeFromRAM, {
    if (length(selected) == 1) {  
      removeFromRamFun(fileID = selected)
    }  
    # Linux seems to not free up memory after delete
    gc()  
  })
  
  observeEvent(input$removeAllFromRAM, {
     filesToRemove <- subset(sampleList, RAM & (!deleted), ID, drop = TRUE)
     remProg <- Progress$new(session, 0, length(filesToRemove))
     for (fileNum in filesToRemove) {
       removeFromRamFun(fileID = fileNum)
       remProg$inc(1, "Removing files from memory")
     }
     gc()
     remProg$close()
  })
  
  observeEvent(input$reloadToRAM, {
    if (length(selected) == 1) {
      if (!sampleList[selected, "RAM"] && !sampleList[selected, "deleted"]) { 
        withProgress(message = 'Loading', value = 1, {
          fid <- selectedR()
          datenList[[fid]] <<- xcms::xcmsRaw(datenList[[fid]]@filepath@.Data, includeMSn = TRUE)
          sampleList[fid, "RAM"] <<- TRUE
          sampleListR(sampleList)
        })
      }
    }  
  })
  
  observeEvent(input$reloadAllToRAM, {
    tryCatch(
      error = function(cnd) 
        showNotification(paste("Bad input:", conditionMessage(cnd)), type = "error"),
      {
        if (input$reloadChoice == "") {
          filesToLoad <- subset(sampleList, (!RAM) & (!deleted), ID, drop = TRUE)
        } else {
          filesToLoad <- as.numeric(eval(parse(text = input$reloadChoice)))
        }
        loadProg <- Progress$new(session, 0, length(filesToLoad))
        loadProg$set(message = "Loading files to memory")
        
        datenList[filesToLoad] <<- parallel::mclapply(filesToLoad, function(i) {
          xcms::xcmsRaw(datenList[[i]]@filepath@.Data, includeMSn = TRUE)
        }, mc.preschedule = FALSE, mc.cores = input$numcores)
        
        sampleList[filesToLoad, "RAM"] <<- TRUE
        sampleListR(sampleList)
        loadProg$close()
      })
  })
  
  #### deleteSample ####
  observeEvent(input$deleteSample, {
    sel <- selectedR()
    if (is.null(sel))
      showNotification("First select sample", type = "error")
    req(!is.null(sel))
    # set elements of deleted samples to NULL
    datenList[sel] <<- list(NULL)
    peaklist[sel] <<- list(NULL)
    peakPickSettings[sel] <<- list(NULL)
    sampleList[sel, "deleted"] <<- T
    sampleList[sel, "RAM"] <<- F
    sampleListR(sampleList)
    
    selected <<- NULL
    selectedR(NULL)
    
    if (!is.null(grouped)) {
      # Remove columns and rows of grouped table
      grouped <<- grouped[, -grep(sprintf("_%i$", sel), colnames(grouped))]
      emptyRows <- apply(grouped[, grep("^Int_", colnames(grouped))], 1, function(x) all(x == 0))
      grouped <<- grouped[!emptyRows, ]
      updateReactive$v <- updateReactive$v + 1  # update aligR
      if (exists("clusteringData"))
        rm(clusteringData, pos = globalenv())
      FillAlignmentTable()
    }

  })
  
  
  # TIC ####
  
  observeEvent(input$TICplot_brush, {
    TICplot_ranges$x <- c(input$TICplot_brush$xmin,input$TICplot_brush$xmax)
    TICplot_ranges$y <- c(input$TICplot_brush$ymin,input$TICplot_brush$ymax)
  })
  
  observeEvent(input$TICplot_click, {
    TICplot_RT$x <- input$TICplot_click$x
    position <- which(abs(datenList[[selected]]@scantime-TICplot_RT$x*60)==min(abs(datenList[[selected]]@scantime-TICplot_RT$x*60)))
    spektrum <<- xcms::getScan(datenList[[selected]],position)
    MSplot_plotten()
  })
  
  observeEvent(input$TICplot_dblclick, {
    TICplot_ranges$x <- NULL
    TICplot_ranges$y <- NULL
  })
  
  observeEvent(input$MSplot_click, {
    massesAboveInt <<- c(0,spektrum[spektrum[,2]>input$MSplot_click$y,1])
    masse$y <- input$MSplot_click$y
    masse$mz <- massesAboveInt[which.min(abs(massesAboveInt-input$MSplot_click$x))]
  })
  
  observeEvent(input$MSplot_dblclick, {
    MSplot_ranges$x <- c(min(spektrum[,1]),max(spektrum[,1]))
    MSplot_ranges$y <- NULL
  })
  
  observeEvent(input$MSplot_brush, {
    MSplot_ranges$x <- c(input$MSplot_brush$xmin,input$MSplot_brush$xmax)
    MSplot_ranges$y <- c(input$MSplot_brush$ymin,input$MSplot_brush$ymax)
  })
  
  
  
  #### XIC ####
  
  observeEvent(input$XICplot_brush, {
    XICplot_ranges$x <- c(input$XICplot_brush$xmin,input$XICplot_brush$xmax)
    XICplot_ranges$y <- c(input$XICplot_brush$ymin,input$XICplot_brush$ymax)
  })
  
  observeEvent(input$XICplot_dblclick, {
    XICplot_ranges$x <- NULL
    XICplot_ranges$y <- NULL
  })
  
  
  
  #### Peak Picking ####
  # when Button "Pick Peaks" is pressed
  observeEvent(input$PickPeaks, { 
   
    overallp <- shiny::Progress$new(session, 0, 1)
    overallp$set(value = 0, message = "Overall progress", detail = "Setting up")
    
    mz_min <- isolate(input$PeakPick_massrange[1])
    mz_max <- isolate(input$PeakPick_massrange[2])
    mz_step <- isolate(input$PeakPick_mzstep)
    rt_min <- isolate(input$PeakPick_RTrange[1])*60
    rt_max <- isolate(input$PeakPick_RTrange[2])*60
    PeakWidth_min <- isolate(input$PeakPick_PeakWidthRange[1])
    PeakWidth_max <- isolate(input$PeakPick_PeakWidthRange[2])
    peak_NoiseScans <- isolate(input$PeakPick_NoiseScans)
    sn <- isolate(input$PeakPick_SN)
    int_threshold <- isolate(input$PeakPick_IntTresh)
    maxNumPeaks <- isolate(input$maxPeaks)
    precursormzTol <- isolate(input$PPMS2Fit_ppm)
    componentization_ppm <- isolate(input$componentization_ppm)
    componentization_rt_tol <- isolate(input$componentization_rt_tol)
    componentization_rt_tol_l <- isolate(input$componentization_rt_tol_l)
    componentization_rt_tol_r <- isolate(input$componentization_rt_tol_r)
    componentization_rt_tol_sum <- isolate(input$componentization_rt_tol_sum)
    componentization_dynamic_tolerance <- isolate(input$Componentization_Dynamic_Tolerance)
    
    settings <- list()
    settings$massrange <- c(mz_min,mz_max)
    settings$mz_step <- mz_step
    settings$rtrange <- c(rt_min, rt_max)
    settings$peakwidth <- c(PeakWidth_min, PeakWidth_max)
    settings$NoiseScans <- peak_NoiseScans
    settings$sn <- sn
    settings$int_threshold <- int_threshold
    settings$precursormzTol <- precursormzTol
    settings$ppm <- componentization_ppm
    settings$RT_Tol <- componentization_rt_tol
    settings$PPTableRow <- 1
    settings$orderType <- list(list(0,"asc"))
    settings$TableLength <- 10
    settings$displayStart <- 0
    settings$maxPeaks <- maxNumPeaks
    settings$componentization_dynamic_tolerance <- componentization_dynamic_tolerance
    
    # PP single sample ####
    if (input$PPRadioSingleAll == "one" && is.numeric(selected) && !sampleList[selected, "deleted"]) {
    
      removeSampleAgain <- FALSE
      # if sample is not in RAM we need to load it first - afterwards it will be removed again
      if (!sampleList[selected, "RAM"]) { 
        removeSampleAgain <- TRUE
        overallp$set(detail = "Loading file")
        overallp$inc()
        datenList[[selected]] <<- xcms::xcmsRaw(datenList[[selected]]@filepath@.Data, includeMSn = TRUE)
      }
      
      overallp$set(0.2, detail = "Peak-Picking")
      peaklist[[selected]] <<- ntsworkflow::FindPeaks_BfG(daten = datenList[[selected]], 
                      mz_min, mz_max, mz_step, rt_min, rt_max, sn, 
                      int_threshold, peak_NoiseScans, precursormzTol, 
                      PeakWidth_min,PeakWidth_max,maxNumPeaks
                      )
        
      
      peakPickSettings[[selected]] <<- settings

      overallp$set(0.7, detail = "Componentization")
      if (nrow(peaklist[[selected]]) > 0) {#run
        peaklist[[selected]] <<- componentization_BfG(
          peaklist[[selected]], 
          daten = datenList[[selected]], 
          ppm = componentization_ppm, 
          Grenzwert_RT = componentization_rt_tol, 
          Grenzwert_FWHM_left = componentization_rt_tol_l, 
          Grenzwert_FWHM_right = componentization_rt_tol_r,
          Summe_all = componentization_rt_tol_sum,
          adjust_tolerance = componentization_dynamic_tolerance)
        
        peaklist[[selected]]$RealPeak <<- TRUE
      }
      
      # if sample was not in RAM before, remove it again:
      overallp$set(0.9, detail = "Finishing")
      if (removeSampleAgain == TRUE) {
        datenList[[selected]]@env$intensity <<- 0
        datenList[[selected]]@env$mz <<- 0
        datenList[[selected]]@env$profile <<- 0
        datenList[[selected]]@env$msnIntensity <<- 0
        datenList[[selected]]@env$msnMz <<- 0
      }
      
      adjust_tolerance <- input$Componentization_Dynamic_Tolerance
      log_file <<- create_log_file()
      
      overallp$set(1, detail = "Done")
      
      # PP multiple samples ####
    } else if (nrow(sampleList) > 0 && any(!sampleList[, "deleted"])) {
      
      removeSampleAgain <- logical(nrow(sampleList))
      # all samples get the same settings
      peakPickSettings[sampleList$ID] <<- list(settings)
      
      overallp$set(0.2, detail = "Peak-picking")
      
      process_pp <- function(datenx, isDeleted, inRam) {
        if (isDeleted)
          return(NULL)
        if (!inRam) 
          datenx <- xcms::xcmsRaw(datenx@filepath@.Data, includeMSn = TRUE)
        
        ntsworkflow::FindPeaks_BfG(daten = datenx, 
                                          mz_min, 
                                          mz_max, 
                                          mz_step,
                                          rt_min,
                                          rt_max,
                                          sn,
                                          int_threshold,
                                          peak_NoiseScans,
                                          precursormzTol,
                                          PeakWidth_min,
                                          PeakWidth_max,
                                          maxNumPeaks
        )
      }
      peaklist <<- parallel::mcMap(process_pp, datenList, sampleList$deleted, 
                         sampleList$RAM, mc.preschedule = FALSE,
                         mc.cores = input$numcores)
      
      
      # componentization for all peaklists at the same time (parallel)
      overallp$set(0.7, detail = "Componentization")
      process_componentization <- function(peaklistx, datenx) {
        
        if (is.null(peaklistx) || nrow(peaklistx) == 0)
          return(peaklistx)
        
        newpl <- ntsworkflow::componentization_BfG(
          peaklistx, 
          daten = datenx, 
          ppm = componentization_ppm, 
          Grenzwert_RT = componentization_rt_tol, 
          Grenzwert_FWHM_left = componentization_rt_tol_l, 
          Grenzwert_FWHM_right = componentization_rt_tol_r,
          Summe_all = componentization_rt_tol_sum,
          adjust_tolerance = componentization_dynamic_tolerance)
        
        newpl$RealPeak <- TRUE
        newpl
      }
      peaklist <<- parallel::mcMap(process_componentization, peaklist, 
                                        datenList, mc.preschedule = FALSE,
                                        mc.cores = 1)#input$numcores)
      
      overallp$set(0.9, detail = "Finishing")
      
      adjust_tolerance <- componentization_dynamic_tolerance
      log_file <<- create_log_file()
      
      overallp$set(1, detail = "Done")
      
    } else {
      showNotification("No samples loaded", type = "warning")
    }
    
    overallp$close()
    
    req(is.numeric(selected) && length(peaklist) >= selected)
    FillPeakPickTable()
  })
  
  
  
  
  #### PP table ####
  
  #### componentsTable ####
  
  componentsTable <- reactive({
    req(selectedR())
    req(length(peaklist) >= selectedR())
    req(!is.null(peaklist[[selectedR()]]))
    req(nrow(peaklist[[selectedR()]]) > 0)
    # get only groupLeader rows
    o <- subset(peaklist[[selectedR()]], groupleader, c(mz, RT, Intensity, Area, SN, RealPeak,
                                                   groupleader, Gruppe))
    o$RT <- o$RT / 60
    o
  })
  
  ppRow <- reactiveVal()
  # When one row in the Peak Picking results table is selected
  observeEvent(input$PeakPickingTable_cell_clicked, { 
    zeile <- input$PeakPickingTable_rows_selected 
    if (is.numeric(selected)) {
      if (length(peakPickSettings) >= selected) {
        if (!is.null(peakPickSettings[[selected]]$massrange)) {
          if (!is.null(zeile)) {
            if ((length(zeile)==1) && sampleList[selected, "RAM"]) {
              output$PeakPickXIC <- renderPlot(PeakPickXIC_plotten(zeile, selectedR()))
              PeakPickMSPlotten(zeile)
              output$PeakPickMS2 <- renderPlot(PeakPickMS2Plotten(zeile))
            }
            peakPickSettings[[selected]]$PPTableRow <<- zeile
            ppRow(zeile)
            # select row in component table
            # what is groupleader of selected row?
            # get group
            g <- peaklist[[selectedR()]][zeile, "Gruppe"]
            # get row of component table
            ct <- componentsTable()
            r <- which(ct$Gruppe == g)
            componentsTableProxy %>% 
              selectRows(r) %>% 
              selectPage(r %/% 10 + 1)
           
          }  
        }  
      }
    }  
    
  })
  
  
  output$ComponentsTable <- renderDT(
     
    datatable(componentsTable(),
              colnames = c("m/z","RT","Intensity","Area","S/N","RealPeak","groupleader", "CompID"),
              rownames = FALSE,
              selection = 'single',
              filter = list(position = "bottom", plain = TRUE),
              width = "100%",
              options = list(
                dom = 'ltip',
                autowidth = TRUE,
                stateSave = TRUE,
                stateDuration = -1,
                pageLength = 10,
                columnDefs = list(list(visible = FALSE, targets = 5:7))
              )
    ) %>% formatRound('mz', 4) %>% formatRound('RT', 2) %>% 
      formatRound(c('Intensity', 'SN', 'Area'), 0) 
  ) 
  componentsTableProxy <- dataTableProxy("ComponentsTable")
  
  # When one row in the component table results table is selected
  observe({ 
    # which row of PP table selected?
    input$ComponentsTable_cell_clicked
    zeileComp <- input$ComponentsTable_rows_selected
    g <- componentsTable()[zeileComp, "Gruppe"]
    zeile <- which(peaklist[[selectedR()]]$Gruppe == g & peaklist[[selectedR()]]$groupleader)
   
    req(is.numeric(zeile) && length(zeile) == 1)
    if (is.numeric(selectedR())) {
      if (length(peakPickSettings) >= selectedR()) {
        if (!is.null(peakPickSettings[[selectedR()]]$massrange)) {
          if (!is.null(zeile)) {
            if (sampleList[selectedR(), "RAM"]) {
              output$ComponentsXIC <- renderPlot(ComponentXIC_plotten(zeile))
            }
            peakPickSettings[[selectedR()]]$PPTableRow <<- zeile
            ComponentXIC_ranges$y <<- c(0,max(peaklist[[selectedR()]]$Intensity[peaklist[[selectedR()]]$Gruppe == peaklist[[selectedR()]]$Gruppe[zeile]]))
            selectRows(dataTableProxy("PeakPickTable"), zeile)
            KomponentenZeilen <- which(peaklist[[selectedR()]]$Gruppe == peaklist[[selectedR()]]$Gruppe[zeile])
            KomponentenZeilen <- KomponentenZeilen[order(peaklist[[selectedR()]]$Intensity[KomponentenZeilen])]
            farben <- rainbow(length(KomponentenZeilen))
            ComponentType <- NULL
            for (i in 1:length(KomponentenZeilen)) {
              ComponentType[i] <- ""
              
              if (peaklist[[selectedR()]]$C13[KomponentenZeilen[i]] > 0)  
                ComponentType[i] <- paste("13C of",as.character(round(peaklist[[selectedR()]]$mz[peaklist[[selectedR()]]$C13[KomponentenZeilen[i]]],4)))
              
              if (peaklist[[selectedR()]]$NaAddukt[KomponentenZeilen[i]] > 0)  
                ComponentType[i] <- paste("[Na] of",as.character(round(peaklist[[selectedR()]]$mz[peaklist[[selectedR()]]$NaAddukt[KomponentenZeilen[i]]],4)))
              
              if (peaklist[[selectedR()]]$NH4Addukt[KomponentenZeilen[i]] > 0)  
                ComponentType[i] <- paste("[NH4] of",as.character(round(peaklist[[selectedR()]]$mz[peaklist[[selectedR()]]$NH4Addukt[KomponentenZeilen[i]]],4)))
              
              if (peaklist[[selectedR()]]$KAddukt[KomponentenZeilen[i]] > 0)  
                ComponentType[i] <- paste("[K] of",as.character(round(peaklist[[selectedR()]]$mz[peaklist[[selectedR()]]$KAddukt[KomponentenZeilen[i]]],4)))
              
              if (peaklist[[selectedR()]]$Cl1[KomponentenZeilen[i]] > 0)  
                ComponentType[i] <- paste("37Cl of",as.character(round(peaklist[[selectedR()]]$mz[peaklist[[selectedR()]]$Cl1[KomponentenZeilen[i]]],4)))
              
              if (peaklist[[selectedR()]]$Cl2[KomponentenZeilen[i]] > 0)  
                ComponentType[i] <- paste("37Cl of",as.character(round(peaklist[[selectedR()]]$mz[peaklist[[selectedR()]]$Cl2[KomponentenZeilen[i]]],4)))
              
              if (peaklist[[selectedR()]]$Cl3[KomponentenZeilen[i]] > 0)  
                ComponentType[i] <- paste("37Cl of",as.character(round(peaklist[[selectedR()]]$mz[peaklist[[selectedR()]]$Cl3[KomponentenZeilen[i]]],4)))
              
              if (peaklist[[selectedR()]]$Cl4[KomponentenZeilen[i]] > 0)  
                ComponentType[i] <- paste("37Cl of",as.character(round(peaklist[[selectedR()]]$mz[peaklist[[selectedR()]]$Cl4[KomponentenZeilen[i]]],4)))
            }
            outputtable <- as.data.frame(cbind(
              farben,
              "",
              round(peaklist[[selectedR()]]$mz[KomponentenZeilen],4),
              round(peaklist[[selectedR()]]$RT[KomponentenZeilen]/60,1),
              round(round(peaklist[[selectedR()]]$Intensity[KomponentenZeilen])),
              ComponentType
            ))
            colnames(outputtable) <- c("Farben","Anzeige","mz","RT","Intensity","Component")
            outputtable$mz <- as.numeric(as.character(outputtable$mz))
            outputtable$RT <- as.numeric(as.character(outputtable$RT))
            outputtable$Intensity <- as.numeric(as.character(outputtable$Intensity))
            output$ComponentsTable2 <- DT::renderDataTable(
              DT::datatable(outputtable,
                            rownames = FALSE,
                            selection = list(mode = 'single'),
                            colnames = c("Farben","","mz","RT","Intensity","Component"),
                            options = list(
                              dom = 'tip',
                              autowidth = TRUE,
                              columnDefs = list(list(visible=FALSE, targets = 0),
                                                list(type = "character", width = "5%",className="dt-right", targets = 1),
                                                list(type = "num", width = "20%",className="dt-right", targets = 2),
                                                list(type = "num", width = "15%",className="dt-right", targets = 3),
                                                list(type = "num", width = "20%",className="dt-right", targets = 4),
                                                list(type = "character", width = "40%",className="dt-right", targets = 5)
                              ),
                              order = list(list(2,"asc"))
                            )
                            
              )
              %>% formatRound(columns = "mz",digits = 4)
              %>% formatRound(columns = "RT",digits = 1)
              %>% formatRound(columns = "Intensity",digits = 0)
              %>% formatStyle(columns = "Anzeige", valueColumns = "Farben", target = "cell", backgroundColor = styleEqual(farben,farben))
            )
          }  
        }  
      }
    }  
    
  })
  
  
  observeEvent(input$lastkeypresscode, {
     
    if (length(selected) == 1 &&selected!= 0) {
      if (length(peakPickSettings) >= selected) {
        zeilen <- isolate(input$PeakPickingTable_rows_all)
        angezeigteZeilen <- isolate(input$PeakPickingTable_rows_current)
        TableState <- isolate(input$PeakPickingTable_state)
        aktuelleZeile <- which(zeilen == peakPickSettings[[selected]]$PPTableRow)
        if (length(zeilen) > 0) {
          if (input$lastkeypresscode[1] == 40) { #arrow down
            if (length(aktuelleZeile) == 0) {
              peakPickSettings[[selected]]$PPTableRow <<- zeilen[1]
            } else {
              if (aktuelleZeile < length(zeilen)) {
                peakPickSettings[[selected]]$PPTableRow <<- zeilen[aktuelleZeile+1]
                if ((aktuelleZeile %% TableState$length) == 0) selectPage(dataTableProxy("PeakPickingTable"), aktuelleZeile/TableState$length+1)
              }
            }
            selectRows(dataTableProxy("PeakPickingTable"), peakPickSettings[[selected]]$PPTableRow)
            if (sampleList[selected, "RAM"]) {
              output$PeakPickXIC <- renderPlot(PeakPickXIC_plotten(peakPickSettings[[selected]]$PPTableRow, selectedR()))
              PeakPickMSPlotten(peakPickSettings[[selected]]$PPTableRow)
              output$PeakPickMS2 <- renderPlot(PeakPickMS2Plotten(peakPickSettings[[selected]]$PPTableRow))
            }
          }
          if (input$lastkeypresscode[1] == 38) { #arrow up
            if (length(aktuelleZeile) == 0) {
              peakPickSettings[[selected]]$PPTableRow <<- zeilen[1]
            } else {
              if (aktuelleZeile > 1) {
                peakPickSettings[[selected]]$PPTableRow <<- zeilen[aktuelleZeile-1]
                if ((aktuelleZeile %% TableState$length) == 1) selectPage(dataTableProxy("PeakPickingTable"), (aktuelleZeile-1)/TableState$length)
              }
            }
            selectRows(dataTableProxy("PeakPickingTable"), peakPickSettings[[selected]]$PPTableRow)
            if (sampleList[selected, "RAM"]) {
              output$PeakPickXIC <- renderPlot(PeakPickXIC_plotten(peakPickSettings[[selected]]$PPTableRow, selected))
              PeakPickMSPlotten(peakPickSettings[[selected]]$PPTableRow)
              output$PeakPickMS2 <- renderPlot(PeakPickMS2Plotten(peakPickSettings[[selected]]$PPTableRow))
            }
          }
          if (input$lastkeypresscode[1] == 33) { #page up
            if (aktuelleZeile < TableState$length) {
              peakPickSettings[[selected]]$PPTableRow <<- zeilen[1]
            } else {
              peakPickSettings[[selected]]$PPTableRow <<- zeilen[aktuelleZeile-TableState$length]
              selectPage(dataTableProxy("PeakPickingTable"), (aktuelleZeile-1) %/% TableState$length)
            }
            selectRows(dataTableProxy("PeakPickingTable"), peakPickSettings[[selected]]$PPTableRow)
            if (sampleList[selected, "RAM"]) {
              output$PeakPickXIC <- renderPlot(PeakPickXIC_plotten(peakPickSettings[[selected]]$PPTableRow, selectedR()))
              PeakPickMSPlotten(peakPickSettings[[selected]]$PPTableRow)
              output$PeakPickMS2 <- renderPlot(PeakPickMS2Plotten(peakPickSettings[[selected]]$PPTableRow))
            }
          }
          if (input$lastkeypresscode[1] == 34) { #page down
            if (aktuelleZeile > (length(zeilen)-TableState$length)) {
              peakPickSettings[[selected]]$PPTableRow <<- zeilen[length(zeilen)]
              selectPage(dataTableProxy("PeakPickingTable"), (peakPickSettings[[selected]]$PPTableRow %/% TableState$length) + 1)
            } else {
              peakPickSettings[[selected]]$PPTableRow <<- zeilen[aktuelleZeile+TableState$length]
              selectPage(dataTableProxy("PeakPickingTable"), (aktuelleZeile %/% TableState$length) + 2)
            }
            selectRows(dataTableProxy("PeakPickingTable"), peakPickSettings[[selected]]$PPTableRow)
            if (sampleList[selected, "RAM"]) {
              output$PeakPickXIC <- renderPlot(PeakPickXIC_plotten(peakPickSettings[[selected]]$PPTableRow, selectedR()))
              PeakPickMSPlotten(peakPickSettings[[selected]]$PPTableRow)
              output$PeakPickMS2 <- renderPlot(PeakPickMS2Plotten(peakPickSettings[[selected]]$PPTableRow))
              
            }
          }
          if (input$lastkeypresscode[1] == 46) { #delete
            peaklist[[selected]]$RealPeak[peakPickSettings[[selected]]$PPTableRow] <<- FALSE
            FillPeakPickTable(orderType = TableState$order, displayStart = which(zeilen == angezeigteZeilen[1])-1, TableLength = TableState$length)
            
          }
          if (input$lastkeypresscode[1] == 45) { #insert
            peaklist[[selected]]$RealPeak[peakPickSettings[[selected]]$PPTableRow] <<- TRUE
            FillPeakPickTable(orderType = TableState$order, displayStart = which(zeilen == angezeigteZeilen[1])-1, TableLength = TableState$length)
          }
        }
      }
    }
  })
  
  
  observeEvent(input$PeakPickXIC_brush, {#zooming - change the reactive values according to the brushed area
    PPXIC_ranges$x <- c(input$PeakPickXIC_brush$xmin,input$PeakPickXIC_brush$xmax)
    PPXIC_ranges$y <- c(input$PeakPickXIC_brush$ymin,input$PeakPickXIC_brush$ymax)
  })
  
  observeEvent(input$PeakPickXIC_dblclick, {#reset zooming
    PPXIC_ranges$x <- NULL
    PPXIC_ranges$y <- NULL
  })
  
  observeEvent(input$ComponentsXIC_brush, {#zooming - change the reactive values according to the brushed area
    ComponentXIC_ranges$x <- c(input$ComponentsXIC_brush$xmin,input$ComponentsXIC_brush$xmax)
    ComponentXIC_ranges$y <- c(input$ComponentsXIC_brush$ymin,input$ComponentsXIC_brush$ymax)
  })
  
  observeEvent(input$ComponentsXIC_dblclick, {#reset zooming
    ComponentXIC_ranges$x <- NULL
    ComponentXIC_ranges$y <<- c(0,max(peaklist[[selected]]$Intensity[peaklist[[selected]]$Gruppe == peaklist[[selected]]$Gruppe[peakPickSettings[[selected]]$PPTableRow]]))
  })
  
  observeEvent(input$PeakPickMS2_click, {#copy the recent MS2 spectrum into a table and print it in a separate window (bsModal)
    output$PeakPickMS2Table <- DT::renderDataTable(
      DT::datatable(ms2spektrumR(),
                    colnames = c("Fragment m/z","Intensity"),
                    rownames = FALSE,
                    selection = 'single'
      )
      %>% formatRound(columns = "mz",digits = 4)
      %>% formatSignif(columns = "intensity",digits = 3)
    )
    toggleModal(session,"PeakPickMS2popup",toggle="toggle")
  })
  
  observeEvent(input$PeakPick_GetChlorinated,{
    if (!is.null(selected)) {
      if (length(peaklist)  >= selected) {
        if (!is.null(peaklist[[selected]])) {
          if (nrow(peaklist[[selected]]) > 0) { 
            IsotopeType <- matrix(0, ncol = 1, nrow=nrow(peaklist[[selected]]))
            
            anzahlCl1Isotope <- table(peaklist[[selected]]$Cl1)
            IsotopeType[as.numeric(names(which(anzahlCl1Isotope == 1)))] <- 1
            IsotopeType[as.numeric(names(which(anzahlCl1Isotope == 2)))] <- 2
            
            anzahlCl2Isotope <- table(peaklist[[selected]]$Cl2)
            IsotopeType[as.numeric(names(which((anzahlCl2Isotope == 1) | (anzahlCl2Isotope == 2))))] <- 1
            IsotopeType[as.numeric(names(which(anzahlCl2Isotope == 3)))] <- 2
            
            anzahlCl3Isotope <- table(peaklist[[selected]]$Cl3)
            IsotopeType[as.numeric(names(which((anzahlCl3Isotope == 1) | (anzahlCl3Isotope == 2))))] <- 1
            IsotopeType[as.numeric(names(which(anzahlCl3Isotope == 3)))] <- 2
            
            anzahlCl4Isotope <- table(peaklist[[selected]]$Cl4)
            IsotopeType[as.numeric(names(which((anzahlCl4Isotope == 1) | (anzahlCl4Isotope == 2) | (anzahlCl4Isotope == 3))))] <- 1
            IsotopeType[as.numeric(names(which(anzahlCl4Isotope == 4)))] <- 2
            
            FillPeakPickTable(
              filterMaske = list(
                NULL,
                NULL,
                NULL,
                NULL,
                NULL,
                list(search = "1 ... 1"),
                NULL,
                NULL,
                NULL
              ),
              HighlightType = IsotopeType,
              orderType = peakPickSettings[[selected]]$order,
              TableLength = peakPickSettings[[selected]]$length
            )
          }
        }
      }
    }
  }
  )
  
  observeEvent(input$PeakPick_ShowAll,{
    if (!is.null(selected)) {
      if (length(peaklist)  >= selected) {
        if (!is.null(peaklist[[selected]])) {
          if (nrow(peaklist[[selected]]) > 0) {
            FillPeakPickTable(filterMaske = list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL))
          }
        }
      }
    }   
  }
  )
  
  observeEvent(input$PeakPick_HideFalsePositives,{
    if (!is.null(selected)) {
      if (length(peaklist)  >= selected) {
        if (!is.null(peaklist[[selected]])) {
          if (nrow(peaklist[[selected]]) > 0) {
            FillPeakPickTable(filterMaske = list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,list(search = "1 ... 1")))
          }
        }
      }
    }   
  }
  )
  
  
  observeEvent(input$PeakPick_GetBrominated,{
    if (!is.null(selected)) {
      if (length(peaklist)  >= selected) {
        if (!is.null(peaklist[[selected]])) {
          if (nrow(peaklist[[selected]]) > 0) { 
            IsotopeType <- matrix(0, ncol = 1, nrow=nrow(peaklist[[selected]]))
            
            anzahlBr1Isotope <- table(peaklist[[selected]]$Br1)
            IsotopeType[as.numeric(names(which(anzahlBr1Isotope == 1)))] <- 1
            IsotopeType[as.numeric(names(which(anzahlBr1Isotope == 2)))] <- 2
            
            anzahlBr2Isotope <- table(peaklist[[selected]]$Br2)
            IsotopeType[as.numeric(names(which((anzahlBr2Isotope == 1) | (anzahlBr2Isotope == 2))))] <- 1
            IsotopeType[as.numeric(names(which(anzahlBr2Isotope == 3)))] <- 2
            
            anzahlBr3Isotope <- table(peaklist[[selected]]$Br3)
            IsotopeType[as.numeric(names(which((anzahlBr3Isotope == 1) | (anzahlBr3Isotope == 2) | (anzahlBr3Isotope == 3))))] <- 1
            IsotopeType[as.numeric(names(which(anzahlBr3Isotope == 4)))] <- 2
            
            anzahlBr4Isotope <- table(peaklist[[selected]]$Br4)
            IsotopeType[as.numeric(names(which((anzahlBr4Isotope == 1) | (anzahlBr4Isotope == 2) | (anzahlBr4Isotope == 3) | (anzahlBr4Isotope == 4))))] <- 1
            IsotopeType[as.numeric(names(which(anzahlBr4Isotope == 5)))] <- 2
            
            FillPeakPickTable(filterMaske = list(NULL,NULL,NULL,NULL,NULL,NULL,list(search = "1 ... 1"),NULL,NULL),HighlightType = IsotopeType)
          }
        }
      }
    }
  })
  
  # Randomize peak picking table ####
  observeEvent(input$peakpick_random, {
    req(!is.null(selected))
    req(length(peaklist) >= selected)
    req(!is.null(peaklist[[selected]]))
    req(nrow(peaklist[[selected]]) > 0)
    
    FillPeakPickTable(filterMaske = list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL), 
                      orderType = list(list(9,"asc")), TableLength = 100)
  })
  
  shinyFileSave(input, 'PickPeaksTableExport', roots=home, session=session)
  observeEvent(input$PickPeaksTableExport,{
     
    req(is.list(input$PickPeaksTableExport))
    speicherort <- as.character(parseSavePath(home, input$PickPeaksTableExport)$datapath)
    
    if ((length(speicherort) > 0) && !is.null(selected)) {
      if (length(peaklist)  >= selected) {
        if (nrow(peaklist[[selected]]) > 0) { 
          # saveto <- NULL
          # tryCatch({saveto <- file.choose(new = TRUE)}, error = function(ex) {})
          exportTabelle <- cbind(seq(nrow(peaklist[[selected]])),peaklist[[selected]])
          colnames(exportTabelle)[1] <- "ID"
          # mz round 4 d.p.
          exportTabelle$mz <- round(exportTabelle$mz, 4)
          # RT numbers in min and 2 d.p.
          for (co in c("RT", "LeftendRT", "FWHM_left", "FWHM_right", "RightendRT"))
            exportTabelle[, co] <- round(exportTabelle[, co] / 60, 2)
          # Intensity etc. round 3 d.p.
          for (co in c("Intensity", "XIC_Intensity", "Area", "NoiseDeviation",  "Baseline", "SN"))
            exportTabelle[, co] <- round(exportTabelle[, co], 3)
          if (!is.null(speicherort)) 
            write.csv(exportTabelle, file = speicherort, row.names = FALSE)  
        }
      }
    }
  }
  )
  
  # Event handler Alignment
  alignPeaksEvent <- reactiveVal()  # this is done to allow programatic starting of alignment
  observeEvent(input$AlignPeaks, {
    if (is.null(alignPeaksEvent()))
      alignPeaksEvent(1) else alignPeaksEvent(alignPeaksEvent() + 1)
  })
  
  # GenForm ####
  observeEvent(input$genformGo, {
    
    # compile genform? use 
    # $ g++ main.cpp ms*.cpp -o genform 
    # in the genform source folder
    # The program was then placed in ~/bin
    # R/genformR_kj.R was added to ntsworkflow (code from E. Schymanski)
    
    # create MS1 and MS2 files in current working dir
    prog <- shiny::Progress$new(session, 0, 1)
    prog$set(value = 0, message = "Preparing")
    
    if (is.null(selectedR()) || is.null(ppRow())) {
      showNotification("First select a peak with an MS2scan in the Featurelist (see Data)", type = "error")
      prog$close()
    }
    req(is.numeric(selectedR()) && is.numeric(ppRow()))
    
    searchMz <- peaklist[[selectedR()]][ppRow(), "mz"]
    ms1scan <- peaklist[[selectedR()]][ppRow(), "Scan"]
    ms2scan <- peaklist[[selectedR()]][ppRow(), "MS2scan"]
    
    if (ms2scan == 0) {
      showNotification("You must select a peak with an MS2 scan for GenForm to work", type = "error")
      prog$close()
    }
    req(ms2scan != 0)
    
    # check that file is in ram
    if (!subset(sampleList, ID == selectedR(), RAM, drop = TRUE)) {
      showNotification("The sample must be loaded", type = "error")
      prog$close()
    }
    req(subset(sampleList, ID == selectedR(), RAM, drop = TRUE))
    
    gfms1 <- xcms::getScan(datenList[[selectedR()]], ms1scan)
    gfms2 <- xcms::getMsnScan(datenList[[selectedR()]], ms2scan)
    write.table(gfms1, file = file.path(globalwd, "GenFormMS1.txt"), row.names = F, col.names = F)
    write.table(gfms2, file = file.path(globalwd, "GenFormMS2.txt"), row.names = F, col.names = F)
    
    prog$set(value = 0.5, message = "Running GenForm")
    tryCatch(
      error = function(cnd) { 
        showNotification(paste("GenForm did not accept that. Report:", conditionMessage(cnd)), type = "error")
        prog$close()
        },
      {
        # function copied from package GenFormR by Emma Schymanski, edited for use on CentOS
        # GenForm itself from M. Meringer et al.
        ntsworkflow::RunGenForm(ms_file = file.path(globalwd, "GenFormMS1.txt"),
                                mz = searchMz,
                                msms_file = file.path(globalwd, "GenFormMS2.txt"), 
                                ppm = input$genformMS1tol,
                                acc = input$genformMS2tol,
                                elements = input$genformElements,
                                ff = input$genformFf,
                                ion = input$genformIon, 
                                GenFormDir = GENFORM_LOC, 
                                ResultsDir = globalwd)
      })
    prog$set(value = 0.7, message = "Collecting results")
    
    # read in result
    req(file.exists(file.path(globalwd, "GenFormMS2_GenForm.txt")))
    gf_res <- read.table(file.path(globalwd, "GenFormMS2_GenForm.txt"), header = F,
                         col.names = c("formula", "dbe", "dev_ppm", "MS1_fit", "MS2_fit", "combi_fit"))
    
    # display output
    
    # DT output
    output$genFormResTab <- DT::renderDataTable(
      DT::datatable(gf_res)
    )
    prog$close()
  })
  
  # Alignment ####
  
  # When button "Align Peaks" is pressed
  observeEvent(alignPeaksEvent(), {
    # write aligment table into a global variable:
     
    # first correct peaks lists for parallel 
    
    # to avoid problems with deleted samples, if you click on alignment, sampleList is reduced
    # to only non-deleted samples, and the sample IDs are reassigned. In new versions of alignment,
    # it should be able to deal with deleted samples.
    
    progress <- shiny::Progress$new()
    progress$set(value = .5, message = "Preparing...")
    #browser()
    # reset sample list
    sampleList <<- sampleList[!sampleList$deleted, ]
    sampleList$ID <<- seq_len(nrow(sampleList))
    sampleListR(sampleList)  
    # remove NULL values from datenlist, headerlist, peaklist, peaklist settings
    datenList <<- compact(datenList)
    
    peaklist <<- compact(peaklist)
    peakPickSettings <<- compact(peakPickSettings)
    
    if (length(peaklist) != nrow(sampleList))
      showNotification("Not all samples were picked", type = "error")
    req(length(peaklist) == nrow(sampleList))
    
    # new_peaklist will add peak IDs to all peaks
    peaklist <<- ntsworkflow::new_peaklist(peaklist, datenList)
    
    # reorder all peaklists by mass
    for (i in 1:length(peaklist)) {
      peaklist[[i]] <<- peaklist[[i]][order(peaklist[[i]]$mz),]
    }
    
    # reset annotation
    if (!is.null(annotationTable)) {
      annotationTable <<- NULL
    }
    
    progress$set(value = 1, message = "Aligning features...")
    
    # logging
    ppm_dev <- input$Alignment_deltamz
    DeltaRT <- isolate(input$Alignment_deltaRT)
    log_file <<- create_log_file()
    
    # if the user only wants to align the group leaders, then first reduce the peaklists
    if (input$alignLeaders) {
      peaklist2 <- lapply(peaklist, function(x) x[x$groupleader, ])
    } else {
      peaklist2 <- peaklist
    }

    # ready for alignment
    # For more than 10 samples and not just group leaders try parallel alignment, 
    # for anything less it is not worth it!
    # Be aware this does produce small differences in the result, some peaks 
    # are grouped differently, reason is unknown
    if (!input$alignLeaders && input$parallelAlign) {
      result <- ntsworkflow::alignment_BfG_cpp_par(
        peaklists = peaklist2, 
        ppm_dev = isolate(input$Alignment_deltamz), 
        DeltaRT = isolate(input$Alignment_deltaRT), 
        mz_dev_unit = isolate(input$aligMzTolType),
        mDa_split = 100,
        numSplits = isolate(input$numcores))
    } else {
      result <- ntsworkflow::alignment_BfG_cpp(peaklist2, isolate(input$Alignment_deltamz), 
                                  isolate(input$Alignment_deltaRT), isolate(input$aligMzTolType))
    }
    # fill out the Gruppe column (summarised componentisation)
    result <- ntsworkflow::summarize_groups(result)
    
    # add id column to alignment table!
    oldnames <- colnames(result)
    result <- cbind(result, seq_len(nrow(result)))
    colnames(result) <- c(oldnames, "alignmentID")
    
    # write result to the global environment it is called "grouped" for historical reasons
    grouped <<- result 
    
    # clean up
    progress$close()
    updateReactive$v <- updateReactive$v + 1  # update aligR
    if (exists("clusteringData")) 
      rm(clusteringData, pos = globalenv())  # remove clustering data
    if (exists("groupedWithAnnotation"))
      rm(groupedWithAnnotation, pos = globalenv())  # remove annotation data
    annotationTable <<- NULL
    FillAlignmentTable()
  })
  
  # Blank correction ####
  observeEvent(input$BlankCorrection, {#when button "Blank Correction" is pressed
    if (is.null(grouped)) {
      showNotification("first align the peak lists")
      return(NULL)
    }
    # There must be at least one blank
    blanksPresent <- any(sampleListR()$sampleType == "Blank")
    
    if (!blanksPresent)
      showNotification("No blank samples", type = "error")
    req(blanksPresent)
    if (!is.null(annotationTable)) {
      annotationTable <<- NULL
      clearAnnotationTable()
    }
    progress <- shiny::Progress$new()
    progress$set(message = "Removing blanks...")
    
    nrow_grouped_before_blankcorrection <- nrow(grouped)
    grouped <<- ntsworkflow::blankCorrection(grouped, sampleListR(), intensityFactor = input$blankCorrFac,
                                deleteGrouped = input$deleteGrouped)
    # logging
    intensityFactor <- input$blankCorrFac
    deleteGrouped <- input$deleteGrouped
    log_file <<- create_log_file()
    
    progress$close()
    
    # Need to invalidate the aligR reactive so that it using only the blank corrected table
    updateReactive$v <- updateReactive$v + 1
    if (exists("clusteringData"))
      rm(clusteringData, pos = globalenv())
    FillAlignmentTable()
  })
  
  # Normalize ####
  observeEvent(input$Normalize, {
     
    if (is.null(input$AlignmentTable_rows_selected)) {
      showNotification("First select a row from the table")
      return(NULL)
    }
    # record peak IDs used to do the normalization
    pidColsNorm <- grep("^PeakID", colnames(grouped))
    pidNorm <- as.numeric(grouped[input$AlignmentTable_rows_selected, pidColsNorm, drop = T])
    stopifnot(length(pidNorm) == nrow(subset(sampleList, !deleted)))
    stopifnot(all(pidNorm != 0), all(!is.na(pidNorm)))
    sampleList[!sampleList$deleted, "normalizePeakID"] <<- pidNorm
    sampleListR(sampleList)
    
    intCols_normalize <- grep("^Int", colnames(grouped))
    intensities <- as.numeric(grouped[input$AlignmentTable_rows_selected, intCols_normalize])
    if (any(intensities == 0)) {
      showNotification("Samples with intensity 0 were not changed")
      intCols_normalize <- intCols_normalize[intensities != 0]
      intensities <- intensities[intensities != 0]
    }
    
    grouped[, intCols_normalize] <<- sweep(grouped[, intCols_normalize, drop = F], 2, intensities, "/")
    
    # Logging 
    mean_mz_chosen_feature <- as.numeric(grouped[input$AlignmentTable_rows_selected, grep("^mean_mz", colnames(grouped))])
    mean_RT_chosen_feature <- as.numeric(grouped[input$AlignmentTable_rows_selected, grep("^mean_RT", colnames(grouped))])
    row_number <- input$AlignmentTable_rows_selected
    log_file <<- create_log_file()
    
    FillAlignmentTable()
  })
  
  # Plot alignment results #### 
  observeEvent(input$AlignmentTable_cell_clicked, {
    AlignmentTableZeile <<- input$AlignmentTable_rows_selected
    alignmentIdSelected(grouped[AlignmentTableZeile, "alignmentID"])
    if (length(AlignmentTableZeile)==1) {
      output$AlignmentXICs <- renderPlot(alignmentXICsPlotten(alignmentIdSelected()))
      output$AlignmentTrend <- renderPlot(alignmentTrendPlotten(AlignmentTableZeile, tabelle = grouped))
      output$AlignmentMS2 <- renderPlot(alignmentMS2Plotten(alignmentIdSelected()))
    }
  })
  

  
  shinyFileSave(input, 'groupedExport', roots=home, session=session)
  observeEvent(input$groupedExport, {
    req(is.list(input$groupedExport))  # so that shinyFiles does not crash

    saveto <- as.character(parseSavePath(home, input$groupedExport)$datapath)
    if (!is.null(grouped) && !is.null(saveto)) {
      # RT seconds to minutes and round 2 d.p.
      export <- grouped
      rtCols <- grep("RT", colnames(export))
      mzCols <- grep("mz", colnames(export))
      intCols <- grep("Int", colnames(export))
      for (col in rtCols)
        export[, col] <- round(export[, col] / 60, 2)
      # round mz to 4 d.p.
      for (col in mzCols)
        export[, col] <- round(export[, col], 4)
      # round int to 3 d.p.
      for (col in intCols)
        export[, col] <- round(export[, col], 3)
      
      export <- tibble::as_tibble(export)
      # add samples names
      
      # get sample names from table
      intColNames <- colnames(grouped)[grep("^Int_", colnames(grouped))]
      sampsInTable <- as.numeric(stringr::str_match(intColNames, "_(\\d+)$")[,2])
      
      for (sId in sampsInTable) {
        # get column number of first column for sample
         
        firstCol <- grep(paste0("^PeakID_", sId, "$"), colnames(export))
        # insert sample name
        export <- tibble::add_column(export, sampleName = basename(sampleList[sId, "File"]), .before = firstCol)
        colnames(export)[firstCol] <- paste0("sampleName_", sId)
      }
      write.csv(export, file = saveto, row.names = FALSE)
    }  
  })
  

  # Alignment table filters ####
  
  observeEvent(input$removeRare, {
    p <- Progress$new()
     
    p$set(message = "Filtering...")
    # filter table
    tryCatch(
      error = function(cnd) 
        showNotification(paste("Bad input:", conditionMessage(cnd)), type = "error"),
      {
        nrow_grouped_before_filter <- nrow(grouped)
         
        filteredAlig <- removeRare(grouped, input$minDetections)
        # first fill alig, since if there is a crash, grouped will not be overwritten
        FillAlignmentTable(filteredAlig)
        # update global and reactive variables
        grouped <<- filteredAlig
        updateReactive$v <- updateReactive$v + 1  ## invalidate AligR reactive
        if (exists("clusteringData"))
          rm(clusteringData, pos = globalenv())
        
        # Logging
        minDetections <- input$minDetections
        log_file <<- create_log_file()
      }
    )
    p$close()
  })
  
  observeEvent(input$onlyGroupLeaders, {
    p <- Progress$new()
    p$set(message = "Filtering...")
    tryCatch(
      error = function(cnd) 
        showNotification(paste("Bad input:", conditionMessage(cnd)), type = "error"),
      {
        stopifnot(eval(parse(text = input$leaderSamples)) %in% sampleListR()$ID)
        nrow_grouped_before_filter <- nrow(grouped)
         
        filteredAlig <- onlyGroupLeaders(grouped, peaklist, eval(parse(text = input$leaderSamples)))
        FillAlignmentTable(filteredAlig)
        grouped <<- filteredAlig
        updateReactive$v <- updateReactive$v + 1  ## invalidate AligR reactive
        if (exists("clusteringData"))
          rm(clusteringData, pos = globalenv())
        
        # Logging
        leader_first_sample <- min(eval(parse(text = input$leaderSamples)))
        leader_last_sample <- max(eval(parse(text = input$leaderSamples)))
        log_file <<- create_log_file()
      }
    )
    p$close()
  })
  
  observeEvent(input$keepReps, {
    p <- Progress$new()
    p$set(message = "Filtering...")
    tryCatch(
      error = function(cnd) 
        showNotification(paste("Bad input:", conditionMessage(cnd)), type = "error"),
      {
        stopifnot(eval(parse(text = input$repSamples)) %in% sampleListR()$ID)
        nrow_grouped_before_filter <- nrow(grouped)
        filteredAlig <- keepReps(grouped, eval(parse(text = input$repSamples)), input$repNum,
                             input$repLeast)
        FillAlignmentTable(filteredAlig)
        grouped <<- filteredAlig
        updateReactive$v <- updateReactive$v + 1  ## invalidate AligR reactive
        if (exists("clusteringData"))
          rm(clusteringData, pos = globalenv())
        
        # Logging
        rep_first_sample <- min(eval(parse(text = input$repSamples)))
        rep_last_sample <- max(eval(parse(text = input$repSamples)))
        rep_number_replicates <- input$repNum
        in_at_least_samples <- input$repLeast
        log_file <<- create_log_file()
      }
    )
    p$close()
  })
  
  observeEvent(input$aveReps, {
    p <- Progress$new()
    p$set(message = "Filtering...")
    tryCatch(
      error = function(cnd) 
        showNotification(paste("Bad input:", conditionMessage(cnd)), type = "error"),
      {
        stopifnot(eval(parse(text = input$repAveSamples)) %in% sampleListR()$ID)
        filteredAlig <- averageReplicates(grouped, eval(parse(text = input$repAveSamples)), input$repAveNum)
        FillAlignmentTable(filteredAlig)
        grouped <<- filteredAlig
        updateReactive$v <- updateReactive$v + 1  ## invalidate AligR reactive
        if (exists("clusteringData"))
          rm(clusteringData, pos = globalenv())
        
        # Logging
        repAve_first_sample <- min(eval(parse(text = input$repAveSamples)))
        repAve_last_sample <- max(eval(parse(text = input$repAveSamples)))
        repAve_number_replicates <- input$repAveNum
        log_file <<- create_log_file()
       
      }
    )
    p$close()
  })
  
  observeEvent(input$consecGo, {
    p <- Progress$new()
    p$set(message = "Filtering...")
    tryCatch(
      error = function(cnd)
        showNotification(paste("Bad input:", conditionMessage(cnd)), type = "error"),
      {
        stopifnot(eval(parse(text = input$consecSamples)) %in% sampleListR()$ID)
        nrow_grouped_before_filter <- nrow(grouped)
        filteredAlig <- keepConsecutive(grouped, eval(parse(text = input$consecSamples)), input$consecNum)
        FillAlignmentTable(filteredAlig)
        grouped <<- filteredAlig
        updateReactive$v <- updateReactive$v + 1  ## invalidate AligR reactive
        if (exists("clusteringData"))
          rm(clusteringData, pos = globalenv())
        
        #create log_file
        consec_first_sample <- min(eval(parse(text = input$consecSamples)))
        consec_last_sample <- max(eval(parse(text = input$consecSamples)))
        number_consecutives <- input$consecNum
        log_file <<- create_log_file()
      }
    )
    p$close()
  })
  
  #### Alignment componentisation 2 ####
  
  observeEvent(input$componen2go, {
    p <- Progress$new()
    p$set(message = "2nd stage componentisation...")
    if (nrow(subset(sampleListR(), !deleted)) < 5) {
      showNotification("Session must have more than five data files", type = "error")
      p$close()
    }
    req(nrow(subset(sampleListR(), !deleted)) >= 5)
    aligTemp <- ntsworkflow::alig_componentisation(
      altable = grouped,
      rttols = isolate(input$componen2rttol), 
      fracComponMatch = isolate(input$componen2fracShape), 
      mztol = isolate(input$componen2mztol) / 1000, 
      pearsonCorr = isolate(input$componen2corr), 
      pol = isolate(input$componen2pol)
    )
    p$set(message = "Complete, plotting results...")
    output$componen2text <- renderText({
      sprintf(
      "DBSCAN clustering for %i aligned features.
       The clustering contains %i cluster(s) and %i noise points.",
       nrow(aligTemp), 
       length(table(aligTemp[, "Gruppe"]))-1, 
       table(aligTemp[, "Gruppe"])["0"]
      )
    })
    output$componen2table <- DT::renderDataTable({
      DT::datatable(
        as.data.frame(table(aligTemp[, "Gruppe"])),
        colnames = c("Group ID", "Size"),
        rownames = FALSE,
      )
    })
    output$componen2hist <- renderPlot({
      clusterSize <- unname(table(aligTemp[, "Gruppe"])[-1])
      hist(clusterSize, breaks = seq(2, max(clusterSize)+1, 2), 
           main = "Histogram of cluster sizes", xlab = "Number of features per cluster")
    })
    grouped <<- aligTemp
    FillAlignmentTable()
    p$close()
  })
  
  #### Alignment Overview ####
  
  # build alignment table for overview, highest int and clustering tabs 
  # including calculation for rel intensity
  # range and max intensity and prepare calculation for clustering
  updateReactive <- reactiveValues(v=0)
  aligR <- reactive({
     
    updateReactive$v  # for invalidating aligR
    p <- Progress$new(session, 1, 2)
    p$set(m = "Processing")
    alig <- as.data.frame(grouped)
    alig$row <- seq_len(nrow(alig))
    intCols <- grep("^Int_", colnames(alig))
    alig[, intCols] <- lapply(alig[, intCols, drop = F], signif, digits = 4)
    alig$maxInt <- apply(alig[, intCols, drop = F], 1, max)
    p$inc(1)
    # relative intensity range as minimum use lower quartile 
    # (since there will always be some zeros!)
    # the choice of lower quartile was arbitrary and may need to be changed
    alig$intRange <- apply(alig[, intCols, drop = F], 1, function(x) {
      1 - quantile(x / max(x), .25, na.rm = TRUE, names = FALSE)
    })
    alig$intRange <- round(alig$intRange, 3)
    alig$mean_RT <- alig$mean_RT / 60
    if (!is.null(annotationTable)) {
      annot <- annotationTable[, c("alignmentID", "name")]
      # ignore duplicates
      annot <- annot[!duplicated(annot$alignmentID), ]
      alig <- merge(alig, annot, by = "alignmentID", all.x = TRUE)
      alig <- alig[order(alig$row), ]
    }
    # select only necessary columns
    ic <- grep("^Int_", colnames(alig))
    nc <- grep("name", colnames(alig))
    rc <- grep("row", colnames(alig))
    mc <- grep("maxInt", colnames(alig))
    Rc <- grep("intRange", colnames(alig))
    idc <- grep("alignmentID", colnames(alig))
    mzc <- grep("mean_mz", colnames(alig))
    rtc <- grep("mean_RT", colnames(alig))
    cols <- if (length(nc) == 1) 
      c(rc, mzc, rtc, ic, mc, Rc, nc, idc) else c(rc, mzc, rtc, ic, mc, Rc, idc)
    alig <- alig[, cols]
     
    p$close()
    alig
  })
  overviewRanges <- reactiveValues(x = NULL, y = NULL)
  output$overview <- renderPlot({
    ggplot(aligR(), aes(mean_RT, mean_mz, size = maxInt)) +  geom_point(alpha = .5) +
      ylab("Mean m/z (Da)") + xlab("Mean retention time (min)") + theme_bw() +
      coord_cartesian(xlim = overviewRanges$x, ylim = overviewRanges$y)
  })
  observeEvent(input$overview_dblclick, {
    brush <- input$overview_brush
    if (!is.null(brush)) {
      overviewRanges$x <- c(brush$xmin, brush$xmax)
      overviewRanges$y <- c(brush$ymin, brush$ymax)
    } else {
      overviewRanges$x <- NULL
      overviewRanges$y <- NULL
    }
  })
  
  # keep a store of the mouse position to draw chart
  overviewRowStore <<- NULL
  overviewRow <- reactive({
    r <- which(nearPoints(aligR(), input$overviewHover, allRows = TRUE, maxpoints = 1)$selected_)
    if (is.integer(r) && length(r) == 1) {
      overviewRowStore <<- r
    } else {
      r <- overviewRowStore
    }
    r
  })
  output$intTrendOverview <- renderPlot({
    # get row which is selected
    selRow <- overviewRow()
    req(selRow)
    tl <- sprintf("m/z = %.4f Da, RT = %.2f min", aligR()$mean_mz[selRow], aligR()$mean_RT[selRow])
    if ("name" %in% colnames(aligR()) && !is.na(aligR()$name[selRow]))
      tl <- paste0(tl, ", ", aligR()$name[selRow])
    alignmentTrendPlotten(selRow, titel = tl, tabelle = aligR())
  })
  
  # for similar trends store rows of table
  overviewSelectedIds <- reactiveVal()
  observeEvent(input$overview_click, {
    req(input$overview_brush)
    brush <- input$overview_brush
    alig <- aligR()
    alig <- alig[alig$mean_RT >= brush$xmin & alig$mean_RT <= brush$xmax, ]
    alig <- alig[alig$mean_mz >= brush$ymin & alig$mean_mz <= brush$ymax, ]
    overviewSelectedIds(alig$alignmentID)
    # columns not needed
    alig$intRange <- NULL
    alig$maxInt <- NULL
    alig$row <- NULL
    alig$alignmentID <- NULL
    alig$MS2Fit <- NULL
    coln <- colnames(alig)
     
    
    output$overviewSelected <- DT::renderDataTable(
      DT::datatable(alig,
                    colnames = coln,
                    rownames = FALSE,
                    selection = 'single',
                    width = "100%",
                    options = list(scrollX = TRUE)
      )
      %>% formatRound(columns = "mean_mz",digits = 4)
      %>% formatRound(columns = "mean_RT",digits = 1)
    )
  })
  

  #### Alignment highest intensities ####
  intRangePlRanges <- reactiveValues(x = NULL, y = NULL)
  output$intRangePl <- renderPlot({
    ggplot(aligR(), aes(intRange, maxInt)) + geom_point(alpha = .5) + scale_y_log10() +
      ylab("Max. intensity (counts)") + xlab("Rel. intensity range (-)") + theme_bw() +
      coord_cartesian(xlim = intRangePlRanges$x, ylim = intRangePlRanges$y)
  })
  observeEvent(input$intRangePl_dblclick, {
    brush <- input$intRangePl_brush
    if (!is.null(brush)) {
      intRangePlRanges$x <- c(brush$xmin, brush$xmax)
      intRangePlRanges$y <- c(brush$ymin, brush$ymax)
    } else {
      intRangePlRanges$x <- NULL
      intRangePlRanges$y <- NULL
    }
  })
  
  highIntRowStore <<- NULL
  highIntRow <- reactive({
    r <- which(nearPoints(aligR(), input$highIntHover, allRows = TRUE, maxpoints = 1)$selected_)
    if (is.integer(r) && length(r) == 1) {
      highIntRowStore <<- r
    } else {
      r <- highIntRowStore
    }
    r
  })
  output$intTrendHighInt <- renderPlot({
    selRow <- highIntRow()
    req(selRow)
    tl <- sprintf("m/z = %.4f Da, RT = %.2f min", aligR()$mean_mz[selRow], aligR()$mean_RT[selRow])
    if ("name" %in% colnames(aligR()) && !is.na(aligR()$name[selRow]))
      tl <- paste0(tl, ", ", aligR()$name[selRow])
    alignmentTrendPlotten(selRow, titel = tl, tabelle = aligR())
  })
  
  # for similar trends store rows of table
  intRangeSelectedIds <- reactiveVal()
  observeEvent(input$intRangePl_click, {
    req(input$intRangePl_brush)
    brush <- input$intRangePl_brush
    alig <- aligR()
    alig <- alig[alig$intRange >= brush$xmin & alig$intRange <= brush$xmax, ]
    alig <- alig[alig$maxInt >= brush$ymin & alig$maxInt <= brush$ymax, ]
    
    intRangeSelectedIds(alig$alignmentID)
    # remove rows
    alig$row <- NULL
    alig$MS2Fit <- NULL
    alig$alignmentID <- NULL
    coln <- colnames(alig)
    
    output$intRangeSelected <- DT::renderDataTable(
      DT::datatable(alig,
                    colnames = coln,
                    rownames = FALSE,
                    selection = 'single',
                    width = "100%",
                    options = list(scrollX = TRUE)
      )
      %>% formatRound(columns = "mean_mz",digits = 4)
      %>% formatRound(columns = "mean_RT",digits = 1)
    )
  })
  
  #### Alignment cluster ####
  
  hcData <- reactive({
     
    updateReactive$v
    req(aligR())
    # check that more than one files are present
    req(ncol(aligR()) >= 7)
    # use clustering data from environment if available, skip calculation
    if (exists("clusteringData"))
      return(clusteringData)
     
    prb <- Progress$new(session, 1, 4)
    prb$set(m = "Clustering rows...")
    intens <- aligR()[, grep("^Int_", colnames(aligR()))]
    sumStat <- apply(intens, 1, max)
    prb$inc(1)
    intens <- sweep(intens, 1, sumStat, "/")
    intens <- as.matrix(intens)
    prb$inc(1)
    distM <- parallelDist::parDist(intens, method = "dtw", threads = input$numcores)
    prb$inc(1)
    hc <- hclust(distM, method = "average")
    # if annotation has happened, need to add labels to hc
    if ("name" %in% colnames(aligR())) {
      alig <- aligR()
      alig$name[is.na(alig$name)] <- ""
      lab <- mapply(paste, hc$order, alig$name, MoreArgs = list(sep= "-"))
      lab <- substr(lab, 1, 12)
      hc$labels <- lab
    }
    # Copy to global environment so that the data can still be used after closing/crash
    clusteringData <<- hc  
    prb$close()
    hc
  })
  
  intDendroRanges <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$dendroDblClick, {
    if (!is.null(input$dendroBrush)) {
    brush <- input$dendroBrush
    xmin <- round(brush$xmin)
    if (xmin < 1)
      xmin <- 1
    xmax <- round(brush$xmax)
    if (xmax == -1)
      xmax <- xmin
    intDendroRanges$x <- c(xmin, xmax)
    intDendroRanges$y <- c(brush$ymin, brush$ymax)
    } else {
      intDendroRanges$x <- NULL
      intDendroRanges$y <- NULL
    }
  })
  output$intDendro <- renderPlot({
     
    req(hcData())
    dc <- as.dendrogram(hcData())  # convert to dendrogram for better plotting
    plot(dc, leaflab = "perpendicular", type = "rectangle", horiz = FALSE,
         main = "Alignment table cluster dendrogram", xlim = intDendroRanges$x,
         ylim = intDendroRanges$y)
  })

  
  # denroLocStore stores the old value of dendro hover so that the plot 
  # remains even if you move the cursor off the dendrogram
  dendroLocStore <<- NULL
  dendroLocX <- reactive({
     
    if (is.null(input$dendHover)) {
      return(dendroLocStore)
    } else {
      dendroLocStore <<- input$dendHover$x
      return(input$dendHover$x)
    }
  })
 
  output$intTrendDendro <- renderPlot({
    req(dendroLocX())
    aligRow <- hcData()$order[round(dendroLocX())]
    mmz <- aligR()[aligRow, "mean_mz"]
    mrt <- aligR()[aligRow, "mean_RT"]
    tl <- sprintf("m/z: %.4f rt: %.2f", mmz, mrt)
    if ("name" %in% colnames(aligR()) && !is.na(aligR()$name[aligRow]))
      tl <- paste0(tl, ", ", aligR()$name[aligRow])
    alignmentTrendPlotten(aligRow, tl, tabelle = aligR())
  })
  
  # for similar trends
  dendroSelectedIds <- reactiveVal()
  observeEvent(input$dendroClick, {
    req(input$dendroBrush)
    brush <- input$dendroBrush
    xmin <- round(brush$xmin)
    if (xmin < 1)
      xmin <- 1
    xmax <- round(brush$xmax)
    if (xmax == -1)
      xmax <- xmin
    alig <- aligR()
    alig <- alig[hcData()$order[xmin:xmax], ]
    dendroSelectedIds(alig$alignmentID)

    # columns not needed
    alig$intRange <- NULL
    alig$maxInt <- NULL
    alig$MS2Fit <- NULL
    alig$row <- NULL
    
    coln <- colnames(alig)
    output$dendroSelected <- DT::renderDataTable(
      DT::datatable(alig,
                    colnames = coln,
                    rownames = FALSE,
                    selection = 'single',
                    width = "100%",
                    options = list(scrollX = TRUE)
      )
      %>% formatRound(columns = "mean_mz",digits = 4)
      %>% formatRound(columns = "mean_RT",digits = 1)
    )
  })

  #### Alignment similar trends ####

  observeEvent(input$overviewSelected_cell_clicked, {
    req(input$overviewSelected_rows_selected)
     
    x <- overviewSelectedIds()[input$overviewSelected_rows_selected]
    alignmentIdSelected(x)
  }, ignoreInit = TRUE)

  observeEvent(input$intRangeSelected_cell_clicked, {
    req(input$intRangeSelected_rows_selected)
    x <- intRangeSelectedIds()[input$intRangeSelected_rows_selected]
    alignmentIdSelected(x)
  }, ignoreInit = TRUE)

  observeEvent(input$dendroSelected_cell_clicked, {
    req(input$dendroSelected_rows_selected)
    x <- dendroSelectedIds()[input$dendroSelected_rows_selected]
    alignmentIdSelected(x)
  }, ignoreInit = TRUE)
  
  observe({
    input$AlignmentTable_cell_clicked
    input$overviewSelected_cell_clicked
    input$intRangeSelected_cell_clicked
    input$dendroSelected_cell_clicked
    req(!is.na(input$trendFac))
    req(input$trendFac >= 0.3)  # otherwise will never stop computations
    req(alignmentIdSelected())
    FillSimilarTrendsTable(alignmentIdSelected())
    SimilarTrends_plotten(alignmentIdSelected())  # otherwise old plot remains
  })
  
  observeEvent(input$SimilarTrendsTable_cell_clicked, {
    req(input$SimilarTrendsTable_rows_selected)
    zeile <- input$SimilarTrendsTable_rows_selected
    SimilarTrends_plotten(alignmentIdSelected(), highlight = similarTrends()[zeile])
  })
  
  output$whichSimTren <- renderText({
    req(input$simTrenHover)
    req(simTrenInten())
     
    r <- nearPoints(simTrenInten(), input$simTrenHover, xvar = "samp", yvar = "int", maxpoints = 1)
    if ("name" %in% colnames(aligR())) {
      i <- aligR()[aligR()$alignmentID == r$alignmentID, c("mean_mz", "mean_RT", "name")]
      return(sprintf("m/z = %.4f Da, RT = %.2f min., %s", 
                     i$mean_mz, i$mean_RT, i$name))
    } else {
      i <- aligR()[aligR()$alignmentID == r$alignmentID, c("mean_mz", "mean_RT")]
      return(sprintf("m/z = %.4f Da, RT = %.2f min.", i$mean_mz, i$mean_RT))
    }
  })
  
  # Annotation ####
  
  observeEvent(input$dbHelp, {
    showModal(modalDialog(
      title = "Database files and formats",
      p(
        "There are several options for databases. The system will choose the
        correct annotation method based on the format of the database file, or, 
        in the case of no file, an adduct annotation is done."
      ),
      p(
        "1: MS2 'Spektrendatenbank' SQLite database, E.g. MS2_db_v9.db. To add standards to 
        DB, ask G2."
      ),
      p(
        "2: Custom database in yaml format for MS2 search (incl. fragments and neutral-loss search) 
        for fomating examples ask around G2. E.g. example_db.yaml"
      ),
      p(
        "3: Custom database in csv format for only doing MS1 (m/z, RT) suspect 
        screening (MS2 is ignored). Comma separated with headers: 
        'name', 'mz', 'rt' in any order. RT in minutes and optional.
        Warning: High number of errors expected. Estimate 
        false positives and false negatives.
        For examples ask around G2. E.g. example_sus_list.csv"
      ),
      p(
        "4: Labelling database in SQLite format. Ask G2."
      ),
      p(
        "5: No database file."
      ),

      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  shinyFileChoose(input, 'ms2Db', roots=home3, session=session, filetypes=c('db', 'yaml', 'csv'))
  
  output$dbLoc <- renderText(as.character(parseFilePaths(home, input$ms2Db)$datapath))
  
  # populate source select list and chrom method list
  observe({
    dbPath <- as.character(parseFilePaths(home, input$ms2Db)$datapath)
    req(length(dbPath) > 0)
    if (grepl("\\.db$", dbPath)) {
      con <- DBI::dbConnect(RSQLite::SQLite(), dbPath)
      if ("experimentGroup" %in% DBI::dbListTables(con)) {
        updateSelectInput(
          session, 
          "annotExpSource", 
          choices = DBI::dbReadTable(con, "experimentGroup")$name,
          selected = DBI::dbReadTable(con, "experimentGroup")$name[1]
        )
      }
      if ("retention_time" %in% DBI::dbListTables(con)) {
        updateSelectInput(
          session, 
          "annotChromMeth", 
          choices = DBI::dbReadTable(con, "retention_time")$chrom_method,
          selected = DBI::dbReadTable(con, "retention_time")$chrom_method[1]
        )
      }
      DBI::dbDisconnect(con)
    }
  })
  
  observeEvent(input$annotGo, {
    dbPath <- as.character(parseFilePaths(home, input$ms2Db)$datapath)
    if (is.null(grouped))
      showNotification("First align the samples", type = "error")
    validate(need(!is.null(grouped), "First align the samples"))
    progress <- shiny::Progress$new()
    progress$set(message = "Annotating alignment table...")
    
    # if no db chosen, annotate with peaklist component info (Cl and Br)
    if (length(dbPath) == 0) {
      showNotification("No database chosen, annotating just with components")
      progress$set(detail = "Searching Peaklists...")
      
      annotationTableNew <- ntsworkflow::annotate_grouped_components(
        sampleListLocal = sampleList,
        peakListList = peaklist,
        alignmentTable = grouped,
        numcores = input$numcores
      )
    } else if (grepl("\\.csv$", dbPath)) {  
      # csv library: just do m/z-rt screening in alignment table
      progress$set(detail = "Searching average m/z, rt")
      annotationTableNew <- ntsworkflow::annotate_grouped_mz_rt(
        grouped, dbPath, input$annotMzTolmDa / 1000, input$annotRtTolM
      )
      
      
    } else if (grepl("\\.db$", dbPath) || grepl("\\.ya?ml$", dbPath)) {
      # yaml or SQLite db: do MS2 searching
      
      # first check
      
      progress$set(detail = "Searching MS2...")
      # make sure the source and retention time is selected
      expSource <- if (is.null(input$annotExpSource)) {
        showNotification("No source chosen or none available, setting BfG as standard")
        "BfG"
      } else {
        input$annotExpSource
      }
      chromMethod <- if (is.null(input$annotChromMeth)) {
        showNotification(
          "No retention times chosen or none available, setting BfG method as standard"
        )
        "dx.doi.org/10.1016/j.chroma.2015.11.014"
      } else {
        input$annotChromMeth
      }
      
      annotationTableNew <- ntsworkflow::annotate_grouped(  # instrument: default settings
        sampleListLocal = sampleList,
        peakListList = peaklist,
        alignmentTable = grouped,
        db_path = dbPath,
        threshold_score = input$annotThresh,
        mztolu = input$annotMzTolmDa / 1000,
        rttol = input$annotRtTolM,
        polarity = input$annotPol,
        CE = input$annotCE,
        CES = input$annotCES,
        mztolu_ms2 = input$annotDpWind / 1000,
        rtoffset = input$annotRtOffset,
        intCutData = input$annotIntCut,
        numcores = input$numcores,
        datenListLocal = datenList,
        expGroups = expSource,
        chrom_method = chromMethod
      )
    } else {
      showNotification("Unknown database file", type = "error", duration = NULL)
    }
    
    if (is.null(annotationTableNew)) {
      showNotification("No compounds from the database were found in the alignment table.", 
                       duration = NULL, id = "nothingFound")
      progress$close()
    }
    req(!is.null(annotationTableNew))
    
    if (input$annotAppend) {
      progress$set(detail = "Appending annotations...")
      
      # need to have unified column names
      allColNm <- union(colnames(annotationTable), colnames(annotationTableNew))
      newColAnT <- setdiff(allColNm, colnames(annotationTable))
      newColAnTN <- setdiff(allColNm, colnames(annotationTableNew))

      for (cn in newColAnT)
        annotationTable[, cn] <<- vector(unlist(dplyr::summarise_all(annotationTableNew[, cn], class)), nrow(annotationTable))
      for (cnn in newColAnTN)
        annotationTableNew[, cnn] <- vector(unlist(dplyr::summarise_all(annotationTable[, cnn], class)), nrow(annotationTableNew))
      
      #browser()
      stopifnot(setequal(colnames(annotationTable), colnames(annotationTableNew)))
      annotationTable <<- rbind(annotationTable, annotationTableNew)
    } else {
      annotationTable <<- annotationTableNew
    }
    
    progress$set(detail = "Annotation complete, cleaning up...")

    # build annotation table
    FillAnnotationTable()
    # need to invalidate the aligR reactive so that it includes the annotations
    updateReactive$v <- updateReactive$v + 1
    if (exists("clusteringData"))
      rm(clusteringData, pos = globalenv())
    #create log_file
    threshold_score <- input$annotThresh
    mztolu <- input$annotMzTolmDa / 1000
    rttol <- input$annotRtTolM
    polarity <- input$annotPol
    CE <- input$annotCE
    CES <- input$annotCES
    mztolu_ms2 <- input$annotDpWind / 1000
    rtoffset <- input$annotRtOffset
    log_file <<- create_log_file()
    
    FillAlignmentTable()  # to include annotations
    gc()
   
    progress$close()
  })
  
  alignmentIdSelected <- reactiveVal(NULL)
  
  shinyFileSave(input, 'annotationTableExport', roots=home, session=session)
  observeEvent(input$annotationTableExport,{
    req(is.list(input$annotationTableExport))
     
    saveto <- as.character(parseSavePath(home, input$annotationTableExport)$datapath)
    if (!is.null(annotationTable) && !is.null(saveto)) {
      # merge annotation table with grouped
      grouped2 <- as.data.frame(grouped)
   
      newTable <- merge(grouped2, annotationTable, all.x = TRUE, by = "alignmentID")
      rtCols <- grep("RT", colnames(newTable))
      mzCols <- grep("mz", colnames(newTable))
      intCols <- grep("Int", colnames(newTable))
      for (col in rtCols)
        newTable[, col] <- round(newTable[, col] / 60, 2)
      # round mz to 4 d.p.
      for (col in mzCols)
        newTable[, col] <- round(newTable[, col], 4)
      # round int to 3 d.p.
      for (col in intCols)
        newTable[, col] <- round(newTable[, col], 3)
      
      
      newTable <- tibble::as_tibble(newTable)
       
      # add samples names
      # get sample names from table
      intColNames <- colnames(grouped2)[grep("^Int_", colnames(grouped2))]
      sampsInTable <- as.numeric(stringr::str_match(intColNames, "_(\\d+)$")[,2])
      
      for (sId in sampsInTable) {
        # get column number of first column for sample
         
        firstCol <- grep(paste0("^PeakID_", sId, "$"), colnames(newTable))
        # insert sample name
        newTable <- tibble::add_column(newTable, sampleName = basename(sampleList[sId, "File"]), .before = firstCol)
        colnames(newTable)[firstCol] <- paste0("sampleName_", sId)
      }
      write.csv(newTable, file = saveto, row.names = FALSE)  
    }  
  })
  
  observeEvent(input$AlignmentAnnotTable_rows_selected, {
    cr <- input$AlignmentAnnotTable_rows_selected
    alignmentIdSelected(groupedWithAnnotation[cr, "alignmentID"])
    if (groupedWithAnnotation[input$AlignmentAnnotTable_rows_selected, "multHits"])
      fillMultiHitsTable() else output$multiHitsTable <- NULL
  })
  # plots for showing annotations
  output$annotationXIC <- renderPlot({
    validate(need(input$AlignmentAnnotTable_rows_selected, "Load DB, click 'Annotate', select row"))
    alignmentXICsPlotten(alignmentIdSelected())
  })
  output$annotationMS1 <- renderPlot({
    validate(need(input$AlignmentAnnotTable_rows_selected, ""))
    validate(
      need(
        !is.na(groupedWithAnnotation[input$AlignmentAnnotTable_rows_selected, "name"]), ""
      )
    )
    annotationMS1Plotten(input$AlignmentAnnotTable_rows_selected)
  })
  output$annotationMS2 <- renderPlot({
    validate(need(input$AlignmentAnnotTable_rows_selected, ""))
    validate(
      need(
        !is.na(groupedWithAnnotation[input$AlignmentAnnotTable_rows_selected, "name"]), ""
      )
    )
    annotationMS2Plotten(input$AlignmentAnnotTable_rows_selected)
  })
  output$annotationTrend <- renderPlot({
    validate(need(input$AlignmentAnnotTable_rows_selected, ""))
    alignmentTrendPlotten(input$AlignmentAnnotTable_rows_selected, tabelle = grouped)
  })
  
  
  # Save request batch ####
  # This has to be placed at the very end so that other reactive events are done
  # before saving
  saveLater <- reactiveVal()
  observeEvent(saveLater(), {
    save_data(file.path(globalwd, paste0(input$saveBatchName, ".RDS")))
  })
  
  # Exit ####
  session$onSessionEnded(function() {
    if (exists("cl")) {
      parallel::stopCluster(cl)
      rm(cl)
    }
  
    # remove files created by GenForm
    suppressWarnings({ # if the files were not created, gives a bunch of warnings
      file.remove(file.path(globalwd, "GenFormMS1.txt"))
      file.remove(file.path(globalwd, "GenFormMS2.txt"))
      file.remove(file.path(globalwd, "GenFormMS2_GenForm.txt"))
      file.remove(file.path(globalwd, "GenFormMS2_GenForm_log.txt"))
      file.remove(file.path(globalwd, "GenFormMS2_GenForm_an.txt"))
      file.remove(file.path(globalwd, "GenFormMS2_GenForm_loss.txt"))
      file.remove(file.path(globalwd, "GenFormMS2_GenForm_cleanMSMS.txt"))
    })
    
    stopApp()
  }) 
  
}

shinyApp(ui = ui, server = server)


