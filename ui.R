library(shiny)
library(shinyBS)
library(shinydashboard)
library(markdown)
names0 <- c('nTreatment', paste('effect', 1:5, sep = ''), 'effectType',
            paste('adherence', 1:5, sep = ''), 'burninSampleSize', 'interimSampleSize', 'supThreshold',
            'futThreshold', 'maxSampleSize', paste('power', 1:4, sep = '') , 
            paste('typeIerror', 1:4, sep = ''), 'ESS', 'minSS', 'firstQuartileSS',
            'medianSS', 'thirdQuartileSS', 'maxSS', 'ESavedSS')


shinyUI(fluidPage(#theme="bootstrap.css",

  # Application title
  titlePanel(title=div(img(src="43186754.png", width = 150), "HECT - Highly Efficient Clinical Trial simulator")),
  # navbarPage('',
  #                  dashboardPage(
  #                    dbHeader,
  #                    dashboardSidebar(),
  #                    dashboardBody()
  #                  ),
  tabsetPanel(
    tabPanel("Trial Simulation",
             sidebarLayout(
               sidebarPanel(
                 div(style="display: inline-block;vertical-align:top; ", 
                     numericInput("nt", "Number of treatments:", value = 3, min = 2, max = 10)),
                 div(style="", radioButtons("efftype", "Primary outcome type", 
                                            c('Continuous (standardized: variance = 1)' = "absolute",
                                              "Proportion (between 0 and 1)" = "rate"),
                                            selected = "absolute")),
                 # checkboxInput("goodout", "Bigger is better"),
                 # bsTooltip("goodout", "Check this box if a larger outcome value is associated with a positive effect.",
                 #           "right", options = list(container = "body")),
                 div(style="display: inline-block;vertical-align:top; ", uiOutput("effboxes")),
                 
                 div(style="", uiOutput("adhere")),
                 
                 checkboxInput("so", "Add secondary outcome"),
                 conditionalPanel(
                   condition = "input.so == true",
                   radioButtons("efftypeso", "Secondary outcome type", 
                                c('Continuous (standardized: variance = 1)' = "absolute",
                                  #"OR" = "OR",
                                  "Proportion (between 0 and 1)" = "rate"),
                                selected = "absolute"),
                   div(style="display: inline-block;vertical-align:top; ", uiOutput("effboxesso"))
                 ),
                 
                 checkboxInput("platf", "Platform design", value = FALSE),
                 bsTooltip("platf", "When this box is checked the trial starts without the last arm. The last arm is added when at least one arm is dropped.",
                           "right", options = list(container = "body")),
                 
                 radioButtons("compCon", "Comparison:",
                              choices = c("Compare all arms simultaneously" = 'FALSE',
                                "Compare arms against reference treatment" = 'TRUE')),
                 conditionalPanel(
                   condition = "input.compCon == 'TRUE'",
                   div(style="display: inline-block;vertical-align:top; ",
                       numericInput("MID", "Minimum important difference (MID):",
                                    value = 0, min = 0, max = 1, step = 0.01)),
                   bsTooltip("MID", "The minimum important difference is restricted to be between 0 and 1. The sign is automatically determined based on the direction of effect.",
                             "right", options = list(container = "body")),
                   div(style="",sliderInput("uppfut", "Futility probability threshold:", step = .005,
                                            min = 0, max = 1, value = 0.95)),
                   bsTooltip("uppfut", "If the probability that a treatment effect being smaller than MID, for any of the treatment arms, falls below this threshold, the arm will be stopped (i.e. no more patients are assined to it) .",
                             "right", options = list(container = "body")),
                   
                   div(style="",sliderInput("upthreshT",
                                            "Superiority probability threshold:", step = .005,
                                            min = 0,
                                            max = 1,
                                            value = 0.975)),
                   bsTooltip("upthreshT", "If the probability that a treatment is better than the reference exceeds this threshold, at any point after burn-in and before the maximum sample size is reached, the trial will be terminated.",
                             "right", options = list(container = "body"))
                 ),
                 
                 conditionalPanel( # "Compare all arms simultaneously" = 'FALSE'
                   condition = "input.compCon == 'FALSE'",
                   checkboxInput("adapt", "Employ response-adaptive randomization", value = TRUE),
                   bsTooltip("adapt", "Based on patient responses accrued, allocate more patients to the better treatment(s).",
                             "right", options = list(container = "body")),
                   
                   div(style="",sliderInput("lowthresh", "Stop arm for futility at:", step = .005,
                                            min = 0, max = 1, value = 0.01)),
                   bsTooltip("lowthresh", "If the probability that a treatment is better than the rest, for any of the treatment arms, falls below this threshold, the arm will be stopped (i.e. no more patients are assined to it) .",
                             "right", options = list(container = "body")),
                   
                   div(style="",sliderInput("upthreshF",
                                            "Terminate trial for superiority at:", step = .005,
                                            min = 0,
                                            max = 1,
                                            value = 0.975)),
                   bsTooltip("upthreshF", "If the probability that a treatment is better than the rest exceeds this threshold, at any point after burn-in and before the maximum sample size is reached, the trial will be terminated.",
                             "right", options = list(container = "body"))
                 
                 ),
                 
                 div(style="display: inline-block;",
                     numericInput("burnin", "Burn-in sample size:",
                                  width = "270px", value = 100, min = 10)),
                 bsTooltip("burnin", "This is the number of patients that need to be enrolled into the trial before adaptation begins.",
                           "right", options = list(container = "body")),
                 
                 div(style="display: inline-block;",
                     numericInput("batchsize", "Sample size between two interim looks:",
                                  width = "270px", value = 100, min = 1)),
                 bsTooltip("batchsize", "This is the number of patients enrolled into the trial between two adaptations. If this number is small adaptations are performed more frequently resulting in more computation and increasing simulation time.",
                           "right", options = list(container = "body")),
                 
                 div(style="display: inline-block;",
                     numericInput("max", "Maximum total sample size:",
                                  width = "270px", value = 500, min = 10)),
                 bsTooltip("max", "If the total number of patients enrolled into the trial reaches this number, the trial will be terminated.",
                           "right", options = list(container = "body")),
                 
                 div(style="display: inline-block;",
                     numericInput("ec", "Enrollment cost per patient ($):",
                                  width = "270px", value = 0, min = 0))
               ),
               
               mainPanel(
                 
                 tabsetPanel( 
                   tabPanel(# START 1st row tabs
                     
                         fluidRow(align = "left", style = "width:1220px; margin-left:10px",
                                  h3("Trial design properties")
                         ),
                         br(),
                         fluidRow(
                           column(2, align = "middle",
                                  numericInput("M", label = "Number of simulations", value = 100),
                                  bsTooltip("M", "Number of Monte Carlo simulations to estimate the power.",
                                            "right", options = list(container = "body")) 
                           ),
                           column(2, align = "middle",
                                  br(),
                                  checkboxInput("tlimit", "Set run time limit")
                           ),
                           column(2, align = "middle",
                                  conditionalPanel(
                                    condition = "input.tlimit == true",
                                    numericInput("Tpower", "Maximum run time (sec)", value = NULL),
                                    bsTooltip("Tpower", "The simulation is terminated if run time exceeds this limit.",
                                              "right", options = list(container = "body")))
                           ),
                           column(5, align = "right", style = "margin-right:80px", 
                                  br(),
                                  div(actionButton("button2", "Run", #qqq
                                                   style="padding:4px; font-size:100%;
                                              color:#FFF; background-color: #0095ff; border-color:#07c"))
                           )
                           
                         )
                         
                   ),  # END 1st row tabs
                     
                     
                     tabsetPanel( # START 2nd row tabs
                       tabPanel("Power and type I error",
                                br(),
                                DT::dataTableOutput("power")
                                # ,
                                # DT::dataTableOutput("alpha")
                       ),
                       tabPanel("Sample size distribution/cost evaluation", 
                                plotOutput("ssdist"), DT::dataTableOutput("sssummary"),
                                htmlOutput('sssave'), plotOutput("cdist"), DT::dataTableOutput("csummary"),
                                htmlOutput('csave')),
                      
                       tabPanel("Saved simulation configurations and results",
                                br(),
                                actionButton("submit", "Save recent results"),
                                bsTooltip("submit", "This button saves results to the application shared data folder. To see the most recent results in the table, reload after saving. The application data folder is cleared on a daily basis and can be modified by other users. To save the table contents to your machine use the download button below.",
                                          "right", options = list(container = "body")),
                                actionButton("load", "Load saved results"),
                                bsTooltip("load", "Load simulation results from the application shared data folder.",
                                          "right", options = list(container = "body")),
                                br(),
                                br(),
                                selectInput('column', 'Columns', choices = names0, 
                                            selected = names0[c(1,13:18, 22, 26)],
                                            multiple = T),
                                actionButton("deleteRows", "Delete rows"), 
                                DT::dataTableOutput('responses'), 
                                tabPanel("Files list", DT::dataTableOutput("tbl")),
                                downloadButton("downloadData", "Download"),
                                bsTooltip("downloadData", "Use this button to download simulation data.",
                                          "right", options = list(container = "body")),
                                actionButton("clearData", "Clear data"),
                                bsTooltip("clearData", "By clicking this button all the data in the application shared data folder is cleared.",
                                          "right", options = list(container = "body"))
                                #tags$div(tags$label("Upload", class="btn btn-primary",
                                #                    tags$input(id = "fileIn", webkitdirectory = TRUE, type = "file", style="display: none;", onchange="pressed()")))
                                #, shinyDirButton('directory', 'Folder select', 'Please select a folder')
                       ),
                       tabPanel("Comparison with RCT",
                                checkboxInput("mc", "Correct for multiple comparison", value = TRUE),
                                bsTooltip("mc", "When this box is checked multiple comparison correction is made to control the overall type I error at 5%.",
                                          "right", options = list(container = "body")),
                                radioButtons("fix", "Fixed for RCT:",
                                             choices = c("Sample size" = 'ss',
                                                         "Power" = 'pow')),
                                bsTooltip("fix", "If sample size is selected, the sample size for RCT is fixed at maximum sample size and the designs are compared with respect to power. Otherwise, the power is fixed at 80% and the required sample size for RCT to acheive 80% is compared to the expected sample size for HECT.",
                                          "right", options = list(container = "body")),
                                br(),
                                # div(actionButton("compRCT", "Run", 
                                #                  style="padding:4px; font-size:100%;
                                #                  color:#FFF; background-color: #0095ff; border-color:#07c")),
                                #textOutput('test'),
                                plotOutput("HECTvsBRCTplot"),
                                br()
                                #,
                                #DT::dataTableOutput('RCT')
                       )
                       
                     ) # END 2nd row tabs
                     
                   ), # END top panel 
                
                 ##############################################################################
                 # Single trial simulation ####################################################
                 ##############################################################################
                 
                   br(),
                   br(),
                   
                   tabsetPanel(
                   tabPanel(
                     
                     fluidRow(align = "left", style = "width:1220px; margin-left:10px",
                              h3("Single trial simulation")
                     ),
                     fluidRow(
                       column(6, align = "middle",
                              div("TEXT FOR SPACING", style = "color:#FFF")),
                       
                       column(5, align = "right", style = "margin-right:80px", 
                              div(actionButton("button", "Run", #qqq
                                               style="padding:4px; font-size:100%;
                                               color:#FFF; background-color: #0095ff; border-color:#07c")))
                              
                        # other colour option
                         # color: #39739d; background-color: #E1ECF4; border-color: #96bdd9
                     ),
                     
                            tabsetPanel(
                              tabPanel("Design",
                                       bsAlert('alert0'), plotOutput("designPlot")),
                              tabPanel("Data",
                                       bsAlert('alert0'), plotOutput("dataPlot"),
                                       DT::dataTableOutput("datatable")),
                              tabPanel("Probability of superiority",
                                       plotOutput("supPlot"),
                                       DT::dataTableOutput("psuptable")),
                              tabPanel("Estimates",
                                       plotOutput("estPlot"),
                                       DT::dataTableOutput("estable")),
                              tabPanel("Posterior",
                                       plotOutput("postPlot")),
                              tabPanel("Secondary analysis",
                                       
                                       # wellPanel(
                                       #   conditionalPanel(
                                       #
                                       #     condition = "input.so != true",
                                       #     tags$div(tags$b(
                                       #       HTML(paste(tags$span(style="color:darkred", "Oops ... You did not include a secondary outcome!")))),id="message")
                                       #   ))
                                       bsAlert('alert1')
                                       , bsAlert('alertso'),
                                       plotOutput('secData'),
                                       plotOutput('secEst'),
                                       DT::dataTableOutput("sec.estable")
                              )
                              
                              
                              
                            )))
               )
             )
    ),
    tabPanel("Sample Size Calculation",
             sidebarLayout(
               sidebarPanel(

                 div(style="display: inline-block;vertical-align:top; ", numericInput("ntss", "Number of treatments:",
                                                                                      value = 3)),
                 div(style="", radioButtons("efftypess", "Primary outcome type", c('Continuous' = "absolute",
                                                                            #"OR" = "OR",
                                                                            "Proportion (between 0 and 1)" = "rate"),
                                            selected = "absolute")),
                 div(style="display: inline-block;vertical-align:top; ", uiOutput("effboxesss")),
                 div(style="display: inline-block;vertical-align:top; ", uiOutput("varboxesss")),
                 div(style="", sliderInput("power", "Required power:",
                                                                                      min = 0, max = 1, value = 0.80)),
                 div(style="", sliderInput("alpha", "Required type I error rate:",
                                                                                      min = 0, max = 1, value = 0.05)),
                 div(style="", sliderInput("dropout", "Drop-out rate:",
                                                                                      min = 0, max = 1, value = 0.2))

               ),
               mainPanel(
                 br(),
                 actionButton("button1", "Calculate", style="padding:4px; font-size:100%;
                                  color:#FFF; background-color: #0095ff; border-color:#07c"),
                 bsAlert('alertss'),
                 htmlOutput('ss')

               )
             )
    ),
    tabPanel("User Manual",
             fluidRow(uiOutput('about'))

    )
  )
)
)
