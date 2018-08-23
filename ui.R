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
  titlePanel(title=div(img(src="MTEK LOGO - v0.4.png", width = 200), "HECT - Highly Efficient Clinical Trial simulator")),
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
                     numericInput("nt", "Number of treatments:", value = 3)),
                 div(style="", radioButtons("efftype", "Primary outcome type", 
                                            c('Continuous' = "absolute",
                                              "Proportion (between 0 and 1)" = "rate"),
                                            selected = "absolute")),
                 
                 div(style="display: inline-block;vertical-align:top; ", uiOutput("effboxes")),
                 
                 div(style="", uiOutput("adhere")),
                 
                 checkboxInput("so", "Add secondary outcome"),
                 conditionalPanel(
                   condition = "input.so == true",
                   radioButtons("efftypeso", "Secondary outcome type", 
                                c('Continuous' = "absolute",
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
                                    value = 0)),
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
                                  width = "270px", value = 30)),
                 bsTooltip("burnin", "This is the number of patients that need to be enrolled into the trial before adaptation begins.",
                           "right", options = list(container = "body")),
                 
                 div(style="display: inline-block;",
                     numericInput("batchsize", "Sample size between two interim looks:",
                                  width = "270px", value = 20)),
                 bsTooltip("batchsize", "This is the number of patients enrolled into the trial between two adaptations. If this number is small adaptations are performed more frequently resulting in more computation and increasing simulation time.",
                           "right", options = list(container = "body")),
                 
                 div(style="display: inline-block;",
                     numericInput("max", "Maximum total sample size:",
                                  width = "270px", value = 500)),
                 bsTooltip("max", "If the total number of patients enrolled into the trial reaches this number, the trial will be terminated.",
                           "right", options = list(container = "body")),
                 
                 div(style="display: inline-block;",
                     numericInput("ec", "Enrollment cost per patient ($):",
                                  width = "270px", value = 0))
               ),
               
               mainPanel(
                 
                 tabsetPanel( 
                   tabPanel(# START 1st row tabs
                     
                         fluidRow(align = "left", style = "width:1220px; margin-left:10px",
                                  h1("Trial design properties")
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
                                    numericInput("Tpower", "Maximum run time (s)", value = NULL),
                                    bsTooltip("Tpower", "The simulation is terminated if run time exceeds this limit.",
                                              "right", options = list(container = "body")))
                           ),
                           column(5, align = "right",
                                  br(),
                                  div(actionButton("button2", "Run", 
                                                   style="padding:4px; font-size:100%;
                                              color:#FFF; background-color: #0095ff; border-color:#07c"))
                           )
                           
                         )
                         
                   ),  # END 1st row tabs
                     
                     
                     tabsetPanel( # START 2nd row tabs
                       tabPanel("Power",
                                br(),
                                DT::dataTableOutput("power")
                       ),
                       tabPanel("Sample size distribution/cost evaluation", 
                                plotOutput("ssdist"), DT::dataTableOutput("sssummary"),
                                htmlOutput('sssave'), plotOutput("cdist"), DT::dataTableOutput("csummary"),
                                htmlOutput('csave')),
                       tabPanel("Type I error rate",
                                
                                # wellPanel(tags$div(tags$em(
                                #   HTML(paste(tags$span(style="color:darkred", "Note: Changing the effect sizes and adherence rates does not affect the type I error rate as it is computed under the null hypothesis (i.e., none of the treatments are effective).")))
                                # ))),
                                # numericInput("Malpha", "Number of simulations", value = 100),
                                # bsTooltip("Malpha", "Number of Monte Carlo simulations to estimate the type I error rate.",
                                #           "right", options = list(container = "body")),
                                # checkboxInput("tlimita", "Set a limit for run time"),
                                # conditionalPanel(
                                #   condition = "input.tlimita == true",
                                #   numericInput("Talpha", "Maximum permitted run time in seconds", value = NULL),
                                #   bsTooltip("Talpha", "The simulation is terminated if run time exceeds this limit.",
                                #             "right", options = list(container = "body"))),
                                
                                #actionButton("button01", "Compute"),
                                br(),
                                DT::dataTableOutput("alpha")
                       ),
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
                                br(),
                                div(actionButton("compRCT", "Run", 
                                                 style="padding:4px; font-size:100%;
                                                 color:#FFF; background-color: #0095ff; border-color:#07c")),
                                
                                plotOutput("BRATvsBRCTplot"),
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
                              h1("Single trial simulation")
                     ),
                     fluidRow(align = "right", style = "margin-right:108px",
                              actionButton("button", "Run", 
                                           style="padding:4px; font-size:100% ;margin-bottom:19px;
                                          color:#FFF; background-color: #0095ff; border-color:#07c")
                              
                     # div(style="display:inline-block; text-align: middle; width:280px",
                     #     # fluidRow(align = "middle",
                     #     #          div("Single trial simulation \n", 
                     #     #              style = "font-size:120%; font-weight:bold; margin-bottom:10px")
                     #     # ),
                     #     fluidRow(align = "right",
                     #              actionButton("button", "Run", 
                     #                           style="padding:4px; font-size:100% ;margin-bottom:19px;
                     #                      color:#FFF; background-color: #0095ff; border-color:#07c")
                     #              )
                        # other colour option
                         # color: #39739d; background-color: #E1ECF4; border-color: #96bdd9
                     ),
                     
                     # div(style="display:inline-block; margin-left:25px; margin-right:25px;",
                     #     actionButton("button", "Run single trial simulation", style="vertical-align: middle; padding:4px; font-size:100%")),
                     
                            tabsetPanel(
                              tabPanel("Design",
                                       bsAlert('alert0'), plotOutput("designPlot")),
                              tabPanel("Data",
                                       bsAlert('alert0'), plotOutput("dataPlot"),
                                       DT::dataTableOutput("datatable")),
                              tabPanel("Probability of Superiority",
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
                 div(style="", radioButtons("efftypess", "Primary outcome Type", c('Continuous' = "absolute",
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
    tabPanel("About",
             fluidRow(uiOutput('about'))

    )
  )
)
)
