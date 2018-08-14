library(shiny)
library(shinyBS)
library(shinydashboard)
library(markdown)
names0 <- c('nTreatment', paste('effect', 1:5, sep = ''), 'effectType',
            paste('adherence', 1:5, sep = ''), 'burninSampleSize', 'interimSampleSize', 'supThreshold',
            'futThreshold', 'maxSampleSize', 'power', 'typeIerror', 'ESS', 'minSS', 'firstQuartileSS',
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
                 
                 div(style="display: inline-block;vertical-align:top; ", numericInput("nt", "Number of treatments:",
                                                                                      value = 3)),
                 div(style="", radioButtons("efftype", "Primary outcome type", c('Continuous' = "absolute",
                                                                                 #"OR" = "OR",
                                                                                 "Proportion (between 0 and 1)" = "rate"),#, "Count (greater than 0)" = "count"),
                                            selected = "absolute")),
                 checkboxInput("con", "First arm is control", value = TRUE),
                 checkboxInput("ai", "Add last arm later", value = FALSE),
                 conditionalPanel(
                   condition = "input.ai == true",
                   div(style="display: inline-block;vertical-align:top; ", numericInput("ail", "Add last arm at interim look:",
                                                                                        value = 0))
                 ),
                 div(style="display: inline-block;vertical-align:top; ", uiOutput("effboxes")),
                 checkboxInput("so", "Add secondary outcome"),
                 conditionalPanel(
                   condition = "input.so == true",
                   radioButtons("efftypeso", "Secondary outcome type", c('Continuous' = "absolute",
                                                                         #"OR" = "OR",
                                                                         "Proportion (between 0 and 1)" = "rate"),#, "Count (greater than 0)" = "count"),
                                selected = "absolute"),
                   div(style="display: inline-block;vertical-align:top; ", uiOutput("effboxesso"))
                 ),
                 #div(style="display: inline-block;vertical-align:top; ", uiOutput("ai")),
                 div(style="", uiOutput("adhere")),
                 checkboxInput("adapt", "Adapt randomization probabilities", value = TRUE),
                 div(style="display: inline-block;",numericInput("burnin", "Burn-in sample size:",
                                                                 value = 30)),
                 
                 bsTooltip("burnin", "This is the number of patients that need to be enrolled into the trial before adaptation begins.",
                           "right", options = list(container = "body")),
                 div(style="display: inline-block;",numericInput("batchsize", "Sample size between two interim looks:",
                                                                 value = 20)),
                 bsTooltip("batchsize", "This is the number of patients enrolled into the trial between two adaptations. If this number is small adaptations are performed more frequently resulting in more computation and increasing simulation time.",
                           "right", options = list(container = "body")),
                 div(style="",sliderInput("upthresh",
                                          "Terminate trial for superiority at:", step = .005,
                                          min = 0,
                                          max = 1,
                                          value = 0.975)),
                 bsTooltip("upthresh", "If the probability that a treatments is better than the rest exceeds this threshold, at any point after burn-in and before the maximum sample size is reached, the trial will be terminated.",
                           "right", options = list(container = "body")),
                 div(style="",sliderInput("lowthresh", "Stop arm for futility at:",
                                          min = 0, max = 1, value = 0.01)),
                 bsTooltip("lowthresh", "If the probability that a treatment is better than the rest, for any of the treatment arms, falls below this threshold, the arm will be stopped (i.e. no more patients are assined to it) .",
                           "right", options = list(container = "body")),
                 div(style="display: inline-block;",numericInput("max", "Maximum total sample size:",
                                                                 value = 500)),
                 bsTooltip("max", "If the total number of patients enrolled into the trial reaches this number, the trial will be terminated.",
                           "right", options = list(container = "body")),
                 div(style="display: inline-block;",numericInput("ec", "Enrollment cost per patient ($):",
                                                                 value = 0))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Single trial simulation",
                            tabsetPanel(
                              tabPanel("Data",
                                       actionButton("button", "Run"), bsAlert('alert0'), plotOutput("dataPlot"),
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
                 ,
                 tabsetPanel(
                   tabPanel("Trial design properties",
                            tabsetPanel(
                              tabPanel("Power",
                                       numericInput("M", "Number of simulations", value = 100),
                                       bsTooltip("M", "Number of Monte Carlo simulations to estimate the power.",
                                                 "right", options = list(container = "body")),
                                       checkboxInput("tlimitp", "Set a limit for run time"),
                                       conditionalPanel(
                                         condition = "input.tlimitp == true",
                                         numericInput("Tpower", "Maximum permitted run time in seconds", value = NULL),
                                         bsTooltip("Tpower", "The simulation is terminated if run time exceeds this limit.",
                                                   "right", options = list(container = "body"))),
                                       
                                       actionButton("button0", "Compute"),
                                       
                                       wellPanel(conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                  tags$div("Computing power ...",id="loadmessage"))), htmlOutput('power')),
                              tabPanel("Sample size distribution/cost evaluation", plotOutput("ssdist"), DT::dataTableOutput("sssummary"),
                                       htmlOutput('sssave'), plotOutput("cdist"), DT::dataTableOutput("csummary"),
                                       htmlOutput('csave')),
                              tabPanel("Type I error rate",
                                       
                                       wellPanel(tags$div(tags$em(
                                         HTML(paste(tags$span(style="color:darkred", "Note: Changing the effect sizes and adherence rates does not affect the type I error rate as it is computed under the null hypothesis (i.e., none of the treatments are effective).")))
                                       ))),
                                       numericInput("Malpha", "Number of simulations", value = 100),
                                       bsTooltip("Malpha", "Number of Monte Carlo simulations to estimate the type I error rate.",
                                                 "right", options = list(container = "body")),
                                       checkboxInput("tlimita", "Set a limit for run time"),
                                       conditionalPanel(
                                         condition = "input.tlimita == true",
                                         numericInput("Talpha", "Maximum permitted run time in seconds", value = NULL),
                                         bsTooltip("Talpha", "The simulation is terminated if run time exceeds this limit.",
                                                   "right", options = list(container = "body"))),
                                       
                                       actionButton("button01", "Compute"),
                                       wellPanel(
                                         conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                          tags$div("Computing type I error rate ...",id="loadmessage"))), htmlOutput('alpha')),
                              tabPanel("Saved simulation configurations and results",
                                       
                                       actionButton("submit", "Save recent results"),
                                       bsTooltip("submit", "This button saves results to the application shared data folder. To see the most recent results in the table, reload after saving. The application data folder is cleared on a daily basis and can be modified by other users. To save the table contents to your machine use the download button below.",
                                                 "right", options = list(container = "body")),
                                       actionButton("load", "Load saved results"),
                                       bsTooltip("load", "Load simulation results from the application shared data folder.",
                                                 "right", options = list(container = "body")),
                                       selectInput('column', 'Columns', choices = names0, selected = names0[c(1, 13:20)],
                                                   multiple = T),
                                       actionButton("deleteRows", "Delete rows"), DT::dataTableOutput('responses'), tabPanel("Files list", DT::dataTableOutput("tbl")),
                                       downloadButton("downloadData", "Download"),
                                       bsTooltip("downloadData", "Use this button to download simulation data.",
                                                 "right", options = list(container = "body")),
                                       actionButton("clearData", "Clear data"),
                                       bsTooltip("clearData", "By clicking this button all the data in the application shared data folder is cleared.",
                                                 "right", options = list(container = "body"))
                                       #tags$div(tags$label("Upload", class="btn btn-primary",
                                       #                    tags$input(id = "fileIn", webkitdirectory = TRUE, type = "file", style="display: none;", onchange="pressed()")))
                                       #, shinyDirButton('directory', 'Folder select', 'Please select a folder')
                              )
                            )),
                   tabPanel("Compare with RCT",
                            actionButton("compRCT", "Generate comparison table"), wellPanel(
                              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                               tags$div("Just a moment, please! Fetching results ...",id="loadmessage"))), 
                            DT::dataTableOutput('RCT'))
                 )
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
                 actionButton("button1", "Calculate"),
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
