
library(abind)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(DT)
library(R.utils)
source('sim_funs.R')
source('plots.R')
source('tables.R')
source('dsl.R')

t_now <- Sys.time()
load('Time/time_last')
if (difftime(t_now, t0, units = 'days') > 1) {
  clearData()
  t0 = Sys.time()
  save(t0, file = 'Time/time_last')
}

# The palette with grey:
cbPalette <<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <<- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

shinyServer(function(input, output, session) {

  # volumes <- c(Home = '~')
  # shinyFileChoose(input, 'file', roots=volumes, session=session)
  # shinyDirChoose(input, 'directory', roots=volumes, session=session, filetypes = c('', 'txt'))
  # shinyFileSave(input, 'save', roots=volumes, session=session)
  # output$filepaths <- renderPrint({parseFilePaths(volumes, input$file)})
  # output$directorypath <- renderPrint({parseDirPath(volumes, input$directory)})
  # output$savefile <- renderPrint({parseSavePath(volumes, input$save)})
  #outputDir <<- parseDirPath(volumes, input$directory)

  observe({
    x = as.integer(input$nt)
    updateSelectInput(session, "treatment",
                      label = "Treatment:",
                      choices = c(paste('treatment', 1:x))
    )
    
  })

  # Workaround for hover box not showing up properly for superiority bar given diffrerent options

  dataInput = eventReactive(input$button, {
    if(isValid_num0() | input$efftype == 'absolute'){
      closeAlert(session, 'Alert0')
      progress <- shiny::Progress$new()
      progress$set(message = "Simulating trial", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())

      # Create a callback function to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 50
        }
        progress$set(value = value, detail = detail)
      }
      eff0 = paste('eff', 1:input$nt, sep = '')
      eff = c()
      if (input$efftype == "absolute") {
        for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
      }
      if (input$efftype == "count") {
        for (i in 1:input$nt) eff[i] = log(as.numeric(input[[eff0[i]]]))
      }
      if (input$efftype == "rate") {
        for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
      }
      adh0 = paste('adh', 1:input$nt, sep = '')
      adh = c()
      for (i in 1:input$nt) adh[i] = as.numeric(input[[adh0[i]]])
      if (is.null(input$adapt)) adp = F else adp = input$adapt
      if (is.null(input$MID)) mid = 0 else mid = input$MID
      
      if (input$compCon == "TRUE") {
        upthreshold <- input$upthreshT
      } else {
        upthreshold <- input$upthreshF
      }
      
      RAR_sim(nt = input$nt, theta0 = eff, nb = input$batchsize, maxN = input$max, N = 1000,
              upper = upthreshold, uppfut = input$uppfut, lower = input$lowthresh,
              burn = input$burnin, response.type = input$efftype, conjugate_prior = T, padhere = adh, adapt = adp,
              platf = input$platf, compCon = input$compCon, MID = mid)
    } else {
      createAlert(session, "alert0", "Alert0", title = "Invalid Entry",
                  content = "Proportions must be between 0 and 1.", append = FALSE)

    }

  })

  SA = eventReactive(input$button, {
    trial = dataInput()
    if (input$so == T) {
      eff0 = paste('effso', 1:input$nt, sep = '')
      eff = c()
      if (input$efftypeso == "absolute") {
        for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
      }
      if (input$efftypeso == "rate") {
        for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
      }
      AE_sim(trial, nt = input$nt, theta_ae = eff, response.type = input$efftypeso)
    } else return()

  })
  
  power_alpha_calc = eventReactive(input$button2, {

    adh0 = paste('adh', 1:input$nt, sep = '')
    adh = c()
    for (i in 1:input$nt) adh[i] = as.numeric(input[[adh0[i]]])
    eff0 = paste('eff', 1:input$nt, sep = '')
    eff = c()
    for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
    if (is.null(input$adapt)) adp = F else adp = input$adapt
    if (is.null(input$MID)) mid = 0 else mid = input$MID
    
    if (input$compCon == "TRUE") {
      upthreshold <- input$upthreshT
    } else {
      upthreshold <- input$upthreshF
    }
    
    withProgress(message = 'Estimating power and type I error', value = 0, max = 2*input$M, {
    
    power_out <- withTimeout({power_compute(nt = input$nt, theta0 = eff, 
                                            nb = input$batchsize, maxN = input$max, N = 1000,
                               upper = upthreshold, uppfut = input$uppfut, lower = input$lowthresh,
                               burn = input$burnin, response.type = input$efftype, 
                               conjugate_prior = T, padhere = adh, adapt = adp,
                               platf = input$platf, compCon = input$compCon, MID = mid, M = input$M)},
                timeout = ifelse(!is.null(input$Tpower), input$Tpower, Inf), onTimeout = 'silent')
    
    alpha_out <- withTimeout({alpha_compute(nt = input$nt, theta0 = eff, 
                                            nb = input$batchsize, maxN = input$max, N = 1000,
                                            upper = upthreshold, uppfut = input$uppfut, lower = input$lowthresh,
                                            burn = input$burnin, response.type = input$efftype, 
                                            conjugate_prior = T, padhere = adh, adapt = adp,
                                            platf = input$platf, compCon = input$compCon, MID = mid,
                               M = input$M)},
                timeout = ifelse(!is.null(input$Tpower), input$Tpower, 1000000), onTimeout = 'silent')
    
    }) # END withProgress    
    
    return(list(power_out = power_out, 
                alpha_out = alpha_out))
    
  })
  
  power_calc_RCT = eventReactive(input$compRCT, {
    adh0 = paste('adh', 1:input$nt, sep = '')
    adh = c()
    for (i in 1:input$nt) adh[i] = as.numeric(input[[adh0[i]]])
    eff0 = paste('eff', 1:input$nt, sep = '')
    eff = c()
    for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
    
    if (input$compCon == "TRUE") {
      upthreshold <- input$upthreshT
    } else {
      upthreshold <- input$upthreshF
    }
    
    withTimeout({power_compute_RCT(nt = input$nt, theta0 = eff, maxN = input$max, N = 1000,
                               upper = upthreshold,
                               response.type = input$efftype, conjugate_prior = T,
                               padhere = adh , compCon = input$compCon, M = input$M)}, 
                timeout = ifelse(!is.null(input$Tpower), input$Tpower, Inf), onTimeout = 'silent')
    
  })
  
  alpha_calc_RCT = eventReactive(input$compRCT, {
    
    eff0 = paste('eff', 1:input$nt, sep = '')
    eff = c()
    for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
    
    if (input$compCon == "TRUE") {
      upthreshold <- input$upthreshT
    } else {
      upthreshold <- input$upthreshF
    }
    
    withTimeout({alpha_compute_RCT(nt = input$nt, theta0 = eff, maxN = input$max, N = 1000,
                               upper = upthreshold, 
                               response.type = input$efftype, conjugate_prior = T, compCon = input$compCon, 
                               M = input$M)}, 
                timeout = ifelse(!is.null(input$Tpower), input$Tpower, 1000000), onTimeout = 'silent')
  })


  isValid_num0 <- eventReactive(input$button, {
    eff0 = paste('eff', 1:input$nt, sep = '')
    eff = c()
    for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
    sum(eff>=0) == input$nt && sum(eff<=1) == input$nt && sum(!is.na(eff)) == input$nt
  })

  isValid_num <- eventReactive(input$button1,{
    eff0 = paste('effss', 1:input$ntss, sep = '')
    eff = c()
    for (i in 1:input$ntss) eff[i] = as.numeric(input[[eff0[i]]])
    v0 = paste('var', 1:input$ntss, sep = '')
    var0 = c()
    for (i in 1:input$ntss) var0[i] = as.numeric(input[[v0[i]]])
    sum(eff>=0) == input$ntss && sum(eff<=1) == input$ntss && sum(!is.na(eff)) == input$ntss && sum(var0>=0) == input$ntss
  })

  isValid_num1 <- eventReactive(input$button,{
    eff0 = paste('effso', 1:input$nt, sep = '')
    eff = c()
    for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
    sum(eff>=0) == input$nt && sum(eff<=1) == input$nt && sum(!is.na(eff)) == input$nt
  })

  isValid_na0 <- eventReactive(input$button, {
    eff0 = paste('eff', 1:input$nt, sep = '')
    eff = c()
    for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
    sum(!is.na(eff)) == input$nt
  })

  isValid_na <- eventReactive(input$button1,{
    eff0 = paste('effss', 1:input$ntss, sep = '')
    eff = c()
    for (i in 1:input$ntss) eff[i] = as.numeric(input[[eff0[i]]])
    v0 = paste('var', 1:input$ntss, sep = '')
    var0 = c()
    for (i in 1:input$ntss) var0[i] = as.numeric(input[[v0[i]]])
    sum(!is.na(var0)) == input$ntss && sum(!is.na(eff)) == input$ntss
  })

  isValid_na1 <- eventReactive(input$button,{
    eff0 = paste('effso', 1:input$nt, sep = '')
    eff = c()
    for (i in 1:input$nt) eff[i] = as.numeric(input[[eff0[i]]])
    sum(!is.na(eff)) == input$nt
  })

  sample_size = eventReactive(input$button1, {
    if(isValid_num() | (input$efftypess == 'absolute' && isValid_na())){
      eff0 = paste('effss', 1:input$ntss, sep = '')
      eff = c()
      for (i in 1:input$ntss) eff[i] = as.numeric(input[[eff0[i]]])
      v0 = paste('var', 1:input$ntss, sep = '')
      var0 = c()
      for (i in 1:input$ntss) var0[i] = as.numeric(input[[v0[i]]])

      sampsize(pars = eff, power = input$power, alpha = input$alpha, dropout = input$dropout, type = input$efftypess,
               sigmas = var0)
    }
  })

  output$effboxes <- renderUI({
    numT <- as.integer(input$nt)
    typ <- as.character(input$efftype)
    if (typ == "absolute") {
      typ  = "Effect size"
      step0 = 0.5
      min0 = -Inf
      max0 = Inf
      value0 = 0
    }
    if (typ == "rate") {
      typ  = "Proportion"
      step0 = .05
      min0 = 0
      max0 = 1
      value0 = 0.5
    }
    if (typ == "count") {
      typ  = "Count"
      step0 = 1
      min0 = 0
      max0 = Inf
      value0 = 0
    }
    
      cid = paste('eff', 1, sep = '')
      c1 = paste(typ, "for Treatment 1 (Reference)")
      list(div(style="display: inline-block; ", numericInput(cid, c1, step = step0 , min = min0, max = max0, value = value0, width = '100%')),
      lapply(1:(numT-1), function(i) {
        div(style="", numericInput(paste('eff', i + 1, sep = ''), paste(typ, "for Treatment", i+1),
                                                                             step = step0 , min = min0, max = max0, value = value0, width = '100%'))
        #
        # list(div(style="display: inline-block;vertical-align:top; ", numericInput(paste('eff', i + 1, sep = ''), paste(typ, "for Treatment", i),
        #              step = step0 , min = min0, max = max0, value = value0, width = '100%')),
        #      div(style="display: inline-block;vertical-align:top; ",
        #          checkboxInput(paste("aal", (i+1), sep = ''), "Add later")))
      }))
           
    # When "First arm is control" checkbox was still being used
    # if (input$con == T) {
    #   cid = paste('eff', 1, sep = '')
    #   c1 = paste(typ, "for Control")
    #   list(div(style="display: inline-block; ", numericInput(cid, c1, step = step0 , min = min0, max = max0, value = value0, width = '100%')),
    #   lapply(1:(numT-1), function(i) {
    #     div(style="", numericInput(paste('eff', i + 1, sep = ''), paste(typ, "for Treatment", i),
    #                                                                          step = step0 , min = min0, max = max0, value = value0, width = '100%'))
    #     # 
    #     # list(div(style="display: inline-block;vertical-align:top; ", numericInput(paste('eff', i + 1, sep = ''), paste(typ, "for Treatment", i),
    #     #              step = step0 , min = min0, max = max0, value = value0, width = '100%')),
    #     #      div(style="display: inline-block;vertical-align:top; ", 
    #     #          checkboxInput(paste("aal", (i+1), sep = ''), "Add later")))
    #   }))
    # } else {
    #   lapply(1:numT, function(i) {
    #     div(style="", numericInput(paste('eff', i, sep = ''), paste(typ, "for Treatment", i),
    #                                                                          step = step0 , min = min0, max = max0, value = value0, width = '100%'))
    #     # 
    #     # list(div(style="display: inline-block;vertical-align:top; ", numericInput(paste('eff', i + 1, sep = ''), paste(typ, "for Treatment", i),
    #     #                                                                           step = step0 , min = min0, max = max0, value = value0, width = '100%')),
    #     #      div(style="display: inline-block;vertical-align:top; ", 
    #     #          checkboxInput(paste("aal", i, sep = ''), "Add later"))
    #     #      )
    #   })
    # }
  })
  
  # output$ai <- renderUI({
  #   numT <- as.integer(input$nt)
  #   aal0 <- paste("aal", 1:numT, sep = '')
  #   if (input$con == T) {
  #   lapply(X = 2:numT, FUN = function(i){
  #           if (input[[aal0[i]]]) 
  #             numericInput(paste('ai', i,sep = ''), 'at interim', value = 1)
  #   })
  #   } else {
  #     lapply(X = 1:numT, FUN = function(i){
  #       if (input[[aal0[i]]]) 
  #         numericInput(paste('ai', i,sep = ''), 'at interim', value = 1)
  #   }
  #     )
  #   }
  # })
  
  
  # observeEvent(input[[paste("aal", 2, sep = '')]], {
  #   insertUI(
  #     selector = "#aal",
  #     where = "afterEnd",
  #     ui = numericInput('ai', 'at interim', value = 1)
  #   )
  # })

  output$effboxesso <- renderUI({
    numT <- as.integer(input$nt)
    typ <- as.character(input$efftypeso)
    if (typ == "absolute") {
      typ  = "Effect size"
      step0 = 0.5
      min0 = -Inf
      max0 = Inf
      value0 = 0
    }
    if (typ == "rate") {
      typ  = "Proportion"
      step0 = 0.05
      min0 = 0
      max0 = 1
      value0 = 0.5
    }
    if (typ == "count") {
      typ  = "Count"
      step0 = 1
      min0 = 0
      max0 = Inf
      value0 = 0
    }
    lapply(1:numT, function(i) {
      numericInput(paste('effso', i, sep = ''), paste(typ, "for Treatment", i),
                   step = step0, min = min0, max = max0, value = value0, width = '100%')
    })
  })



  output$adhere <- renderUI({
    numT <- as.integer(input$nt)
    lapply(1:numT, function(i) {
      sliderInput(paste('adh', i, sep = ''), paste("Probability of adherence to Treatment", i),
                  min = 0, max = 1, value = 1)
    })
  })


  output$effboxesss <- renderUI({
    numT <- as.integer(input$ntss)
    typ <- as.character(input$efftypess)
    if (typ == "absolute") {
      typ  = "Effect size"
      step0 = 0.5
      min0 = -Inf
      max0 = Inf
      value0 = 0

    }
    if (typ == "rate") {
      typ  = "Proportion"
      step0 = .05
      min0 = 0
      max0 = 1
      value0 = 0.5
    }
    lapply(1:numT, function(i) {
      numericInput(paste('effss', i, sep = ''), paste(typ, "for Treatment", i),
                   step = step0, min = min0, max = max0, value = value0)
    })
  })


  output$varboxesss <- renderUI({
    numT <- as.integer(input$ntss)
    if (as.character(input$efftypess) == "absolute") {
      lapply(1:numT, function(i) {
        numericInput(paste('var', i, sep = ''), paste("Variance", "for Treatment", i),
                     step = .5, min = 0, max = Inf, value = 1)
      })
    }
  })

  output$ss <- renderUI({
    if(isValid_num() | (input$efftypess == 'absolute' && isValid_na())){
      closeAlert(session, 'Alertss')
      ss = sample_size()
      HTML(paste('<br/>'), paste("<pre>",'Sample size per arm:', "<font color=\"#FF0000\"><b>", ss$n, "</b></font>"),
           paste('<br/>'), paste('Total sample size:', "<font color=\"#FF0000\"><b>", ss$N, "</b></font>"))
    } else {
      createAlert(session, "alertss", "Alertss", title = "Invalid Entry",
                  content = "Proportions must be between 0 and 1.", append = FALSE)
    }
  })

  output$supPlot <- renderPlot({
    trial = dataInput()
    if (input$compCon == "TRUE") {
      upthreshold <- input$upthreshT
    } else {
      upthreshold <- input$upthreshF
    }
    psup_plot(trial, upthreshold)
  })
  
  output$designPlot <- renderPlot({
    if(isValid_num0() | (input$efftype == 'absolute' && isValid_na0())){
      closeAlert(session, 'Alert0')
      trial = dataInput()
      designPlot(trial)
    } else {
      createAlert(session, "alert0", "Alert0", title = "Invalid Entry",
                  content = "Proportions must be between 0 and 1.", append = FALSE)
    }
  })

  output$dataPlot <- renderPlot({
    if(isValid_num0() | (input$efftype == 'absolute' && isValid_na0())){
      closeAlert(session, 'Alert0')
      trial = dataInput()
      data_plot(trial)
    } else {
      createAlert(session, "alert0", "Alert0", title = "Invalid Entry",
                  content = "Proportions must be between 0 and 1.", append = FALSE)
    }
  })

  output$secData <- renderPlot({
    if (input$so == T) {
      closeAlert(session, 'Alert1')
      if(isValid_num1() | (input$efftypesso == 'absolute' && isValid_na1())){
        closeAlert(session, 'Alertso')
        trial = SA()
        data_plot(trial)
      } else {
        createAlert(session, "alertso", "Alertso", title = "Invalid Entry",
                    content = "Proportions must be between 0 and 1.", append = FALSE)
      }
    } else {
      closeAlert(session, 'Alert1')
      createAlert(session, "alert1", "Alert1", title = "Oops",
                  content = "You did not include a secondary outcome.", append = FALSE)
    }

  })

  output$postPlot <- renderPlot({
    if(isValid_num0() | (input$efftype == 'absolute' & isValid_na0())){
      closeAlert(session, 'Alert0')
      trial = dataInput()
      post_plot(trial, type = input$efftype)
    }

  })

  output$estPlot <- renderPlot({
    if(isValid_num0() | (input$efftype == 'absolute' & isValid_na0())){
      closeAlert(session, 'Alert0')
      trial = dataInput()
      estPlot(trial)
    }
  })

  output$secEst <- renderPlot({
    if(isValid_num1() | (input$efftypeso == 'absolute' & isValid_na1())){
      closeAlert(session, 'Alertso')
      trial = SA()
      if (input$so == T){
        estPlot(trial)
      } else return()
    }
  })
  
  
  output$power <- DT::renderDataTable({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    
    if (is.null(p0)) {
      d0 <- data.frame(power = "Calculation terminated: run time exceeded the maximum permitted time by user.")
      colnames(d0) = "Error:"
      row.names(d0) = ""
      return(datatable(d0, rownames = T, escape = FALSE,
                options = list(bLengthChange=0,                       # show/hide records per page dropdown
                               bFilter=0,                             # global search box on/off
                               bInfo=0,
                               bPaginate=0,
                               bSort=0)))
        }
    
    if (!is.null(p0)) {
      pow0 = c()
      for (i in 1:length(p0$power)) {
        if (p0$power[i]>=0.8) 
          pow0[i] = paste("<font color=\"#1b7c45\"><b>",round(p0$power[i], 2), "</b></font>") else
            pow0[i] = paste("<font color=\"#FF0000\"><b>",round(p0$power[i], 2), "</b></font>")
      }
      d0 = t(data.frame(power = pow0))
      if (input$compCon == T) {
        colnames(d0) = paste('Treatment', 2:input$nt)
        row.names(d0) = "<b>Power</b>"
        datatable(d0, rownames = T, escape = FALSE,
                  caption = paste('Power to detect the effect of all treatments against the reference treatment:'),
                  options = list(bLengthChange=0,                       # show/hide records per page dropdown
                                 bFilter=0,                             # global search box on/off
                                 bInfo=0,
                                 bPaginate=0,
                                 bSort=0))
      } else {
        colnames(d0) = paste(' ')
        row.names(d0) = "<b>Power</b>"
        datatable(d0, rownames = T, escape = FALSE,
                  caption = paste('Power to detect the superior treatment:'), 
                  options = list(bLengthChange=0,                       # show/hide records per page dropdown
                                 bFilter=0,                             # global search box on/off
                                 bInfo=0,
                                 bPaginate=0,
                                 bSort=0))
      }
      
    }
    
  })
  
  
  output$alpha <- DT::renderDataTable({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$alpha_out #power_calc()
    
    if (is.null(p0)) {
      d0 <- data.frame(power = "Calculation terminated: run time exceeded the maximum permitted time by user.")
      colnames(d0) = "Error:"
      row.names(d0) = ""
      return(datatable(d0, rownames = T, escape = FALSE,
                       options = list(bLengthChange=0,                       # show/hide records per page dropdown
                                      bFilter=0,                             # global search box on/off
                                      bInfo=0,
                                      bPaginate=0,
                                      bSort=0)))
    }
    
    if (!is.null(p0)) {
      a0 = c()
      for (i in 1:length(p0$alpha)) {
        if (p0$alpha[i]<=0.05) 
          a0[i] = paste("<font color=\"#1b7c45\"><b>",round(p0$alpha[i], 2), "</b></font>") else
            a0[i] = paste("<font color=\"#FF0000\"><b>",round(p0$alpha[i], 2), "</b></font>")
      }
      d0 = t(data.frame(alpha = a0))
      if (input$compCon == T) {
        colnames(d0) = paste('Treatment', 2:input$nt)
        row.names(d0) = "<b>Type I error rate</b>"
        datatable(d0, rownames = T, escape = FALSE,
                  caption = paste('Probability of detecting a false positive effect for each treatment against the reference treatment:'),
                  options = list(bLengthChange=0,                       # show/hide records per page dropdown
                                 bFilter=0,                             # global search box on/off
                                 bInfo=0,
                                 bPaginate=0,
                                 bSort=0))
      } else {
        colnames(d0) = paste(' ')
        row.names(d0) = "<b>Type I error rate</b>"
        datatable(d0, rownames = T, escape = FALSE,
                  caption = paste('Probability of false positive:'),
                  options = list(bLengthChange=0,                       # show/hide records per page dropdown
                                 bFilter=0,                             # global search box on/off
                                 bInfo=0,
                                 bPaginate=0,
                                 bSort=0))
      }
      
    }
  })

  # output$alpha <- renderUI({
  #   power_alpha0 <- power_alpha_calc()
  #   p0 = power_alpha0$alpha_out #alpha_calc()
  #   if (is.null(p0)) {
  #     HTML(paste('<br/>'), paste("<pre>","<font color=\"#FF0000\"><b>",'Power calculation terminated: run time exceeded the maximum permitted time by user.', "</b></font>"))
  #   } else {
  #     if (round(p0$alpha, 2) <= 0.05) {
  #       HTML(paste('<br/>'), paste("<pre>",'Estimated type I error rate:', "<font color=\"#1b7c45\"><b>",round(p0$alpha, 2), "</b></font>"))
  #     }
  #     else HTML(paste('<br/>'), paste("<pre>",'Estimated type I error rate:', "<font color=\"#FF0000\"><b>",round(p0$alpha, 2), "</b></font>"))
  #   }
  # 
  # })

  output$ssdist <- renderPlot({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    
    if (is.null(p0)) HTML(paste('<br/>'), paste("<pre>","<font color=\"#FF0000\"><b>",'Power calculation terminated: run time exceeded the maximum permitted time by user.', "</b></font>"))
    else {
      d0 = as.data.frame(p0$Nt)
      names(d0) = 'sample.size'
      ggplot(d0, aes(x = sample.size)) + geom_histogram(fill = '#33FFFF', alpha = .5, color = 'grey') +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              strip.text.y = element_text(size = 12, face = "bold"))
    }
    
  })

  output$sssummary <- DT::renderDataTable({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    
    if (!is.null(p0)) {
      d0 = t(data.frame(quantile(p0$Nt, c(0,.25,.5,.75,1))))
      row.names(d0) = 'sample.size'
      colnames(d0) = paste(colnames(d0), 'quantile')
      datatable(d0, rownames = T)
      }

  })

  output$sssave <- renderUI({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    
    if (!is.null(p0)) {
      d00 = mean(input$max - p0$Nt)
      HTML(paste('<br/>'), paste("<pre>",'Expected saved sample size (difference between maximum sample size and sample size at trial termination) ', "<font color=\"#FF0000\"><b>", round(d00), "</b></font>"))
    }
  })
  
  output$cdist <- renderPlot({
    
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    
    if (is.null(p0)) HTML(paste('<br/>'), paste("<pre>","<font color=\"#FF0000\"><b>",'Power calculation terminated: run time exceeded the maximum permitted time by user.', "</b></font>"))
    else {
      d0 = as.data.frame(p0$Nt*input$ec)
      names(d0) = 'cost'
      ggplot(d0, aes(x = cost)) + geom_histogram(fill = 'purple', alpha = .5, color = 'grey') +
        xlab('cost($)') +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              strip.text.y = element_text(size = 12, face = "bold"))
    }
    
  })
  
  output$csummary <- DT::renderDataTable({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()

    if (!is.null(p0)) {
      d0 = t(data.frame(quantile(p0$Nt * input$ec, c(0,.25,.5,.75,1))))
      row.names(d0) = 'cost'
      colnames(d0) = paste(colnames(d0), 'quantile')
      datatable(d0, rownames = T)
    }
    
  })
  
  output$csave <- renderUI({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    
    if (!is.null(p0)) {
      d00 = mean(input$max * input$ec - p0$Nt * input$ec)
      HTML(paste('<br/>'), paste("<pre>",'Expected savings in dollar (difference between maximum cost and expected cost of an adaptive trial) ', "<font color=\"#FF0000\"><b>", d00, "</b></font>"))
    }
  })

  output$estable <- DT::renderDataTable({
    if(isValid_num0() | (input$efftype == 'absolute' & isValid_na0())){
      closeAlert(session, 'Alert0')
      trial = dataInput()
      datatable(est(trial, 0.95, type = input$efftype), rownames = T)
    }
  })

  output$sec.estable <- DT::renderDataTable({
    if(isValid_num1() | (input$efftypeso == 'absolute' & isValid_na1())){
      closeAlert(session, 'Alertso')
      trial = SA()
      if (input$so == T) {
        datatable(est(trial, 0.95, type = input$efftypeso), rownames = T)
      } else return()
    }
  })

  output$psuptable <- DT::renderDataTable({
    if(isValid_num0() | (input$efftype == 'absolute' & isValid_na0())){
      closeAlert(session, 'Alert0')
      trial = dataInput()
      datatable(psup(trial), rownames = T)
    }
  })

  output$datatable <- DT::renderDataTable({
    if(isValid_num0() | (input$efftype == 'absolute' & isValid_na0())){
      closeAlert(session, 'Alert0')
      trial = dataInput()
      datatable(datat(trial), rownames = T)
    }
  })


  formData <- reactive({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    p00 = power_alpha0$alpha_out #alpha_calc()
    n0 = length(p0$power)
    pow = rep(NA, 4)
    pow[1:n0] = p0$power
    alpha = ifelse(!is.null(p00), p00$alpha, NA)
    n1 = length(alpha)
    alph = rep(NA, 4)
    alph[1:n1] = alpha
    if (!is.null(p0)) {
      d0 = round(quantile(p0$Nt, c(0,.25,.5,.75,1)))
      d00 = round(mean(input$max - p0$Nt))
    } else {
      d0 = NA
      d00 = NA
    }
    data <- c(sapply(fields, function(x) ifelse(is.null(input[[x]]), NA, input[[x]])), c(pow, alph, round(mean(p0$Nt)), round(d0), round(d00)))
    names(data) = names0
    data
  })

  # When the Submit button is clicked, save the form data
  observeEvent(input$submit, {
    saveData(formData())
  })

  values <- reactiveValues(df = NULL, dd = NULL, xx = NULL)
  df <- observeEvent(input$load, {
    values$df <- loadData(input$column)
  })

  observeEvent(input$deleteRows,{
    if (!is.null(input$responses_rows_selected)) {
      values$df <- values$df[-as.numeric(input$responses_rows_selected),]
    }
  })

  # Show the previous responses
  # (update with current response when Submit is clicked)
  output$responses <- DT::renderDataTable({
    input$submit
    values$df
    #loadData(input$column)
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("sim_data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(values$df, file)
    }
  )

  observeEvent(input$clearData,{
    clearData()
  })

  output$test <- renderText({
    eff0 = paste('effss', 1:input$ntss, sep = '')
    eff = c()
    for (i in 1:input$ntss) eff[i] = as.numeric(input[[eff0[i]]])
    v0 = paste('var', 1:input$ntss, sep = '')
    var0 = c()
    for (i in 1:input$ntss) var0[i] = as.numeric(input[[v0[i]]])
    eff
  })
  
  output$BRATvsBRCTplot <- renderPlot({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    a0 = power_alpha0$alpha_out #alpha_calc()
    if(!is.null(a0)) al0 = a0$alpha else al0 =  NULL
    if(!is.null(a0)) pow0 = p0$power else pow0 =  NULL
    d0 = ifelse(!is.null(p0), mean(p0$Nt * input$ec), NULL)
    al1 = alpha_calc_RCT()
    pow1 = power_calc_RCT()
    d1 = input$max * input$ec
    BRATvsBRCTPlot(a0 = al0, a1 = al1,
                   p0 = pow0, p1 = pow1, 
                   c0 = d0, c1 = d1)
    
  })
  
  output$test <- renderText({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    a0 = power_alpha0$alpha_out #alpha_calc()
    if(!is.null(a0)) al0 = a0$alpha else al0 =  NULL
    if(!is.null(a0)) pow0 = p0$power else pow0 =  NULL
    d0 = ifelse(!is.null(p0), mean(p0$Nt * input$ec), NULL)
    al1 = alpha_calc_RCT()
    pow1 = power_calc_RCT()
    d1 = input$max * input$ec
    pow0
    
  })
  
  
  
    
      #' @param a0 alpha for BRAT
      #' @param a1 alpha for BRCT
      #' @param p0 power for BRAT
      #' @param p1 power for BRCT
      #' @param c0 cost for BRAT
      #' @param c1 cost for BRCT
  
  output$RCT <- DT::renderDataTable({
    power_alpha0 <- power_alpha_calc()
    p0 = power_alpha0$power_out #power_calc()
    a0 = power_alpha0$alpha_out #alpha_calc()
    al0 = ifelse(!is.null(a0), a0$alpha, NULL)
    pow0 = ifelse(!is.null(p0), p0$power, NULL)
    d0 = ifelse(!is.null(p0), mean(p0$Nt * input$ec), NULL)
    al1 = alpha_calc_RCT()
    pow1 = power_calc_RCT()
    d1 = input$max * input$ec
    dt = data.frame(c(al0, al1), c(pow0, pow1), c(d0, d1))
    row.names(dt) = c('BRAT', 'BRCT')
    colnames(dt) = c('type I error rate', 'power', 'cost($)')
    datatable(dt, rownames = T)
  })

  output$about = renderUI({
    withMathJax(includeHTML('rmd0.html'))
  })


  addPopover(session, "ssdist", "Sample size", content = paste0("<p>Sample size at time of trial termination for the simulated trials.</p>"), trigger = 'click')
  addPopover(session, "ssdist", "Cost", content = paste0("<p>Cost distribution over the simulated trials.</p>"), trigger = 'click')
  addPopover(session, "dataPlot", "Data", content = paste0("<p>Primary outcome observations for every patient and each treatment arm.</p>"), trigger = 'click')
  addPopover(session, "estPlot", "Estimates", content = paste0("<p>Point estimates (posterior mean) and 95% credible intervals for the treatment effects on the primary outcome.</p>"), trigger = 'click')
  addPopover(session, "secData", "Data", content = paste0("<p>Secondary outcome observations for every patient and each treatment arm.</p>"), trigger = 'click')
  addPopover(session, "secEst", "Estimates", content = paste0("<p>Point estimates (posterior mean) and 95% credible intervals for the treatment effects on the secondary outcome.</p>"), trigger = 'click')
  addPopover(session, "postPlot", "Posterior", content = paste0("<p>Posterior density of treatment effects on the primary outcome at five equally spaced interim looks. As more data is collected the posterior distribution becomes more focused. If a treatment arm is stopped early due to futility, the corresponding posterior does not evolve further.</p>"), trigger = 'click')
  addPopover(session, "supPlot", "Probabilities of superiority", content = paste0("<p>Probabilities of superiority of each treatment in terms of the effect on primary outcome. As evidence of better performance is observed for a treatment arm the corresponding curve ascends. The grey horizonthal line shows the upper threshold that if reached by any of the curves the trial is stopped.</p>"), trigger = 'click')
})

