#Hdir <- getwd()
outputDir <- "Data"

names0 <- c('nTreatment', paste('effect', 1:5, sep = ''), 'effectType',
            paste('adherence', 1:5, sep = ''), 'burninSampleSize', 'interimSampleSize', 'supThreshold',
            'futThreshold', 'maxSampleSize', 'power', 'typeIerror', 'ESS', 'minSS', 'firstQuartileSS',
            'medianSS', 'thirdQuartileSS', 'maxSS', 'ESavedSS')
fields <- c('nt', paste('eff', 1:5, sep = ''), 'efftype', paste('adh', 1:5, sep = ''), 'burnin',
           'batchsize', 'upthresh', 'lowthresh', 'max')


saveData <- function(data) {
  data <- t(data)
  # Create a unique file name
  fileName <- sprintf("%s_%s.csv", as.integer(Sys.time()), digest::digest(data))
  # Write the file to the local system
  write.csv(
    x = data,
    file = file.path(outputDir, fileName),
    row.names = FALSE, quote = TRUE
  )
}

loadData <- function(column) {
  # Read all the files into a list
  files <- list.files(outputDir, full.names = TRUE)
  data <- lapply(files, read.csv, stringsAsFactors = FALSE)
  # Concatenate all data together into one data.frame
  data <- do.call(rbind, data)
  data[,column]
}

clearData <- function() {
  files <- list.files(outputDir, full.names = TRUE)
  unlink(files)
}
