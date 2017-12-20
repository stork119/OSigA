###
### prepare data librarr ###
### ###

source("R/optimisation/initialise_optimisation.R")

poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"

poster.path.list$output.dir <- "resources/output/poster/"
dir.create(poster.path.list$output.dir)


#### getPathsList ####
getPathsList <- function(date = "2017-12-18", 
                         normalisation){

  poster.path.list <- list()
  poster.path.list$input.dir  <-  paste("resources/input/poster/")
  poster.path.list$input.rds.dir <- paste(poster.path.list$input.dir, date, "rds/", sep = "/")
  poster.path.list$rds.path <- paste(poster.path.list$input.rds.dir,
                                     paste("data_",
                                           normalisation, 
                                           ".RDS",
                                           sep = ""),
                                     sep = "/")
  poster.path.list$output.dir <- paste("resources/output/poster", date, normalisation, sep = "/")
  dir.create(path = poster.path.list$output.dir,
             showWarnings = FALSE,
             recursive = TRUE)
  return(poster.path.list)
}

#### getColumnNames ####
getColumnNames <- function(data){
  col_well <- "well.name"
  
  if("time" %in% colnames(data)){
    col_time <- "time"
  } else if ("time.2.1" %in% colnames(data)){
    col_time <- "time.2.1"
  } else {
    col_time <- "time.1.1"
  }
  
  if("stimulation" %in% colnames(data)){
    col_stimulation <- "stimulation"
  } else if ("stimulation.2.1" %in% colnames(data)){
    col_stimulation <- "stimulation.2.1"
  } else {
    col_stimulation <- "stimulation.1.1"
  }
  
  
  if("Intensity_MeanIntensity_Alexa" %in% colnames(data)){
    col_response <- "Intensity_MeanIntensity_Alexa"
  } else if("Intensity_MeanIntensity_Alexa488" %in% colnames(data)){
    col_response <- "Intensity_MeanIntensity_Alexa488"
  } else {
    col_response <- "Intensity_MeanIntensity_Alexa555"
  }
  
  col_priming <- "priming.1.1"
  
  return(
    list(well = col_well,
         time = col_time,
         stimulation = col_stimulation,
         response = col_response,
         priming = col_priming
    )
  )
}
