### ###
### experiments to RDS
### ###

source("R/scripts/2017-10-05-KarolinaPoster/readCamcor.R")
source("R/scripts/2017-10-05-KarolinaPoster/readMapplate.R")

#### reading data ####
poster.data.list <- list()
poster.path.list$experiments <- 
  list(
    list(id = "2016-01-26-KA10-C"),
    list(id = "2016-01-26-KA10-IFNB"),
    list(id = "2016-01-26-KA10-IFNG"),
    list(id = "2016-01-28-KA11-C"),
    list(id = "2016-01-28-KA11-IFNB"),
    list(id = "2016-01-28-KA11-IFNG"),
    list(id = "2017-07-18-KS25", path = "2017-07-18-KS25/488"),
    list(id = "2017-07-18-KS26", path = "2017-07-18-KS26/488"),
    list(id = "2017-07-25-KS27", path = "2017-07-25-KS27/488"),
    list(id = "2017-07-25-KS28", path = "2017-07-25-KS28/488"),
    list(id = "2017-07-19-KZ63", path = "2017-07-19-KZ63/488"),
    list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", output_id = "2017-08-25-KZ73-Nuclei"),
    list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", output_id = "2017-09-15-KZ78-Nuclei"),
    list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", file = "Cells555.csv", output_id = "2017-08-25-KZ73-Cells"),
    list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", file = "Cells555.csv", output_id = "2017-09-15-KZ78-Cells"),
    list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", file = "Cytoplasm555.csv", output_id = "2017-08-25-KZ73-Cytoplasm"),
    list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", file = "Cytoplasm555.csv", output_id = "2017-09-15-KZ78-Cytoplasm")
  )

#### ####
normalisation.list <- c("raw", "ffc", "bck", "pbs", "pbsbck")

#### ####
no_cores <- 6
registerDoParallel(no_cores)
foreach(experiment = poster.path.list$experiments) %dopar% {
  # if(is.null(experiment[["file"]])){
  #   experiment$file <- "ShrinkedNuclei.csv"
  # }
  if(is.null(experiment[["path"]])){
    experiment$path <- experiment$id
  }
  path <- paste(poster.path.list$input.dir,
                experiment$path,
      #          "ffc/data_quantify/",
      #          experiment$file,
                sep = "/")
  
  experiment$mapplate <- readMapplate(
    input.path = poster.path.list$input.dir,
    experiment.path = experiment$path
    )
  
  experiment$camcor <- readCamcor(
    input.path = poster.path.list$input.dir,
    experiment.path = experiment$path
  )
  
  experiment.path.list <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  experiment.path.list <- experiment.path.list[experiment.path.list %in% normalisation.list]
  
  f <- foreach(normalistion = experiment.path.list) %do% {
    experiment[[normalistion]] <- list()
    input.path <- paste(poster.path.list$input.dir, experiment$path, normalistion, "data_quantify", sep = "/")
    quantification.list <- list.files(input.path)
    f <- foreach(quantification = quantification.list) %do% {
      tryCatch({
        experiment[[normalistion]][[strsplit(quantification, ".csv")[[1]]]] <-
          read.table(
            file = paste(input.path,
                         quantification,
                         sep = "/"),
            header = TRUE,
            sep = ","
          ) %>% data.table()
        if(ncol(experiment[[normalistion]][[strsplit(quantification, ".csv")[[1]]]]) == 1){
          experiment[[normalistion]][[strsplit(quantification, ".csv")[[1]]]] <-
            read.table(file = 
                         paste(input.path,
                               quantification,
                               sep = "/"),
                       header = TRUE,
                       sep = "\t")  %>% 
            data.table()
        }
      }, error = function(e){ print(e) })
      return()
    }
    return()
  }
  output.path <- paste(poster.path.list$output.dir, 
                       "rds", 
                       sep = "/")
  dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(object = experiment, 
          file = paste(output.path, 
                       paste(experiment$id, ".RDS", sep = ""),
                       sep = "/"))
}
stopImplicitCluster()


#### ####
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