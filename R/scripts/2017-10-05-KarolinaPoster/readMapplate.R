### ###
### readMapplate
### ###

readMapplate <- function(
  input.path = poster.path.list$input.dir,
  experiment.path = experiment$path,
  ...
){

  path.mapplate <-  paste(input.path,
                          experiment$path,
                          "metadata",
                          #          "ffc/data_quantify/",
                          #          experiment$file,
                          sep = "/")
  
  mapplate.files <- list.files(
    path = path.mapplate,
    include.dirs = FALSE,
    full.names = FALSE,
    pattern = ".csv")
  
  mapplate <- list()
  foreach(mapplate.file = strsplit(mapplate.files, ".csv")) %do% {
    if(mapplate.file[[1]] == "args_ind"){
      mapplate.input.name  <- mapplate.file[[1]]
      mapplate.output.name <- "args_tag"
    } else if(mapplate.file[[1]] == "args_names"){
      mapplate.input.name  <- mapplate.file[[1]]
      mapplate.output.name <- "args_id"
    } else {
      mapplate.input.name  <- mapplate.file[[1]]
      mapplate.output.name <- mapplate.file[[1]]
    }
    mapplate[[mapplate.output.name]] <-
      read.table(file = 
                   paste(path.mapplate,
                         paste(mapplate.input.name, "csv", sep = "."),
                         sep = "/"),
                 header = FALSE,
                 sep = ",")
    if(ncol(mapplate[[mapplate.output.name]]) == 1){
      mapplate[[mapplate.output.name]] <-
        read.table(file = 
                     paste(path.mapplate,
                           paste(mapplate.input.name, "csv", sep = "."),
                           sep = "/"),
                   header = FALSE,
                   sep = "\t")
    }
  }
  
  mapplate.dirs <- list.dirs(
    path = path.mapplate,
    recursive = FALSE,
    full.names = FALSE)
  
  foreach(mapplate.dir = mapplate.dirs) %do% {
    mapplate.files.list <- strsplit(list.files(paste(path.mapplate, mapplate.dir, sep = "/")), ".csv")
    foreach(mapplate.file = mapplate.files.list) %do% {
      mapplate[[paste(mapplate.dir, mapplate.file[[1]], sep = ".")]] <-
        read.table(file = 
                     paste(path.mapplate,
                           mapplate.dir,
                           paste(mapplate.file[[1]], "csv", sep = "."),
                           sep = "/"),
                   header = FALSE,
                   sep = ",")
      if(ncol(mapplate[[paste(mapplate.dir, mapplate.file[[1]], sep = ".")]]) == 1){
        mapplate[[paste(mapplate.dir, mapplate.file[[1]], sep = ".")]] <-
          read.table(file = 
                       paste(path.mapplate,
                             mapplate.dir,
                             paste(mapplate.file[[1]], "csv", sep = "."),
                             sep = "/"),
                     header = FALSE,
                     sep = "\t")
      }
    }
  }
  return(mapplate)
}