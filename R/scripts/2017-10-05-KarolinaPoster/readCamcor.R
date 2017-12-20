### ### 
### read CAMCROR
### ###


readCamcor <- function(
  input.path = poster.path.list$input.dir,
  experiment.path = experiment$path,
  ... 
){
  path.camcor <- paste(input.path,
                       experiment$path,
                       "data_camcor",
                       sep = "/")
  path.camcor.list <- list.files(path.camcor)
  path.camcor <- paste(path.camcor,
                       path.camcor.list[[1]],
                       sep = "/")
  path.camcor.list <- list.files(path.camcor)
  path.camcor <- paste(path.camcor,
                       path.camcor.list[[1]],
                       sep = "/")
  camcor <- list()
  camcor$camcor <-  read.table(file = 
                                 paste(path.camcor,
                                       "camcor.csv",
                                       sep = "/"),
                               header = TRUE,
                               sep = ",")
  camcor$camcor.images <-  read.table(file = 
                                 paste(path.camcor,
                                       "camcor-images.csv",
                                       sep = "/"),
                               header = TRUE,
                               sep = ",")
  return(camcor)
}