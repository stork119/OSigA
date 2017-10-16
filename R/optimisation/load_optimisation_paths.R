### ###
### load optimisation directories
### ###

LoadOptimisationPaths <- function(
  path.output,
  id,
  id.list = NULL
){
  path.list <- list() 
  
  if(!is.null(id.list)){
    id <- paste(id.list, sep = "/", collapse = "/")
    path.list$id  <- paste(id.list, sep = "-", collapse = "-")
  } else {
    path.list$id          <- id
  }
  path.list$optimisation          <- paste(path.output, "optimisation", id, "/", sep = "/")
  path.list$optimisation.data     <- paste(path.list$optimisation, "data/", sep = "/")
  path.list$optimisation.results  <- paste(path.list$optimisation, "results/", sep = "/")
  path.list$optimisation.analysis <- paste(path.list$optimisation, "analysis/", sep = "/")
  path.list$optimisation.conditions     <- paste(path.list$optimisation, "conditions/", sep = "/")
  
  sapply(path.list, dir.create, recursive = TRUE, showWarnings = FALSE)   
  
  return(path.list)
}
