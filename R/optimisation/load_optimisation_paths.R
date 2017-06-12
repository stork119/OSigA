### ###
### load optimisation directories
### ###

LoadOptimisationPaths <- function(
  path.output,
  id
){
  path.list <- list() 
  path.list$optimisation          <- paste(path.output, "optimisation", id, "/", sep = "/")
  path.list$optimisation.data     <- paste(path.list$optimisation, "data/", sep = "/")
  path.list$optimisation.results  <- paste(path.list$optimisation, "results/", sep = "/")
  path.list$optimisation.analysis <- paste(path.list$optimisation, "analysis/", sep = "/")
  
  sapply(path.list, dir.create, recursive = TRUE, showWarnings = FALSE)   
  
  return(path.list)
}
