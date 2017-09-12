### ###
### initialise optimisation directoreis
### ###

#### sources ####
source("R/libraries.R")
source("R/sources.R")
source("R/initialise.R")

source("R/logging/logging-initnialisation.R")

filename.optimisation       <- "initialise_optimisation.R"
path.optimisation  <- "R/optimisation/"
files.optimisation <- list.files(path = path.optimisation,
                             pattern = ".R$", 
                             full.names = FALSE, 
                             recursive = TRUE)
files.optimisation <- files.optimisation[filename.optimisation != files.optimisation]

comment.optimisation <- sapply(
  paste(path.optimisation, files.optimisation, sep = "/"),
  source,
  echo = FALSE,
  verbose = FALSE)

rm(filename.optimisation,path.optimisation,files.optimisation,comment.optimisation)

source("R/computations/parallel_computing.R")

#### ####
path.list <- list()
path.list$id <- "2017-09-08-density_ver0"
path.list <- 
  append(
    path.list,
    LoadOptimisationPaths(
      path.output = "resources/output/",
      id = path.list$id
    ))
  

path.optimisation <- path.list$optimisation
path.optimisation.data <- path.list$optimisation.data
path.optimisation.analysis <- path.list$optimisation.analysis
path.optimisation.results <- path.list$optimisation.results