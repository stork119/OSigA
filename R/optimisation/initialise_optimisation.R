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

source("R/computations/parallel_computing_cv.R")
source("R/computations/load_conditions_cv.R")
source("R/computations/computations_summary_library.R")

#### ####
path.list <- 
    LoadOptimisationPaths(
      path.output = "resources/output/",
      id = "2017-10-17-2"
    )
  

# path.optimisation <- path.list$optimisation
# path.optimisation.data <- path.list$optimisation.data
# path.optimisation.analysis <- path.list$optimisation.analysis
# path.optimisation.results <- path.list$optimisation.results