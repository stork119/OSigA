### ###
### 2017-07-11-sigmapoints-script
### ###

setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")
source("R/sum_data_initialise.R")

gplot.list <- list()
attach(LoadOptimisationConditions(path.optimisation = path.list$optimisation,
                                  path.optimisation.data = path.list$optimisation.data))


sigmapoints <- LoadSigmapointsConditions()