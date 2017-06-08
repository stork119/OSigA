###
### parallel_computing_run
###
setwd("~/Documents/modelling/")
source("R/parallel_computing.R")

# path.optimisation <- paste(path.output, "cmaes/normalized/2017-02-04-5/", sep = "/")

# data.exp.grouped <-  read.table(
#   file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
#   sep = ",",
#   header = TRUE)

path.optimisation <- paste(path.output, "optimisation/2017-06-08-2/", sep = "/")
data.model.list <- LoadInitialModels(path.optimisation = path.optimisation)
no_cores <- 16

run_parallel_computations(path.optimisation = path.optimisation,
                          # data.exp.grouped = data.exp.grouped,
                          no_cores = no_cores,
                          maxit.tmp   =  Inf,
                          # fun.optimisation = pureCMAES,
                          # optimisation.res.par = "xmin",
                          data.model.list = data.model.list)
