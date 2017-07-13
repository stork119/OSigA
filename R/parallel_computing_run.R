###
### parallel_computing_run
###
setwd("~/Documents/modelling/")

# path.optimisation <- paste(path.output, "cmaes/normalized/2017-02-04-5/", sep = "/")

# data.exp.grouped <-  read.table(
#   file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
#   sep = ",",
#   header = TRUE)

source(file = "R/optimisation/initialise_optimisation.R")

# run_parallel_computations(path.optimisation = path.list$optimisation,
#                           # data.exp.grouped = data.exp.grouped,
#                           no_cores = 12,
#                           maxit.tmp   =  Inf,
#                           # fun.optimisation = pureCMAES,
#                           # optimisation.res.par = "xmin",
#                           data.model.list = list(),
#                             # LoadInitialModels(path.optimisation = path.list$optimisation,
#                             #                   path.optimisation.data = path.list$optimisation.data),
#                           fun_modify_input = function(...){PrepareModelArguments.pSTAT_extrinsic(priming_constant = 3.4, ...)}
#                           )


run_parallel_computations(path.optimisation = path.list$optimisation,
                          # data.exp.grouped = data.exp.grouped,
                          no_cores = 12,
#                          maxit.tmp   =   1,
                          # fun.optimisation = pureCMAES,
                          # optimisation.res.par = "xmin",
                          data.model.list = list(),
                          # LoadInitialModels(path.optimisation = path.list$optimisation,
                          #                   path.optimisation.data = path.list$optimisation.data),
                          sigmapoints = LoadSigmapointsConditions(path.optimisation = path.list$optimisation),
                          fun_modify_input = PrepareModelArguments.ut,
                          optimisation_procedure = optimisation_ut,
 #                         par.list.ids.part = 1
)
