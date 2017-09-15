###
### parallel_computing_run
###
#setwd("~/Documents/Projects/Modelling/modelling/")
setwd("~/Documents//modelling/")

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


# run_parallel_computations(path.optimisation = path.list$optimisation,
#                           # data.exp.grouped = data.exp.grouped,
#                           no_cores = 12,
#                           #maxit.tmp   =   1,
#                           #par.list.ids.part = 1,
#                           # fun.optimisation = pureCMAES,
#                           # optimisation.res.par = "xmin",
#                           data.model.list = list(),
#                           # LoadInitialModels(path.optimisation = path.list$optimisation,
#                           #                   path.optimisation.data = path.list$optimisation.data),
#                           sigmapoints = LoadSigmapointsConditions(path.optimisation = path.list$optimisation),
#                           model.computations = list(raw = TRUE, priming = TRUE),
#                           fun_modify_input = PrepareModelArguments.ut,
#                           fun_modify_parameters = PrepareModelParameters.ut,
#                           optimisation_procedure = optimisation_ut
# )


# basicConfig()
# addHandler(writeToConsole)
# setLevel(loglevels[["DEBUG"]], getHandler('basic.stdout'))

# run_parallel_computations(path.list = path.list,
#                           # data.exp.grouped = data.exp.grouped,
#                           no_cores = 1,
#                           maxit.tmp   =   1,
#                           par.list.ids.part = 1,
#                           # fun.optimisation = pureCMAES,
#                           # optimisation.res.par = "xmin",
#                           data.model.list = list(),
#                           # LoadInitialModels(path.optimisation = path.list$optimisation,
#                           #                   path.optimisation.data = path.list$optimisation.data),
#                           sigmapoints = LoadSigmapointsConditions(path.optimisation = path.list$optimisation),
#                           model.computations = list(raw = TRUE, priming = TRUE),
#                           fun_modify_input      = PrepareModelArguments.ut.multiple,
#                           fun_modify_parameters = PrepareModelParameters.ut,
#                           optimisation_procedure = optimisation_ut,
#                           fun_parameters_penalty =  NULL#fun_parameters_penalty_sigmapoints
# )


run_parallel_computations_cv(path.list = path.list,
                          # data.exp.grouped = data.exp.grouped,
                          fun_model_ode = rmainmean,
                          no_cores = 16,
                          #maxit.tmp   =   1,
                          #par.list.ids.part = 1,
                          #computations.ids = 1,
                          # fun.optimisation = pureCMAES,
                          # optimisation.res.par = "xmin"warni,
                          data.model.list = list(),
                          # LoadInitialModels(path.optimisation = path.list$optimisation,
                          #                   path.optimisation.data = path.list$optimisation.data),
                          sigmapoints = LoadSigmapointsConditions(path.optimisation = path.list$optimisation),
                          model.computations = list(raw = TRUE, priming = TRUE),
                          fun_modify_input      = PrepareModelArguments.ut.multiple,
                          fun_modify_parameters = PrepareModelParameters.ut,
                          optimisation_procedure = optimisation_ut,
                          fun_parameters_penalty =  NULL#fun_parameters_penalty_sigmapoints
)

run_summary(path.list = path.list)