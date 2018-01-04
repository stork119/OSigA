###
### parallel_computing_run
###
#setwd("~/Documents/Projects/Modelling/modelling/")
setwd("~/Documents//modelling/")
source(file = "R/optimisation/initialise_optimisation.R")

# basicConfig()
# addHandler(writeToConsole)
# setLevel(loglevels[["DEBUG"]], getHandler('basic.stdout'))

run_parallel_computations_cv(path.list = path.list,
                          # data.exp.grouped = data.exp.grouped,
                          fun_model_ode = rmainmean,
                          no_cores = 16,
                          stopfitness = -56000,
                           # maxit.tmp   =   1,
                           # par.list.ids.part = 1,
                           # computations.ids = 1,
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

fun.computations$run_summary(path.list = path.list, data.list = data.list)

#### testing ####
# fun_model_ode = rmainmean
# no_cores = 1
# maxit.tmp   =   1
# par.list.ids.part = 1
# computations.ids = 1
# stopfitness = -10000000
# fun.optimisation = cma_es
# optimisation.res.par = "par"
# optimisation_procedure = optimisation_ut
# fun_modify_input      = PrepareModelArguments.ut.multiple
# fun_modify_parameters = PrepareModelParameters.ut
# optimisation_procedure = optimisation_ut
# fun_parameters_penalty =  NULL#
# sigmapoints = LoadSigmapointsConditions(path.optimisation = path.list$optimisation)