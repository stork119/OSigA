###
### parallel_computing_run
###
#setwd("~/Documents/Projects/Modelling/modelling/")
setwd("~/Documents//modelling/")

source(file = "R/optimisation/initialise_optimisation.R")


path.list <-
  LoadOptimisationPaths(
    path.output = "resources/output/",
    id = "2017-10-17-2"
  )

conditions.computations.list <- 
  read.table(file = paste(path.list$optimisation, "conditions_computations_list.csv", sep = ""),
            header = TRUE,
            sep = ",",
            stringsAsFactors = FALSE)

parameters.list <- (conditions.computations.list%>%
  dplyr::filter(state != "finished") %>%
  dplyr::arrange(state))$parameters

for( par in parameters.list[2]){
  conditions.computations.list[which(conditions.computations.list$parameters == par),]$state <- "computations"
  path.list.par <-
    LoadOptimisationPaths(
      path.output = paste("resources/output/", sep = "/"),
      id.list = list(path.list$id, par)
    )
  write.table(x = conditions.computations.list, 
              file = paste(path.list$optimisation, "conditions_computations_list.csv", sep = ""),
              col.names = TRUE,
              row.names = FALSE,
              sep = ",")
  run_parallel_computations_cv(path.list = path.list.par,
                               # data.exp.grouped = data.exp.grouped,
                               fun_model_ode = rmainmean,
                               no_cores = 18,
                               stopfitness = -56000,
                               # no_cores= 1,
                               #  maxit.tmp   =   1,
                               par.list.ids.part = 1:10,
                               #  computations.ids = 1,
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
  
  fun.computations$run_summary(path.list = path.list.par, data.list = data.list)
  conditions.computations.list[which(conditions.computations.list$parameters == par),]$state <- "finished"
  write.table(x = conditions.computations.list, 
              file = paste(path.list$optimisation, "conditions_computations_list.csv", sep = ""),
              col.names = TRUE,
              row.names = FALSE,
              sep = ",")
}
