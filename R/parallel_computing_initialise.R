### ###
### parallel_computing_initialise
### ###
source("R/parallel_computing.R")

path.optimisation <- paste(path.output, "cmaes/mvn/2017-01-28-4-stm/", sep = "/")
dir.create(path.optimisation, recursive = TRUE)

parameters.factor <- par.def
parameters.base <- rep(x = 10, times = length(parameters.factor))
par.lower <- rep(x = -2, times = length(parameters.factor))
par.upper <- rep(x = 2, times = length(parameters.factor))
lhs.res <- randomLHS(1000, length(par.def))

write.table(x = lhs.res,
            file = paste(path.optimisation, "parameters_list.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

parameters.conditions <- data.table(
  factor = par.def,
  base   = parameters.base,
  lower  = par.lower,
  upper  = par.upper)

write.table(x = parameters.conditions,
             file = paste(path.optimisation, "parameters_conditions.csv", sep = ""),
             sep = ",",
             row.names = FALSE,
             col.names = TRUE)

stimulation.list <- unique(data.exp.grouped$stimulation)[c(2,4,6,8)]

write.table(x = matrix(stimulation.list, ncol = 1),
            file = paste(path.optimisation, "stimulation_list.txt", sep ="/"),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)


write.table(x = data.table(maxit = 1000,
                           fun.optimisation.likelihood = "fun.likelihood.mvn.mean",
                           fun_run_model = "run_model_mean"),
            file = paste(path.optimisation, "optimisation_conditions.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
