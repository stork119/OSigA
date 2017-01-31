### ###
### parallel_computing_initialise
### ###
source("R/parallel_computing.R")

path.optimisation <- paste(path.output, "cmaes/mvn/2017-01-31/", sep = "/")
dir.create(path.optimisation, recursive = TRUE)

parameters.filename <- paste(path.optimisation, "parameters_conditions.csv", sep = "")
if(file.exists(parameters.filename)){
parameters.conditions <- read.table(
  file = paste(path.optimisation, "parameters_conditions.csv", sep = ""),
  sep = ",",
  header = TRUE)
parameters.base <- parameters.conditions$base
parameters.factor <- parameters.conditions$factor
par.lower <- parameters.conditions$lower
par.upper <- parameters.conditions$upper
} else {
parameters.factor <- par.def
parameters.base <- rep(x = 10, times = length(parameters.factor))
par.lower <- rep(x = -2, times = length(parameters.factor))
par.upper <- rep(x = 2, times = length(parameters.factor))
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
}

lhs.res <- randomLHS(10000, length(par.def))

write.table(x = lhs.res,
            file = paste(path.optimisation, "parameters_list.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

stimulation.list <- unique(data.exp$stimulation)[c(4,7)]

write.table(x = matrix(stimulation.list, ncol = 1),
            file = paste(path.optimisation, "stimulation_list.txt", sep ="/"),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)


write.table(x = data.table(maxit = 1000,
                           fun.optimisation.likelihood = "fun.likelihood.mvn.sd_const",
                           fun_run_model = "run_model_mean"),
            file = paste(path.optimisation, "optimisation_conditions.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)


#### ####
data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time) %>% filter(stimulation %in% stimulation.list)
data.exp.grouped <- get_equal_data(data.exp.grouped)
data.exp.grouped <- data.exp.grouped %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))
write.table(x = data.exp.grouped,
            file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)