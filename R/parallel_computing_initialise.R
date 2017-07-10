### ###
### parallel_computing_initialise
### ###
setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")

parameters.filename <- paste(path.list$optimisation, "parameters_conditions.csv", sep = "")
if(file.exists(parameters.filename)){
parameters.conditions <- read.table(
  file = parameters.filename,
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
            file = parameters.filename,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
}
par.optimised   <- which(par.lower != par.upper)

lhs.res <- randomLHS(1000, length(par.optimised))

write.table(x = rbind(matrix(0.5 + 0*par.optimised, nrow = 1), lhs.res),
            file = paste(path.list$optimisation, "parameters_list.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

stimulation.list <- (data.list$data.exp %>%
                       dplyr::filter(stimulation != 0) %>%
                       dplyr::distinct(stimulation))$stimulation

write.table(x = matrix(stimulation.list, ncol = 1),
            file = paste(path.list$optimisation, "stimulation_list.txt", sep ="/"),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)


optimisation.filename <- paste(path.list$optimisation, "optimisation_conditions.csv", sep = "")
if(!file.exists(optimisation.filename)){
  write.table(x = data.table(maxit = 1000,
                             fun.optimisation.likelihood = "fun.likelihood.lmvn.data.norm",
                             fun_run_model = "run_model"),
              file = optimisation.filename,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
}

### prepare data ###
data.exp.grouped <- data.list$data.exp.norm
write.table(x = data.exp.grouped,
            file = paste(path.list$optimisation, "data_exp_grouped.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)