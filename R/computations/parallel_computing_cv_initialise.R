### ###
### parallel_computing_initialise
### ###
#setwd("~/Documents/modelling/")
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

#### optimisation.conditions ####

optimisation.conditions <- 
  list(
    maxit = 250,
    data.size = 1000,
    data.opt.size = 250,
    cross_validation_num = 25,
    optimisation = "optimisation_ut",
    fun.optimisation.likelihood = "fun.likelihood.list.sd",
    fun_run_model = "run_model_ut"
  )
optimisation.filename <- paste(path.list$optimisation, "optimisation_conditions.csv", sep = "")
if(!file.exists(optimisation.filename)){
  write.table(x = data.table(
    maxit = optimisation.conditions$maxit,
    data.size = optimisation.conditions$data.size,
    data.opt.size = optimisation.conditions$data.opt.size,
    cross_validation_num = 
      optimisation.conditions$cross_validation_num,
    fun.optimisation.likelihood =
      optimisation.conditions$fun.optimisation.likelihood,
    fun_run_model = optimisation.conditions$fun_run_model),
              file = optimisation.filename,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
} else {
  optimisation.conditions <- LoadOptimisationConditions(optimisation.filename)
}

#### create computation data ####
data.list$data.exp.norm <- 
  get_equal_data(
    data = data.list$data.exp,
    sample_size = optimisation.conditions$data.size)

path <- paste(path.list$optimisation.conditions, 1, sep = "/")
dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
write.table(
  x = data.list$data.exp.norm,
  file = paste(path, "data_exp_grouped.csv", sep = "/"),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE)


write.table(
  x = data.list$data.exp.norm,
  file = paste(path.list$optimisation, "data_exp_grouped.csv", sep = "/"),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE)

data.opt.list <- list()
data.opt.list[[1]] <- data.list$data.exp.norm
for(i in 2:(optimisation.conditions$cross_validation_num+1)){
  data.opt.list[[i]] <- 
    get_equal_data(
      data = data.list$data.exp.norm,
      sample_size = optimisation.conditions$data.opt.size)
  path <- paste(path.list$optimisation.conditions, i, sep = "/")
  dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
  write.table(
    x = data.opt.list[[i]],
    file = paste(path, "data_exp_grouped.csv", sep = "/"),
    sep = ",",
    row.names = FALSE,
    col.names = TRUE)
}