### ###
### parallel computing initialise set of computations
### ###
setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")
replace.cond <- FALSE
path.list <-
  LoadOptimisationPaths(
    path.output = "resources/output/",
    id = "2017-10-17-2"
  )
lhs.num <- 10
no_cores <- 33
#parameters.hypothesis <- c(6, 10,4, 2, 3, 7, 8, 9)
parameters.hypothesis <- c(1, 2, 4, 0)
conditions.computations.list <- 
  data.table(parameters = parameters.hypothesis, 
           state = "not started",
           started.time  = NA,
           finished.time = NA)
write.table(x = conditions.computations.list, 
            file = paste(path.list$optimisation, "conditions_computations_list.csv", sep = ""),
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

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
  optimisation.conditions <- LoadOptimisationConditions(optimisation.filename = optimisation.filename)
}
optimisation.conditions <- read.table(file = optimisation.filename, header = TRUE, sep = ",")
sigmapoints.conditions <- read.table(file = paste(path.list$optimisation, "sigmapoints_conditions.csv", sep = ""),
                                 header = TRUE, sep = ",")
sigmapoints.parameters.conditions <- read.table(file = paste(path.list$optimisation, "sigmapoints_parameters_conditions.csv", sep = ""),
                                 header = TRUE, sep = ",")
variables <- read.table(file = paste(path.list$optimisation, "variables.csv", sep = ""),
                                     header = FALSE, sep = ",")
variables.priming <- read.table(file = paste(path.list$optimisation, "variables-priming.csv", sep = ""),
                        header = TRUE, sep = ",")
registerDoParallel(no_cores)
foreach(par = parameters.hypothesis) %dopar% {

  path.list.par <-
    LoadOptimisationPaths(
      path.output = "resources/output/",
      id = paste(path.list$id, par, sep = "/") 
  )
  parameters.conditions.par <- parameters.conditions
  if(par != 0){
    parameters.conditions.par <- rbind(parameters.conditions.par,
                                       parameters.conditions.par[par,] %>% 
                                         dplyr::mutate(parameters = 0))
    parameters.conditions.par[par,] <- parameters.conditions.par[par,] %>% 
      dplyr::mutate(parameters.priming = 0)
  }
  
  write.table(x = parameters.conditions.par,
              file = paste(path.list.par$optimisation, "parameters_conditions.csv", sep = ""),
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  
  par.optimised   <- which(parameters.conditions.par$lower != parameters.conditions.par$upper)
  lhs.res <- randomLHS(lhs.num, length(par.optimised))

  if(replace.cond | !file.exists(paste(path.list.par$optimisation, "parameters_list.csv", sep = "/"))){
    write.table(x = rbind(matrix(0.5 + 0*par.optimised, nrow = 1), lhs.res),
              file = paste(path.list.par$optimisation, "parameters_list.csv", sep = ""),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  }
  stimulation.list <- (data.list$data.exp %>%
                       dplyr::filter(stimulation != 0) %>%
                       dplyr::distinct(stimulation))$stimulation

  write.table(x = matrix(stimulation.list, ncol = 1),
            file = paste(path.list.par$optimisation, "stimulation_list.txt", sep ="/"),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

  optimisation.filename <- paste(path.list.par$optimisation, "optimisation_conditions.csv", sep = "")
  if(!file.exists(optimisation.filename)){
    write.table(x = optimisation.table,
                file = optimisation.filename,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
  }
  write.table(x = sigmapoints.conditions,
              file = paste(path.list.par$optimisation, "sigmapoints_conditions.csv", sep = ""),
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)

  write.table(x = sigmapoints.parameters.conditions,
              file = paste(path.list.par$optimisation, "sigmapoints_parameters_conditions.csv", sep = ""),
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  write.table(x = variables.priming,
              file = paste(path.list.par$optimisation, "variables-priming.csv", sep = ""),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  write.table(x = variables,
              file = paste(path.list.par$optimisation, "variables.csv", sep = ""),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)

  path <- paste(path.list.par$optimisation.conditions, 1, sep = "/")
  if(replace.cond | !file.exists(paste(path, "data_exp_grouped.csv", sep = "/"))){
    data.list$data.exp.norm <- 
      get_equal_data(
        data = data.list$data.exp,
        sample_size = optimisation.conditions$data.size)
    dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
    write.table(
      x = data.list$data.exp.norm,
      file = paste(path, "data_exp_grouped.csv", sep = "/"),
      sep = ",",
      row.names = FALSE,
      col.names = TRUE)
    
    write.table(
      x = data.list$data.exp.norm,
      file = paste(path.list.par$optimisation, "data_exp_grouped.csv", sep = "/"),
      sep = ",",
      row.names = FALSE,
      col.names = TRUE)
  }
  data.opt.list <- list()
  data.opt.list[[1]] <- data.list$data.exp.norm
  for(i in 2:(optimisation.conditions$cross_validation_num+1)){
    path <- paste(path.list.par$optimisation.conditions, i, sep = "/")
    if(replace.cond | !file.exists(paste(path, "data_exp_grouped.csv", sep = "/"))){
      data.opt.list[[i]] <- 
        get_equal_data(
          data = data.list$data.exp.norm,
          sample_size = optimisation.conditions$data.opt.size)
      dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
      write.table(
        x = data.opt.list[[i]],
        file = paste(path, "data_exp_grouped.csv", sep = "/"),
        sep = ",",
        row.names = FALSE,
        col.names = TRUE)
    }
  }
  return()
}
stopImplicitCluster()