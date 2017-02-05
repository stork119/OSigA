### ###
### parallel computing ###
### ###
source("R/libraries.R")
source("R/initialise.R")
#### ####
#### ####
#### ####

#### ####
run_parallel_computations <- function(path.optimisation,
                                      path.optimisation.data = paste(path.optimisation, "data", sep = "/"),
                                      data.exp.grouped,
                                      no_cores = 1,
                                      maxit.tmp    =  Inf,
                                      fun.optimisation = pureCMAES,
                                      optimisation.res.par = "xmin"){
                                      
                                      
  ### initialization ###
  dir.create(path.optimisation.data, showWarnings = FALSE, recursive = TRUE)
  optimisation.conditions <- read.table(
    file = paste(path.optimisation, "optimisation_conditions.csv", sep = ""),
    sep = ",",
    header = TRUE, stringsAsFactors = FALSE)
  fun.optimisation.likelihood <- get(optimisation.conditions$fun.optimisation.likelihood)
  fun_run_model <-  get(optimisation.conditions$fun_run_model)
  maxit <- min(optimisation.conditions$maxit, maxit.tmp)
  
  parameters.conditions <- read.table(
      file = paste(path.optimisation, "parameters_conditions.csv", sep = ""),
      sep = ",",
      header = TRUE)
  parameters.base <- parameters.conditions$base
  parameters.factor <- parameters.conditions$factor
  par.lower <- parameters.conditions$lower
  par.upper <- parameters.conditions$upper
  par.optimised   <- which(par.lower != par.upper)
  
  stimulation.list <- scan(paste(path.optimisation, "stimulation_list.txt", sep ="/"))
  data.exp.grouped.optimisation <- data.exp.grouped %>% filter(stimulation %in% stimulation.list)

  lhs.res <- read.table(file = paste(path.optimisation, "parameters_list.csv", sep = ""),
                      sep = ",",
                      header = FALSE)
  par.list <- lapply(1:nrow(lhs.res), function(i){(par.upper[par.optimised] - par.lower[par.optimised])*lhs.res[i,] + par.lower[par.optimised]})
  
  
  ids <- list.dirs(path.optimisation.data, full.names = FALSE)
  ids <- as.numeric(ids[which(!is.na(as.numeric(ids)))])
  par.list.ids <- 1:length(par.list)
  if(length(ids) > 0){
    par.list.ids <- par.list.ids[-which(ids %in% par.list.ids)]
  }
#### ####
registerDoParallel(no_cores)
test <- foreach(i = par.list.ids, .combine = list, .multicombine = TRUE ) %dopar%
{
  par <- as.numeric(par.list[[i]])
  optimisation.res <- do.call(
    fun.optimisation,
    list(par = par, 
        fun = optimisation,
        lower = par.lower[par.optimised],
        upper = par.upper[par.optimised],
        stopeval = maxit,
        fun_run_model = fun_run_model,
        variables = variables,
        variables.priming = variables.priming,
        parameters.base = parameters.base,
        parameters.factor = parameters.factor,
        tmesh = tmesh,
        tmesh.list = tmesh.list,
        stimulation.list = stimulation.list,
        background = background,
        data.exp.grouped = data.exp.grouped.optimisation,
        fun.likelihood = fun.optimisation.likelihood,
        par.optimised = par.optimised))
  
  parameters <- parameters.factor
  parameters[par.optimised] <- parameters.factor[par.optimised]*(parameters.base[par.optimised])^(optimisation.res[[optimisation.res.par]])
  par.exp.opt <- rep(0, times = length(parameters)) 
  par.exp.opt[par.optimised] <-  optimisation.res[[optimisation.res.par]]
  
  model.simulation <- do.call(run_model,
                              list(parameters = parameters,
                                variables = variables,
                                variables.priming = variables.priming,
                                tmesh = tmesh,
                                tmesh.list = tmesh.list,
                                stimulation.list = stimulation.list,
                                background = background))
  
  error <- model.simulation$error
  if(model.simulation$error){
    model.simulation <- do.call(run_model_mean,
                                list(parameters = parameters,
                                     variables = variables,
                                     variables.priming = variables.priming,
                                     tmesh = tmesh,
                                     tmesh.list = tmesh.list,
                                     stimulation.list = stimulation.list,
                                     background = background))
    
  }
  
  result <- sapply(fun.likelihood.list, 
                   function(fun.likelihood){
                     sum( likelihood( 
                       fun.likelihood = fun.likelihood,
                       data.model = model.simulation$data.model,
                       data.exp.grouped = data.exp.grouped.optimisation))
                     })
  print(parameters)
  print(result)
  
  path.optimisation.i <- paste(path.optimisation.data, i, sep = "/")
  dir.create(path.optimisation.i, recursive = TRUE, showWarnings = FALSE)
  
  save_results(path.opt = path.optimisation.i,
               data.model.opt = model.simulation$data.model,
               par.opt = parameters,
               par.exp.opt = par.exp.opt,
               optimisation.opt = result,
               res.list = data.model.list,
               data.exp.grouped = data.exp.grouped.optimisation,
               error = error,
               variables = variables,
               variables.priming = variables.priming,
               grid.ncol = ceiling(length(stimulation.list)/2))
  
  return(list(par = parameters))
}
stopImplicitCluster()
}

#### run computations ####