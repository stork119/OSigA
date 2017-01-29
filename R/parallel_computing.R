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
                                      data.exp.grouped,
                                      no_cores = 1,
                                      maxit.tmp    =  Inf,
                                      fun.optimisation = pureCMAES,
                                      optimisation.res.par = "xmin"){
                                      
                                      
  ### initialization ###
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

  stimulation.list <- scan(paste(path.optimisation, "stimulation_list.txt", sep ="/"))
  data.exp.grouped.optimisation <- data.exp.grouped %>% filter(stimulation %in% stimulation.list)

  lhs.res <- read.table(file = paste(path.optimisation, "parameters_list.csv", sep = ""),
                      sep = ",",
                      header = FALSE)
par.list <- lapply(1:nrow(lhs.res), function(i){(par.upper - par.lower)*lhs.res[i,] + par.lower})


ids <- list.dirs(path.optimisation, full.names = FALSE)
ids <- as.numeric(ids[which(!is.na(as.numeric(ids)))])
previous.computations <- max(c(ids,0)) + 1
#### ####
registerDoParallel(no_cores)
test <- foreach(i = previous.computations:length(par.list), .combine = list, .multicombine = TRUE ) %dopar%
{
  par <- as.numeric(par.list[[i]])
  optimisation.res <- do.call(
    fun.optimisation,
    list(par = par, 
        fun = optimisation,
        lower = par.lower,
        upper = par.upper,
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
        fun.likelihood = fun.optimisation.likelihood))
  
  parameters <- parameters.factor*(parameters.base)^(optimisation.res[[optimisation.res.par]])
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
                                list(parameters = optimisation.res[[optimisation.res.par]],
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
  
  path.optimisation.i <- paste(path.optimisation, i, sep = "/")
  dir.create(path.optimisation.i, recursive = TRUE, showWarnings = FALSE)
  
  save_results(path.opt = path.optimisation.i,
               data.model.opt = model.simulation$data.model,
               par.opt = parameters,
               par.exp.opt = optimisation.res[[optimisation.res.par]],
               optimisation.opt = result,
               res.list = data.model.list,
               data.exp.grouped = data.exp.grouped.optimisation,
               error = error,
               variables = variables,
               variables.priming = variables.priming,
               grid.ncol = length(stimulation.list))
  
  return(list(par = optimisation.res[[optimisation.res.par]]))
}
stopImplicitCluster()
}

#### run computations ####