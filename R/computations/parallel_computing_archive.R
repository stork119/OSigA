### ###
### parallel computing ###
### ###


#### ####
no_cores <- 6

#### ####
maxit <- 1000

#### ####
fun.optimisation.likelihood <- fun.likelihood.list[[1]]
fun_run_model <- run_model_mean
fun.optimisation <- pureCMAES
path.optimisation <- paste(path.output, "cmaes/mvn/2017-01-27/", sep = "/")
dir.create(path.optimisation, recursive = TRUE)

ids <- list.dirs(path.analysis, full.names = FALSE)
ids <- as.numeric(ids[which(!is.na(as.numeric(ids)))])
previous.computations <- max(ids)

#### ####
registerDoParallel(no_cores)
test <- foreach(i = 1:length(par.list), .combine = list, .multicombine = TRUE ) %dopar%
{
  par <- par.list[[i]]
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
         data.exp.grouped = data.exp.grouped,
         fun.likelihood = fun.optimisation.likelihood))
  
  parameters <- parameters.factor*(parameters.base)^optimisation.res$xmin
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
                                list(parameters = optimisation.res$xmin,
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
                       data.exp.grouped = data.exp.grouped))
                   })
  print(parameters)
  print(result)
  
  path.optimisation.i <- paste(path.optimisation, i + previous.computations, sep = "/")
  dir.create(path.optimisation.i, recursive = TRUE, showWarnings = FALSE)
  
  save_results(path.opt = path.optimisation.i,
               data.model.opt = model.simulation$data.model,
               par.opt = parameters,
               optimisation.opt = result,
               res.list = data.model.list,
               data.exp.grouped = data.exp.grouped,
               error = error,
               variables = variables,
               variables.priming = variables.priming)
  
  return(list(par = optimisation.res$xmin))
}
stopImplicitCluster()
