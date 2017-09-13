### ###
### Unscented transformation model
### optimisation procedure
### ###

#### optimisation ####
optimisation_ut <- function(par,
                            fun_run_model = run_model_ut,
                            parameters.conditions,
                            parameters.base,
                            parameters.factor,
                            variables, 
                            variables.priming, 
                            tmesh, 
                            tmesh.list,
                            stimulation.list,
                            background,
                            data.exp.grouped,
                            data.exp.summarise,
                            return.model = FALSE,
                            fun.likelihood,
                            par.optimised = rep(1, times = length(par)),
                            sigmapoints,
                            fun_parameters_penalty = NULL,
                            ...)
{
  
  flog.debug("optimisation_ut.R optimisation_ut", name="logger.optimisation")
  
  model.simulation <- do.call(fun_run_model,
                              list(par = par,
                                    parameters.base = parameters.base,
                                    parameters.factor = parameters.factor,
                                    variables = variables, 
                                    variables.priming = variables.priming, 
                                    tmesh = tmesh, 
                                    tmesh.list = tmesh.list,
                                    stimulation.list = stimulation.list,
                                    background = background,
                                    par.optimised = par.optimised,
                                    sigmapoints = sigmapoints, 
                                    parameters.conditions = parameters.conditions,
                                    ...))
  
  if(model.simulation$error){
    return(Inf)
  }
  
  result <- sum(
    likelihood(fun.likelihood = fun.likelihood,
               data.model = model.simulation$data.model,
               data.exp.grouped = data.exp.grouped,
               data.exp.summarise = data.exp.summarise))
  if(!is.null(fun_parameters_penalty)){
    result <- result + nrow(data.exp.summarise)*fun_parameters_penalty(par, ...)
  }
  
  if(return.model){
    return(list( error = FALSE, 
                 #optimisation = result,
                 data.model = model.simulation$data.model, 
                 data.model.sigmapoints = model.simulation$data.model.sigmapoints,
                 data.model.ut = model.simulation$data.model.ut,
                 arguments.list = model.simulation$arguments.list
    ))
  }
  return(result)
}
