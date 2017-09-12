### ### ###
### jakstat_estimation
### ### ###

#### optimisation ####
optimisation <- function(par,
                         fun_run_model = run_model,
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
                         fun_modify_input = 
                           function(parameters,
                                    parameters.priming = parameters,
                                    variables,
                                    variables.priming) {
                             return(list(parameters = parameters,
                                         parameters.priming = parameters.priming,
                                         variables = variables,
                                         variables.priming = variables.priming
                             ))},
                         ...)
{
  parameters <- parameters.factor
  parameters[par.optimised] <- parameters.factor[par.optimised]*(parameters.base[par.optimised])^par
  
  input <- fun_modify_input(parameters = parameters,
                            variables = variables,
                            variables.priming = variables.priming)
  parameters <- input$parameters
  variables  <- input$variables
  variables.priming <- input$variables.priming
  
  model.simulation <- do.call(fun_run_model,
                              list(
                                parameters = parameters,
                                variables = variables,
                                variables.priming = variables.priming,
                                tmesh = tmesh,
                                tmesh.list = tmesh.list,
                                stimulation.list = stimulation.list,
                                background = background,
                                parameters.conditions = parameters.conditions,
                                ...))
  
  
  if(model.simulation$error){
 #   print("kupa")
    return(Inf)
  }
#  print("hurrra")

  result <- sum(
    likelihood(fun.likelihood = fun.likelihood,
      data.model = model.simulation$data.model, 
      data.exp.grouped = data.exp.grouped,
      data.exp.summarise = data.exp.summarise))
  print(result)
  if(return.model){
    return(c(model.simulation, list(optimisation = result)))
  }
  return(result)
}