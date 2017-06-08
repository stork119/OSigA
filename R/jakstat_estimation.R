### ### ###
### jakstat_estimation
### ### ###

#### model  ####

run_model <- function(parameters, 
                      variables,
                      variables.priming,
                      tmesh,
                      tmesh.list,
                      stimulation.list,
                      background,
                      time_interval = 100,
                      time_computation = 1000*60*5){
  data.model <- data.table(time = numeric(),
                           m = numeric(),
                           sd = numeric(),
                           priming = numeric(), 
                           stimulation = numeric())
  # print(parameters)
  for(stm in stimulation.list){
    res <- rmain(parameters = parameters, 
                 variables = variables, 
                 stm = stm, 
                 tmesh = tmesh, 
                 time_interval = time_interval, 
                 time_computation = time_computation)
    res.priming <- rmain(parameters = parameters, 
                         variables = variables.priming, 
                         stm = stm, 
                         tmesh = tmesh, 
                         time_interval = time_interval, 
                         time_computation = time_computation)
    if(res$success & res.priming$success){
      for(tmesh.i in tmesh.list){
        data.model <- rbind(data.model,
                            data.table(time = c(tmesh[tmesh.i], tmesh[tmesh.i]),
                                       m = c(res$output[[tmesh.i]][14],
                                             res.priming$output[[tmesh.i]][14]),
                                       sd =  c(res$output[[tmesh.i]][31],
                                               res.priming$output[[tmesh.i]][31]),
                                       priming = c(0, 1000), 
                                       stimulation = c(stm, stm))
        )
      }
    } else {
      return(list(error = TRUE))
    }
  }
  data.model <- normalization(data.model, background = background)
  data.model <- lmvn(data.model)
  return(list(error = FALSE, data.model = data.model))
}


run_model_mean <- function(parameters, 
                           variables,
                           variables.priming,
                           tmesh,
                           tmesh.list,
                           stimulation.list,
                           background,
                           time_interval = 100,
                           time_computation = 1000*60*5,
                           ...){
  
  data.model <- data.table(time = numeric(),
                           m = numeric(),
                           sd = numeric(),
                           priming = numeric(), 
                           stimulation = numeric())
  # print(parameters)
  for(stm in stimulation.list){
    res <- rmainmean(parameters = parameters, 
                     variables = variables, 
                     stm = stm, 
                     tmesh = tmesh, 
                     time_interval = time_interval, 
                     time_computation = time_computation)
    res.priming <- rmainmean(parameters = parameters, 
                             variables = variables.priming, 
                             stm = stm, 
                             tmesh = tmesh, 
                             time_interval = time_interval, 
                             time_computation = time_computation)
    if(res$success & res.priming$success){
      for(tmesh.i in tmesh.list){
        data.model <- rbind(data.model,
                            data.table(time = c(tmesh[tmesh.i], tmesh[tmesh.i]),
                                       m = c(res$output[[tmesh.i]][14],
                                             res.priming$output[[tmesh.i]][14]),
                                       sd =  c(0,0),
                                       priming = c(0, 1000), 
                                       stimulation = c(stm, stm))
        )
      }
    } else {
      return(list(error = TRUE))
    }
  }
  data.model <- normalization(data.model, background = background)
  data.model <- lmvn(data.model)
  return(list(error = FALSE, data.model = data.model))
}
#### ####




#### ####

likelihood <- function(data.model,
                       data.exp.grouped,
                       data.exp.summarise,
                       fun.likelihood,
                       ...
){
    sapply(1:nrow(data.model),
           function(data.model.i){
             data.model.tmp <- data.model[data.model.i,]
             return((data.exp.grouped %>%
                       filter(priming == data.model.tmp$priming,
                              time == data.model.tmp$time,
                              stimulation == data.model.tmp$stimulation) %>%
                       mutate(likelihood = 
                                do.call(fun.likelihood,list(logintensity = logintensity, 
                                                            intensity = intensity, 
                                                            data.model.tmp = data.model.tmp,
                                                            data.exp.summarise = data.exp.summarise))) %>%
                       summarise(likelihood.sum = 
                                   sum(likelihood)))$likelihood.sum)
           }
    )
}

#### optimisation ####
optimisation <- function(fun_run_model = run_model,
                         par,
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
                         ...
                         
){
  parameters <- parameters.factor
  parameters[par.optimised] <- parameters.factor[par.optimised]*(parameters.base[par.optimised])^par
  model.simulation <- do.call(fun_run_model,
                              list(
                                parameters = parameters,
                                variables = variables,
                                variables.priming = variables.priming,
                                tmesh = tmesh,
                                tmesh.list = tmesh.list,
                                stimulation.list = stimulation.list,
                                background = background))
  
  
  if(model.simulation$error){
    return(Inf)
  }
  result <- sum(
    likelihood(fun.likelihood = fun.likelihood,
      data.model = model.simulation$data.model, 
      data.exp.grouped = data.exp.grouped,
      data.exp.summarise = data.exp.summarise))
  # print(result)
  if(return.model){
    return(c(model.simulation, list(optimisation = result)))
  }
  return(result)
}