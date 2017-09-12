run_model <- function(parameters, 
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