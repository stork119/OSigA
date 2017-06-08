### ###
### model_execution
### ###
source("R/jakstat_model.R")
#### ####
normalization_simulation <- function(data.model,
                          m.scale = 400,
                          sd.scale = m.scale^2,
                          background,
                          epsilon = 1){
  data.model$m.norm <- data.model$m/m.scale + background
  data.model$sd.norm <- data.model$sd/sd.scale
  data.model$sd.norm <- sapply(data.model$sd.norm, function(sd.norm){ifelse(sd.norm < epsilon, epsilon, sd.norm)})
  data.model$time <- data.model$time
  return(data.model)
}

simulate_model <- function(parameters, 
                           variables,
                           variables.priming,
                           tmesh,
                           tmesh.list,
                           stimulation.list,
                           background,
                           time_interval = 100,
                           time_computation = 1000*60*5,
                           ...){
  
  data.model <- data.table(
    time = numeric(),
    m = numeric(),
    sd = numeric(),
    priming = numeric(), 
    stimulation = numeric()
  )
  data.trajectory <- data.table(
    time = numeric(),
    var = numeric(),
    m = numeric(),
    priming = numeric(), 
    stimulation = numeric()
  )
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
        
        data.trajectory <- rbind(data.trajectory,
                            data.table(
                              time = rep(x = tmesh[tmesh.i],
                                         times= 2*length(res$output[[tmesh.i]])),
                              m = c(res$output[[tmesh.i]],
                                    res.priming$output[[tmesh.i]]),
                              priming = c(rep(x = 0, times= length(res$output[[tmesh.i]])),
                                          rep(x = 1000, times= length(res$output[[tmesh.i]]))), 
                              stimulation = rep(x = stm, times= 2*length(res$output[[tmesh.i]])),
                              var  = rep(x = 1:length(res$output[[tmesh.i]]), times = 2) 
                            )
        )
      }
    } else {
      return(list(error = TRUE))
    }
  }
  data.model <- normalization_simulation(data.model, background = background)
  data.model <- lmvn(data.model)
  return(list(error = FALSE, data.model = data.model, data.trajectory = data.trajectory))
}
