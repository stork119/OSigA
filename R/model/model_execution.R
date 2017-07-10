### ###
### model_execution
### ###

#### ####
normalization_trajectory <- function(data.trajectory,
                                     variables.indices,
                                     m.scale = 400,
                                     sd.scale = m.scale^2,
                                     background,
                                     epsilon = 1){
  
  data.trajectory.means <- data.trajectory %>% 
    dplyr::filter(var %in%  variables.indices$means)
  data.trajectory.vars  <- data.trajectory %>% 
    dplyr::filter(var %in%  variables.indices$vars$ind) %>% 
    dplyr::left_join(variables.indices$vars, by = c("var" = "ind"))
  data.trajectory.cov  <- data.trajectory %>% 
    dplyr::filter(var %in%  variables.indices$cov$ind) %>% 
    dplyr::left_join(variables.indices$cov, by = c("var" = "ind"))
  
  data.trajectory.means$m.norm <- data.trajectory.means$m/m.scale + background
  data.trajectory.vars$m.norm  <- data.trajectory.vars$m/(m.scale^2)# - 2*background*data.trajectory.means$m/m.scale
  data.trajectory.cov$m.norm  <- data.trajectory.cov$m/(m.scale^2)# - 2*background*data.trajectory.means$m/m.scale
  
  data.trajectory.list <- list()
  data.trajectory.list$means <- data.trajectory.means
  data.trajectory.list$vars  <- data.trajectory.vars
  data.trajectory.list$cov   <- data.trajectory.cov
  
  return(data.trajectory.list)
}

normalization_simulation <- function(data.model,
                          m.scale = 400,
                          sd.scale = m.scale^2,
                          background,
                          epsilon = 1){
  data.model$m.norm <- data.model$m/m.scale + background
  data.model$sd.norm <- data.model$sd/(m.scale^2)
  #data.model$sd.norm <- sapply(data.model$sd.norm, function(sd.norm){ifelse(sd.norm < epsilon, epsilon, sd.norm)})
  data.model$time <- data.model$time
  return(data.model)
}

simulate_model <- function(fun_run_model = rmainmean,
                           parameters, 
                           parameters.priming = parameters,
                           variables,
                           variables.priming,
                           tmesh,
                           tmesh.list,
                           stimulation.list,
                           background,
                           time_interval = 100,
                           time_computation = 1000*60*5,
                           tmesh.list.tmp = NULL,
                           ...){
  
  if(is.null(tmesh.list.tmp)){
    tmesh.list.tmp <- tmesh.list 
  }
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
    res <- fun_run_model(parameters = parameters, 
                     variables = variables, 
                     stm = stm, 
                     tmesh = tmesh, 
                     time_interval = time_interval, 
                     time_computation = time_computation)
    res.priming <- fun_run_model(parameters = parameters.priming, 
                             variables = variables.priming, 
                             stm = stm, 
                             tmesh = tmesh, 
                             time_interval = time_interval, 
                             time_computation = time_computation)
    
    
    if(res$success & res.priming$success){
      for(tmesh.i in tmesh.list.tmp){
        data.model <- rbind(data.model,
                            data.table(time = c(tmesh[tmesh.i], tmesh[tmesh.i]),
                                       m = c(res$output[[tmesh.i]][14],
                                             res.priming$output[[tmesh.i]][14]),
                                       sd =  c(ifelse(length(res$output[[tmesh.i]]) > 17, res$output[[tmesh.i]][31], 0),
                                               ifelse(length(res.priming$output[[tmesh.i]]) > 17, res.priming$output[[tmesh.i]][31], 0)),
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
