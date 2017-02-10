### ### ###
### jakstat_estimation
### ### ###

#### LMVN ####
mean.lmvn <- function(m, sd){
  ifelse(m != 0.0, (log(m) - (1/2)*log((sd/(m^2)) + 1)), 0.0)
}

sd.lmvn <- function(m, sd){
  ifelse(m != 0.0, log((sd/(m^2)) + 1), 0.0)
}

lmvn <- function(data){
  data %>%
    mutate(mean.lmvn = mean.lmvn(m.norm, sd.norm),
           sd.lmvn = sd.lmvn(m.norm, sd.norm))
}

#### data normalization ####
normalization <- function(data.model,
                          m.scale = 400,
                          sd.scale = m.scale^2,
                          background,
                          epsilon = 1){
  data.model$m.norm <- data.model$m/m.scale + background
  data.model$sd.norm <- data.model$sd/sd.scale
  data.model$sd.norm <- sapply(data.model$sd.norm, function(sd.norm){ifelse(sd.norm < epsilon, epsilon, sd.norm)})
  return(data.model)
}

#### likelihood functions ####
fun.likelihood.lmvn <- function(logintensity, intensity, data.model.tmp, ...){
  return(log(sqrt(data.model.tmp$sd.lmvn)) +  ((logintensity - data.model.tmp$mean.lmvn)^2)/data.model.tmp$sd.lmvn)
}

fun.likelihood.mvn <- function(logintensity, intensity, data.model.tmp, ...){
  return(log(sqrt(data.model.tmp$sd.norm)) + (((intensity - data.model.tmp$m.norm)^2)/data.model.tmp$sd.norm))
}


fun.likelihood.mvn.mean <- function(logintensity, intensity, data.model.tmp, ...){
  return((intensity - data.model.tmp$m.norm)^2)
}

fun.likelihood.mvn.sd_const <- function(logintensity, intensity, data.model.tmp, intensity.sd, ...){
  return((((intensity - data.model.tmp$m.norm)^2)/intensity.sd))
}


fun.likelihood.lmvn.data <- function(logintensity, intensity, data.model.tmp, intensity.sd, ...){
  nu <- mean.lmvn(intensity, intensity.sd)
  sd <- sd.lmvn(intensity, intensity.sd)
  return((((nu - log(data.model.tmp$m.norm))^2)/sd))
}


fun.normalised <- function(logintensity, intensity, data.model.tmp, ...){
  return(((intensity - data.model.tmp$m.norm)^2)/intensity^2)
}

fun.normalised.by_priming <- function(logintensity, intensity, data.model.tmp, data, ...){
  return(((intensity - data.model.tmp$m.norm)^2)/ as.numeric( 
    data %>% 
      filter(priming == data.model.tmp$priming) %>%
      summarise(intensity = mean(intensity)))^2)
}

fun.likelihood.list <- list(fun.likelihood.mvn.mean,
                            fun.likelihood.mvn,
                            fun.likelihood.lmvn, 
                            fun.likelihood.mvn.sd_const,
                            fun.likelihood.lmvn.data,
                            fun.normalised,
                            fun.normalised.by_priming)

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
                           parameters.base,
                           parameters.index,
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

likelihood <- function(data.model,
                       data.exp.grouped,
                       fun.likelihood
){
  sapply(1:nrow(data.model),
         function(data.model.i){
           data.model.tmp <- data.model[data.model.i,]
           return((data.exp.grouped %>%
                     filter(priming == data.model.tmp$priming,
                            time == data.model.tmp$time,
                            stimulation == data.model.tmp$stimulation) %>%
                     mutate(likelihood = 
                              do.call(fun.likelihood,list(logintensity, 
                                                          intensity, 
                                                          data.model.tmp = data.model.tmp,
                                                          intensity_sd,
                                                          data = data.exp.grouped))) %>%
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
                         return.model = FALSE,
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
    likelihood(
      data.model = model.simulation$data.model, 
      data.exp.grouped = data.exp.grouped,
      ...))
  # print(result)
  if(return.model){
    return(c(model.simulation, list(optimisation = result)))
  }
  return(result)
}