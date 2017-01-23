### ### ###
### jakstat_estimation
### ### ###

mean.lmvn <- function(m, sd){
  if(m != 0.0){
    return(log(m) - log(sd/(m^2)) +1)
  } else {
    return(0.0)
  }
}

sd.lmvn <- function(m, sd){
  if(m != 0.0){
    return(log((sd/(m^2)) + 1 ))
  } else {
    return(0.0)
  }
}


normalization <- function(data.model,
                          m.scale = 400,
                          sd.scale = 400,
                          background,
                          epsilon = 1){
  data.model$m.norm <- data.model$m/m.scale + background
  data.model$sd.norm <- data.model$sd/sd.scale
  data.model$sd.norm <- sapply(data.model$sd.norm, function(sd.norm){ifelse(sd.norm < epsilon, epsilon, sd.norm)})
  return(data.model)
}