### ###
### likelihood functions
### ###

# fun.likelihood.normal <- function(m,
#                                   sd, 
#                                   model.m,
#                                   model.sd,
#                                   ...){
#   
# }

ComputeLikelihood.lmvn <- function(m,
                                sd, 
                                X,
                                ... ){
  nu    <- lmvn.mean(m, sd)
  sigma <- lmvn.sd(m, sd)
  return((((nu - log(X))^2)/sigma))
}

ComputeLikelihood.lmvn.nusigma <- function(nu, # lmvn.mean(m, sd)
                                   sigma, # lmvn.sd(m, sd)
                                   X,
                                ... ){
  return((((nu - log(X))^2)/sigma))
}

ComputeLikelihood.lmvn.rse <- function(nu, # lmvn.mean(m, sd)
                                    sigma, # lmvn.sd(m, sd)
                                    X,
                                    ... ){
  return((((nu - log(X))^2)/nu^{2}))
}


ComputeLikelihood.lmvn.bias <- function(nu, # lmvn.mean(m, sd)
                                       sigma, # lmvn.sd(m, sd)
                                       X,
                                       ... ){
  return(abs(nu - log(X)))
}


fun.likelihood.list.sd_data <- function(logintensity = logintensity, 
                                        data.model.tmp = data.model.tmp,
                                        data.exp.summarise = data.exp.summarise,
                                        ...){
  intensity.sd <- (data.exp.summarise %>% dplyr::filter(
    stimulation == data.model.tmp$stimulation,
    priming == data.model.tmp$priming,
    time == data.model.tmp$time))$sd.norm
  
  nu <- mean.lmvn(data.model.tmp$m.norm, intensity.sd)
  sd <- sd.lmvn(data.model.tmp$m.norm, intensity.sd)
  return(((nu - logintensity)^2)/sd)
}

fun.likelihood.list <- list(
  sd_data = fun.likelihood.list.sd_data 
)

# ComputeLikelihood.lmvn.bias <- function(nu, # lmvn.mean(m, sd)
#                                         sigma, # lmvn.sd(m, sd)
#                                         X,
#                                         ... ){
#   return(abs(nu - log(X)))
# }

#### likelihood functions ####
# fun.likelihood.lmvn <- function(logintensity, intensity, data.model.tmp, ...){
#   return(log(sqrt(data.model.tmp$sd.lmvn)) +  ((logintensity - data.model.tmp$mean.lmvn)^2)/data.model.tmp$sd.lmvn)
# }
# 
# fun.likelihood.mvn <- function(logintensity, intensity, data.model.tmp, ...){
#   return(log(sqrt(data.model.tmp$sd.norm)) + (((intensity - data.model.tmp$m.norm)^2)/data.model.tmp$sd.norm))
# }
# 
# 
# fun.likelihood.mvn.mean <- function(logintensity, intensity, data.model.tmp, ...){
#   return((intensity - data.model.tmp$m.norm)^2)
# }
# 
# fun.likelihood.mvn.sd_const <- function(logintensity, intensity, data.model.tmp, intensity.sd, ...){
#   return((((intensity - data.model.tmp$m.norm)^2)/intensity.sd))
# }
# 
# 
# fun.likelihood.lmvn.data <- function(logintensity, intensity, data.model.tmp, intensity.sd, ...){
#   nu <- mean.lmvn(intensity, intensity.sd)
#   sd <- sd.lmvn(intensity, intensity.sd)
#   return((((nu - log(data.model.tmp$m.norm))^2)/sd))
# }
# 
# 
# fun.normalised <- function(logintensity, intensity, data.model.tmp, ...){
#   return(((intensity - data.model.tmp$m.norm)^2)/intensity^2)
# }
# 
# fun.normalised.by_priming <- function(logintensity, intensity, data.model.tmp, data, ...){
#   return(((intensity - data.model.tmp$m.norm)^2)/ as.numeric( 
#     data %>% 
#       filter(priming == data.model.tmp$priming) %>%
#       summarise(intensity = mean(intensity)))^2)
# }
# 
# fun.likelihood.list <- list(fun.likelihood.mvn.mean,
#                             fun.likelihood.mvn,
#                             fun.likelihood.lmvn, 
#                             fun.likelihood.mvn.sd_const,
#                             fun.likelihood.lmvn.data,
#                             fun.normalised,
#                             fun.normalised.by_priming)
