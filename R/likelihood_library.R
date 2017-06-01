### ###
### likelihood library 
### ###

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
# fun.likelihood.lmvn.data.norm <- function(logintensity, intensity, data.model.tmp, intensity.sd, ...){
#   nu <- mean.lmvn(intensity, intensity.sd)
#   sd <- sd.lmvn(intensity, intensity.sd)
#   norm <- ((log(intensity) - nu)^2)/sd
#   return((((nu - log(data.model.tmp$m.norm))^2)/sd)/norm)
# }
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
# fun.likelihood.list <- list(fun.likelihood.mvn.mean = fun.likelihood.mvn.mean,
#                             fun.likelihood.mvn = fun.likelihood.mvn,
#                             fun.likelihood.lmvn = fun.likelihood.lmvn, 
#                             fun.likelihood.mvn.sd_const = fun.likelihood.mvn.sd_const,
#                             fun.likelihood.lmvn.data = fun.likelihood.lmvn.data,
#                             fun.likelihood.lmvn.data.norm = fun.likelihood.lmvn.data.norm,
#                             fun.likelihood.normalised = fun.normalised,
#                             fun.likelihood.normalised.by_priming = fun.normalised.by_priming)
