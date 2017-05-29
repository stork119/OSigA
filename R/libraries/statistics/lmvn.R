### ###
### LOG - MVN distribution
### ###

#### mean.lmvn ####
lmvn.mean <- function(m, sd){
  ifelse(m != 0.0, (log(m) - (1/2)*log((sd/(m^2)) + 1)), 0.0)
}

#### sd.lmvn ####
lmvn.sd <- function(m, sd){
  ifelse(m != 0.0, log((sd/(m^2)) + 1), 0.0)
}

#### lmvn ####
### function computes mean and sd of lognormal distribution
lmvn <- function(data){
  data %>%
    mutate(mean.lmvn = 
             mean.lmvn(m.norm, sd.norm),
           sd.lmvn = 
             sd.lmvn(m.norm, sd.norm))
}