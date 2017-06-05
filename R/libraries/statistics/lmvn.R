### ###
### LOG - MVN distribution
### ###

#### mean.lmvn ####
lmvn.mean <- function(m.norm, sd.norm){
  ifelse(m.norm != 0.0, (log(m.norm) - (1/2)*log((sd.norm/(m.norm^2)) + 1)), 0.0)
}

#### sd.norm.lmvn ####
lmvn.sd <- function(m.norm, sd.norm){
  ifelse(m.norm != 0.0, log((sd.norm/(m.norm^2)) + 1), 0.0)
}

### ###
#### LMVN ####
mean.lmvn <- function(m.norm, sd.norm){
  ifelse(m.norm != 0.0, (log(m.norm) - (1/2)*log((sd.norm/(m.norm^2)) + 1)), 0.0)
}

sd.lmvn <- function(m.norm, sd.norm){
  ifelse(m.norm != 0.0, log((sd.norm/(m.norm^2)) + 1), 0.0)
}

lmvn <- function(data){
  data %>%
    mutate(mean.lmvn = mean.lmvn(m.norm, sd.norm),
           sd.lmvn = sd.lmvn(m.norm, sd.norm))
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


