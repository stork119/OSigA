### ###
### lmvn
### ###

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