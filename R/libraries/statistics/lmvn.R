### ###
### LOG - MVN distribution
### ###

#### mean.lmvn ####
lmvn.mean <- function(m.norm, sd.norm){
  ifelse(m.norm != 0.0, (log(m.norm) - (1/2)*log((sd.norm/(m.norm^2)) + 1)), 0.0)
}

lmvn.mean.vector <- function(m.norm, sd.norm){
  if(length(m.norm) != length(sd.norm)){
    stop("lmvn.mean.vector: vectors m.norm and sd.norm should have the same length")
  }
  sapply(1:length(m.norm), 
         function(i){
           lmvn.mean(m.norm = m.norm[i],
                     sd.norm = sd.norm[i])
         })
}

#### sd.norm.lmvn ####
lmvn.sd <- function(m.norm, sd.norm){
  ifelse(m.norm != 0.0, log((sd.norm/(m.norm^2)) + 1), 0.0)
}

lmvn.sd.vector <- function(m.norm, sd.norm){
  if(length(m.norm) != length(sd.norm)){
    stop("lmvn.sd.vector: vectors m.norm and sd.norm should have the same length")
  }
  sapply(1:length(m.norm), 
         function(i){
           lmvn.sd(m.norm = m.norm[i],
                     sd.norm = sd.norm[i])
         })
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


