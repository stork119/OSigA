### ###
### IRF model
### ###
#install.packages("entropy")
library(entropy)
library(FNN)
#### GetParametersRanges.irf1 ####
# hn_pstat hn_irf theta1 theta2 theta3 theta4 theta5 scale_pstat scale_irf sd_pstat sd_irf
GetParametersRanges.irf1 <- 
  function(
    ymin =  min(irfmodel.data.list$irfsum$logresponse),
    ymax =  max(irfmodel.data.list$irfsum$logresponse),
    scale.factor = ymax - ymin,
    scale.min = c(-4),
    scale.max = c(4),
    scale.base = c(2),
    sd.factor = mean(data.raw.sum$irf.sd[!is.na(data.raw.sum$irf.sd)]),
    sd.min = c(-4),
    sd.max = c(-4),
    sd.base = c(2),
    theta.factor = c(1,1,1),
    theta.min = c(-4,-4,-4),
    theta.max = c(4,4,4),
    theta.base = c(2,2,2),
    hn.factor = c(1),
    hn.min = c(-4),
    hn.max = c(4),
    hn.base = c(2),
    ...
  ){
    # hn_pstat theta1 theta2 theta3 scale_pstat sd_pstat
    ranges <- list()  
    ranges$max <- c(hn.max, theta.max, scale.max, sd.max)
    ranges$min <- c(hn.min, theta.min, scale.min, sd.min)
    ranges$factor <- c(hn.factor, theta.factor, scale.factor, sd.factor)
    ranges$base <- c(hn.base, theta.base, scale.base, sd.base)
    ranges$opt <- which(ranges$min < ranges$max)
    ranges$par <- 0*ranges.opt
    return(ranges)
  }
#### GetParametersList.irf1 ####
# hn_pstat hn_irf theta1 theta2 theta3 theta4 theta5 scale_pstat scale_irf sd_pstat sd_irf
GetParametersList.irf1 <- 
  function(
    params,
    ...
  ){
    # hn_pstat theta1 theta2 theta3 scale_pstat sd_pstat
    params.list <- list()  
    params.list$hn <- params[1]
    params.list$theta <- params[2:4]
    params.list$scale <- params[5]
    params.list$sd <- params[6]
    return(params.list)
  }
#### model_fun.irf1 ####
model_fun.irf1 <- function(
  ps1,
  ymax = c(1),
  ymin = c(0),
  theta = c(1, 1, 1),
  stm = 0,
  hn = 1,
  scale = ymax - ymin,
  ...){
  
  x <- ((theta[1]*ps1)^hn[1])/((theta[1]*ps1)^hn[1] + theta[2]^hn[1]) + theta[3] 
  y <- scale[1]*x
    return(list(x = x, y = y, sd.y = sd, sd.x = sd/scale[1]))
}

#### model_fun.df.irf1 ####
model_fun.df.irf1 <- function(
  ...){
  return((model_fun.irf1(...))$y)
}

#### simulateModel.irf ####
simulateModel.irf <- function(
  ranges,
  par = ranges$par,
  data.sample,
  nsimulations = 1000,
  ...
){
  data.model <- 
    data.sample  %>%
    data.table() %>% 
    dplyr::mutate(
      response.irf1.mean = 
        model_fun.df.irf1(
          theta = params$theta,
          hn    = params$hn,
          scale = params$scale,
          sd = params$sd,
          ps1 = response.pstat
        )) %>%
    dplyr::mutate(
      response.irf1.sd = params$sd
    ) %>%
    data.table()
    
  sample.list <- list()
  #dodac parallelize
  sample.list <- foreach(i = 1:nrow(data.model)) %do% {
    return(
      data.frame(
        response.irf = 
          exp(
            rnorm(n = nsimulations, 
                  mean = data.model[i,]$response.irf1.mean, 
                  sd = data.model[i,]$response.irf1.sd)),
        stimulation = data.model[i,]$stimulation
      )
    )
  }
  data.sample.new <- do.call(rbind, sample.list) %>%
    data.table()
  return(data.sample.new)
}

#### optimiset.fun.irf.1 ####
optimise.fun.irf1 <- function(
  par,
  stimulations,
  data.raw.list,
  data.sample.ps1,
  ranges.irf1,
  ...
){
  tryCatch({
    data.sample.irf1 <-  
      simulateModel.irf(
        ranges = ranges.irf1,
        par = par,
        data.sample = data.sample.ps1,
        nsimulations = nsimulations)
    
    stimulation.list <- (data.sample.irf1 %>% dplyr::distinct(stimulation))$stimulation
    kl.list <- foreach(stm = stimulation.list) %do% {
      kl.div <- KL.divergence(
        X = (irfmodel.data.list$irf %>% 
               dplyr::filter(stimulation == stm))$logresponse,
        Y = (data.sample.irf1 %>% 
          dplyr::filter(stimulation == stm))$response.irf,
        k = k)[k]
      return(kl.div)
    }
    
    kl <- do.call(sum, kl.list)
    return(kl)
  }, error = function(e){print(e)})
  return(Inf)
}