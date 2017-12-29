### ###
### IRF1 model ps1
### ###

#### ps1 ####
#### GetParametersRanges.ps1 ####
# hn_pstat hn_irf theta1 theta2 theta3 theta4 theta5 scale_pstat scale_irf sd_pstat sd_irf
GetParametersRanges.ps1 <- 
  function(
    ymin =  min(irfmodel.data.list$pSTATsum$logresponse),
    ymax =  max(irfmodel.data.list$pSTATsum$logresponse),
    scale.factor = ymax - ymin,
    scale.min = c(-4),
    scale.max = c(-4),
    scale.base = c(2),
    sd.factor = mean(data.raw.sum$pstat.sd[!is.na(data.raw.sum$pstat.sd)]),
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
  ranges$par <- 0*ranges$opt
  return(ranges)
  }

#### GetParametersList.ps1 ####
# hn_pstat hn_irf theta1 theta2 theta3 theta4 theta5 scale_pstat scale_irf sd_pstat sd_irf
GetParametersList.ps1 <- 
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

#### model_fun.ps1 ####
model_fun.ps1 <- function(
  ymax = c(1),
  ymin = c(0),
  theta = c(1, 1, 1),
  stm = 0,
  hn = 1,
  scale = ymax - ymin,
  sd = c(1),
  ...){
  
  x <- c(0)
  x[1] <- (theta[2] * (stm^hn[1])/(stm^hn[1] + theta[1]^hn[1]))+  theta[3]
  #x[2] <- (x[1]^hn[2])/(x[1]^hn[2] + theta[4]^hn[2]) + theta[5] 
  
  y <- c(0)
  y[1] <- scale[1]*x[1]
  #y[2] <- scale[2]*x[2]
  
  return(list(x = x, y = y, sd.y = sd, sd.x = sd/scale[1], sd = sd))
}
#### model model_fun_stm.ps1 ####
model_fun_stm.ps1 <- function(
  stimulations,
  ...
){
  y.list <- list()
  foreach( stm = stimulations ) %do% {
    result <-  model_fun.ps1(
      stm  = stm,
      ...)
    y.list[[as.character(stm)]] <-
    data.frame( 
        pstat = result$y[1],
        pstat.sd = result$sd.y[1],
        pstat.model.sd = result$sd.x[1],
        pstat.model = result$x[1],
        #irf   = result$y[2],
        #irf.sd = sd[2],
        #irf.model   = result$x[2],
        stimulation = stm,
        type = "model")
  }
  data.model <- do.call(rbind, y.list)
}

#### model model_fun_stm_params.ps1 ####
model_fun_stm_params.ps1 <- 
  function(params,
           ...){
    params.list <- GetParametersList.ps1(params)
    data.model <- model_fun_stm.ps1(
      hn = params.list$hn,
      theta = params.list$theta, 
      scale = params.list$scale,
      sd    = params.list$sd,
      ...
    )
    return(data.model)
}

# #### optimisation function ####
# optimise.fun.ps1 <- function(par,
#                      stimulations,
#                      data.raw.list,
#                      ranges.factor,
#                      ranges.base,
#                      ranges.opt,
#                      ...
#                      ){
#   params <- ranges.factor
#   params[ranges.opt] <- ranges.factor[ranges.opt]*ranges.base[ranges.opt]^par
#   
#   data.model <- model_fun_stm_params.ps1(
#     stimulations = stimulations,
#     params = params
#   )
#   likelihood.list <- foreach( data.i = 1:length(data.raw.list) ) %do% {
#       data  <- data.raw.list[[data.i]]
#       normalise <- (data %>%
#                       dplyr::mutate(normalise = (logresponse^2)/(logresponse.sd^2)) %>%
#                       dplyr::summarise(normalise = mean(normalise)))$normalise
#       likelihood <-
#         (data %>% 
#         left_join(
#           by =  "stimulation",
#           ((data.model %>% 
#               dplyr::mutate_(
#                 "logmodel" = data.model.colnames[data.i],
#                 "logmodel.sd" = paste(data.model.colnames[data.i], "sd", sep = ".")
#                 ))[, 
#                                                            c("logmodel", 
#                                                              "logmodel.sd",
#                                                              "stimulation")])) %>%
#         dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(2*(logmodel.sd^2))  + log(logmodel.sd)) %>%
#         #dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(logresponse.sd^2)) %>%
#         #  dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(logresponse^2)) %>%
#         dplyr::summarise(likelihood = sum(likelihood)))$likelihood
#       return(likelihood/normalise)         
#   }
#   likelihood <- do.call(what = sum, args = likelihood.list)
#   # print(likelihood)
#   # print("\n")
#   # print(par)
#   return(as.numeric(likelihood))
# }

#### simulateModel ####
simulateModel.ps1 <- function(
  ranges,
  par = ranges$par,
  stimulations,
  nsimulations = 1000,
  ...
){
  params <- ranges$factor
  params[ranges$opt] <- ranges$factor[ranges$opt]*ranges$base[ranges$opt]^par
  data.model <- model_fun_stm_params.ps1(
    stimulations = stimulations,
    params = params
  )
  sample.list <- list()
  sample.list <- foreach(i = 1:nrow(data.model)) %do% {
    return(
      data.frame(
        response.pstat = 
          exp(
            rnorm(n = nsimulations, 
                  mean = data.model[i,]$pstat.model, 
                  sd = data.model[i,]$pstat.model.sd)),
          stimulation = data.model[i,]$stimulation
        )
      )
  }
  data.sample <- do.call(rbind, sample.list) %>% data.table()
  return(data.sample)
}
