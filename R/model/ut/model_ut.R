### ###
### Unscented transformation model
### running model
### ###
##

#### fun_parameters_penalty ####
fun_parameters_penalty_sigmapoints <- function(par, 
                                   parameters.conditions,
                                   ...){
  a <- parameters.conditions$lower[which(parameters.conditions$lower != parameters.conditions$upper)]
  b <- parameters.conditions$upper[which(parameters.conditions$lower != parameters.conditions$upper)]
  base <- parameters.conditions$base[which(parameters.conditions$lower != parameters.conditions$upper)]
  return(sum(((base^par - base^0)^2)/(base^(2*b) - base^(2*a))))
}
#### PrepareModelParameters.ut ####
PrepareModelParameters.ut <-
  function(parameters.model,
           priming_constant = 3.4,
           ...)
    {
    flog.debug("model_ut.R PrepareModelParameters.ut", name="logger.optimisation")
    
   # parameters.model[15] <- priming_constant*parameters.model[11]
    parameters.model[11] <- parameters.model[17]*parameters.model[11]
    parameters.model[12] <- (parameters.model[17]^2)*parameters.model[12]
    parameters.model[15] <- parameters.model[17]*parameters.model[15]
    parameters.model[16] <- (parameters.model[17]^2)*parameters.model[16]
    parameters.model[13] <- parameters.model[18]*parameters.model[13]
    parameters.model[14] <- (parameters.model[18]^2)*parameters.model[13]
    return(parameters.model)
}

#### PrepareModelArguments.ut.multiple ####
PrepareModelArguments.ut.multiple <-
  function(parameters,
           parameters.priming = parameters,
           variables,
           variables.priming,
           priming_constant = 3.4,
           parameters.conditions,
           p1p6_constant = 1.635657e-07,
           ...) {
    
    flog.debug("model_ut.R PrepareModelArguments.ut.multiple", 
               name = "logger.optimisation")
    
    parameters.model <- rep(0, times = length(which(parameters.conditions$parameters != 0)))
    parameters.priming.model <- rep(0, times = length(which(parameters.conditions$parameters.priming != 0)))
    
    parameters.model[parameters.conditions$parameters[which(parameters.conditions$parameters != 0)]] <- 
      parameters[which(parameters.conditions$parameters != 0)]
    
    parameters.priming.model[parameters.conditions$parameters.priming[which(parameters.conditions$parameters.priming != 0)]] <- 
      parameters[which(parameters.conditions$parameters.priming != 0)]
    
    # parameters.model[1] <- parameters.model[6]*p1p6_constant
    # parameters.priming.model[1] <- parameters.priming.model[6]*p1p6_constant
    # 
    return(list(parameters = parameters.model,
                parameters.priming = parameters.priming.model,
                variables = variables,
                variables.priming = variables.priming
    ))}

#### AggreagateSimulationData ####
AggreagateSimulationData <- function(
  data.model,
  data.trajectory,
  res,
  tmesh,
  tmesh.list.tmp,priming, 
  stm,
  # model.observable.mean = 14,
  # model.observable.sd   = 31,
  # model.variables.num   = 17, 
  model.observable.mean = 4,
  model.observable.sd   = 12,
  model.variables.num   = 8 
  ){
  for(tmesh.i in tmesh.list.tmp){
    data.model <- rbind(data.model,
                        data.table(
                          time = 
                            c(tmesh[tmesh.i]),
                          m    = 
                            res$output[[tmesh.i]][model.observable.mean],
                          sd   = 
                            c(ifelse(
                              length(res$output[[tmesh.i]]) > model.variables.num,
                              res$output[[tmesh.i]][model.observable.sd],
                              0)),
                          priming = priming, 
                          stimulation = stm
                        )
    )
    
    data.trajectory <- rbind(data.trajectory,
                             data.table(
                               time = rep(x = tmesh[tmesh.i],
                                          times= length(res$output[[tmesh.i]])),
                               m = res$output[[tmesh.i]],
                               priming = rep(x = priming, times= length(res$output[[tmesh.i]])),
                               stimulation = rep(x = stm, times= length(res$output[[tmesh.i]])),
                               var  = 1:length(res$output[[tmesh.i]]) 
                             )
    )
  }
  return(list(data.trajectory = data.trajectory, data.model = data.model)) 
}
  

#### simulate_model_ut ####
simulate_model_ut <- function(
  fun_model_ode = rmain,
  parameters.model, 
  parameters.priming.model = parameters.model,
  variables,
  variables.priming,
  tmesh,
  tmesh.list,
  tmesh.list.tmp = NULL,
  stimulation.list,
  background,
  time_interval = 100,
  time_computation = 1000*60*5,
  model.computations = list(raw = TRUE, priming = TRUE),
  #stm,
  ...
  ){
#  print(model.computations)
  if(is.null(tmesh.list.tmp)){
    tmesh.list.tmp <- tmesh.list 
  }
  data.model <- data.table(
    time = numeric(),
    m = numeric(),
    sd = numeric(),
    priming = numeric(), 
    stimulation = numeric()
  )
  data.trajectory <- data.table(
    time = numeric(),
    var = numeric(),
    m = numeric(),
    priming = numeric(), 
    stimulation = numeric()
  )
  for(stm in stimulation.list){
    # print(stm)
    if(model.computations$raw){
      # print(parameters.model)
      # print(variables)
      res <- do.call(
        fun_model_ode,
        list(parameters = parameters.model, 
             variables = variables, 
             stm = stm, 
             tmesh = tmesh, 
             time_interval = time_interval, 
             time_computation = time_computation)
      )
      # print(res$success)
      if(res$success){
        res$aggregate <- AggreagateSimulationData(
          data.model = data.model,
          data.trajectory = data.trajectory,
          res = res,
          tmesh = tmesh,
          tmesh.list.tmp = tmesh.list.tmp,
          priming = 0, 
          stm = stm)
          data.model <- res$aggregate$data.model
          data.trajectory <- res$aggregate$data.trajectory
      } else {
        return(list(error = TRUE))
      }
    }
    if(model.computations$priming){
      res.priming <- do.call(
        fun_model_ode,
        list(parameters = parameters.priming.model, 
             variables = variables.priming, 
             stm = stm, 
             tmesh = tmesh, 
             time_interval = time_interval, 
             time_computation = time_computation))
      # print(res$success)
      if(res.priming$success){
        res$aggregate <- AggreagateSimulationData(
          data.model = data.model,
          data.trajectory = data.trajectory,
          res = res.priming,
          tmesh = tmesh,
          tmesh.list.tmp = tmesh.list.tmp,
          priming = 1000, 
          stm = stm)
        data.model <- res$aggregate$data.model
        data.trajectory <- res$aggregate$data.trajectory
      } else {
        return(list(error = TRUE))
      }
    }
  }
  data.model <- normalization_simulation(data.model, background = background)
  data.model <- lmvn(data.model)
  return(list(error = FALSE, 
              data.model = data.model,
              data.trajectory = data.trajectory))
}

#### run_model_ut ####
run_model_ut <- function(
  parameters,
  par,
  parameters.base,
  parameters.factor,
  parameters.conditions,
  variables, 
  variables.priming, 
  tmesh, 
  tmesh.list,
  stimulation.list,
  background,
  fun.likelihood,
  par.optimised = rep(1, times = length(par)),
  fun_modify_input = PrepareModelArguments.ut.multiple,
  sigmapoints,
  ...){
  
  
  # flog.debug("run_model_ut %s", parameters.conditions)
  ### run 
  par.all <- rep(0, times = length(parameters.factor))
  par.all[par.optimised] <- par
  
  arguments.list  <-  GetSigmapointsArguments(
    par.all = par.all,
    sigmapoints.parameters.conditions = sigmapoints$parameters.conditions,
    sigmapoints.conditions = sigmapoints$conditions,
    parameters.conditions = parameters.conditions,
    parameters.factor = parameters.factor,
    parameters.base = parameters.base,
    variables = variables,
    variables.priming = variables.priming,
    fun_modify_input = fun_modify_input
    #)
    ,...)
  
  sigmapoints$weights <- 
    GetSigmapointsWeigths(alpha = sigmapoints$conditions$alpha, 
                          kappa = sigmapoints$conditions$kappa,
                          beta = sigmapoints$conditions$beta,
                          D = nrow(sigmapoints$parameters.conditions))
  ###
  data.model.list <- list()
  data.trajectory.list <- list()
  for(argument.i in 1:length(arguments.list)){
    # print(argument.i)
    parameters <- arguments.list[[argument.i]]$parameters
    parameters.priming <- arguments.list[[argument.i]]$parameters.priming
    variables  <- arguments.list[[argument.i]]$variables
    variables.priming <- arguments.list[[argument.i]]$variables.priming
    

    model.simulation<- do.call(simulate_model_ut,
                               list(
                                 parameters.model = parameters,
                                 parameters.priming.model = parameters.priming,
                                 variables = variables,
                                 variables.priming = variables.priming,
                                 tmesh = tmesh,
                                 tmesh.list = tmesh.list,
                                 stimulation.list = stimulation.list,
                                 background = background
                                 #fun_model_ode =fun_model_ode))
                                 ,...))
    if(model.simulation$error){
      return(list(error = TRUE))
    } 
    
    data.model.list[[argument.i]] <- model.simulation$data.model
    data.model.list[[argument.i]]$sigmapoint <- argument.i
    data.trajectory.list[[argument.i]] <- model.simulation$data.trajectory
    data.trajectory.list[[argument.i]]$sigmapoint <- argument.i
    data.model.list[[argument.i]] <- model.simulation$data.model
    data.model.list[[argument.i]]$sigmapoint <- argument.i
     
  }
  
  data.model.sigmapoints <- do.call(
    rbind,
    data.model.list) %>% 
    data.table()
  data.trajectory <- do.call(
    rbind,
    data.trajectory.list) %>% 
    data.table()
  data.model.ut <- 
    data.model.sigmapoints %>% 
    dplyr::mutate(mean.lmvn.ut = 
                    sigmapoints[["weights"]][["weights.means"]][sigmapoint]*mean.lmvn,
                  sd_intrinsic.lmvn.ut = 
                    sigmapoints[["weights"]][["weights.means"]][sigmapoint]*sd.lmvn,
                  m.norm.ut = 
                    sigmapoints[["weights"]][["weights.means"]][sigmapoint]*m.norm,
                  sd_intrinsic.norm.ut = 
                    sigmapoints[["weights"]][["weights.means"]][sigmapoint]*sd.norm) %>%
    dplyr::group_by(time, priming, stimulation) %>% 
    dplyr::summarise(mean.lmvn.ut = sum(mean.lmvn.ut),
                     sd_intrinsic.lmvn.ut = sum(sd_intrinsic.lmvn.ut),
                     m.norm.ut = sum(m.norm.ut),
                     sd_intrinsic.norm.ut = sum(sd_intrinsic.norm.ut)) %>%
    data.table()
  
  data.model.ut <-
    data.model.sigmapoints %>% 
    left_join(data.model.ut, by = c("time", "priming", "stimulation")) %>% 
    data.table() %>% 
    dplyr::mutate(sd_extrinsic.lmvn.ut = 
                    sigmapoints[["weights"]][["weights.variance"]][sigmapoint]*((mean.lmvn - mean.lmvn.ut)^2), 
                  sd_extrinsic.norm.ut = 
                    sigmapoints[["weights"]][["weights.variance"]][sigmapoint]*((m.norm - m.norm.ut)^2)) %>%
    dplyr::group_by(time, priming, stimulation) %>% 
    dplyr::summarise(mean.lmvn.ut = mean(mean.lmvn.ut),
                     sd_intrinsic.lmvn.ut = mean(sd_intrinsic.lmvn.ut),
                     sd_extrinsic.lmvn.ut = sum(sd_extrinsic.lmvn.ut),
                     m.norm.ut = mean(m.norm.ut),
                     sd_intrinsic.norm.ut = mean(sd_intrinsic.norm.ut),
                     sd_extrinsic.norm.ut = sum(sd_extrinsic.norm.ut)) %>%
    dplyr::mutate(sd.lmvn.ut = sd_intrinsic.lmvn.ut + sd_extrinsic.lmvn.ut, 
                  sd.norm.ut = sd_intrinsic.norm.ut + sd_extrinsic.norm.ut) %>%
    data.table()
  
  
  data.model <- 
    data.model.ut %>% 
    dplyr::mutate(mean.lmvn = mean.lmvn.ut,
                  sd.lmvn = sd.lmvn.ut,
                  m.norm = m.norm.ut,
                  sd.norm = sd.norm.ut) %>%
    dplyr::select(time, priming, stimulation, m.norm, sd.norm, mean.lmvn, sd.lmvn)
  
  return(list( error = FALSE, 
               data.model = data.model, 
               data.model.sigmapoints = data.model.sigmapoints,
               data.model.ut = data.model.ut,
               arguments.list = arguments.list,
               data.trajectory = data.trajectory
  ))
}