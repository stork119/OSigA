### ###
### Unscented transform model 
### ###

#### PrepareModelArguments.ut ####
PrepareModelParameters.ut <-
  function(parameters.model,
           priming_constant = 3.4,
           ...)
    {
    parameters.model[11] <- parameters.model[17]*parameters.model[11]
    parameters.model[12] <- (parameters.model[17]^2)*parameters.model[12]
    parameters.model[15] <- parameters.model[17]*parameters.model[15]
    parameters.model[16] <- (parameters.model[17]^2)*parameters.model[16]
    parameters.model[13] <- parameters.model[18]*parameters.model[13]
    parameters.model[14] <- (parameters.model[18]^2)*parameters.model[13]
    return(parameters.model)
}

PrepareModelArguments.ut <-
  function(parameters,
           parameters.priming = parameters,
           variables,
           variables.priming,
           priming_constant = 3.4,
           ...) {
    
    parameters<- parameters[1:10]
    parameters.priming<- parameters[1:10]
    
    return(list(parameters = parameters,
                parameters.priming = parameters.priming,
                variables = variables,
                variables.priming = variables.priming
    ))}


PrepareModelArguments.ut.multiple <-
  function(parameters,
           parameters.priming = parameters,
           variables,
           variables.priming,
           parameters.conditions,
           priming_constant = 3.4,
           ...) {
    
    # print(parameters.conditions)
    parameters.model <- rep(0, times = length(which(parameters.conditions$parameters != 0)))
    parameters.priming.model <- rep(0, times = length(which(parameters.conditions$parameters.priming != 0)))
    
    parameters.model[parameters.conditions$parameters[which(parameters.conditions$parameters != 0)]] <- 
      parameters[which(parameters.conditions$parameters != 0)]
    
    parameters.priming.model[parameters.conditions$parameters.priming[which(parameters.conditions$parameters.priming != 0)]] <- 
      parameters[which(parameters.conditions$parameters.priming != 0)]
    
    return(list(parameters = parameters.model,
                parameters.priming = parameters.priming.model,
                variables = variables,
                variables.priming = variables.priming
    ))}

#### LoadSigmapointsConditions ####
LoadSigmapointsConditions <- function(path.optimisation){
  path.sigmapoints <-  paste(path.optimisation, "sigmapoints_conditions.csv", sep = "")
  if(file.exists(path.sigmapoints)){
    sigmapoints.conditions <- read.table(file = path.sigmapoints, 
                                         header = TRUE, 
                                         sep = ",")
  }
  path.sigmapoints.parameters <-  paste(path.optimisation, "sigmapoints_parameters_conditions.csv", sep = "")
  if(file.exists(path.sigmapoints.parameters)){
    sigmapoints.parameters.conditions <- read.table(file = path.sigmapoints.parameters, 
                                         header = TRUE, 
                                         sep = ",")
  }
  #optimisation.procedure <- optimisation_ut
  fun_modify_input <- PrepareModelArguments.ut
  fun_run_model <- run_model_ut
  return(list(conditions = sigmapoints.conditions,
              parameters.conditions = sigmapoints.parameters.conditions,
          #    optimisation.procedure = optimisation.procedure,
          fun_run_model = fun_run_model
              ))
  
}
####  GetWeigths ####
GetWeigths <- function(
  alpha = 0.7,
  kappa = 0,
  beta = 2,
  D # number of sigma points
){
  weights.means <- c()
  weights.variance <- c()
  weights.means[1]    <- 1 - D/((alpha^2)*(D + kappa))
  weights.variance[1] <- weights.means[1] + 1 - (alpha^2) + beta 
  weights.means[2:(2*D + 1)] <- 1/ (2*(alpha^2)*(D+kappa))
  weights.variance[2:(2*D + 1)] <- weights.means[2:(2*D + 1)] 
  return(
    list(weights.means = weights.means, 
         weights.variance =weights.variance))
}
#### AggreagateSimulationData ####
AggreagateSimulationData <- function(
  data.model,
  data.trajectory,
  res,
  tmesh,
  tmesh.list.tmp,priming, 
  stm
  ){
  for(tmesh.i in tmesh.list.tmp){
    data.model <- rbind(data.model,
                        data.table(time = c(tmesh[tmesh.i]),
                                   m = res$output[[tmesh.i]][14],
                                   sd =  c(ifelse(length(res$output[[tmesh.i]]) > 17, res$output[[tmesh.i]][31], 0)),
                                   priming = priming, 
                                   stimulation = stm)
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
  fun_run_model = rmain,
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
  stm,
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
   # print(parameters.model)
   # print(parameters.priming.model)
  for(stm in stimulation.list){
    if(model.computations$raw){
      res <- rmain(parameters = parameters.model, 
                         variables = variables, 
                         stm = stm, 
                         tmesh = tmesh, 
                         time_interval = time_interval, 
                         time_computation = time_computation)
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
      res.priming <- rmain(parameters = parameters.priming.model, 
                                 variables = variables.priming, 
                                 stm = stm, 
                                 tmesh = tmesh, 
                                 time_interval = time_interval, 
                                 time_computation = time_computation)
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

#### GetSigmapoints ####
GetSigmapoints <- function(
  sigmapoints.parameters,
  sigmapoints.parameters.sd,
  alpha = 0.7,
  kappa = 0,
  beta = 2
  
){
  D <- length(sigmapoints.parameters)
  
  sigmapoints.parameters.lmvn.mean <- lmvn.mean.vector(m.norm = sigmapoints.parameters,
                                                sd.norm = sigmapoints.parameters.sd)
  sigmapoints.parameters.lmvn.sd   <- lmvn.sd.vector(m.norm = sigmapoints.parameters, 
                                              sd.norm = sigmapoints.parameters.sd)
  
  sigmapoints.list <- list()
  sigmapoints.list[[1]] <- exp(sigmapoints.parameters.lmvn.mean)
  for(i in 1:length(sigmapoints.parameters.lmvn.mean)){
    sigmapoints.parameters.lmvn.sd_tmp <- sigmapoints.parameters.lmvn.sd
    sigmapoints.parameters.lmvn.sd_tmp[-i] <- 0
    sigmapoints.list[[i+1]] <- exp(sigmapoints.parameters.lmvn.mean + alpha*(sqrt(D + kappa))*sqrt(sigmapoints.parameters.lmvn.sd_tmp))
    sigmapoints.list[[D + i + 1]] <- exp(sigmapoints.parameters.lmvn.mean - alpha*(sqrt(D + kappa))*sqrt(sigmapoints.parameters.lmvn.sd_tmp))
  }
  return(sigmapoints.list)
}

#### ut.fun_sigmapoints ####
ut.fun_sigmapoints <- 
  function(sigmapoints.parameters.conditions,
           sigmapoints.conditions,
           par.all,
           parameters.factor,
           parameters.base,
           variables,
           variables.priming,
           fun_modify_input,
           fun_modify_parameters = function(parameters.model){return(parameters.model)},
           ...
           ){
    
    # print(fun_modify_input)
    
    # sigmapoints.par <- par.all
    # sigmapoints.par[sigmapoints.parameters.conditions$id.par] <- sigmapoints.list[[i]]
    # 
    parameters <- parameters.factor*(parameters.base)^par.all
    parameters <- fun_modify_parameters(parameters.model = parameters)
    
    sigmapoints.list <- 
      GetSigmapoints(sigmapoints.parameters = parameters[sigmapoints.parameters.conditions$id.par],
                      sigmapoints.parameters.sd = parameters[sigmapoints.parameters.conditions$id.par.sd],
                     alpha = sigmapoints.conditions$alpha,
                     beta = sigmapoints.conditions$beta,
                     kappa = sigmapoints.conditions$kappa)
    arguments.list <- list()
    for(i in 1:length(sigmapoints.list)){
      sigmapoints.parameters <- parameters 
      sigmapoints.parameters[sigmapoints.parameters.conditions$id.par] <- sigmapoints.list[[i]]
      sigmapoints.variables <- variables
      sigmapoints.variables.priming <- variables.priming
      sigmapoints.variables[ 
        sigmapoints.parameters.conditions$variables[
          which(sigmapoints.parameters.conditions$variables != 0)]] <- 
        sigmapoints.parameters[
            sigmapoints.parameters.conditions$id.par[
              which(sigmapoints.parameters.conditions$variables != 0)]]
      sigmapoints.variables.priming[ 
        sigmapoints.parameters.conditions$variables.priming[
          which(sigmapoints.parameters.conditions$variables.priming != 0)]] <- 
        sigmapoints.parameters[
            sigmapoints.parameters.conditions$id.par[
              which(sigmapoints.parameters.conditions$variables.priming != 0)]]
      
      arguments.list[[i]] <- fun_modify_input(parameters = sigmapoints.parameters,
                                variables = sigmapoints.variables,
                                variables.priming = sigmapoints.variables.priming,
                                sigmapoints.parameters.conditions = sigmapoints.parameters.conditions,
                                ...)
    }
    return(arguments.list)
}

#### run_model_ut ####
run_model_ut <- function(
  parameters,
  par,
  parameters.base,
  parameters.factor,
  variables, 
  variables.priming, 
  tmesh, 
  tmesh.list,
  stimulation.list,
  background,
  fun.likelihood,
  par.optimised = rep(1, times = length(par)),
  fun_modify_input = PrepareModelArguments.ut,
  sigmapoints,
  ...){
  ### run 
  #par <- as.numeric(par.list[[11]])
  par.all <- rep(0, times = length(parameters.factor))
  par.all[par.optimised] <- par
  
  arguments.list  <-  ut.fun_sigmapoints(
    par.all = par.all,
    sigmapoints.parameters.conditions = sigmapoints$parameters.conditions,
    sigmapoints.conditions = sigmapoints$conditions,
    parameters.factor = parameters.factor,
    parameters.base = parameters.base,
    variables = variables,
    variables.priming = variables.priming,
    fun_modify_input = fun_modify_input,
    ...)

  
  
  sigmapoints$weights <- GetWeigths(alpha = sigmapoints$conditions$alpha, 
                                    kappa = sigmapoints$conditions$kappa,
                                    beta = sigmapoints$conditions$beta,
                                    D = nrow(sigmapoints$parameters.conditions))
  ###
  data.model.list <- list()
  data.trajectory.list <- list()
  for(argument.i in 1:length(arguments.list)){
    
    input <- arguments.list[[argument.i]]
    parameters <- input$parameters
    parameters.priming <- input$parameters.priming
    variables  <- input$variables
    variables.priming <- input$variables.priming
    

    model.simulation<- do.call(simulate_model_ut,
                               list(
                                 parameters.model = parameters,
                                 parameters.priming.model = parameters.priming,
                                 variables = variables,
                                 variables.priming = variables.priming,
                                 tmesh = tmesh,
                                 tmesh.list = tmesh.list,
                                 stimulation.list = stimulation.list,
                                 background = background,
                                 ...))
    data.trajectory.list[[argument.i]] <- model.simulation$data.trajectory
    data.trajectory.list[[argument.i]]$sigmapoint <- argument.i
    data.model.list[[argument.i]] <- model.simulation$data.model
    data.model.list[[argument.i]]$sigmapoint <- argument.i
    
    if(model.simulation$error){
      return(list(error = TRUE))
    } 
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
               data.trajectory = data.trajectory,
               arguments.list = arguments.list
  ))
}

#### optimisation ####
optimisation_ut <- function(par,
                            fun_run_model = run_model_ut,
                            parameters.base,
                            parameters.factor,
                            variables, 
                            variables.priming, 
                            tmesh, 
                            tmesh.list,
                            stimulation.list,
                            background,
                            data.exp.grouped,
                            data.exp.summarise,
                            return.model = FALSE,
                            fun.likelihood,
                            par.optimised = rep(1, times = length(par)),
                            fun_modify_input = PrepareModelArguments.ut,
                            sigmapoints,
                            ...)
{
  
  model.simulation <- fun_run_model(par = par,
                                    parameters.base = parameters.base,
                                    parameters.factor = parameters.factor,
                                    variables = variables, 
                                    variables.priming = variables.priming, 
                                    tmesh = tmesh, 
                                    tmesh.list = tmesh.list,
                                    stimulation.list = stimulation.list,
                                    background = background,
                                    par.optimised = par.optimised,
                                    fun_modify_input = fun_modify_input,
                                    sigmapoints = sigmapoints, 
                                    ...)
            
  if(model.simulation$error){
    return(Inf)
  }
  
  result <- sum(
    likelihood(fun.likelihood = fun.likelihood,
               data.model = model.simulation$data.model,
               data.exp.grouped = data.exp.grouped,
               data.exp.summarise = data.exp.summarise))
  #print(paste(c(result, par), sep = " "))
  
  if(return.model){
    return(list( error = FALSE, 
                 #optimisation = result,
                 data.model = model.simulation$data.model, 
                 data.model.sigmapoints = model.simulation$data.model.sigmapoints,
                 data.model.ut = model.simulation$data.model.ut,
                 arguments.list = model.simulation$arguments.list
    ))
  }
  return(result)
}
