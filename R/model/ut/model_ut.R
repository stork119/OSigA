### ###
### Unscented transform model 
### ###

#### PrepareModelArguments.ut ####
PrepareModelArguments.ut <-
  function(parameters,
           parameters.priming = parameters,
           variables,
           variables.priming,
           priming_constant = 2.6,
           ...) {
    
    variables.priming[1] <- priming_constant*variables[1]
    
    parameters<- parameters[1:10]
    parameters.priming<- parameters[1:10]
    
    return(list(parameters = parameters,
                parameters.priming = parameters.priming,
                variables = variables,
                variables.priming = variables.priming
    ))}

#### LoadSigmapointsConditions ####
LoadSigmapointsConditions <- function(){
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
              fun_modify_input = PrepareModelArguments.ut,
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

#### simulate_model_ut ####
simulate_model_ut <- function(
  fun_run_model = rmain,
  parameters, 
  parameters.priming = parameters,
                      variables,
                      variables.priming,
                      tmesh,
                      tmesh.list,
                      tmesh.list.tmp = NULL,
                      stimulation.list,
                      background,
                      time_interval = 100,
                      time_computation = 1000*60*5,
                      stm,
                      ...){
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
  # print(parameters)
  for(stm in stimulation.list){
    res <- rmain(parameters = parameters, 
                         variables = variables, 
                         stm = stm, 
                         tmesh = tmesh, 
                         time_interval = time_interval, 
                         time_computation = time_computation)
    res.priming <- rmain(parameters = parameters.priming, 
                                 variables = variables.priming, 
                                 stm = stm, 
                                 tmesh = tmesh, 
                                 time_interval = time_interval, 
                                 time_computation = time_computation)
    
    
    if(res$success & res.priming$success){
      for(tmesh.i in tmesh.list.tmp){
        data.model <- rbind(data.model,
                            data.table(time = c(tmesh[tmesh.i], tmesh[tmesh.i]),
                                       m = c(res$output[[tmesh.i]][14],
                                             res.priming$output[[tmesh.i]][14]),
                                       sd =  c(ifelse(length(res$output[[tmesh.i]]) > 17, res$output[[tmesh.i]][31], 0),
                                               ifelse(length(res.priming$output[[tmesh.i]]) > 17, res.priming$output[[tmesh.i]][31], 0)),
                                       priming = c(0, 1000), 
                                       stimulation = c(stm, stm))
        )
        
        data.trajectory <- rbind(data.trajectory,
                                 data.table(
                                   time = rep(x = tmesh[tmesh.i],
                                              times= 2*length(res$output[[tmesh.i]])),
                                   m = c(res$output[[tmesh.i]],
                                         res.priming$output[[tmesh.i]]),
                                   priming = c(rep(x = 0, times= length(res$output[[tmesh.i]])),
                                               rep(x = 1000, times= length(res$output[[tmesh.i]]))), 
                                   stimulation = rep(x = stm, times= 2*length(res$output[[tmesh.i]])),
                                   var  = rep(x = 1:length(res$output[[tmesh.i]]), times = 2) 
                                 )
        )
      }
    } else {
      return(list(error = TRUE))
    }
  }
  data.model <- normalization_simulation(data.model, background = background)
  data.model <- lmvn(data.model)
  return(list(error = FALSE, data.model = data.model, data.trajectory = data.trajectory))
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
           fun_modify_input = fun_modify_input,
           ...
           ){
    
    # sigmapoints.par <- par.all
    # sigmapoints.par[sigmapoints.parameters.conditions$id.par] <- sigmapoints.list[[i]]
    # 
    parameters <- parameters.factor*(parameters.base)^par.all
    
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
        sigmapoints.parameters[sigmapoints.parameters.conditions$id.par]
      sigmapoints.variables.priming[ 
        sigmapoints.parameters.conditions$variables.priming[
          which(sigmapoints.parameters.conditions$variables.priming != 0)]] <- 
        sigmapoints.parameters[sigmapoints.parameters.conditions$id.par]
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
    fun_modify_input = sigmapoints$fun_modify_input)
  
  
  
  sigmapoints$weights <- GetWeigths(alpha = sigmapoints$conditions$alpha, 
                                    kappa = sigmapoints$conditions$kappa,
                                    beta = sigmapoints$conditions$beta,
                                    D = nrow(sigmapoints$parameters.conditions))
  ###
  data.model.list <- list()
  for(argument.i in 1:length(arguments.list)){
    
    input <- arguments.list[[argument.i]]
    parameters <- input$parameters
    variables  <- input$variables
    variables.priming <- input$variables.priming
    
    model.simulation<- do.call(simulate_model_ut,
                               list(
                                 parameters = parameters,
                                 variables = variables,
                                 variables.priming = variables.priming,
                                 tmesh = tmesh,
                                 tmesh.list = tmesh.list,
                                 stimulation.list = stimulation.list,
                                 background = background))
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
    dplyr::select(time, priming, stimulation, m.norm, sd.norm, mean.lmvn, sd.norm)
  
  return(list( error = FALSE, 
               data.model = data.model, 
               data.model.sigmapoints = data.model.sigmapoints,
               data.model.ut = data.model.ut,
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
                                    sigmapoints =sigmapoints,
                                    ...)
            
  if(model.simulation$error){
    return(Inf)
  }
  
  result <- sum(
    likelihood(fun.likelihood = fun.likelihood,
               data.model = model.simulation$data.model,
               data.exp.grouped = data.exp.grouped,
               data.exp.summarise = data.exp.summarise))
  print(result)
  
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
