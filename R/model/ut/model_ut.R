### ###
### Unscented transform model 
### ###

#### ####
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

#### ####
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
  optimisation.procedure <- optimisation_ut
  fun_modify_input <- PrepareModelArguments.ut
  
  
}
#### ####
GetSigmapoints <- function(
  sigmapoints.parameters,
  sigmapoints.parameters.sd,
  alpha = 0.7,
  kappa = 0,
  beta = 2
  
){
  D <- length(sigmapoints.parameters)
  sigmapoints.list <- list()
  sigmapoints.list[[1]] <- sigmapoints.parameters
  for(i in 1:length(sigmapoints.parameters)){
    sigmapoints.parameters.sd_tmp <- sigmapoints.parameters.sd
    sigmapoints.parameters.sd_tmp[-i] <- 0
    sigmapoints.list[[i+1]] <- sigmapoints.parameters + alpha*(sqrt(D + kappa))*sigmapoints.parameters.sd_tmp
    sigmapoints.list[[D + i + 1]] <- sigmapoints.parameters - alpha*(sqrt(D + kappa))*sigmapoints.parameters.sd_tmp
  }
  return(sigmapoints.list)
}

#### ####
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
  weights.means[2:(2*D + 1)] <- weights.variance[2:(2*D + 1)] 
  return(
    list(weights.means = weights.means, 
         weights.variance =weights.variance))
}

#### ####
run_model_ut <- function(parameters, 
                      variables,
                      variables.priming,
                      tmesh,
                      tmesh.list,
                      stimulation.list,
                      background,
                      time_interval = 100,
                      time_computation = 1000*60*5){
  data.model <- data.table(time = numeric(),
                           m = numeric(),
                           sd = numeric(),
                           priming = numeric(), 
                           stimulation = numeric())
  # print(parameters)
  for(stm in stimulation.list){
    res <- rmain(parameters = parameters, 
                 variables = variables, 
                 stm = stm, 
                 tmesh = tmesh, 
                 time_interval = time_interval, 
                 time_computation = time_computation)
    res.priming <- rmain(parameters = parameters, 
                         variables = variables.priming, 
                         stm = stm, 
                         tmesh = tmesh, 
                         time_interval = time_interval, 
                         time_computation = time_computation)
    if(res$success & res.priming$success){
      for(tmesh.i in tmesh.list){
        data.model <- rbind(data.model,
                            data.table(time = c(tmesh[tmesh.i], tmesh[tmesh.i]),
                                       m = c(res$output[[tmesh.i]][14],
                                             res.priming$output[[tmesh.i]][14]),
                                       sd =  c(res$output[[tmesh.i]][31],
                                               res.priming$output[[tmesh.i]][31]),
                                       priming = c(0, 1000), 
                                       stimulation = c(stm, stm))
        )
      }
    } else {
      return(list(error = TRUE))
    }
  }
  data.model <- normalization(data.model, background = background)
  data.model <- lmvn(data.model)
  return(list(error = FALSE, data.model = data.model))
}
#### ####
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
    sigmapoints.list <- 
      GetSigmapoints(sigmapoints.parameters = par.all[sigmapoints.parameters.conditions$id.par],
                      sigmapoints.parameters.sd = par.all[sigmapoints.parameters.conditions$id.par.sd],
                     alpha = sigmapoints.conditions$alpha,
                     beta = sigmapoints.conditions$beta,
                     kappa = sigmapoints.conditions$kappa)
    arguments.list <- list()
    for(i in 1:length(sigmapoints.list)){
      sigmapoints.par <- par.all
      sigmapoints.par[sigmapoints.parameters.conditions$id.par] <- sigmapoints.list[[i]]
      
      parameters <- parameters.factor*(parameters.base)^sigmapoints.par
      variables[ 
        sigmapoints.parameters.conditions$variables[
          which(sigmapoints.parameters.conditions$variables != 0)]] <- 
        parameters[sigmapoints.parameters.conditions$id.par]
      variables.priming[ 
        sigmapoints.parameters.conditions$variables.priming[
          which(sigmapoints.parameters.conditions$variables.priming != 0)]] <- 
        parameters[sigmapoints.parameters.conditions$id.par]
      arguments.list[[i]] <- fun_modify_input(parameters = parameters.factor*(parameters.base)^sigmapoints.par,
                                variables = variables,
                                variables.priming = variables.priming,
                                sigmapoints.parameters.conditions = sigmapoints.parameters.conditions,
                                ...)
    }
    return(arguments.list)
}

#### optimisation ####
optimisation_ut <- function(par,
                         fun_run_model = run_model,
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
                         ut.fun_sigmapoints,
                         sigmapoints.parameters.conditions,
                         sigmapoints.conditions,
                         ...)
{
  
  
  
  ### run 
  par.all <- rep(0, times = length(parameters.factor))
  par.all[par.optimised] <- par
  arguments.list  <-  ut.fun_sigmapoints(
    par.all = par.all,
    sigmapoints.parameters.conditions = sigmapoints.parameters.conditions,
    sigmapoints.conditions = sigmapoints.conditions,
    parameters.factor = parameters.factor,
    parameters.base = parameters.base,
    variables = variables,
    variables.priming = variables.priming,
    fun_modify_input = fun_modify_input)
  
  
  ###
  model.simulation.list <- list()
  for(argument.i in 1:length(arguments.list)){
  
    input <- arguments.list[[argument.i]]
    parameters <- input$parameters
    variables  <- input$variables
    variables.priming <- input$variables.priming
    
    model.simulation.list[[argument.i]] <- do.call(fun_run_model,
                                list(
                                  parameters = parameters,
                                  variables = variables,
                                  variables.priming = variables.priming,
                                  tmesh = tmesh,
                                  tmesh.list = tmesh.list,
                                  stimulation.list = stimulation.list,
                                  background = background))
    
    
    if(model.simulation.list[[argument.i]]$error){
      
      return(Inf)
    } else {
      model.simulation.list[[argument.i]]$data.model$ 
    }
    #  print("hurrra")
    
    result <- sum(
      likelihood(fun.likelihood = fun.likelihood,
                 data.model = model.simulation$data.model, 
                 data.exp.grouped = data.exp.grouped,
                 data.exp.summarise = data.exp.summarise))
    print(result)
    
  }
  if(return.model){
    return(c(model.simulation, list(optimisation = result)))
  }
  return(result)
}