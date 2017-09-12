### ###
### unscented transformation functions 
### ###

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
  fun_run_model <- run_model_ut
  return(list(conditions = sigmapoints.conditions,
              parameters.conditions = sigmapoints.parameters.conditions,
              #    optimisation.procedure = optimisation.procedure,
              fun_run_model = fun_run_model
  ))
}

####  GetSigmapointsWeigths ####
GetSigmapointsWeigths <- function(
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

#### GetSigmapoints ####
GetSigmapoints <- function(
  sigmapoints.parameters,
  sigmapoints.parameters.sd,
  alpha = 0.7,
  kappa = 0,
  beta = 2
  
){
  ### ###
  ### function returns values of sigmapoints basing on sigmapoints conditions
  ### ###
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

#### GetSigmapointsArguments ####
GetSigmapointsArguments <- 
  function(sigmapoints.parameters.conditions,
           sigmapoints.conditions,
           parameters.conditions,
           par.all,
           parameters.factor,
           parameters.base,
           variables,
           variables.priming,
           fun_modify_input,
           fun_modify_parameters = function(parameters.model){return(parameters.model)},
           ...
  ){
    ### ### ###
    ### function returns list of argumnets for each of the sigmapoint
    ### output : 
    ### list(
    ###   parameters,
    ###   parameters-priming,
    ###   variables,
    ###   variables-priming)
    ### ### ###
    
    flog.debug("model_ut.R GetSigmapointsArguments", name="logger.optimisation")
    
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
                                              parameters.conditions = parameters.conditions,
                                              ...)
    }
    return(arguments.list)
  }
