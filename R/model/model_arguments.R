### ###
###
### ###

#### ####

PrepareModelArguments.default <- function(parameters, variables, variables.priming){
  return(list(parameters = parameters, variables = variables, variables.priming = variables.priming))
}


PrepareModelArguments.variables_mean <- function(parameters, variables, variables.priming){
  
  variables[1:17]         <- variables[1:17]*parameters[11:27]
  variables.priming[1:17] <- variables.priming[1:17]*parameters[11:27]
  return(list(parameters = parameters, variables = variables, variables.priming = variables.priming))
}

#### ####
PrepareModelArguments.pSTAT <-
  function(parameters,
           parameters.priming = parameters,
           variables,
           variables.priming,
           mm_constant = 339.3174,
           priming_constant = 2.6,
           ...) {
    
    variables[1] <- parameters[6]/mm_constant
    variables.priming[1] <- priming_constant*variables[1]
    return(list(parameters = parameters,
                parameters.priming = parameters.priming,
                variables = variables,
                variables.priming = variables.priming
    ))}


#### ####
PrepareModelArguments.pSTAT_extrinsic <-
  function(parameters,
           parameters.priming = parameters,
           variables,
           variables.priming,
           mm_constant = 339.3174,
           priming_constant = 2.6,
           ...) {
    
    variables[1] <- parameters[6]/mm_constant
    variables.priming[1] <- priming_constant*variables[1]
    
    variables[18] <- parameters[11]
    variables.priming[18] <- parameters[11]
    variables[17+15] <- parameters[12]
    variables.priming[17+15] <- parameters[12]
    
    parameters<- parameters[1:10]
    parameters.priming<- parameters[1:10]
    
    return(list(parameters = parameters,
                parameters.priming = parameters.priming,
                variables = variables,
                variables.priming = variables.priming
    ))}