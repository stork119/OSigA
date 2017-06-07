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