### ###
### optimisation load conditions
### ###

LoadOptimisationConditions <- function(
  optimisation.path,
  optimisation.filename = paste(optimisation.path, "optimisation_conditions.csv", sep = "")){
  optimisation.conditions.dt <- read.table(
    stringsAsFactors = FALSE,
    file = optimisation.filename,
    header = TRUE,
    sep = ",")
  optimisation.conditions <- list()
  optimisation.conditions$maxit <- optimisation.conditions.dt$maxit  
  optimisation.conditions$data.size <- optimisation.conditions.dt$data.size 
  optimisation.conditions$data.opt.size <- optimisation.conditions.dt$data.opt.size
  optimisation.conditions$cross_validation_num <- optimisation.conditions.dt$cross_validation_num
  optimisation.conditions$optimisation <- optimisation.conditions.dt$optimisation
  optimisation.conditions$fun.optimisation.likelihood <- optimisation.conditions.dt$fun.optimisation.likelihood
  optimisation.conditions$fun_run_model <- optimisation.conditions.dt$fun_run_model
  return(optimisation.conditions)
}


# LoadParametersConditions <- function(){
#   
# }

#### ####
InitiateOptimisation <- function(path.list,
                                 maxit.tmp = Inf,
                                 ...){
  
  flog.debug("load_conditions.R InitiateOptimisation", name="logger.optimisation")
  
  variables         <- scan(paste(path.list$optimisation, "variables.csv", sep = "/"))
  variables.priming <- scan(paste(path.list$optimisation, "variables-priming.csv", sep = "/"))
  
  optimisation.conditions <-
    LoadOptimisationConditions(
      optimisation.path = path.list$optimisation)
  
  if(maxit.tmp != Inf){
    optimisation.conditions$maxit <- maxit.tmp
  }
  
  parameters.conditions <- read.table(
    file = paste(path.list$optimisation, "parameters_conditions.csv", sep = ""),
    sep = ",",
    header = TRUE)
  parameters.base <- parameters.conditions$base
  parameters.factor <- parameters.conditions$factor
  par.lower <- parameters.conditions$lower
  par.upper <- parameters.conditions$upper
  par.optimised   <- which(par.lower != par.upper)
  
  stimulation.list <- scan(paste(path.list$optimisation, 
                                 "stimulation_list.txt", sep ="/"))
  
  data.list <- list()
  data.list$data.exp.grouped <- read.table(
    file = paste(path.list$optimisation, "data_exp_grouped.csv", sep = ""),
    sep = ",",
    header = TRUE)
  data.list$data.exp.grouped.optimisation <-
    data.list$data.exp.grouped %>% 
    filter(stimulation %in% stimulation.list)
  
  data.list$data.exp.summarise.optimisation <- 
    data.list$data.exp.grouped.optimisation %>%
    dplyr::group_by(priming, stimulation, time) %>% 
    dplyr::summarise(m.norm = mean(intensity),
                     mean.lmvn = mean(logintensity),
                     sd.norm = var(intensity),
                     sd.lmvn = var(logintensity))
  data.opt.list <- list()
  data.opt.summary.list <- list()
  if(optimisation.conditions$cross_validation_num > 0 ){
    for(i in 1:optimisation.conditions$cross_validation_num){
      path <- paste(path.list$optimisation.conditions, i, sep = "/")
      data.opt.list[[i]] <- 
        read.table(
          file = paste(path, "data_exp_grouped.csv", sep = "/"),
          sep = ",",
          header = TRUE
        ) %>% 
        filter(stimulation %in% stimulation.list)
      data.opt.summary.list[[i]] <- 
        data.opt.list[[i]] %>%
        dplyr::group_by(priming, stimulation, time) %>% 
        dplyr::summarise(m.norm = mean(intensity),
                         mean.lmvn = mean(logintensity),
                         sd.norm = var(intensity),
                         sd.lmvn = var(logintensity))
    }
  } else {
    data.opt.list[[1]] <- data.list$data.exp.grouped.optimisation
  }
  
  lhs.res <- read.table(file = paste(path.list$optimisation, "parameters_list.csv", sep = ""),
                        sep = ",",
                        header = FALSE)
  par.list <- lapply(1:nrow(lhs.res), 
                     function(i){(par.upper[par.optimised] - par.lower[par.optimised])*lhs.res[i,] + par.lower[par.optimised]})
  
  computations.list <- expand.grid(
    data.id = 1:length(data.opt.list),
    par.id =  1:length(par.list))
  
  computations.list$remove <- 0
  
  par.ids <- list.dirs(path.list$optimisation.data, full.names = FALSE, recursive = FALSE)
  par.ids <- as.numeric(par.ids[which(!is.na(as.numeric(par.ids)))])
  
  for(id in par.ids){
    data.ids <- list.dirs(paste(path.list$optimisation.data, id, sep = "/"),
                          full.names = FALSE)
    data.ids <- as.numeric(data.ids[which(!is.na(as.numeric(data.ids)))])
    computations.list <- computations.list %>% 
      dplyr::filter(!(par.id == id & data.id %in% data.ids))
  }
  
  results.list <- list(
    variables = variables,
    variables.priming = variables.priming,
    optimisation.conditions = optimisation.conditions,
    parameters.conditions = parameters.conditions,
    parameters.base = parameters.base,
    parameters.factor = parameters.factor,
    par.lower = par.lower,
    par.upper = par.upper,
    par.optimised = par.optimised,
    stimulation.list = stimulation.list,
    data.list = data.list,
    data.opt.list = data.opt.list,
    data.opt.summary.list = data.opt.summary.list,
    par.list = par.list,
    computations.list =computations.list
  )
  
  return(results.list)
}