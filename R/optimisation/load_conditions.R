### ###
### optimisation load conditions
### ###


LoadOptimisationConditions <- function(path.optimisation,
                                       path.optimisation.data,
                                       maxit.tmp = Inf,
                                       ...){

  optimisation.conditions <- read.table(
    file = paste(path.optimisation, "optimisation_conditions.csv", sep = ""),
    sep = ",",
    header = TRUE, stringsAsFactors = FALSE)
  fun.optimisation.likelihood <- get(optimisation.conditions$fun.optimisation.likelihood)
  fun_run_model <-  get(optimisation.conditions$fun_run_model)
  maxit <- min(optimisation.conditions$maxit, maxit.tmp)
  
  parameters.conditions <- read.table(
    file = paste(path.optimisation, "parameters_conditions.csv", sep = ""),
    sep = ",",
    header = TRUE)
  parameters.base <- parameters.conditions$base
  parameters.factor <- parameters.conditions$factor
  par.lower <- parameters.conditions$lower
  par.upper <- parameters.conditions$upper
  par.optimised   <- which(par.lower != par.upper)
  
  stimulation.list <- scan(paste(path.optimisation, "stimulation_list.txt", sep ="/"))
  data.exp.grouped.optimisation <- data.exp.grouped %>% filter(stimulation %in% stimulation.list)
  
  data.exp.summarise.optimisation <- data.exp.grouped.optimisation %>%
    dplyr::group_by(priming, stimulation, time) %>% 
    dplyr::summarise(m.norm = mean(intensity),
                     mean.lmvn = mean(logintensity),
                     sd.norm = var(intensity),
                     sd.lmvn = var(logintensity))
  
  lhs.res <- read.table(file = paste(path.optimisation, "parameters_list.csv", sep = ""),
                        sep = ",",
                        header = FALSE)
  par.list <- lapply(1:nrow(lhs.res), 
                     function(i){(par.upper[par.optimised] - par.lower[par.optimised])*lhs.res[i,] + par.lower[par.optimised]})
  
  
  ids <- list.dirs(path.optimisation.data, full.names = FALSE)
  ids <- as.numeric(ids[which(!is.na(as.numeric(ids)))])
  par.list.ids <- 1:length(par.list)
  if(length(ids) > 0){
    par.list.ids <- par.list.ids[-which(ids %in% par.list.ids)]
  }
  
  return(list(optimisation.conditions = optimisation.conditions,
              fun.optimisation.likelihood = fun.optimisation.likelihood,
              fun_run_model = fun_run_model,
              maxit = maxit, 
              parameters.conditions = parameters.conditions,
              parameters.base = parameters.base,
              parameters.factor = parameters.factor,
              par.lower = par.lower,
              par.upper = par.upper,
              par.optimised = par.optimised,
              stimulation.list = stimulation.list,
              data.exp.grouped.optimisation = data.exp.grouped.optimisation,
              data.exp.summarise.optimisation = data.exp.summarise.optimisation,
              lhs.res = lhs.res,
              par.list = par.list,
              par.list.ids = par.list.ids
              ))
}