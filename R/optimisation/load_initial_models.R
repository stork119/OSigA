### ###
###
### ###


LoadInitialModels <- function(path.optimisation,
                              path.optimisation.data = paste(path.optimisation, "data", sep = "/")){
  
  attach(LoadOptimisationConditions(path.optimisation = path.optimisation,
                                    path.optimisation.data = path.optimisation.data))
  #### default ####
  variables <- rep(0.0, times = 629)
  variables.priming <- rep(0.0, times = 629)
  variables[1:17] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
  variables.priming[1:17] <- scan(file = paste(path.parameters, "var-priming.txt", sep = ""))
  variables[44:52] <- scan(file = paste(path.parameters, "var-extrinsic.txt", sep = ""))
  variables.priming[44:52] <- variables[44:52]
  variables[27:43] <- varscale*(variables[1:17]^2)
  variables.priming[27:43] <- varscale*(variables.priming[1:17]^2)
  
  
  model.simulation.def <- do.call(run_model,
                                  list(parameters = par.def,
                                       variables = variables,
                                       variables.priming = variables.priming,
                                       tmesh = tmesh,
                                       tmesh.list = tmesh.list,
                                       stimulation.list = stimulation.list,
                                       background = background))
  model.simulation.def$data.model$type <- "single"
  model.simulation.def$data.model$likelihood  <- 
    likelihood(data.model =model.simulation.def$data.model,
               data.exp.grouped = data.exp.grouped,
               data.exp.summarise =  data.exp.summarise.optimisation,
               fun.likelihood = fun.likelihood.list.sd_data)
  
  data.model.likelihood <- model.simulation.def$data.model %>% filter(stimulation %in% stimulation.list)
  result.likelihood.list <-
    lapply(fun.likelihood.list, 
           function(fun.likelihood){
             likelihood( 
               fun.likelihood = fun.likelihood,
               data.model  = data.model.likelihood,
               data.exp.summarise = data.exp.summarise.optimisation,
               data.exp.grouped = data.exp.grouped.optimisation)
           })
  result <- sapply(result.likelihood.list, sum)
  
  for(likelihood.name in names(result.likelihood.list)){
    data.model.likelihood[,likelihood.name] <- result.likelihood.list[[likelihood.name]]
  }
  
  
  path.single <- paste(path.optimisation.data, "single", sep = "/")
  dir.create(path.single, recursive = TRUE, showWarnings = FALSE)
  
  write.table(x = data.model.likelihood, file = paste(path.single, "data_model_likelihood.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = TRUE)
  
  
  save_results(path.opt = path.single,
               variables = variables,
               variables.priming = variables.priming,
               data.model.opt = model.simulation.def$data.model,
               optimisation.opt = result,
               optimisation.opt.colnames = TRUE,
               par.opt = par.def, 
               res.list = list(),
               data.exp.grouped = data.list$data.exp.norm,
               grid.ncol = 4,
               grid.nrow = 2)
  
  
  #### double ####
  
  variables <- rep(0.0, times = 629)
  variables.priming <- rep(0.0, times = 629)
  variables[1:17] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
  variables.priming[1:17] <- scan(file = paste(path.parameters, "var-receptors-priming.txt", sep = ""))
  variables[44:52] <- scan(file = paste(path.parameters, "var-extrinsic.txt", sep = ""))
  variables.priming[44:52] <- variables[44:52]
  variables[27:43] <- varscale*(variables[1:17]^2)
  variables.priming[27:43] <- varscale*(variables.priming[1:17]^2)
  
  
  model.simulation.double <- do.call(run_model,
                                     list(parameters = par.def,
                                          variables = variables,
                                          variables.priming = variables.priming,
                                          tmesh = tmesh,
                                          tmesh.list = tmesh.list,
                                          stimulation.list = stimulation.list,
                                          background = background))
  
  model.simulation.double$data.model$type <- "double"
  model.simulation.double$data.model$likelihood  <- 
    likelihood( data.model =model.simulation.double$data.model,
                data.exp.grouped = data.exp.grouped.optimisation,
                data.exp.summarise =  data.exp.summarise.optimisation,
                fun.likelihood = fun.likelihood.list.sd_data)
  
  data.model.double.likelihood <- model.simulation.double$data.model %>% filter(stimulation %in% stimulation.list)
  result.double.likelihood.list <-
    lapply(fun.likelihood.list, 
           function(fun.likelihood){
             likelihood( 
               fun.likelihood = fun.likelihood,
               data.model = data.model.double.likelihood,
               data.exp.summarise = data.exp.summarise.optimisation,
               data.exp.grouped = data.exp.grouped.optimisation)
           })
  result.double <- sapply(result.double.likelihood.list, sum)
  
  for(likelihood.name in names(result.double.likelihood.list)){
    data.model.double.likelihood[,likelihood.name] <- result.double.likelihood.list[[likelihood.name]]
  }
  
  path.receptors <- paste(path.optimisation.data, "receptors", sep = "/")
  dir.create(path.receptors, recursive = TRUE, showWarnings = FALSE)
  
  write.table(x = data.model.double.likelihood, 
              file = paste(path.receptors, "data_model_likelihood.csv", sep = "/"), 
              sep = ",", 
              row.names = FALSE, 
              col.names = TRUE)
  save_results(path.opt = path.receptors,
               variables = variables,
               variables.priming = variables.priming,
               data.model.opt = model.simulation.double$data.model,
               optimisation.opt = result.double,
               optimisation.opt.colnames = TRUE,
               par.opt = par.def, 
               res.list = list(),
               data.exp.grouped = data.list$data.exp.norm,
               grid.ncol = 4,
               grid.nrow = 2)
  
  #### ####
  data.model.list <- list(single = model.simulation.def$data.model, double = model.simulation.double$data.model)
  return(data.model.list)
}


#### ####
InitiOptimisationTable <- function(path.optimisation,
                              path.optimisation.data = paste(path.optimisation, "data", sep = "/")){
  path.single <- paste(path.optimisation.data, "single", sep = "/")
  optimisation.table <- read_optimisation(path = path.single,
                                          id = "single",
                                          names = names(fun.likelihood.list))
  path.receptors <- paste(path.optimisation.data, "receptors", sep = "/")
  optimisation.table <- rbind(optimisation.table,
                              read_optimisation(
                                path = path.receptors,
                                id = "receptors",
                                names(fun.likelihood.list)))
  return(optimisation.table)
}
