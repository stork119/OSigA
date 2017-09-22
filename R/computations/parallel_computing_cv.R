### ###
### parallel computing ###
### ###

#### ####
#### ####
#### ####
# write.table(
#   x = data.frame(name = name,
#                  mse = results[["cmaes"]]$fmin, 
#                  maxit = maxit,
#                  stopfitness = stopfitness),
#   file = paste(output.path, paste("optimisation_", name, ".csv", sep = ""), sep = "/"),
#   sep = ",", 
#   row.names = FALSE, 
#   col.names = TRUE
# )
#### parallel computations_cv ####
run_parallel_computations_cv <- function(
  path.list = list("id" = NULL,
                   "optimisation" = NULL,
                   "optimisation.data" = NULL),
  no_cores = 18,
  stopfitness = -60000,
  #fun.optimisation = pureCMAES,
  #optimisation.res.par = "xmin"
  fun.optimisation = cma_es,
  optimisation.res.par = "par",
  data.model.list,
  fun_modify_input,
  optimisation_procedure = optimisation_ut,
  computations.ids = NULL,
  ...
){
  
  #TODO get path to logfile$path outside
  logfile <- list()
  logfile$name <- paste("optimisation", path.list$id, Sys.time(), sep = "-")
  logfile$path <- "scripts/"
  logfile$filename <- paste(logfile$path, logfile$name, ".log", sep = "")

  InitLogging(filename = logfile$filename)
  #remove(logfile)
  
  flog.info("run_parallel_computations", name ="logger.optimisation")

  ### initialization ###
  dir.create(path.list$optimisation.data, showWarnings = FALSE, recursive = TRUE)
  flog.info("run_parallel_computations optimisation %s",
            path.list$optimisation.data,
            name ="logger.optimisation")
  
  optimisation.initiate <- InitiateOptimisation(
    path.list = path.list,
    #maxit.tmp = maxit.tmp)
    ...)
  #rm(list = labels(optimisation.conditions.toload))
  #attach(optimisation.conditions.toload)
  
  variables <- optimisation.initiate$variables
  variables.priming <- optimisation.initiate$variables.priming
  optimisation.conditions <- optimisation.initiate$optimisation.conditions
  parameters.conditions <- optimisation.initiate$parameters.conditions
  parameters.base <- optimisation.initiate$parameters.base
  parameters.factor <- optimisation.initiate$parameters.factor
  par.lower <- optimisation.initiate$par.lower
  par.upper <- optimisation.initiate$par.upper
  par.optimised <- optimisation.initiate$par.optimised
  stimulation.list <- optimisation.initiate$stimulation.list
  data.list <- optimisation.initiate$data.list
  data.opt.list <- optimisation.initiate$data.opt.list
  data.opt.summary.list <- optimisation.initiate$data.opt.summary.list
  par.list <- optimisation.initiate$par.list
  computations.list <-optimisation.initiate$computations.list
  
  if(!is.null(computations.ids)){
    computations.list <- computations.list[computations.ids,]
  }
  #### ####
  registerDoParallel(no_cores)
  test <- foreach(computation.i = 1:nrow(computations.list),
                  .combine = list,
                  .multicombine = TRUE ) %dopar%
  {
    tryCatch({
      time.list <- list()
      time.list[[1]] <- Sys.time()
      par <- as.numeric(par.list[[computations.list[computation.i,]$par.id]])
      
      
      flog.info("run_parallel_computations 
                par.id %s
                data.set %s",
                computations.list[computation.i,]$par.id,
                computations.list[computation.i,]$data.id,
                name ="logger.optimisation")

        
      
      optimisation.res <- do.call(
        fun.optimisation,
        list(par = par, 
             # fun = optimisation,
             fn = optimisation_procedure,
             control = list(maxit = optimisation.conditions$maxit,
                            stopfitness = stopfitness,
                            diag.sigma = TRUE,
                            diag.eigen = TRUE,
                            diag.pop = TRUE,
                            diag.value = TRUE),
             lower = par.lower[par.optimised],
             upper = par.upper[par.optimised],
             #stopeval = maxit,
             #stopfitness = stopfitness, 
             fun_run_model = optimisation.conditions$fun_run_model,
             variables = variables,
             variables.priming = variables.priming,
             parameters.conditions = parameters.conditions,
             parameters.base = parameters.base,
             parameters.factor = parameters.factor,
             tmesh = tmesh,
             tmesh.list = tmesh.list,
             stimulation.list = stimulation.list,
             background = background,
             data.exp.grouped = 
               data.opt.list[[computations.list[computation.i,]$data.id]],
               #data.opt.list[[computations.list[computation.i,]$data.id]] %>% dplyr::filter(stimulation != 0, time != 0),
             data.exp.summarise =
               data.opt.summary.list[[computations.list[computation.i,]$data.id]],
               #data.opt.summary.list[[computations.list[computation.i,]$data.id]] %>% dplyr::filter(stimulation != 0, time != 0),
             fun.likelihood = fun.optimisation.likelihood,
             par.optimised = par.optimised,
             fun_modify_input = fun_modify_input,
             #sigmapoints = sigmapoints, fun_model_ode = rmainmean))
             ...))
      
      path.optimisation.i <- paste(path.list$optimisation.data, 
                                   computations.list[computation.i,]$par.id,
                                   computations.list[computation.i,]$data.id,
                                   sep = "/")
      dir.create(path.optimisation.i, recursive = TRUE, showWarnings = FALSE)
      
      par.exp.opt <- optimisation.res[[optimisation.res.par]]
      parameters <- parameters.factor
      parameters[par.optimised] <- parameters.factor[par.optimised]*(parameters.base[par.optimised])^par.exp.opt
      
      input <- fun_modify_input(parameters = parameters,
                                variables = variables,
                                variables.priming = variables.priming,
                                parameters.conditions = parameters.conditions)
      
      parameters <- input$parameters
      variables  <- input$variables
      variables.priming <- input$variables.priming
      
      model.simulation <- do.call(optimisation.conditions$fun_run_model,
                                  list(
                                    parameters = parameters,
                                    variables = variables,
                                    variables.priming = variables.priming,
                                    tmesh = tmesh,
                                    tmesh.list = tmesh.list,
                                    stimulation.list = stimulation.list,
                                    background = background,
                                    parameters.base = parameters.base,
                                    parameters.factor = parameters.factor,
                                    par.optimised = par.optimised,
                                    fun_modify_input = fun_modify_input,
                                    par = par.exp.opt,
                                    parameters.conditions = parameters.conditions,
                                    #sigmapoints = sigmapoints, fun_model_ode = "rmain")),
                                    ...))
      
      error <- model.simulation$error
      # if(model.simulation$error){
      #   model.simulation <- do.call(run_model_mean,
      #                               list(parameters = parameters,
      #                                    variables = variables,
      #                                    variables.priming = variables.priming,
      #                                    tmesh = tmesh,
      #                                    tmesh.list = tmesh.list,
      #                                    stimulation.list = stimulation.list,
      #                                    background = background))
      #   
      # }
      data.exp.opt <- 
        data.list$data.exp.grouped.optimisation
      # data.exp.opt <- 
      #   data.list$data.exp.grouped.optimisation  %>% 
      #   dplyr::filter(stimulation != 0, time != 0)
      data.exp.opt.summarise <- 
        data.list$data.exp.summarise.optimisation
      # data.exp.opt.summarise <- 
      #   data.list$data.exp.summarise.optimisation %>% 
      #   dplyr::filter(stimulation != 0, time != 0)
      result <- sapply(fun.likelihood.list, 
                       function(fun.likelihood){
                         sum( likelihood( 
                           fun.likelihood = fun.likelihood,
                           data.model = model.simulation$data.model,
                           data.exp.grouped = data.exp.opt,
                           data.exp.summarise = data.exp.opt.summarise))
                       })
      flog.info("run_parallel_computations 
                par.id %s
                data.set %s
                likelihood %s",
                computations.list[computation.i,]$par.id,
                computations.list[computation.i,]$data.id,
                paste(result, collapse = " "),
                name ="logger.optimisation")
      
      
      model.simulation$data.model$likelihood  <- 
        likelihood(data.model = model.simulation$data.model,
                   data.exp.grouped = data.exp.opt,
                   data.exp.summarise =  data.exp.opt.summarise,
                   fun.likelihood = fun.optimisation.likelihood)
      
      model.simulation$data.model$type <- "optimised"
      #print(sum(model.simulation$data.model$likelihood))
      gplot <- ggplot(model.simulation$data.model ,
                      mapping = aes(x = time,
                                    y = mean.lmvn,
                                    ymin = mean.lmvn - sqrt(sd.lmvn),
                                    ymax = mean.lmvn + sqrt(sd.lmvn),
                                    group = type,
                                    color = type)) +
        geom_point() +
        geom_line() +
        geom_errorbar() +
        facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
        geom_point(data = data.list$data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
        #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
        geom_errorbar(data = data.list$data.exp.summarise.optimisation %>% mutate(type = "data"),
                      color = "black") +
        ggtitle(
          paste("parameter", computations.list[computation.i,]$par.id,
                "data.set", computations.list[computation.i,]$data.id, 
                collapse = " "))
      
      ggsave(filename = paste(path.optimisation.i, "model_compare_variance.pdf", sep = "/"), 
             plot = gplot,
             width = plot.args$width,
             height =plot.args$height,
             useDingbats = plot.args$useDingbats)
      
      data <- do.call(rbind,append(data.model.list, list(optimised = model.simulation$data.model)))
      gplot <- list()
      gplot[[1]] <- ggplot(data = data,
                           mapping = aes(x = factor(time), y = likelihood, color = type)) +
        geom_point() + 
        do.call(theme_jetka, args = plot.args) + 
        facet_grid(priming ~ stimulation)
      
      gplot[[2]] <- ggplot(data = data %>% dplyr::filter(stimulation != 5),
                           mapping = aes(x = factor(time), y = likelihood, color = type)) +
        geom_point() + 
        do.call(theme_jetka, args = plot.args) + 
        facet_grid(priming ~ stimulation)
      ggsave(filename = paste(path.optimisation.i, "likelihood_comparison.pdf", sep = "/"), 
             plot = marrangeGrob(grobs = gplot, ncol = 1, nrow =1),
             width = plot.args$width,
             height =plot.args$height,
             useDingbats = plot.args$useDingbats)
      
      save_results(path.opt = path.optimisation.i,
                   data.model.opt = model.simulation$data.model,
                   par.opt = parameters,
                   par.exp.opt = par.exp.opt,
                   optimisation.opt = result,
                   optimisation.opt.colnames = names(fun.likelihood.list),
                   res.list = data.model.list,
                   data.exp.grouped = data.list$data.exp.grouped.optimisation,
                   error = error,
                   variables = variables,
                   variables.priming = variables.priming,
                   grid.ncol = ceiling(length(stimulation.list)/2))
      
      parameters.conditions.opt <- parameters.conditions
      parameters.conditions.opt$init <- 0
      parameters.conditions.opt$init[par.optimised] <- par
      parameters.conditions.opt$opt <- 0
      parameters.conditions.opt$opt[par.optimised] <- optimisation.res[[optimisation.res.par]]
      write.table(x = parameters.conditions.opt, 
                  file = paste(path.optimisation.i, "parameters_conditions.csv", sep = "/"),
                  sep = ",", 
                  row.names = FALSE, 
                  col.names = TRUE)
      
      time.list[[2]] <- Sys.time()
      
      write.table(x = difftime(time.list[[2]],
                               time.list[[1]], 
                               units = c("secs"))[[1]], 
                  file = paste(path.optimisation.i, "time.txt", sep = "/"),
                  sep = ",", 
                  row.names = FALSE, 
                  col.names = FALSE)
      
      # write.table(
      #   x = data.frame(mse = optimisation.res$fmin,
      #                  maxit = maxit,
      #                  stopfitness = stopfitness),
      #   file = paste(path.optimisation.i, paste("optimisation_steps", ".csv", sep = ""), sep = "/"),
      #   sep = ",",
      #   row.names = FALSE,
      #   col.names = TRUE
      # )
      saveRDS(object = optimisation.res, file = paste(path.optimisation.i, "optimisation.rds", sep = "/"))
      
      return(list(par = parameters))
      
    }, error =function(e){
      flog.error("run_parallel_computations 
                error %s",
                e,
                name ="logger.optimisation")
      return(list(error = TRUE))
    })
  }
  stopImplicitCluster()
}