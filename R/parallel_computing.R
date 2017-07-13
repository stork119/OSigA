### ###
### parallel computing ###
### ###
# source("R/libraries.R")
# source("R/sources.R")
# source("R/initialise.R")
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
#### ####
run_parallel_computations <- function(path.optimisation,
                                      path.optimisation.data = paste(path.optimisation, "data", sep = "/"),
                                      no_cores = 18,
                                      stopfitness = -10000000,
                                      #fun.optimisation = pureCMAES,
                                      #optimisation.res.par = "xmin"
                                      fun.optimisation = cma_es,
                                      optimisation.res.par = "par",
                                      data.model.list,
                                      fun_modify_input,
                                      optimisation_procedure = optimisation,
                                      par.list.ids.part = NULL,
                                      ...
                                      ){
                                      
                                      
  ### initialization ###
  dir.create(path.optimisation.data, showWarnings = FALSE, recursive = TRUE)
  print(path.optimisation.data)
  optimisation.conditions.toload <- LoadOptimisationConditions(
    path.optimisation = path.optimisation,
    path.optimisation.data = path.optimisation.data,
    #maxit.tmp = maxit.tmp)
    ...)
  rm(list = labels(optimisation.conditions.toload))
  attach(optimisation.conditions.toload)
  if(!is.null(par.list.ids.part)){
    par.list.ids <- par.list.ids.part
  }
#### ####
    registerDoParallel(no_cores)
    test <- foreach(i = par.list.ids, .combine = list, .multicombine = TRUE ) %dopar%
    {
      tryCatch({
        #print(par.list.ids)
        par <- as.numeric(par.list[[i]])
        print(par)
        optimisation.res <- do.call(
          fun.optimisation,
          list(par = par, 
               # fun = optimisation,
               fn = optimisation_procedure,
               control = list(maxit = maxit,
                              stopfitness = stopfitness,
                              diag.sigma = TRUE,
                              diag.eigen = TRUE,
                              diag.pop = TRUE,
                              diag.value = TRUE),
               lower = par.lower[par.optimised],
               upper = par.upper[par.optimised],
               #stopeval = maxit,
               #stopfitness = stopfitness, 
               fun_run_model = fun_run_model,
               variables = variables,
               variables.priming = variables.priming,
               parameters.conditions = parameters.conditions,
               parameters.base = parameters.base,
               parameters.factor = parameters.factor,
               tmesh = tmesh,
               tmesh.list = tmesh.list,
               stimulation.list = stimulation.list,
               background = background,
               data.exp.grouped = data.exp.grouped.optimisation,
               data.exp.summarise = data.exp.summarise.optimisation,
               fun.likelihood = fun.optimisation.likelihood,
               par.optimised = par.optimised,
               fun_modify_input = fun_modify_input,
               #sigmapoints = sigmapoints))
               ...))
        
        path.optimisation.i <- paste(path.optimisation.data, i, sep = "/")
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
        
        model.simulation <- do.call(fun_run_model,
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
                                      #sigmapoints = sigmapoints)),
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
        result <- sapply(fun.likelihood.list, 
                         function(fun.likelihood){
                           sum( likelihood( 
                             fun.likelihood = fun.likelihood,
                             data.model = model.simulation$data.model,
                             data.exp.grouped = data.exp.grouped.optimisation,
                             data.exp.summarise = data.exp.summarise.optimisation))
                         })
        
        print(paste(c(path.optimisation.i, result, par.exp.opt), sep = " "))
        
        model.simulation$data.model$likelihood  <- 
          likelihood(data.model = model.simulation$data.model,
                     data.exp.grouped = data.exp.grouped.optimisation,
                     data.exp.summarise =  data.exp.summarise.optimisation,
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
          geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
          #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
          geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data"),
                        color = "black") +
          ggtitle(paste(i, collapse = " "))
        
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
                     data.exp.grouped = data.exp.grouped.optimisation,
                     error = error,
                     variables = variables,
                     variables.priming = variables.priming,
                     grid.ncol = ceiling(length(stimulation.list)/2))
        
        
        init <- 0*par.optimised 
        init[par.optimised] <- par
        opt <- 0*par.optimised 
        opt[par.optimised] <- optimisation.res[[optimisation.res.par]]
        write.table(x = data.frame(factor = parameters.factor,
                                   base = parameters.base, 
                                   init = init,
                                   opt = opt,
                                   lower = par.lower,
                                   upper = par.upper), 
                    file = paste(path.optimisation.i, "parameters.csv", sep = "/"),
                    sep = ",", 
                    row.names = FALSE, 
                    col.names = TRUE)
        
        
        
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
        print(e)
        return(list(error = TRUE))
      })
    }
  stopImplicitCluster()
}

