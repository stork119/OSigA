### ###
### 2017-07-11-sigmapoints-script
### ###

setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")
source("R/sum_data_initialise.R")
source("R/model/ut/model_ut.R")

gplot.list <- list()
optimisation.conditions.toload <-
  LoadOptimisationConditions(path.optimisation = path.list$optimisation,
                             path.optimisation.data = path.list$optimisation.data,
                             maxit.tmp = Inf)
rm(list = labels(optimisation.conditions.toload))
attach(optimisation.conditions.toload)


sigmapoints <- LoadSigmapointsConditions(path.optimisation = path.list$optimisation)


optimisation.results <- 
  optimisation_ut(par = as.numeric(par.list[[1]]),
                fun_run_model = run_model_ut,
                parameters.base = parameters.base,
                parameters.factor = parameters.factor,
                variables = variables, 
                variables.priming = variables.priming, 
                tmesh = tmesh, 
                tmesh.list = tmesh.list,
                stimulation.list = stimulation.list,
                            background = background,
                            data.exp.grouped = data.exp.grouped.optimisation,
                            data.exp.summarise = data.exp.summarise.optimisation,
                            return.model = TRUE,
                            fun.likelihood = fun.likelihood.list.sd,
                            par.optimised = par.optimised,
                            fun_modify_input = PrepareModelArguments.ut,
                            fun_modify_parameters = PrepareModelParameters.ut,
                            sigmapoints = sigmapoints)





gplot <- ggplot(optimisation.results$data.model.ut %>% mutate(type = "model") ,
                mapping = aes(x = time,
                              y = mean.lmvn.ut,
                              ymin =  mean.lmvn.ut - sqrt(sd.lmvn.ut),
                              ymax =  mean.lmvn.ut + sqrt(sd.lmvn.ut),
                              group = type,
                              color = type)) +
  geom_point() +
  geom_line() +
  geom_errorbar(color = "blue") +
  geom_errorbar(mapping = aes(x = time,
                              y = mean.lmvn.ut,
                              ymin =  mean.lmvn.ut - sqrt(sd_extrinsic.lmvn.ut),
                              ymax =  mean.lmvn.ut + sqrt(sd_extrinsic.lmvn.ut),
                              group = type), color = "red") +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  geom_point(data =  data.exp.summarise.optimisation %>% mutate(mean.lmvn.ut = mean.lmvn,sd.lmvn.ut = sd.lmvn, type = "data")) +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(mean.lmvn.ut = mean.lmvn,sd.lmvn.ut = sd.lmvn, type = "data"),
                color = "black")



gplot <- ggplot(optimisation.results$data.model %>% mutate(type = factor(sigmapoint)) ,
                mapping = aes(x = time,
                              y = mean.lmvn,
                              ymin =  mean.lmvn - sqrt(sd.lmvn),
                              ymax =  mean.lmvn + sqrt(sd.lmvn),
                              group = type,
                              color = type)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  geom_point(data =  data.exp.summarise.optimisation %>% mutate(type = "data")) +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  geom_errorbar(data =  data.exp.summarise.optimisation %>% mutate(type = "data")),
                color = "black")

#### analyse_model_ut ####
analyse_model_ut <- function(variables.model,
                             variables.priming.model,
                             save = TRUE,
                             plot = TRUE,
                             title = "",
                             analyse_name = NULL,
                             fun.likelihood = fun.likelihood.list$sd,
                             sigmapoints,
                             parameters.df,
                             ...){
  if(is.null(analyse_name)){
    analyse_name <- Sys.time()
  }
  print(analyse_name)
  results <- list()
  path <- ""
  if(save){
    path <- paste(path.list$optimisation.analysis, analyse_name, sep = "/")
    dir.create(path,recursive = TRUE, showWarnings = FALSE)
    print(path)
  }
  par.optimised <- which(parameters.df$lower != parameters.df$upper)
  
  model<- run_model_ut(par = parameters.df$par[par.optimised],
                       parameters.base = parameters.df$base,
                       parameters.factor = parameters.df$factor,
                       variables = variables.model, 
                       variables.priming = variables.priming.model, 
                       tmesh = tmesh, 
                       tmesh.list = tmesh.list,
                       stimulation.list = stimulation.list,
                       background = background,
                       par.optimised = par.optimised,
                       fun_modify_input = PrepareModelArguments.ut,
                       sigmapoints = sigmapoints,
                       ...)
  data.model <- model$data.model
  
  data.model$likelihood  <- 
    likelihood(data.model = data.model,# %>% dplyr::filter(time %in% tmesh[tmesh.list]),
               data.exp.grouped = data.exp.grouped.optimisation,
               data.exp.summarise =   data.exp.summarise.optimisation,
               fun.likelihood = fun.likelihood)
  
  optimisation.opt <- sum(data.model$likelihood)
  print(optimisation.opt)
  if(plot){
    results[["compare_log"]] <- 
      ggplot(data.model %>% mutate(type = "model") ,
             mapping = aes(x = time,
                           y = mean.lmvn,
                           ymin =  mean.lmvn - sqrt(sd.lmvn),
                           ymax =  mean.lmvn + sqrt(sd.lmvn),
                           group = factor(type),
                           color = factor(type))) +
      geom_point() +
      geom_line() +
      geom_errorbar() +
      facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
      geom_point(data =  data.exp.summarise.optimisation %>%
                   mutate(type = "data")) +
      geom_errorbar(data = data.exp.summarise.optimisation %>% 
                      mutate(type = "data")) +
      ggtitle(paste(title, "compare log", collapse = ""))
    
    print(results[["compare_log"]])
    
    results[["compare"]]<-
      ggplot(data.model %>% mutate(type = "model") ,
             mapping = aes(x = time,
                           y = m.norm,
                           ymin =  m.norm - sqrt(sd.norm),
                           ymax =  m.norm + sqrt(sd.norm),
                           group = factor(type),
                           color = factor(type))) +
      geom_point() +
      geom_line() +
      geom_errorbar() +
      facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
      geom_point(data =  data.exp.summarise.optimisation %>%
                   mutate(type = "data")) +
      geom_errorbar(data = data.exp.summarise.optimisation %>% 
                      mutate(type = "data")) +
      ggtitle(paste(title, "compare log", collapse = ""))
    
    if(save){
      do.call(what = ggsave,
              args = append(plot.args.ggsave,
                            list(filename = paste(path, "models_compare_log.pdf", sep = "/"),
                                 plot = results[["compare_log"]])))
      print("model_compare_log saved")
      do.call(what = ggsave,
              args = append(plot.args.ggsave,
                            list(filename = paste(path, "models_compare_raw.pdf", sep = "/"),
                                 plot = results[["compare"]])))
      print("model_compare_raw saved")
    }
    
  }
  if(save){
    
    write.table(x = parameters.df, 
                file = paste(path, "parameters-conditions.csv", sep ="/"),
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
    
    write.table(x = matrix(variables.model, ncol = 1), 
                file = paste(path, "variables.csv", sep ="/"),
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)
    write.table(x = matrix(variables.priming.model, ncol = 1), 
                file = paste(path, "variables-priming.csv", sep ="/"),
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)
    
    write.table(x = matrix(optimisation.opt, nrow = 1), 
                file = paste(path, "optimisation.csv", sep ="/"),
                sep = ",",
                row.names = FALSE)
    
    write.table(x = data.model,
                file = paste(path, "data_model.csv", sep ="/"),
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
  }
  return(append(results,
                list( analyse_name = analyse_name,
                    likelihood = optimisation.opt,
                     data.model = data.model,
                     argumeents.list = model$arguments.list
                )))
}
####compare_distribution ####
compare_distribution <- function(data.exp.grouped.optimisation,
                                 data.model,
                                 analyse_name = "opt"){
  data.exp.grouped.optimisation.zscore <-
    data.exp.grouped.optimisation %>% 
    data.table() %>% 
    left_join((data.model %>% data.table()),
              by = c("priming", "stimulation", "time")) %>% 
    dplyr::mutate(zscore = (logintensity - mean.lmvn)/sqrt(sd.lmvn))%>%
    select(priming, stimulation, time, zscore) 
  
  #exp.grid <- data.exp.grouped.optimisation.zscore %>% dplyr::distinct(priming, stimulation, time)
  dt <- data.exp.grouped.optimisation.zscore %>% 
    dplyr::distinct(priming, stimulation) %>% 
    mutate(time = "sample") %>% 
    merge(data.table(zscore  = rnorm(n = 10000)))
  
  min(data.exp.grouped.optimisation.zscore$zscore)
  gplot <- ggplot(data = data.exp.grouped.optimisation.zscore,
                  mapping = aes( x = zscore, group = factor(time), color = factor(time) )) +
    facet_grid(priming ~ stimulation) +
    geom_density() +
    do.call(theme_jetka, args = plot.args)  +
    geom_density(data = dt, color = "black") +
    ylim(c(0,1))
  
  plot.args.ggsave.tmp <- plot.args.ggsave
  plot.args.ggsave.tmp$width <- 48
  do.call(what = ggsave,
          args = append(plot.args.ggsave.tmp,
                        list(filename = paste(path.list$optimisation.analysis, paste(analyse_name, sep = "/"), "density_variance_comp.pdf", sep = "/"),
                             plot = gplot)))
}

#### model simulation ####
sigmapoints <- LoadSigmapointsConditions(path.optimisation = path.list$optimisation)

variables.model <- scan(paste(path.list$optimisation, "variables.csv", sep = "/"))
variables.priming.model <- scan(paste(path.list$optimisation, "variables-priming.csv", sep = "/"))
results.id <- "1"
parameters.df <- read.table(paste(path.list$optimisation.data, results.id, "parameters.csv", sep = "/"), header = TRUE, sep = ",")

parameters.df <-parameters.df %>% 
  dplyr::mutate(par = opt) %>%
  dplyr::select(-c(opt, init)) %>%
  dplyr::mutate(factor = factor*base^par) %>% 
  dplyr::mutate(par = 0)

variables.model[31] <- 0.5*variables.model[31]
variables.priming.model[31] <- variables.model[31]
parameters.df$par[10] <- 2.5
parameters.df$par[12] <- -0.25
parameters.df$par[14] <- 0

res <- analyse_model_ut(variables.model = variables.model,
                 variables.priming.model = variables.priming.model,
                 sigmapoints = sigmapoints,
                 analyse_name = "karol",
                 model.computations = list(raw = TRUE, priming = TRUE),
                 parameters.df = parameters.df,
                 fun_modify_parameters = PrepareModelParameters.ut
                 )
compare_distribution(data.exp.grouped.optimisation = data.exp.grouped.optimisation, data.model = res$data.model, analyse_name = res$analyse_name)


