### ###
### sigmapointsscript library
### ###

#### analyse_model_ut ####
analyse_model_ut <- function(variables.model,
                             variables.priming.model,
                             save = TRUE,
                             plot = TRUE,
                             plot.title = "",
                             analyse_name = NULL,
                             fun.likelihood = fun.likelihood.list$sd,
                             sigmapoints,
                             parameters.df,
                             fun_parameters_penalty = NULL,
                             data.list,
                             data.exp.grouped.optimisation = data.list$data.exp.grouped.optimisation,
                             data.exp.summarise.optimisation = data.list$data.exp.summarise.optimisation,
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
  results$path <- path
  par.optimised <- which(parameters.df$lower != parameters.df$upper)
  
  parameters.conditions$factor <- parameters.df$factor*parameters.df$base^parameters.df$par
  
  model<- run_model_ut(
    par = parameters.df$par[par.optimised],
    parameters.base = parameters.df$base,
    parameters.factor = parameters.df$factor,
    variables = variables.model, 
    variables.priming = variables.priming.model, 
    tmesh = tmesh, 
    tmesh.list = tmesh.list,
    stimulation.list = stimulation.list,
    background = background,
    par.optimised = par.optimised,
    sigmapoints = sigmapoints,
    parameters.conditions = parameters.conditions,
    fun_modify_input = fun_modify_input,
    fun_modify_parameters = fun_modify_parameters,
    fun_model_ode =fun_model_ode)
  #  ...)
  print(model$error)
  
  data.model <- model$data.model
  data.trajectory <- model$data.trajectory
  
  # do.call(what = ggsave,
  #         args = append(plot.args.ggsave,
  #                       list(filename = paste(path, "models_trajectory.pdf", sep = "/"),
  #                            plot = marrangeGrob(grobs = gplot.trajectory.list, ncol = 1, nrow = 1))))
  # 
  data.model$likelihood  <- 
    likelihood(data.model = data.model,# %>% dplyr::filter(time %in% tmesh[tmesh.list]),
               data.exp.grouped = data.exp.grouped.optimisation,
               data.exp.summarise =   data.exp.summarise.optimisation,
               fun.likelihood = fun.likelihood)
  
  optimisation.opt <- sum(data.model$likelihood)
  if(!is.null(fun_parameters_penalty)){
    
    optimisation.opt <- optimisation.opt + nrow(data.exp.summarise.optimisation)*
      fun_parameters_penalty(par = parameters.df$par[par.optimised],
    #                         parameters.conditions = parameters.conditions)
                              ...)
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
  
  print(optimisation.opt)
  if(plot){
    
    results[["compare_log_noise"]] <- 
      ggplot(model$data.model.ut %>% mutate(type = "model") ,
             mapping = aes(x = time,
                           y = mean.lmvn.ut,
                           
                           ymin =  mean.lmvn.ut - sqrt(sd.lmvn.ut),
                           ymax =  mean.lmvn.ut + sqrt(sd.lmvn.ut),
                           group = factor(type),
                           color = factor(type),
                           fill = factor(type))) +
      geom_ribbon(alpha = .5) +
      geom_ribbon(model$data.model.ut %>% mutate(type = "model_intrinsic") ,
                mapping = aes(x = time,
                              ymax = mean.lmvn.ut + sqrt(sd_intrinsic.lmvn.ut),
                              ymin = mean.lmvn.ut - sqrt(sd_intrinsic.lmvn.ut),
                              group = factor(type),
                              color = factor(type),
                              fill = factor(type)),
                alpha = .5) +
      geom_ribbon(model$data.model.ut %>% mutate(type = "model_extrinsic") ,
                  mapping = aes(x = time,
                                ymax = mean.lmvn.ut + sqrt(sd_extrinsic.lmvn.ut),
                                ymin = mean.lmvn.ut - sqrt(sd_extrinsic.lmvn.ut),
                                group = factor(type),
                                color = factor(type),
                                fill = factor(type)),
                  alpha = .5) +
      geom_point() +
      geom_line() +
      facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
      geom_point(data =  data.exp.summarise.optimisation %>%
                   mutate(type = "data"),
                 mapping = aes(x = time,
                               y = mean.lmvn,
                               ymin =  mean.lmvn - sqrt(sd.lmvn),
                               ymax =  mean.lmvn + sqrt(sd.lmvn),
                               group = factor(type),
                               color = factor(type))) +
      geom_errorbar(data = data.exp.summarise.optimisation %>% 
                      mutate(type = "data"),
                    mapping = aes(x = time,
                                  y = mean.lmvn,
                                  ymin =  mean.lmvn - sqrt(sd.lmvn),
                                  ymax =  mean.lmvn + sqrt(sd.lmvn),
                                  group = factor(type),
                                  color = factor(type))) +
      ggtitle(paste(plot.title, "compare log", collapse = ""))
    
    print(results[["compare_log"]])
    
    
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
      ggtitle(paste(plot.title, "compare log", collapse = ""))
    
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
      ggtitle(paste(plot.title, "compare log", collapse = ""))
    
    if(save){
      gplot.trajectory.list <- plot_trajectories(path = path,
                                                 data.trajectory = data.trajectory %>% filter(sigmapoint == 1, var <= 34),
                                                 plot.args = plot.args,
                                                 plot.args.ggsave = plot.args.ggsave,
                                                 save = save)
      do.call(what = ggsave,
              args = append(plot.args.ggsave,
                            list(filename = paste(path, "models_compare_log_noise.pdf", sep = "/"),
                                 plot = results[["compare_log_noise"]])))
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
  return(append(results,
                list( analyse_name = analyse_name,
                      likelihood = optimisation.opt,
                      data.model = data.model,
                      argumeents.list = model$arguments.list,
                      data.trajectory = data.trajectory,
                      gplot.trajectory.list = gplot.trajectory.list
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

