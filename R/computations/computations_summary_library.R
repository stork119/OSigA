### ###
### computations summary library
### ###

fun.computations <- list()
#### local_summary ####
fun.computations$local_summary <- function(
  data.list,
  path.list = list("id" = NULL,
                   "optimisation" = NULL,
                   "optimisation.data" = NULL),
  par.id = 1,
  optimisation.name = "sd",
  ...){
  results.local <- list()
  
  data.ids <- list.dirs(paste(path.list$optimisation.data, par.id, sep = "/"),
                        full.names = FALSE)
  data.ids <- as.numeric(data.ids[which(!is.na(as.numeric(data.ids)))])
  path.par <- paste(path.list$optimisation.data, par.id, sep = "/") 
  likelihood.df <- data.frame(
    id = numeric(), 
    likelihood = numeric()
  )
  parameters.list <- list()
  parameters.exp.list <- list()
  for(data.id in data.ids){
    path <- paste(path.par, data.id, "optimisation.csv", sep = "/")
    if(file.exists(path)){
      likelihood.df <- rbind(
        likelihood.df,
        data.frame(
          id = 
            data.id,
          likelihood = 
            read.table(file = path, 
                       header = TRUE,
                       sep = ",")[,optimisation.name]))
      
      
    }
    # path <- paste(path.par, data.id, "par.txt", sep = "/")
    # if(file.exists(path)){
    #   parameters.list[[as.character(data.id)]] <- 
    #     read.table(
    #       file = path, 
    #       header = FALSE,
    #       sep = ","
    #     )
    # }
    path <- paste(path.par, data.id, "parameters_conditions.csv", sep = "/")
    if(file.exists(path)){
      #parameters.list[[as.character(data.id)]] <- 
      parameters.conditions <-  read.table(
        file = path, 
        header = TRUE,
        sep = ","
      )
      
      parameters.list[[as.character(data.id)]] <- (parameters.conditions %>% 
                                                     dplyr::mutate(par = factor*(base^opt)))$par
      parameters.exp.list[[as.character(data.id)]] <- parameters.conditions$opt
    }
  }
  
  likelihood.fit.df <- 
    likelihood.df %>%
    dplyr::filter(id == 1) %>%
    dplyr::select(likelihood)
  if(nrow(likelihood.fit.df) != 0){
    colnames(likelihood.fit.df) <- optimisation.name
    write.table(file = paste(path.par, "otimisation_fit.csv", sep = "/"), 
                x = likelihood.fit.df,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
  }
  
  likelihood.cv.df <- 
    likelihood.df %>%
    dplyr::filter(id != 1) %>%
    dplyr::arrange(likelihood)
  colnames(likelihood.cv.df) <- c("id", optimisation.name)
  write.table(file = paste(path.par, "otimisation_cv.csv", sep = "/"), 
              x = likelihood.cv.df,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  
  likelihood.cv.comp.df <- 
    likelihood.df %>%
    dplyr::mutate(type = ifelse(id == 1, "fit", "cv")) %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(likelihood = mean(likelihood)) %>%
    dplyr::arrange(type)
  
  colnames(likelihood.cv.df) <- c("type", optimisation.name)
  write.table(file = paste(path.par, "otimisation_comparison.csv", sep = "/"), 
              x = likelihood.cv.comp.df,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE) 
  
  parameters.df <- 
    t(do.call(cbind, parameters.list)) %>%
    data.frame(row.names = NULL) 
  colnames(parameters.df) <- paste("p", 1:ncol(parameters.df), sep = "")
  parameters.df$id <- as.numeric(labels(parameters.list))
  parameters.df <-
    parameters.df %>% 
    dplyr::left_join(likelihood.df, by = "id")
  write.table(file = paste(path.par, "parameters.csv", sep = "/"), 
              x = parameters.df,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  
  parameters.exp.df <- 
    t(do.call(cbind, parameters.exp.list)) %>%
    data.frame(row.names = NULL) 
  colnames(parameters.exp.df) <- paste("p", 1:ncol(parameters.exp.df), sep = "")
  parameters.exp.df$id <- as.numeric(labels(parameters.exp.list))
  parameters.exp.df <-
    parameters.exp.df %>% 
    dplyr::left_join(likelihood.df, by = "id")
  write.table(file = paste(path.par, "parameters_exp.csv", sep = "/"), 
              x = parameters.exp.df,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  
  likelihood.df$par.id <- par.id
  parameters.df$par.id <- par.id
  parameters.exp.df$par.id <- par.id
  
  return(
    list(
      parameters = parameters.df, 
      parameters.exp = parameters.exp.df,
      likelihood = likelihood.df
    )
  )
}

#### global_summary ####
fun.computations$global_summary <- function(
  path.list = list("id" = NULL,
                   "optimisation" = NULL,
                   "optimisation.data" = NULL),
  data.list,
  optimisation.name = "sd",
  ...)
{
  
  optimisation.conditions <- 
    LoadOptimisationConditions(
      optimisation.path = path.list$optimisation)
  
  par.ids <- list.dirs(path.list$optimisation.data, full.names = FALSE, recursive = FALSE)
  par.ids <- as.numeric(par.ids[which(!is.na(as.numeric(par.ids)))])
  
  local_summary.list <- list()
  likelihood.df <- data.table(id = numeric(),
                              likelihood = numeric(),
                              par.id = numeric())
  for(id in par.ids){
    local_summary.list[[as.character(id)]] <-
      fun.computations$local_summary(path.list = path.list,
                        par.id = id,
                        optimisation.name = optimisation.name)
    if(!is.null(local_summary.list[[as.character(id)]]$likelihood) &&
       sum(colnames(likelihood.df) != colnames(local_summary.list[[as.character(id)]]$likelihood)) == 0
    ){
      likelihood.df <- rbind(
        likelihood.df,
        local_summary.list[[as.character(id)]]$likelihood)
    }
  }
  
  write.table(file = paste(path.list$optimisation.results,
                           "otimisation_fit.csv",
                           sep = "/"), 
              x = 
                likelihood.df %>% 
                dplyr::filter(id == 1) %>%
                dplyr::arrange(likelihood) %>%
                dplyr::select(par.id, likelihood),
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  
  likelihood.df.cv <-
    likelihood.df %>% 
    dplyr::group_by(id) %>%
    dplyr::summarise(likelihood = min(likelihood)) %>%
    dplyr::arrange(id)
  
  write.table(file = paste(path.list$optimisation.results,
                           "otimisation_cv.csv",
                           sep = "/"), 
              x = likelihood.df.cv,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  
  likelihood.cv.comp.df <- 
    likelihood.df.cv %>%
    dplyr::mutate(type = ifelse(id == 1, "fit", "cv")) %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(likelihood = mean(likelihood)) %>%
    dplyr::arrange(type)
  
  colnames(likelihood.cv.comp.df) <- c("type", optimisation.name)
  write.table(file = paste(path.list$optimisation.results,
                           "otimisation_comparison.csv", sep = "/"), 
              x = likelihood.cv.comp.df,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE) 
}

#### plot_model ####
fun.computations$plot_model <- function(data,
                       id){
  gplot.list <- list()
  gplot.list[["models.compare_log"]] <- 
    ggplot(data,
           mapping = aes(x = factor(time), 
                         y = mean.lmvn,
                         group = factor(type),
                         color = factor(type),
                         ymin = mean.lmvn - sqrt(sd.lmvn), 
                         ymax = mean.lmvn + sqrt(sd.lmvn))) +
    geom_point() +
    geom_line() +
    geom_errorbar() +
    facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
    geom_point(data = data.list$data.exp.summarise %>% mutate(type = "data"), color = "black") +
    geom_errorbar(data = data.list$data.exp.summarise %>% mutate(type = "data"),
                  color = "black") +
    ggtitle(paste("Compare", id, collapse = " "))
  
  gplot.list[["models.compare"]] <- ggplot(data = data,
                                                               mapping = aes(x = factor(time), 
                                                                             y = m.norm,
                                                                             group = factor(type),
                                                                             color = factor(type),
                                                                             ymin = m.norm - sqrt(sd.norm), 
                                                                             ymax = m.norm + sqrt(sd.norm))) +
    geom_point() +
    geom_line() +
    geom_errorbar() +
    facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
    geom_point(data = data.list$data.exp.summarise %>% mutate(type = "data"), color = "black") +
    geom_errorbar(data = data.list$data.exp.summarise %>% mutate(type = "data"),
                  color = "black") +
    ggtitle(paste("Compare", id, collapse = " "))
  
  return(gplot.list)
}
#### compare_models ####
fun.computations$compare_models <- function(
  path.list = list("id" = NULL,
                   "optimisation" = NULL,
                   "optimisation.data" = NULL),
  data.list,
  no_cores = 12,
  num = 10,
  ...){
  optimisation.table <- read.table(paste(path.list$optimisation.results, "otimisation_fit.csv", sep = "/" ),
                                   header = TRUE,
                                   sep = ",")
  optimisation.table <- optimisation.table[1:min(num, nrow(optimisation.table)),]
  data.id <- "1"
  registerDoParallel(no_cores)
  data.model.list <- 
    foreach( i = 1:length(optimisation.table$par.id)) %dopar% {
      par.id <- optimisation.table$par.id[i]
      try({
        filename <- paste(path.list$optimisation.data, par.id, data.id, "data_model.csv", sep = "/")
        if(file.exists(filename)){
          data.model <- read.table(
            file = filename,
            sep = ",",
            header = TRUE)
          data.model$type <- par.id
          return( data.model )
        }
      })
    }
  stopImplicitCluster()
  data.model <- do.call(rbind, data.model.list)
  
  gplot.list <- list()
  gplot.list[["models.compare_log"]] <- list()
  gplot.list[["models.compare"]] <- list()
  for(id in optimisation.table$par.id[1:min(num, nrow(optimisation.table))]){
    gplot.list.tmp <- fun.computations$plot_model(data = data.model %>% 
                                   dplyr::filter(type == id),
                                 id = id)
    gplot.list[["models.compare_log"]][[as.character(id)]] <- gplot.list.tmp$models.compare_log
    gplot.list[["models.compare"]][[as.character(id)]] <- gplot.list.tmp$models.compare
  }
  do.call(what = ggsave, 
          args = append(plot.args.ggsave, 
                        list(filename = paste(path.list$optimisation.results, "models_compare_log.pdf", sep = ""),
                             plot = marrangeGrob(grobs = gplot.list[["models.compare_log"]], ncol = 1, nrow = 1))))
  
  do.call(what = ggsave, 
          args = append(plot.args.ggsave, 
                        list(filename = paste(path.list$optimisation.results, "models_compare.pdf", sep = ""),
                             plot = marrangeGrob(grobs = gplot.list[["models.compare"]], ncol = 1, nrow = 1))))
}
#### compare model stochasticity ####
fun.computations$compare_models_cv <- function(
  path.list = list("id" = NULL,
                   "optimisation" = NULL,
                   "optimisation.data" = NULL),
  data.list,
  no_cores = 12,
  num = 10,
  ...){
  optimisation.table <- read.table(paste(path.list$optimisation.results, "otimisation_fit.csv", sep = "/" ),
                                   header = TRUE,
                                   sep = ",")
  optimisation.table <- optimisation.table[1:min(num, nrow(optimisation.table)),]
  registerDoParallel(no_cores)
  data.model.list <- 
    foreach( i = 1:length(optimisation.table$par.id)) %dopar% {
      par.id <- optimisation.table$par.id[i]
      try({
        path <- paste(path.list$optimisation.data, par.id, sep = "/")
        data.id.list <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)
        
        data.model.list.tmp <- foreach( j = 1:length(data.id.list)) %dopar% {    
          data.id <- data.id.list[j]
          filename <- paste(path.list$optimisation.data, par.id, data.id, "data_model.csv", sep = "/")
          if(file.exists(filename)){
            data.model <- read.table(
              file = filename,
              sep = ",",
              header = TRUE)
            data.model$data.id <- data.id
            return(data.model)
          }
        }
        data.model <- do.call(rbind, data.model.list.tmp)
        data.model$par.id <- par.id
        return( data.model )
      })
    }
  stopImplicitCluster()
  data.model <- do.call(rbind, data.model.list)
  data.model <- data.model %>% 
    dplyr::mutate(type = paste(par.id, data.id))
  
  gplot.list <- list()
  gplot.list[["models.compare_log"]] <- list()
  gplot.list[["models.compare"]] <- list()
  for(id in optimisation.table$par.id[1:min(num, nrow(optimisation.table))]){
    gplot.list.tmp <- fun.computations$plot_model(data = data.model %>% 
                                   dplyr::filter(par.id == id),
                                 id = id)
    gplot.list[["models.compare_log"]][[as.character(id)]] <- gplot.list.tmp$models.compare_log
    gplot.list[["models.compare"]][[as.character(id)]] <- gplot.list.tmp$models.compare
  }
  do.call(what = ggsave, 
          args = append(plot.args.ggsave, 
                        list(filename = paste(path.list$optimisation.results, "cv_models_compare_log.pdf", sep = ""),
                             plot = marrangeGrob(grobs = gplot.list[["models.compare_log"]], ncol = 1, nrow = 1))))
  
  do.call(what = ggsave, 
          args = append(plot.args.ggsave, 
                        list(filename = paste(path.list$optimisation.results, "cv_models_compare.pdf", sep = ""),
                             plot = marrangeGrob(grobs = gplot.list[["models.compare"]], ncol = 1, nrow = 1))))

} 

#### compare estimated parameters ####
fun.computations$compare_parameters <- function(
  path.list = list("id" = NULL,
                   "optimisation" = NULL,
                   "optimisation.data" = NULL),
  data.list,
  no_cores = 12,
  num = 10,
  ...){
  optimisation.table <- read.table(paste(path.list$optimisation.results, "otimisation_fit.csv", sep = "/" ),
                                   header = TRUE,
                                   sep = ",")
  
  optimisation.initiate <- InitiateOptimisation(
    path.list = path.list)
  parameters.conditions <- optimisation.initiate$parameters.conditions
  parameters.conditions$variable <- paste("p",
                                          1:nrow(parameters.conditions),
                                          sep ="")
  
  path.ids.list <- list.dirs(path.list$optimisation.data, full.names = FALSE, recursive = FALSE)
  registerDoParallel(no_cores)
  parameters.list <- 
    foreach( i = 1:length(path.ids.list)) %dopar% {
      par.id <- path.ids.list[i]
      try({
        filename <- paste(path.list$optimisation.data, par.id, "parameters_exp.csv", sep = "/")
        if(file.exists(filename)){
          parameters.df <- read.table(
            file = filename,
            sep = ",",
            header = TRUE) %>% 
            dplyr::mutate(par.id = par.id,
                          data.id = id) %>%
            dplyr::select(-id)
          return( parameters.df )
        }
      })
    }
  stopImplicitCluster()
  
  parameters.df <- do.call(rbind, parameters.list)
  parameters.df.melt <- 
    parameters.df %>% 
    melt(id.vars = c("likelihood", "par.id", "data.id")) %>%
    dplyr::left_join(parameters.conditions, by = c("variable"))
  parameters.df.melt$likelihood_round <- round(parameters.df.melt$likelihood/5000)    
  
  gplot.list <- list()
  for(par_id in optimisation.table[1:num,]$par.id){
    gplot.list[[as.character(par_id)]] <- 
      ggplot( 
        parameters.df.melt %>% 
         #dplyr::filter(likelihood < -40000) %>% 
          dplyr::filter(lower != upper) %>%
          dplyr::filter(par.id == par_id),
        aes(y = value,
            x = variable,
            group = interaction(par.id, data.id), 
            colour = factor(likelihood_round))) +
      geom_point() +
      geom_line() + 
      geom_point(parameters.conditions %>% dplyr::filter(lower != upper) %>% dplyr::mutate(par.id = 0, data.id = 0) , mapping = aes(x = variable, y = upper), color = "black") + 
      geom_line(parameters.conditions %>% dplyr::filter(lower != upper) %>% dplyr::mutate(par.id = 0, data.id = 0), mapping = aes(x = variable, y = upper), color = "black") +
      geom_point(parameters.conditions %>% dplyr::filter(lower != upper) %>% dplyr::mutate(par.id = 0, data.id = 0), mapping = aes(x = variable, y = lower), color = "black") + 
      geom_line(parameters.conditions %>% dplyr::filter(lower != upper) %>% dplyr::mutate(par.id = 0, data.id = 0), mapping = aes(x = variable, y = lower), color = "black") +
      do.call(theme_jetka, plot.args) +
      ggtitle(paste("par.id", par_id))
  }
  do.call(what = ggsave, 
          args = append(plot.args.ggsave, 
                        list(filename = paste(path.list$optimisation.results, "parameters_par_id.pdf", sep = ""),
                             plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow = 1))))
  
  gplot.list <- list()
  for(data_id in 1:10){
    gplot.list[[as.character(data_id)]] <- 
      gplot <- ggplot( 
    parameters.df.melt %>% 
      dplyr::filter(likelihood < -30000) %>% 
      dplyr::filter(lower != upper) %>%
      dplyr::filter(data.id == data_id) %>%
      dplyr::arrange(-likelihood_round) %>%
      dplyr::mutate(l = -likelihood_round)  ,
    aes(y = value,
        x = variable,
        group = interaction(par.id, data.id), 
        colour = -likelihood_round)) +
    geom_point(size = 1.5) +
    geom_line( aes_string(alpha = "l", size = "l")) + scale_size("line", range = c(0,2) ) +
   # geom_point(parameters.conditions %>% dplyr::filter(lower != upper) %>% dplyr::mutate(par.id = 0, data.id = 0), mapping = aes(x = variable, y = upper), color = "black") + 
    geom_line(parameters.conditions %>% dplyr::filter(lower != upper) %>% dplyr::mutate(par.id = 0, data.id = 0), mapping = aes(x = variable, y = upper), color = "black") +
  #  geom_point(parameters.conditions %>% dplyr::filter(lower != upper)%>% dplyr::mutate(par.id = 0, data.id = 0) , mapping = aes(x = variable, y = lower), color = "black") + 
    geom_line(parameters.conditions %>% dplyr::filter(lower != upper) %>% dplyr::mutate(par.id = 0, data.id = 0), mapping = aes(x = variable, y = lower), color = "black") +
    do.call(theme_jetka, plot.args) +
    scale_colour_gradient(low = "white", high = "black") +
      ggtitle(paste("data.id", data_id))
  }
  do.call(what = ggsave, 
          args = append(plot.args.ggsave, 
                        list(filename = paste(path.list$optimisation.results, "parameters_data_id.pdf", sep = ""),
                             plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow = 1))))
  
  gplot.list <- list()
  optimisation.table$rownum <- 1:nrow(optimisation.table)
  gplot.list[[1]] <- 
    ggplot(
      data = optimisation.table %>% filter(likelihood < 0),
      mapping = aes(x = factor(rownum),
                    y = log(-likelihood))) + 
    geom_point() + 
    do.call(theme_jetka, plot.args) +
    ylab("log(-Likelihood)") +
    xlab("Position") +
    ggtitle("all samples")
  
  gplot.list[[2]] <- 
    ggplot(
      data = optimisation.table[1:min(num, nrow(optimisation.table)),],
      mapping = aes(x = factor(rownum),
                    y = -likelihood)) + 
    geom_point() + 
    do.call(theme_jetka, plot.args) +
    ylab("log(-Likelihood)") +
    xlab("Position") +
    ggtitle(paste("Best", num, "samples")) 
  
  do.call(what = ggsave, 
          args = append(plot.args.ggsave, 
                        list(filename = paste(path.list$optimisation.results, "parameters_estimation.pdf", sep = ""),
                             plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow = 1))))
}
#### parameters_conditions ####
fun.computations$save_best_conditions <- function(
  path.list = list("id" = NULL,
                   "optimisation" = NULL,
                   "optimisation.data" = NULL),
  data.list,
  no_cores = 12,
  num = 10,
  ...){
  optimisation.table <- 
    read.table(paste(path.list$optimisation.results, "otimisation_fit.csv", sep = "/" ),
               header = TRUE,
               sep = ",")
  par_id <- optimisation.table[1,]$par.id
  data_id <- 1
  parameters.conditions <- read.table(
    file = paste(path.list$optimisation.data, par_id, data_id, "parameters_conditions.csv", sep = "/"),
    header = TRUE, sep = ","
  )
  parameters.conditions <- 
    parameters.conditions %>%
    dplyr::mutate(factor = factor*(base^opt)) %>%
    dplyr::select(-init) %>%
    dplyr::select(-opt)
  write.table(
    x = parameters.conditions,
    file = paste(path.list$optimisation.results,
                 "parameters_conditions.csv",
                 sep = "/"),
    row.names = FALSE,
    col.names = TRUE,
    sep = ",")
  
  file.copy(
    from = paste(path.list$optimisation, "optimisation_conditions.csv", sep = "/"),
    to = paste(path.list$optimisation.results,
                  "optimisation_conditions.csv",
                  sep = "/"))
  file.copy(
    from = paste(path.list$optimisation, "sigmapoints_conditions.csv", sep = "/"),
    to = paste(path.list$optimisation.results,
               "sigmapoints_conditions.csv",
               sep = "/"))
  file.copy(
    from = paste(path.list$optimisation, "sigmapoints_parameters_conditions.csv", sep = "/"),
    to = paste(path.list$optimisation.results,
               "sigmapoints_parameters_conditions.csv",
               sep = "/"))
}



#### run summary ####
fun.computations$run_summary <- 
  function(...){
    fun.computations$global_summary(...)
    fun.computations$compare_models(...)
    fun.computations$compare_models_cv(...)
    fun.computations$compare_parameters(...)
    fun.computations$save_best_conditions(...)
  }
#### TODO ####


#### analyse_model_ut ####
# analyse_model_ut <- function(variables.model,
#                              variables.priming.model,
#                              save = TRUE,
#                              plot = TRUE,
#                              plot.title = "",
#                              analyse_name = NULL,
#                              fun.likelihood = fun.likelihood.list$sd,
#                              sigmapoints,
#                              parameters.df,
#                              fun_parameters_penalty = NULL,
#                              data.list,
#                              data.exp.grouped.optimisation = data.list$data.exp.grouped.optimisation,
#                              data.exp.summarise.optimisation = data.list$data.exp.summarise.optimisation,
#                              ...){
#   if(is.null(analyse_name)){
#     analyse_name <- Sys.time()
#   }
#   print(analyse_name)
#   results <- list()
#   path <- ""
#   if(save){
#     path <- paste(path.list$optimisation.analysis, analyse_name, sep = "/")
#     dir.create(path,recursive = TRUE, showWarnings = FALSE)
#     print(path)
#   }
#   results$path <- path
#   par.optimised <- which(parameters.df$lower != parameters.df$upper)
#   
#   model<- run_model_ut(
#     par = parameters.df$par[par.optimised],
#     parameters.base = parameters.df$base,
#     parameters.factor = parameters.df$factor,
#     variables = variables.model, 
#     variables.priming = variables.priming.model, 
#     tmesh = tmesh, 
#     tmesh.list = tmesh.list,
#     stimulation.list = stimulation.list,
#     background = background,
#     par.optimised = par.optimised,
#     sigmapoints = sigmapoints,
#     #parameters.conditions = parameters.conditions,fun_modify_input = fun_modify_input,fun_modify_parameters= fun_modify_parameters)
#     ...)
#   data.model <- model$data.model
#   data.trajectory <- model$data.trajectory
#   
#   # do.call(what = ggsave,
#   #         args = append(plot.args.ggsave,
#   #                       list(filename = paste(path, "models_trajectory.pdf", sep = "/"),
#   #                            plot = marrangeGrob(grobs = gplot.trajectory.list, ncol = 1, nrow = 1))))
#   # 
#   data.model$likelihood  <- 
#     likelihood(data.model = data.model,# %>% dplyr::filter(time %in% tmesh[tmesh.list]),
#                data.exp.grouped = data.exp.grouped.optimisation,
#                data.exp.summarise =   data.exp.summarise.optimisation,
#                fun.likelihood = fun.likelihood)
#   
#   optimisation.opt <- sum(data.model$likelihood)
#   if(!is.null(fun_parameters_penalty)){
#     
#     optimisation.opt <- optimisation.opt + nrow(data.exp.summarise.optimisation)*
#       fun_parameters_penalty(par = parameters.df$par[par.optimised],
#                              #                         parameters.conditions = parameters.conditions)
#                              ...)
#   }
#   
#   print(optimisation.opt)
#   if(plot){
#     
#     results[["compare_log_noise"]] <- 
#       ggplot(model$data.model.ut %>% mutate(type = "model") ,
#              mapping = aes(x = time,
#                            y = mean.lmvn.ut,
#                            
#                            ymin =  mean.lmvn.ut - sqrt(sd.lmvn.ut),
#                            ymax =  mean.lmvn.ut + sqrt(sd.lmvn.ut),
#                            group = factor(type),
#                            color = factor(type),
#                            fill = factor(type))) +
#       geom_ribbon(alpha = .5) +
#       geom_ribbon(model$data.model.ut %>% mutate(type = "model_intrinsic") ,
#                   mapping = aes(x = time,
#                                 ymax = mean.lmvn.ut + sqrt(sd_intrinsic.lmvn.ut),
#                                 ymin = mean.lmvn.ut - sqrt(sd_intrinsic.lmvn.ut),
#                                 group = factor(type),
#                                 color = factor(type),
#                                 fill = factor(type)),
#                   alpha = .5) +
#       geom_ribbon(model$data.model.ut %>% mutate(type = "model_extrinsic") ,
#                   mapping = aes(x = time,
#                                 ymax = mean.lmvn.ut + sqrt(sd_extrinsic.lmvn.ut),
#                                 ymin = mean.lmvn.ut - sqrt(sd_extrinsic.lmvn.ut),
#                                 group = factor(type),
#                                 color = factor(type),
#                                 fill = factor(type)),
#                   alpha = .5) +
#       geom_point() +
#       geom_line() +
#       facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
#       geom_point(data =  data.exp.summarise.optimisation %>%
#                    mutate(type = "data"),
#                  mapping = aes(x = time,
#                                y = mean.lmvn,
#                                ymin =  mean.lmvn - sqrt(sd.lmvn),
#                                ymax =  mean.lmvn + sqrt(sd.lmvn),
#                                group = factor(type),
#                                color = factor(type))) +
#       geom_errorbar(data = data.exp.summarise.optimisation %>% 
#                       mutate(type = "data"),
#                     mapping = aes(x = time,
#                                   y = mean.lmvn,
#                                   ymin =  mean.lmvn - sqrt(sd.lmvn),
#                                   ymax =  mean.lmvn + sqrt(sd.lmvn),
#                                   group = factor(type),
#                                   color = factor(type))) +
#       ggtitle(paste(plot.title, "compare log", collapse = ""))
#     
#     print(results[["compare_log"]])
#     
#     
#     results[["compare_log"]] <- 
#       ggplot(data.model %>% mutate(type = "model") ,
#              mapping = aes(x = time,
#                            y = mean.lmvn,
#                            ymin =  mean.lmvn - sqrt(sd.lmvn),
#                            ymax =  mean.lmvn + sqrt(sd.lmvn),
#                            group = factor(type),
#                            color = factor(type))) +
#       geom_point() +
#       geom_line() +
#       geom_errorbar() +
#       facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
#       geom_point(data =  data.exp.summarise.optimisation %>%
#                    mutate(type = "data")) +
#       geom_errorbar(data = data.exp.summarise.optimisation %>% 
#                       mutate(type = "data")) +
#       ggtitle(paste(plot.title, "compare log", collapse = ""))
#     
#     print(results[["compare_log"]])
#     
#     results[["compare"]]<-
#       ggplot(data.model %>% mutate(type = "model") ,
#              mapping = aes(x = time,
#                            y = m.norm,
#                            ymin =  m.norm - sqrt(sd.norm),
#                            ymax =  m.norm + sqrt(sd.norm),
#                            group = factor(type),
#                            color = factor(type))) +
#       geom_point() +
#       geom_line() +
#       geom_errorbar() +
#       facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
#       geom_point(data =  data.exp.summarise.optimisation %>%
#                    mutate(type = "data")) +
#       geom_errorbar(data = data.exp.summarise.optimisation %>% 
#                       mutate(type = "data")) +
#       ggtitle(paste(plot.title, "compare log", collapse = ""))
#     
#     if(save){
#       gplot.trajectory.list <- plot_trajectories(path = path,
#                                                  data.trajectory = data.trajectory %>% filter(sigmapoint == 1, var <= 34),
#                                                  plot.args = plot.args,
#                                                  plot.args.ggsave = plot.args.ggsave,
#                                                  
#                                                  save = save)
#       do.call(what = ggsave,
#               args = append(plot.args.ggsave,
#                             list(filename = paste(path, "models_compare_log_noise.pdf", sep = "/"),
#                                  plot = results[["compare_log_noise"]])))
#       do.call(what = ggsave,
#               args = append(plot.args.ggsave,
#                             list(filename = paste(path, "models_compare_log.pdf", sep = "/"),
#                                  plot = results[["compare_log"]])))
#       print("model_compare_log saved")
#       do.call(what = ggsave,
#               args = append(plot.args.ggsave,
#                             list(filename = paste(path, "models_compare_raw.pdf", sep = "/"),
#                                  plot = results[["compare"]])))
#       print("model_compare_raw saved")
#     }
#     
#   }
#   if(save){
#     
#     write.table(x = parameters.df, 
#                 file = paste(path, "parameters-conditions.csv", sep ="/"),
#                 sep = ",",
#                 row.names = FALSE,
#                 col.names = TRUE)
#     
#     write.table(x = matrix(variables.model, ncol = 1), 
#                 file = paste(path, "variables.csv", sep ="/"),
#                 sep = ",",
#                 row.names = FALSE,
#                 col.names = FALSE)
#     write.table(x = matrix(variables.priming.model, ncol = 1), 
#                 file = paste(path, "variables-priming.csv", sep ="/"),
#                 sep = ",",
#                 row.names = FALSE,
#                 col.names = FALSE)
#     
#     write.table(x = matrix(optimisation.opt, nrow = 1), 
#                 file = paste(path, "optimisation.csv", sep ="/"),
#                 sep = ",",
#                 row.names = FALSE)
#     
#     write.table(x = data.model,
#                 file = paste(path, "data_model.csv", sep ="/"),
#                 sep = ",",
#                 row.names = FALSE,
#                 col.names = TRUE)
#   }
#   return(append(results,
#                 list( analyse_name = analyse_name,
#                       likelihood = optimisation.opt,
#                       data.model = data.model,
#                       argumeents.list = model$arguments.list,
#                       data.trajectory = data.trajectory,
#                       gplot.trajectory.list = gplot.trajectory.list
#                 )))
# }
# ####compare_distribution ####
# compare_distribution <- function(data.exp.grouped.optimisation,
#                                  data.model,
#                                  analyse_name = "opt"){
#   data.exp.grouped.optimisation.zscore <-
#     data.exp.grouped.optimisation %>% 
#     data.table() %>% 
#     left_join((data.model %>% data.table()),
#               by = c("priming", "stimulation", "time")) %>% 
#     dplyr::mutate(zscore = (logintensity - mean.lmvn)/sqrt(sd.lmvn))%>%
#     select(priming, stimulation, time, zscore) 
#   
#   #exp.grid <- data.exp.grouped.optimisation.zscore %>% dplyr::distinct(priming, stimulation, time)
#   dt <- data.exp.grouped.optimisation.zscore %>% 
#     dplyr::distinct(priming, stimulation) %>% 
#     mutate(time = "sample") %>% 
#     merge(data.table(zscore  = rnorm(n = 10000)))
#   
#   min(data.exp.grouped.optimisation.zscore$zscore)
#   gplot <- ggplot(data = data.exp.grouped.optimisation.zscore,
#                   mapping = aes( x = zscore, group = factor(time), color = factor(time) )) +
#     facet_grid(priming ~ stimulation) +
#     geom_density() +
#     do.call(theme_jetka, args = plot.args)  +
#     geom_density(data = dt, color = "black") +
#     ylim(c(0,1))
#   
#   plot.args.ggsave.tmp <- plot.args.ggsave
#   plot.args.ggsave.tmp$width <- 48
#   do.call(what = ggsave,
#           args = append(plot.args.ggsave.tmp,
#                         list(filename = paste(path.list$optimisation.analysis, paste(analyse_name, sep = "/"), "density_variance_comp.pdf", sep = "/"),
#                              plot = gplot)))
# }
