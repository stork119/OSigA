### ###
### 2017-06-07 Fit variables and parameters
### ###
source("R/optimisation/initialise_optimisation.R")
source("R/model/model_visualisation.R")
source("R/model/model_execution.R")

attach(LoadOptimisationConditions(
  path.optimisation = path.list$optimisation,
  path.optimisation.data = path.list$optimisation.data))
no_cores <- 8
#### analyse_model ####
analyse_model <- function(parameters.model,
                          parameters.model.priming = parameters.model,
                          variables.model,
                          variables.priming.model,
                          save = TRUE,
                          plot = TRUE,
                          title = "",
                          analyse_name = NULL){
  if(is.null(analyse_name)){
    analyse_name <- Sys.time()
  }
  print(analyse_name)
  results <- list()
  if(save){
    path <- paste(path.list$optimisation.analysis, analyse_name, sep = "/")
    dir.create(path,recursive = TRUE, showWarnings = FALSE)
    print(path)
  }
  model <- 
    simulate_model(parameters = parameters.model, 
                   parameters.priming = parameters.model.priming, 
                   variables  = variables.model,
                   variables.priming = variables.priming.model,
                   tmesh = tmesh,
                   tmesh.list = tmesh.list,
                   stimulation.list = stimulation.list,
                   background = background,
                   time_interval = 100,
                   time_computation = 1000*60*5, 
                   tmesh.list.tmp = 1:length(tmesh))
  data.model <- model$data.model
  data.trajectory <- model$data.trajectory
  data.derivatives <- model_trajectory(data.trajectory = data.trajectory)
  
  
  
  data.model$likelihood  <- 
    likelihood(data.model = data.model,# %>% dplyr::filter(time %in% tmesh[tmesh.list]),
               data.exp.grouped = data.exp.grouped.optimisation,
               data.exp.summarise =   data.exp.summarise.optimisation,
               fun.likelihood = fun.likelihood.list$sd_data)
  
  optimisation.opt <- sum(data.model$likelihood)
  print(optimisation.opt)
  if(plot){
    gplot <- ggplot(data.model %>%  
                      dplyr::mutate(type = "model"),
                    mapping = aes(x = time, y = log(m.norm), group = type, color = type)) +
      # geom_point() +
      geom_line() +
      facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
      geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
      #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
      geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data"),
                    mapping = aes(ymin = mean.lmvn - sqrt(sd.lmvn), 
                                  ymax = mean.lmvn + sqrt(sd.lmvn)),
                    color = "black") +
      ggtitle(paste(title, collapse = " "))
    
    print(gplot)
    results[["gplot"]] <- gplot
    
    gplot.raw <- ggplot(data.model %>%  
                          dplyr::mutate(type = "model"),
                        mapping = aes(x = factor(time), y = m.norm, group = type, color = type)) +
      geom_point() +
      geom_line() +
      facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
      geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
      #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
      geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data"),
                    mapping = aes(ymin = m.norm - sqrt(sd.norm), 
                                  ymax = m.norm + sqrt(sd.norm)),
                    color = "black") +
      ggtitle(paste(title, collapse = " "))
    
    results[["gplot.raw"]] <- gplot.raw
    # print(gplot.raw)
    if(save){
      do.call(what = ggsave,
              args = append(plot.args.ggsave,
                            list(filename = paste(path, "models_compare.pdf", sep = "/"),
                                 plot = gplot)))
      print("model_compare saved")
      do.call(what = ggsave,
              args = append(plot.args.ggsave,
                            list(filename = paste(path, "models_compare_raw.pdf", sep = "/"),
                                 plot = gplot.raw)))
      print("model_compare_raw saved")
    }
    gplot.trajectory.list <- plot_trajectories(path = path,
                                               data.trajectory = data.trajectory,
                                               plot.args = plot.args,
                                               plot.args.ggsave = plot.args.ggsave,
                                               save = save)
    results[["gplot.trajectory.list"]] <- gplot.trajectory.list
    print("trajectory saved")
    
    gplot.derivatives.list <- plot_trajectories(path = path,
                                                data.trajectory = data.derivatives,
                                                plot.args = plot.args,
                                                plot.args.ggsave = plot.args.ggsave,
                                                save = save,
                                                filename = "derivatives.pdf")
    results[["gplot.derivatives.list"]] <- gplot.derivatives.list
    print("derivatives saved")
  }
  if(save){
    write.table(x = matrix(parameters.model, ncol = 1), 
                file = paste(path, "parameters.csv", sep ="/"),
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)
    
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
                list(likelihood = optimisation.opt,
                     data.model = data.model,
                     data.trajectory = data.trajectory,
                     data.derivatives = data.derivatives)))
}



#### ####


#dm.list <- list()
variables <- scan(paste(path.list$optimisation, "variables.csv", sep = "/"))
variables.priming <- scan(paste(path.list$optimisation, "variables-priming.csv", sep = "/"))

parameters.model <- parameters.factor
variables.model <- variables
variables.priming.model <- variables.priming

arguments <- PrepareModelArguments.pSTAT(parameters = parameters.model,
                                         variables = variables.model,
                                         variables.priming = variables.priming)


analyse_name <- "receptors"

n <- 1024

cond.analysed <- list()
cond.analysed$parameters <- 0*arguments$parameters != 0
cond.analysed$variables <- 0*arguments$variables != 0
cond.analysed$variables.priming <- 0*arguments$variables != 0


cond.sd <- list()
cond.sd$parameters <- 0*arguments$parameters
cond.sd$variables <- 0*arguments$variables
cond.sd$variables.priming <- 0*arguments$variables


parameters.list <- matrix(rep(arguments$parameters, times = n), nrow = n, byrow = TRUE) %>% data.frame()
variables.list  <- matrix(rep(arguments$variables, times = n), nrow = n, byrow = TRUE) %>% data.frame()
variables.priming.list  <- matrix(rep(arguments$variables.priming, times = n), nrow = n, byrow = TRUE) %>% data.frame()


cond.analysed$variables[15] <- TRUE 
cond.sd$variables[15] <- 0.1

cond.analysed$variables.priming[15] <- TRUE 
cond.sd$variables.priming[15] <- 0.1

for(i in which(cond.analysed$variables)){
  variables.list[,i] <- variables.list[,i]*2^(rnorm(n = n, mean = 0, sd = cond.sd$variables[i]))
}

for(i in which(cond.analysed$variables.priming)){
  variables.priming.list[,i] <- variables.priming.list[,i]*2^(rnorm(n = n, mean = 0, sd = cond.sd$variables.priming[i]))
}

for(i in which(cond.analysed$parameters)){
  parameters.list[,i] <- parameters.list[,i]*2^(rnorm(n = n, mean = 0, sd = cond.sd$parameters[i]))
}


registerDoParallel(cores = no_cores)
dm.list <- foreach(i = 1:n) %dopar% {
  
  arguments$variables <- variables.list[i,]
  arguments$variables.priming <- variables.priming.list[i,]
  arguments$parameters <- parameters.list[i,]
  arguments$parameters.priming <- parameters.list[i,]
  
  # arguments$parameters[p.analysed[1:length(p.analysed)]] <- as.numeric(p.list[i,1:length(p.analysed)])
  # arguments$parameters.priming[p.analysed[1:length(p.analysed)]] <- as.numeric(p.list[i,1:length(p.analysed)])
  #arguments$parameters.priming[7] <- arguments$parameters.priming[7]*1.1
  if(sum(arguments$parameters < 0 | arguments$parameters.priming < 0) == 0 & sum(arguments$variables < 0  | arguments$variables.priming < 0) == 0){
    dm <- analyse_model(parameters.model = as.numeric(arguments$parameters),
                        parameters.model.priming =  as.numeric(arguments$parameters.priming),
                        variables.model  = as.numeric(arguments$variables),
                        variables.priming.model = as.numeric(arguments$variables.priming),
                        plot = FALSE,
                        save = FALSE
                        #analyse_name = "best-model-3"
    )
    #dm$data.model$p <- parameters.model[p.analysed]
    dm$data.model$sample <- i
    return(dm$data.model)
  }
}
stopImplicitCluster()


dm <- do.call(what = rbind, args = dm.list)
gplot.model.list <- list()
gplot.model.list$trajectory_log <- ggplot(dm %>%  
                  dplyr::mutate(type = "model"),
                mapping = aes(x = time, y = log(m.norm), group = interaction(type,sample), color = type)) +
  # geom_point() +
  geom_line() +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data", sample = 0), color = "black") +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data",sample = 0),
                mapping = aes(ymin = mean.lmvn - sqrt(sd.lmvn), 
                              ymax = mean.lmvn + sqrt(sd.lmvn)),
                color = "black") + 
  ylim(c(0.9,1.1)*as.numeric((data.exp.summarise.optimisation %>% 
            dplyr::ungroup() %>% 
            dplyr::summarise(
              min = min(log(m.norm)),
              max = max(log(m.norm))))[,c("min","max")]))

gplot.model.list$errorbars_log <- ggplot(dm %>%  
                  dplyr::mutate(type = "model") %>% 
                  dplyr::group_by(priming,stimulation,time, type) %>%
                  dplyr::summarise(mean.lmvn = mean(log(m.norm)), 
                                   sd.lmvn = var(log(m.norm))),
                mapping = aes(
                  x = time,
                  y = mean.lmvn, 
                  ymin = mean.lmvn - sqrt(sd.lmvn), 
                  ymax = mean.lmvn + sqrt(sd.lmvn),
                  color = type)) +
  #geom_point() +
  #geom_line() +
  geom_errorbar() +  
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data", sample = 0), color = "black") +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data",sample = 0),
                mapping = aes(ymin = mean.lmvn - sqrt(sd.lmvn), 
                              ymax = mean.lmvn + sqrt(sd.lmvn)),
                color = "black") + 
  ylim(c(0.9,1.1)*as.numeric((data.exp.summarise.optimisation %>% 
                                dplyr::ungroup() %>% 
                                dplyr::summarise(
                                  min = min(log(m.norm)),
                                  max = max(log(m.norm))))[,c("min","max")]))


gplot.model.list$trajectory <- ggplot(dm %>%  
                                        dplyr::mutate(type = "model"),
                                      mapping = aes(x = time, y = m.norm, group = interaction(type,sample), color = type)) +
  # geom_point() +
  geom_line() +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data", sample = 0), color = "black") +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data",sample = 0),
                mapping = aes(ymin = m.norm - sqrt(sd.norm), 
                              ymax = m.norm + sqrt(sd.norm)),
                color = "black")# + 
# ylim(c(0.75,1.25)*as.numeric((data.exp.summarise.optimisation %>% 
#                               dplyr::ungroup() %>% 
#                               dplyr::summarise(
#                                 min = min(m.norm),
#                                 max = max(m.norm)))[,c("min","max")]))
# 


gplot.model.list$errorbars <- ggplot(dm %>%  
                                           dplyr::mutate(type = "model") %>% 
                                           dplyr::group_by(priming,stimulation,time, type) %>%
                                           dplyr::summarise(mean.norm = mean(m.norm), 
                                                            sd.norm = var(m.norm)),
                                         mapping = aes(
                                           x = time,
                                           y = mean.norm, 
                                           ymin = mean.norm - sqrt(sd.norm), 
                                           ymax = mean.norm + sqrt(sd.norm),
                                           color = type)) +
  #geom_point() +
  #geom_line() +
  geom_errorbar() +  
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  geom_point(data = data.exp.summarise.optimisation %>% mutate(mean.norm = m.norm, type = "data", sample = 0), color = "black") +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(mean.norm = m.norm, type = "data",sample = 0),
                mapping = aes(ymin = m.norm - sqrt(sd.norm), 
                              ymax = m.norm + sqrt(sd.norm)),
                color = "black") 


path <-  paste(path.list$optimisation.analysis, analyse_name, sep = "/")
dir.create(path = path, showWarnings = FALSE, recursive = TRUE)
write.table(x = dm, file = paste(path, "data_model.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = TRUE)

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path, "model_compare.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = gplot.model.list, ncol = 1, nrow = 1))))


