### ###
### 2017-06-07 Fit variables and parameters
### ###


path.optimisation <- paste(path.output, "optimisation/2017-06-07-variables_mean/", sep = "/")
path.optimisation.data <- paste(path.optimisation, "data/", sep = "/")
path.optimisation.results <- paste(path.optimisation, "results/", sep = "/")
path.optimisation. <- paste(path.optimisation, "results/", sep = "/")

path.optimisation.analysis <- paste(path.optimisation, "analysis/", sep = "/")

attach(LoadOptimisationConditions(
  path.optimisation = path.optimisation,
  path.optimisation.data = path.optimisation.data))

dm.list <- list()

variables.factor <- rep(x = 1, times = length(variables.model)) 
variables.factor[1] <- 1.1
parameters.model <- parameters.factor[1:10]
parameters.model[c(1)] <- 1.95*parameters.model[c(1)]
parameters.model[c(6)] <- 1*parameters.model[c(6)]
parameters.model[c(10)] <- 2*parameters.model[c(10)]
variables.model <- variables.factor*variables[1:17]
variables.priming.model <- variables.factor*variables.priming[1:17]

dm.list[["big"]] <- analyse_model(parameters.model = parameters.model,
              variables.model  = variables.model,
              variables.priming.model = variables.priming.model)
#### ####
analyse_model <- function(parameters.model,variables.model,variables.priming.model){
  path <- paste(path.optimisation.analysis, Sys.time(), sep = "/")
  dir.create(path)
  print(path)
  data.model <- 
    run_model_mean(parameters = parameters.model, 
                   variables  = variables.model,
                   variables.priming = variables.priming.model,
                   tmesh = tmesh,
                   tmesh.list = tmesh.list,
                   stimulation.list = stimulation.list,
                   background = background,
                   time_interval = 100,
                   time_computation = 1000*60*5)
  data.model <- data.model$data.model
  
  data.model$likelihood  <- 
    likelihood(data.model = data.model,
               data.exp.grouped = data.exp.grouped.optimisation,
               data.exp.summarise =   data.exp.summarise.optimisation,
               fun.likelihood = fun.likelihood.list$sd_data)
  
  optimisation.opt <- sum(data.model$likelihood)
  print(optimisation.opt)
  gplot <- ggplot(data.model %>%  
                    dplyr::mutate(type = "model"),
                  mapping = aes(x = factor(time), y = mean.lmvn, group = type, color = type)) +
    geom_point() +
    geom_line() +
    facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
    geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
    #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
    geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data"),
                  mapping = aes(ymin = mean.lmvn - sqrt(sd.lmvn), 
                                ymax = mean.lmvn + sqrt(sd.lmvn)),
                  color = "black") +
    ggtitle(paste("Compare", type.list, collapse = " "))
  
  print(gplot)
  
  do.call(what = ggsave,
          args = append(plot.args.ggsave,
                        list(filename = paste(path, "models_compare.pdf", sep = "/"),
                             plot = gplot)))
  
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
  return(data.model)
}


#### ####

dm.list$small$type <- "small"
dm.list$big$type <- "big"
dm <- do.call(rbind, dm.list)

ggplot(dm,
       mapping = aes(x = factor(time), y = likelihood, color = type)) +
  geom_point() +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  ggtitle(paste("Likelihood", collapse = " "))

dm %>% 
  filter(!(stimulation %in% c(0.01,0.05,1,5)), !(time %in% c(0,90))) %>% 
  group_by(type) %>% 
  summarise(likelihood_sum = sum(likelihood))

ggplot(dm %>% 
           filter(stimulation != 5)
         ,
       mapping = aes(x = factor(time), y = likelihood, color = type)) +
  geom_point() +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  ggtitle(paste("Likelihood", collapse = " "))
