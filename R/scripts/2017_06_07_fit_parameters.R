### ###
### 2017-06-07 Fit variables and parameters
### ###
source("R/optimisation/initialise_optimisation.R")
source("R/model/model_visualisation.R")
source("R/model/model_execution.R")

attach(LoadOptimisationConditions(
  path.optimisation = path.list$optimisation,
  path.optimisation.data = path.list$optimisation.data))
#### analyse_model ####
analyse_model <- function(parameters.model,
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

variables.factor <- rep(x = 1, times = 17) 
variables.factor[15] <- 4

receptors.factor <- 1

hill.factor <- 1.1

parameters.model <- parameters.factor[1:11]
parameters.model[c(1)] <- 0.9*parameters.model[c(1)]
parameters.model[2] <- 1*parameters.model[2]
parameters.model[4] <- 0.9*parameters.model[4]
parameters.model[c(6)] <- hill.factor*0.3*parameters.model[c(6)]

parameters.model[3] <- 0.05*parameters.model[3]
parameters.model[7] <- 25#2*parameters.model[7]
parameters.model[8] <- 8*parameters.model[8]
parameters.model[9] <- 1*parameters.model[9]
parameters.model[c(10)] <- 1*parameters.model[c(10)]
parameters.model[c(7,8,9,10)] <- receptors.factor*parameters.model[c(7,8,9,10)]


parameters.model[c(9)] <- 0.2*parameters.model[c(9)]
parameters.model[c(10)] <- 2.5*parameters.model[c(10)]
parameters.model[c(11)] <- 100000
variables.model <- variables.factor*variables[1:17]
variables.model[1] <- hill.factor*0.9*variables.model[1]# - 10000
variables.priming.model <- variables.factor*variables.priming[1:17]
variables.priming.model[1] <- 2.6*variables.model[1]

parameters.model[c(6)]/variables.model[1]

dm <- analyse_model(parameters.model = parameters.model,
            variables.model  = variables.model,
            variables.priming.model = variables.priming.model,
            plot = TRUE,
            save = TRUE, 
            analyse_name = "best-model-3"
            )
# data.model <- dm.list[[1]]$data.model
# dm.list[[1]]$data.trajectory
#### ####
data.trajectory <- dm$data.trajectory

p1 <- parameters.model[1]
p6 <- parameters.model[6]
p4 <- parameters.model[4]

#library(reshape2)
data.trajectory.wide <- data.trajectory %>% 
  dplyr::mutate(var = paste("y", var, sep = "")) %>% 
  dcast(formula = time + priming + stimulation ~ var, value.var = "m")

data.trajectory.wide <- data.trajectory.wide %>% dplyr::mutate(hill = y1/(p6 + y1))

data.trajectory.wide <- data.trajectory.wide %>% dplyr::mutate(hill_2 = (p1*y16/(2*p4*y13))-1, p6y1 = p6/y1)

data.trajectory.wide <- data.trajectory.wide %>% dplyr::mutate(hill_p6 = ((p1*y16/(2*p4*y13))-1)*y1)

data.trajectory.wide %>% select(priming, stimulation, time, hill_p6) %>% arrange(priming, stimulation, time) %>% dplyr::filter(time ==15)

ggplot(data.trajectory.wide %>% dplyr::filter(time >= 5), aes(x = factor(time), y = hill)) +
  geom_point() + 
  #geom_point(aes(x = factor(time), y = p6y1), color = "red") + 
  facet_grid(priming~stimulation) + 
  do.call(what = theme_jetka, args = plot.args)  +
  ggtitle(paste("Hill function"))


# data.derivatives <- model_trajectory(data.trajectory = data.trajectory)
# g <- plot_trajectories(path = "", data.trajectory = data.derivatives,plot.args = plot.args, plot.args.ggsave = plot.args.ggsave, save = FALSE)
# g[[1]]
#### ####

#### ####
dm$data.trajectory %>%
  dplyr::filter(var == 1, 
                priming == 1000) %>%
  dplyr::mutate(m_per = m/variables.priming.model[1]) %>%
  dplyr::filter(var == 1, 
                stimulation %in% c(0.5,1,5),
                time %in% c(15,30))

dm$data.trajectory %>%
  dplyr::filter(var == 1, 
                priming == 0) %>%
  dplyr::mutate(m_per = m/variables.model[1]) %>%
  dplyr::filter(var == 1, 
                stimulation %in% c(0.5,1,5),
                time %in% c(15,30))

dm.list[[1]]$data.trajectory %>%
  dplyr::filter(var == 1, 
                stimulation %in% c(0.5,1,5),
                time %in% c(15,30)) 


dm.list[[1]]$data.trajectory %>%
  dplyr::filter(var == 16, 
                stimulation %in% c(0.5,1,5),
                time %in% c(15,30)) 


gplot.trajectory.list <- plot_trajectories(path = "",
                                           data.trajectory = dm.list[[1]]$data.trajectory %>% dplyr::filter(var == 1),
                                           plot.args = plot.args,
                                           plot.args.ggsave = plot.args.ggsave,
                                           save = FALSE)

gplot.trajectory.list$y1 + ylim(c(0,variables.priming.model[1]))


ggplot(dm.list[[1]]$data.trajectory %>%
         
         dplyr::filter(var == 1) %>%
         dplyr::group_by(priming, stimulation, time) %>%
         dplyr::summarise(m = sum(m)) %>%
         dplyr::mutate(type = "type"),
       mapping = aes(x = factor(time), y = m, group = type)) +
  geom_point() +
  geom_line() +
  facet_grid(priming ~ stimulation) + 
  do.call(what = theme_jetka, args = append(plot.args, list(theme.text_size = 18))) +
  ggtitle(labels(var.list)[var.i])
#### losowanie danych  ####

no_cores <- 16
registerDoParallel(no_cores)
likelihood.list <- foreach(i = 1:1000) %dopar% {
  data.exp.grouped.optimisation <- get_equal_data(data = data.list$data.exp,
                                            sample_size = 1000)
  data.exp.summarise.optimisation <- 
    data.exp.grouped.optimisation %>% 
    dplyr::group_by(priming,
                    stimulation,
                    time) %>%
    summarise(m.norm = mean(intensity),
              sd.norm   = var(intensity))
  data.exp.summarise.optimisation<-
    data.exp.summarise.optimisation %>%
    dplyr::mutate(mean.lmvn = lmvn.mean(m = m.norm, sd = sd.norm),
                  sd.lmvn = lmvn.sd(m = m.norm, sd = sd.norm))
  
  data.model$likelihood  <- 
    likelihood(data.model = data.model,
               data.exp.grouped = data.exp.grouped.optimisation,
               data.exp.summarise =   data.exp.summarise.optimisation,
               fun.likelihood = fun.likelihood.list$sd_data)
  optimisation.opt <- sum(data.model$likelihood)
  return(optimisation.opt)
}
stopImplicitCluster()

data.sample <- data.frame(likelihood = unlist(likelihood.list), sample = 1:length(likelihood.list))
write.table(x = data.sample,
            file = paste(path.list$optimisation.analysis,
                         "likelihood_data_sampling.csv", sep = "/"),
            col.names = TRUE,
            row.names = FALSE)
            
plot.args.tmp <- plot.args
plot.args.tmp$theme.title_size <- 24
gplot.list[["sampling"]] <- ggplot(data = data.sample, aes(x = likelihood)) + 
  geom_density() +
  do.call(what = theme_jetka, args = plot.args.tmp)  +
  ggtitle("Dependence of models likelihood on experimental data sample")


do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.list$optimisation.analysis, "likelihood_data_sampling.pdf", sep = "/"),
                           plot = gplot.list[["sampling"]])))
#### ####

# dm.list$small$type <- "small"
# dm.list$big$type <- "big"
# dm <- do.call(rbind, dm.list)
# 
# ggplot(dm,
#        mapping = aes(x = factor(time), y = likelihood, color = type)) +
#   geom_point() +
#   facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
#   ggtitle(paste("Likelihood", collapse = " "))
# 
# dm %>% 
#   filter(!(stimulation %in% c(0.01,0.05,1,5)), !(time %in% c(0,90))) %>% 
#   group_by(type) %>% 
#   summarise(likelihood_sum = sum(likelihood))
# 
# ggplot(dm %>% 
#            filter(stimulation != 5)
#          ,
#        mapping = aes(x = factor(time), y = likelihood, color = type)) +
#   geom_point() +
#   facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
#   ggtitle(paste("Likelihood", collapse = " "))

#### Comparison of p7 p8 p9 p10 ####
pars <- randomLHS(1000, 4)
pars.i <- 1
pars
parameters.model <- parameters.factor
variables.model <- variables[1:17]
variables.priming.model <-variables.priming[1:17]
no_cores <- 16
registerDoParallel(no_cores)
likelihood.list <- foreach(pars.i = 1:nrow(pars)) %dopar% {
  pars.factor <- 2^(pars[pars.i,]*2 - 1)
  parameters.model[7] <- pars.factor[1]*parameters.factor[7]
  parameters.model[8] <- pars.factor[1]*parameters.factor[8]
  parameters.model[9] <- pars.factor[1]*parameters.factor[9]
  parameters.model[10] <- pars.factor[1]*parameters.factor[10]
  results <- analyse_model(parameters.model = parameters.model,
                           variables.model  = variables.model,
                           variables.priming.model = variables.priming.model,
                           plot =FALSE, 
                           save = FALSE)
  df <- matrix(c(results$likelihood, pars.factor), nrow = 1)
  return(df)
}
stopImplicitCluster()
df.likelihood <- do.call(rbind,likelihood.list)
colnames(df.likelihood) <- c("likelihood", "p7", "p8", "p9", "p10")
df.likelihood <- df.likelihood %>% data.table() %>% dplyr::arrange(likelihood)

g <- scatterplot_lieklihood(data = df.likelihood, 
                       colnames.list = c("p7", "p8", "p9", "p10"),
                       path.list = path.list, 
                       filename = "likelihood_p7p8p9p10.pdf")
print(g)


##### ####

lhs_variables <- lhs::randomLHS(n = 100, k = 1)
dm.list <- list()
no_cores <- 8
registerDoParallel(no_cores)
dm.list <- foreach(i = 1:100) %dopar% {

  variables.factor <- rep(x = 1, times = 17) 
  variables.factor[1] <- 1
  receptors.factor <- 1
  parameters.model <- parameters.factor[1:10]
  parameters.model[c(1)] <- 0.003*parameters.model[c(1)]
  parameters.model[2] <- 1*parameters.model[2]
  parameters.model[4] <- 0.9*parameters.model[4]
  parameters.model[c(6)] <- 0.0002*parameters.model[c(6)]
  
  parameters.model[7] <- 0.8*parameters.model[7]
  parameters.model[8] <- 2*parameters.model[8]
  parameters.model[9] <- 1*parameters.model[9]
  parameters.model[c(10)] <- 1*parameters.model[c(10)]
  parameters.model[c(7,8,9,10)] <- receptors.factor*parameters.model[c(7,8,9,10)]
  
  
  parameters.model[c(9)] <- 0*parameters.model[c(9)]
  parameters.model[c(10)] <- 0*parameters.model[c(10)]
  
  variables.model <- variables.factor*variables[1:17]
  variables.model[1] <- lhs_variables[i]*variables.model[1] + variables.model[1]# - 10000
  variables.priming.model <- variables.factor*variables.priming[1:17]
  variables.priming.model[1] <- 2.6*variables.model[1]
  
  dm <- analyse_model(parameters.model = parameters.model,
                                variables.model  = variables.model,
                                variables.priming.model = variables.priming.model,
                                plot = TRUE,
                                save = TRUE, 
                                analyse_name = paste("modify_initial_stat1_2017-06-13-", i+100, sep = ""))
  return(dm)
}
stopImplicitCluster()


dm.list.all <- append(dm.list_1, dm.list)

gplot.list[["stat1_initial"]] <- list()
for(i in order(lhs_variables)){
  gplot.list[["stat1_initial"]][[as.character(i)]] <- dm.list.all[[i]]$gplot + ggtitle(lhs_variables[i])
}

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.list$optimisation.analysis, "model_compare_stat1_initial.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = gplot.list[["stat1_initial"]], nrow =1, ncol =1 ))))



###### #####

parameters.model[6]
variables.model[1]


hill_coef <- function(p6,y1){
  return((2.6*(p6+y1))/(p6+(2.6*y1)))
}


hill_grid <- 2^(-1 + lhs::randomLHS(n = 10000, k =2)*2)
hill_grid <- lhs::randomLHS(n = 10000, k =2)

colnames(hill_grid) <- c("p6", "y1")
hill_grid <- hill_grid %>% data.table() 
#hill_grid$p6 <- hill_grid$p6*parameters.model[6]
#hill_grid$y1 <- hill_grid$y1*variables.model[1]
hill_grid <- hill_grid %>% data.table() %>% dplyr::mutate(hill = hill_coef(p6, y1))

ggplot(data = hill_grid, aes(x = log(p6), y = log(y1), color = hill)) + 
  geom_point() + 
  geom_point(data = hill_grid %>% filter(hill > (quantile(hill_grid$hill, probs = 0.95))), aes(x = log(p6), y = log(y1)), color = "red")


ggplot(data = hill_grid, aes(x = p6/y1, y = hill)) + 
  geom_point() + xlim(0,500)
  geom_point(data = hill_grid %>% filter(hill > (quantile(hill_grid$hill, probs = 0.95))), aes(x = log(p6), y = log(y1)), color = "red")