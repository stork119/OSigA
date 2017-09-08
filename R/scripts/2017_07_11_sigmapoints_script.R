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
                            fun_modify_input = PrepareModelArguments.ut.multiple,
                parameters.conditions = parameters.conditions,
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


#### model simulation ####
source("R/scripts/2017_07_11_sigmapoints_library.R")
sigmapoints <- LoadSigmapointsConditions(path.optimisation = path.list$optimisation)

variables.model <- scan(paste(path.list$optimisation, "variables.csv", sep = "/"))
variables.priming.model <- scan(paste(path.list$optimisation, "variables-priming.csv", sep = "/"))
##
results.id <- "108"
parameters.df <- read.table(paste(path.list$optimisation.data, results.id, "parameters_conditions.csv", sep = "/"), header = TRUE, sep = ",")
parameters.df <-parameters.df %>% 
  dplyr::mutate(par = opt) %>%
  dplyr::select(-c(opt, init)) %>%
  dplyr::mutate(factor = factor*base^par) %>% 
  dplyr::mutate(par = 0)

## no results 
results.id <- 0 
parameters.df <- parameters.conditions
parameters.df <-parameters.df %>% 
  dplyr::mutate(par = 0)

parameters.df$par[20] <- 2
# variables.model[31] <- 0.5*variables.model[31]
# variables.priming.model[31] <- variables.model[31]
# parameters.df$par[10] <- 2.5
# parameters.df$par[12] <- -0.25
# parameters.df$par[14] <- 0

# parameters.df$factor[8] <- parameters.df$factor[8]/2
# parameters.df$factor[7] <- parameters.df$factor[7]/5
# parameters.df$factor[18] <- 1#parameters.df$factor[18]*5

res <- analyse_model_ut(parameters = parameters.factor,
   variables.model = variables.model,
                 variables.priming.model = variables.priming.model,
                 sigmapoints = sigmapoints,
                 analyse_name = "karol-penalty",
                 model.computations = list(raw = TRUE, priming = TRUE),
                 parameters.df = parameters.df,
                 fun_modify_parameters = PrepareModelParameters.ut,
                 fun_modify_input = PrepareModelArguments.ut.multiple,
                 parameters.conditions = parameters.conditions,
                 fun_parameters_penalty = fun_parameters_penalty_sigmapoints
                 )
compare_distribution(data.exp.grouped.optimisation = data.exp.grouped.optimisation, data.model = res$data.model, analyse_name = res$analyse_name)

res$argumeents.list
