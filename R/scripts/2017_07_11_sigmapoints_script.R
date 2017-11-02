### ###
### 2017-07-11-sigmapoints-script
### ###

setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")

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
source("R/optimisation/initialise_optimisation.R")
source("R/scripts/2017_07_11_sigmapoints_library.R")

path.list <-
  LoadOptimisationPaths(
    path.output = paste("resources/output/", sep = "/"),
    id.list = list("2017-10-17-2/0")
  )

gplot.list <- list()
optimisation.initiate <- InitiateOptimisation(
  path.list = path.list)

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

sigmapoints <- LoadSigmapointsConditions(path.optimisation = path.list$optimisation)

##
# results.id <- "108"
# parameters.df <- read.table(paste(path.list$optimisation.data, results.id, "parameters_conditions.csv", sep = "/"), header = TRUE, sep = ",")
# parameters.df <-parameters.df %>% 
#   dplyr::mutate(par = opt) %>%
#   dplyr::select(-c(opt, init)) %>%
#   dplyr::mutate(factor = factor*base^par) %>% 
#   dplyr::mutate(par = 0)

## no results 
par.id <- 1 
data.id <- 1
parameters.df <- read.table(paste(path.list$optimisation.results, "parameters_conditions.csv", sep = "/"), header = TRUE, sep = ",")
#parameters.df <- parameters.conditions
#parameters <-scan('') parameters.df$factor*(parameters.df$base^parameters.df$opt)
#parameters.df$factor <- parameters
parameters.df <-parameters.df %>% 
  dplyr::mutate(par = 0)
#
# parameters.df <- parameters.conditions
# parameters.df <-parameters.df %>% 
#   dplyr::mutate(par = 0)
# parameters <- parameters.factor

# variables.model[31] <- 0.5*variables.model[31]
# variables.priming.model[31] <- variables.model[31]
# parameters.df$par[10] <- 2.5
# parameters.df$par[12] <- -0.25
# parameters.df$par[14] <- 0

# parameters.df$factor[8] <- parameters.df$factor[8]/2
# parameters.df$factor[7] <- parameters.df$factor[7]/5
# parameters.df$factor[18] <- 1#parameters.df$factor[18]*5

analyse_name = "best"

res <- analyse_model_ut(
  fun_model_ode = "rmainmean",
  parameters = parameters.df$factor,
  variables.model = variables,
  variables.priming.model = variables.priming,
  sigmapoints = sigmapoints,
  analyse_name = analyse_name,
  model.computations = list(raw = TRUE, priming = TRUE),
  parameters.df = parameters.df,
  fun_modify_parameters = PrepareModelParameters.ut,
  fun_modify_input = PrepareModelArguments.ut.multiple,
  parameters.conditions = parameters.df,
  #plot = FALSE,
  fun_parameters_penalty = NULL,#fun_parameters_penalty_sigmapoints
  data.list = data.list,
  data.exp.grouped.optimisation = data.list$data.exp.grouped.optimisation %>% dplyr::filter(stimulation != 0, time != 0),
  data.exp.summarise.optimisation = data.list$data.exp.summarise.optimisation %>% dplyr::filter(stimulation != 0, time != 0)
  )

compare_distribution(
  data.exp.grouped.optimisation = data.list$data.exp.grouped.optimisation, 
  data.model = res$data.model, 
  analyse_name = res$analyse_name)
#### ####
data.model <- res$data.model  %>% 
  dplyr::filter(stimulation == 1) %>%
  dplyr::mutate(intensity = m.norm)

data.model[5,]$intensity <- 177

gplot <- 
  ggplot(data = data.list$data.exp.grouped %>% 
       dplyr::filter(stimulation == 1), 
       mapping = aes(x = time, 
                     y = intensity, 
                     group = interaction(time, stimulation, priming),
                     fill = factor(priming))) +
  geom_boxplot() +
  geom_line(data = data.model,
            mapping = aes( x = time, 
                           y = intensity, 
                           group = interaction(priming),
                           color = factor(priming))) +
  geom_errorbar(data = data.model,
            mapping = aes( x = time, 
                           y = intensity, 
                           ymin = intensity - sqrt(sd.norm),
                           ymax = intensity + sqrt(sd.norm),
                           
                           group = interaction(priming),
                           color = factor(priming))) +
  do.call(theme_jetka, args = plot.args) +
  facet_grid(priming~.) +
  ylim(c(0,1000))

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.list$optimisation.analysis, paste(analyse_name, sep = "/"), "intensity_stm_1.pdf", sep = "/"),
                           plot = gplot)))

