# ### ###
# ### 2017-07-11-sigmapoints-script
# ### ###
#### model simulation ####
#setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")
source("R/scripts/2017_07_11_sigmapoints_library.R")

path.list <-
  LoadOptimisationPaths(
    path.output = paste("resources/output/", sep = "/"),
    id.list = list("2017-11-02")
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
parameters.df <- parameters.conditions
parameters.df <-parameters.df %>% 
  dplyr::mutate(par = 0)
#parameters.df$factor[4] <- parameters.df$factor[4]/25 
parameters.df$factor[4] <- parameters.df$factor[4]/25 
parameters.df$factor[22] <- parameters.df$factor[22]/10000

analyse_name = "test"

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
  parameters.conditions = parameters.dfmodel,
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

