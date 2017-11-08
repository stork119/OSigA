#### ####
#### model simulation ####

setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")
source("R/scripts/2017_07_11_sigmapoints_library.R")

path.list <- 
  LoadOptimisationPaths(
    path.output = "resources/output/",
    id = "2017-09-08-density_ver0_mean_1000_55"
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
#parameters.df <- read.table(paste(path.list$optimisation.data, results.id, "parameters_conditions.csv", sep = "/"), header = TRUE, sep = ",")
# parameters.df <-parameters.df %>% 
#   dplyr::mutate(par = opt) %>%
#   dplyr::select(-c(opt, init)) %>%
#   dplyr::mutate(factor = factor*base^par) %>% 
#   dplyr::mutate(par = 0)

## no results 
results.id <- 0 

parameters.df <- parameters.conditions
parameters.df <- read.table(paste(path.list$optimisation.results, "parameters_conditions.csv", sep = "/"), header = TRUE, sep = ",")
parameters.df <-parameters.df %>% 
  dplyr::mutate(par = 0)

analyse_name = "2017-11-08"
# variables.model[31] <- 0.5*variables.model[31]
# variables.priming.model[31] <- variables.model[31]
# parameters.df$par[10] <- 2.5
# parameters.df$par[12] <- -0.25
# parameters.df$par[14] <- 0

# parameters.df$factor[8] <- parameters.df$factor[8]/2
# parameters.df$factor[7] <- parameters.df$factor[7]/5
# parameters.df$factor[18] <- 1#parameters.df$factor[18]*5

res <- analyse_model_ut(parameters = parameters.factor,
                        variables.model = variables,
                        variables.priming.model = variables.priming,
                        sigmapoints = sigmapoints,
                        analyse_name = analyse_name,
                        model.computations = list(raw = TRUE, priming = TRUE),
                        parameters.df = parameters.df,
                        fun_modify_parameters = PrepareModelParameters.ut,
                        fun_modify_input = PrepareModelArguments.ut.multiple,
                        parameters.conditions = parameters.conditions,
                        #plot = FALSE,
                        fun_parameters_penalty = NULL,#fun_parameters_penalty_sigmapoints
                        data.list = data.list
)
compare_distribution(
  data.exp.grouped.optimisation = data.list$data.exp.grouped.optimisation, data.model = res$data.model, analyse_name = res$analyse_name)



###

#### ####
data.model <- res$data.model  %>% 
  dplyr::filter(stimulation == 1) %>%
  dplyr::mutate(intensity = m.norm)

data.model[5,]$intensity <- 177

fill.data <- "#567abc"
color.data <- list()
color.data[["0"]] <- "#000000"
color.data[["1000"]] <- "#175929"
color.model <- list()
color.model[["0"]]  <- "#ea5b18"
color.model[["1000"]]  <- "#7c3a8e"

type <- "nonpriming"
type <- "priming"
gplot.list <- list() 
for(type in c(0,1000)){
  gplot.list[[type]] <- 
    ggplot(data = data.list$data.exp.grouped %>% 
             dplyr::filter(stimulation == 1) %>% 
             dplyr::filter(priming == type), 
           mapping = aes(x = time, 
                         y = intensity, 
                         group = interaction(time, stimulation, priming))) +
    geom_boxplot(fill = fill.data, 
                 color = color.data[[as.character(type)]]) +
    geom_line(data = data.model %>% 
                dplyr::filter(priming == type),
              mapping = aes( x = time, 
                             y = intensity, 
                             group = interaction(priming)),
              color = color.model[[as.character(type)]]
              ) +
    geom_errorbar(data = data.model %>% 
                    dplyr::filter(priming == 0),
                  mapping = aes( x = time, 
                                 ymin = intensity - sqrt(sd.norm),
                                 ymax = intensity + sqrt(sd.norm),
                                 group = interaction(priming)),
                  color = color.model[[as.character(type)]]) +
    do.call(theme_jetka, args = plot.args) +
    ylim(c(0,600)) +
    ggtitle(as.character(type))
}
  
do.call(what = ggsave,
        args = append(list( width = 8,
                            height = 8,
                            useDingbats = plot.args$useDingbats),
                      list(filename = paste(path.list$optimisation.analysis, paste(analyse_name, sep = "/"), "compare_model_to_data.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1))))