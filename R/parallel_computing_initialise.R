### ###
### parallel_computing_initialise
### ###
setwd("~/Documents/modelling/")
source("R/parallel_computing.R")

path.optimisation <- paste(path.output, "optimisation/2017-06-07-variables_mean/", sep = "/")
path.optimisation.data <- paste(path.optimisation, "data/", sep = "/")
path.optimisation.results <- paste(path.optimisation, "results/", sep = "/")
#### initialise ####
dir.create(path.optimisation, recursive = TRUE)
dir.create(path.optimisation.data, recursive = TRUE)
dir.create(path.optimisation.results, recursive = TRUE)

parameters.filename <- paste(path.optimisation, "parameters_conditions.csv", sep = "")
if(file.exists(parameters.filename)){
parameters.conditions <- read.table(
  file = paste(path.optimisation, "parameters_conditions.csv", sep = ""),
  sep = ",",
  header = TRUE)
parameters.base <- parameters.conditions$base
parameters.factor <- parameters.conditions$factor
par.lower <- parameters.conditions$lower
par.upper <- parameters.conditions$upper
} else {
parameters.factor <- par.def
parameters.base <- rep(x = 10, times = length(parameters.factor))
par.lower <- rep(x = -2, times = length(parameters.factor))
par.upper <- rep(x = 2, times = length(parameters.factor))
parameters.conditions <- data.table(
  factor = par.def,
  base   = parameters.base,
  lower  = par.lower,
  upper  = par.upper)

write.table(x = parameters.conditions,
            file = paste(path.optimisation, "parameters_conditions.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
}
par.optimised   <- which(par.lower != par.upper)

lhs.res <- randomLHS(1000, length(par.optimised))

write.table(x = rbind(matrix(0.5 + 0*par.optimised, nrow = 1), lhs.res),
            file = paste(path.optimisation, "parameters_list.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

stimulation.list <- (data.list$data.exp %>%
                       dplyr::filter(stimulation != 0) %>%
                       dplyr::distinct(stimulation))$stimulation

write.table(x = matrix(stimulation.list, ncol = 1),
            file = paste(path.optimisation, "stimulation_list.txt", sep ="/"),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)


optimisation.filename <- paste(path.optimisation, "optimisation_conditions.csv", sep = "")
if(!file.exists(optimisation.filename)){
  write.table(x = data.table(maxit = 1000,
                             fun.optimisation.likelihood = "fun.likelihood.lmvn.data.norm",
                             fun_run_model = "run_model_mean"),
              file = optimisation.filename,
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
}

### prepare data ###
data.exp.grouped <- data.list$data.exp.norm
# data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time) %>% filter(stimulation %in% stimulation.list)
# data.exp.grouped <- get_equal_data(data.exp.grouped)
# data.exp.grouped <- data.exp.grouped %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))
# data.exp.grouped <- data.exp.grouped %>% summarise(intensity = mean(intensity), logintensity = mean(logintensity), intensity_sd = mean(intensity))
# data.exp.grouped <- data.exp.grouped %>% filter(stimulation != 5)
write.table(x = data.exp.grouped,
            file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

#### load ####
# optimisation.conditions <- read.table(
#   file = paste(path.optimisation, "optimisation_conditions.csv", sep = ""),
#   sep = ",",
#   header = TRUE, stringsAsFactors = FALSE)
# fun.optimisation.likelihood <- get(optimisation.conditions$fun.optimisation.likelihood)
# fun_run_model <-  get(optimisation.conditions$fun_run_model)
# stimulation.list <- scan(paste(path.optimisation, "stimulation_list.txt", sep ="/"))
# data.exp.grouped <-  read.table(
#   file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
#   sep = ",",
#   header = TRUE)
# #data.exp.grouped <- data.exp.grouped %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))
# # stimulation.list.all <- unique(data.exp$stimulation)[-1]
# # stimulation.list.all <- stimulation.list.all
# 
# data.exp.summarise <- data.exp.grouped %>%
#   dplyr::group_by(priming, stimulation, time) %>% 
#   dplyr::summarise(m.morm = mean(intensity),
#                    mean.lmvn = mean(logintensity),
#                    sd.norm = var(intensity),
#                    sd.lmvn = var(logintensity))
