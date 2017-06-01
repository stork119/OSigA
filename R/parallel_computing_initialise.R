### ###
### parallel_computing_initialise
### ###
setwd("~/Documents/modelling/")
source("R/parallel_computing.R")

path.optimisation <- paste(path.output, "optimisation/2017-06-01/", sep = "/")
path.optimisation.data <- paste(path.optimisation, "data/", sep = "/")
dir.create(path.optimisation, recursive = TRUE)

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


write.table(x = data.table(maxit = 1000,
                           fun.optimisation.likelihood = "fun.likelihood.lmvn.data.norm",
                           fun_run_model = "run_model_mean"),
            file = paste(path.optimisation, "optimisation_conditions.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)


#### prepare data ####
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

data.exp.summarise <- data.exp.grouped %>%
  dplyr::group_by(priming, stimulation, time) %>% 
  dplyr::summarise(m.morm = mean(intensity),
                   mean.lmvn = mean(logintensity),
                   sd.norm = var(intensity),
                   sd.lmvn = var(logintensity))

#### ####
optimisation.conditions <- read.table(
  file = paste(path.optimisation, "optimisation_conditions.csv", sep = ""),
  sep = ",",
  header = TRUE, stringsAsFactors = FALSE)
fun.optimisation.likelihood <- get(optimisation.conditions$fun.optimisation.likelihood)
fun_run_model <-  get(optimisation.conditions$fun_run_model)
stimulation.list <- scan(paste(path.optimisation, "stimulation_list.txt", sep ="/"))
data.exp.grouped <-  read.table(
  file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
  sep = ",",
  header = TRUE)
#data.exp.grouped <- data.exp.grouped %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))
# stimulation.list.all <- unique(data.exp$stimulation)[-1]
# stimulation.list.all <- stimulation.list.all
#### default ####
variables <- rep(0.0, times = 629)
variables.priming <- rep(0.0, times = 629)
variables[1:17] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
variables.priming[1:17] <- scan(file = paste(path.parameters, "var-priming.txt", sep = ""))
variables[44:52] <- scan(file = paste(path.parameters, "var-extrinsic.txt", sep = ""))
variables.priming[44:52] <- variables[44:52]
variables[27:43] <- varscale*(variables[1:17]^2)
variables.priming[27:43] <- varscale*(variables.priming[1:17]^2)


model.simulation.def <- do.call(run_model,
                                list(parameters = par.def,
                                     variables = variables,
                                     variables.priming = variables.priming,
                                     tmesh = tmesh,
                                     tmesh.list = tmesh.list,
                                     stimulation.list = stimulation.list,
                                     background = background))
model.simulation.def$data.model$type <- "single"
model.simulation.def$data.model$likelihood  <- 
  likelihood(data.model =model.simulation.def$data.model,
             data.exp.grouped = data.exp.grouped,
             data.exp.summarise =  data.exp.summarise,
             fun.likelihood = fun.optimisation.likelihood)

data.model.likelihood <- model.simulation.def$data.model %>% filter(stimulation %in% stimulation.list)
result.likelihood.list <-
  lapply(fun.likelihood.list, 
                 function(fun.likelihood){
                    likelihood( 
                     fun.likelihood = fun.likelihood,
                     data.model  = data.model.likelihood,
                     data.exp.summarise = data.exp.summarise,
                     data.exp.grouped = data.exp.grouped)
                 })
result <- sapply(result.likelihood.list, sum)

for(likelihood.name in names(result.likelihood.list)){
  data.model.likelihood[,likelihood.name] <- result.likelihood.list[[likelihood.name]]
}


path.single <- paste(path.optimisation.data, "single", sep = "/")
dir.create(path.single, recursive = TRUE, showWarnings = FALSE)

write.table(x = data.model.likelihood, file = paste(path.single, "data_model_likelihood.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = TRUE)


save_results(path.opt = path.single,
             variables = variables,
             variables.priming = variables.priming,
             data.model.opt = model.simulation.def$data.model,
             optimisation.opt = result,
             optimisation.opt.colnames = names(fun.likelihood.list),
             par.opt = par.def, 
             res.list = list(),
             data.exp.grouped = data.exp.grouped.all,
             grid.ncol = 4,
             grid.nrow = 2)


#### double ####

variables <- rep(0.0, times = 629)
variables.priming <- rep(0.0, times = 629)
variables[1:17] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
variables.priming[1:17] <- scan(file = paste(path.parameters, "var-receptors-priming.txt", sep = ""))
variables[44:52] <- scan(file = paste(path.parameters, "var-extrinsic.txt", sep = ""))
variables.priming[44:52] <- variables[44:52]
variables[27:43] <- varscale*(variables[1:17]^2)
variables.priming[27:43] <- varscale*(variables.priming[1:17]^2)


model.simulation.double <- do.call(run_model,
                                   list(parameters = par.def,
                                        variables = variables,
                                        variables.priming = variables.priming,
                                        tmesh = tmesh,
                                        tmesh.list = tmesh.list,
                                        stimulation.list = stimulation.list,
                                        background = background))

model.simulation.double$data.model$type <- "double"
model.simulation.double$data.model$likelihood  <- 
  likelihood(data.model =model.simulation.double$data.model,
             data.exp.grouped = data.exp.grouped,
             data.exp.summarise =  data.exp.summarise,
             fun.likelihood = fun.optimisation.likelihood)

data.model.double.likelihood <- model.simulation.double$data.model %>% filter(stimulation %in% stimulation.list)
result.double.likelihood.list <-
  lapply(fun.likelihood.list, 
         function(fun.likelihood){
           likelihood( 
             fun.likelihood = fun.likelihood,
             data.model = data.model.double.likelihood,
             data.exp.summarise = data.exp.summarise,
             data.exp.grouped = data.exp.grouped)
         })
result.double <- sapply(result.double.likelihood.list, sum)

for(likelihood.name in names(result.double.likelihood.list)){
  data.model.double.likelihood[,likelihood.name] <- result.double.likelihood.list[[likelihood.name]]
}

path.receptors <- paste(path.optimisation.data, "receptors", sep = "/")
dir.create(path.receptors, recursive = TRUE, showWarnings = FALSE)

write.table(x = data.model.double.likelihood, file = paste(path.receptors, "data_model_likelihood.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = TRUE)
save_results(path.opt = path.receptors,
             variables = variables,
             variables.priming = variables.priming,
             data.model.opt = model.simulation.double$data.model,
             optimisation.opt = result.double,
             optimisation.opt.colnames = names(fun.likelihood.list),
             par.opt = par.def, 
             res.list = list(),
             data.exp.grouped = data.exp.grouped.all,
             grid.ncol = 4,
             grid.nrow = 2)

#### ####
data.model.list <- list(single = model.simulation.def$data.model, double = model.simulation.double$data.model)