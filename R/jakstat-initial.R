### ### ###
### jakstat default 
### ### ###

mesh.exp <- unique(data.exp$time)
data.exp$logintensity <- log(data.exp$intensity)

data.exp.unique <- distinct(data.exp, priming, stimulation, time)
data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time)

path.parameters <- "resources/input/"
par.def <- scan(file = paste(path.parameters, "par.txt", sep = ""))
varscale <- 0.15

tmesh <- seq(from = 0, to = 100, by = 5)
tmesh.list <- which(tmesh %in% unique(data.exp$time))

stimulation.list <- unique(data.exp$stimulation)[-1]
background <- mean(data.exp[data.exp$time == 5,]$intensity)
background.var <- var(data.exp[data.exp$time == 5,]$intensity)


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

result <- sapply(fun.likelihood.list, 
                 function(fun.likelihood){
                   sum( likelihood( 
                     fun.likelihood = fun.likelihood,
                     data.model = model.simulation.def$data.model,
                     data.exp.grouped = data.exp.grouped))
                 })


path.single <- paste(path.optimisation, "single", sep = "/")
dir.create(path.single, recursive = TRUE, showWarnings = FALSE)

save_results(path.opt = path.single,
             data.model.opt = model.simulation.def$data.model,
             optimisation.opt = result,
             par.opt = par.def, 
             res.list = list(),
             data.exp.grouped = data.exp.grouped.all)


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

result <- sapply(fun.likelihood.list, 
                 function(fun.likelihood){
                   sum( likelihood( 
                     fun.likelihood = fun.likelihood,
                     data.model = model.simulation.double$data.model,
                     data.exp.grouped = data.exp.grouped))
                 })


path.receptors <- paste(path.optimisation, "receptors", sep = "/")
dir.create(path.receptors, recursive = TRUE, showWarnings = FALSE)

save_results(path.opt = path.receptors,
             data.model.opt = model.simulation.double$data.model,
             optimisation.opt = result,
             par.opt = par.def, 
             res.list = list(),
             data.exp.grouped = data.exp.grouped.all)
