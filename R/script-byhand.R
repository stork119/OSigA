#### ####
path.optimisation <- paste(path.output, "tests/2017-02-04/", sep = "/")
path.optimisation.data <- paste(path.optimisation, "initial/", sep = "/")


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
data.exp.grouped <- get_equal_data(data.exp.grouped)
data.exp.grouped <- data.exp.grouped %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))
stimulation.list.all <- unique(data.exp$stimulation)[-1]
par <- scan(file = paste(path.parameters, "par.txt", sep = ""))

id <- 6
path.optimisation.data <- paste(path.optimisation, paste("phosphorylation_rate", id, sep = "-"),"/", sep = "/")
dir.create(path.optimisation.data, recursive = TRUE)
par[1] <- 17500
par[6] <- 2500000
par[7] <- 0.0475
par[8] <- 0.005
par[10] <- 0.005
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
                                list(parameters = par,
                                     variables = variables,
                                     variables.priming = variables.priming,
                                     tmesh = tmesh,
                                     tmesh.list = tmesh.list,
                                     stimulation.list = stimulation.list.all,
                                     background = background))

result <- sapply(fun.likelihood.list, 
                 function(fun.likelihood){
                   sum( likelihood( 
                     fun.likelihood = fun.likelihood,
                     data.model = model.simulation.def$data.model %>% filter(stimulation %in% stimulation.list),
                     data.exp.grouped = data.exp.grouped))
                 })


path.single <- paste(path.optimisation.data, "single", sep = "/")
dir.create(path.single, recursive = TRUE, showWarnings = FALSE)

save_results(path.opt = path.single,
             variables = variables,
             variables.priming = variables.priming,
             data.model.opt = model.simulation.def$data.model,
             optimisation.opt = result,
             par.opt = par, 
             res.list = list(),
             data.exp.grouped = data.exp.grouped.all,
             grid.ncol = 4,
             grid.nrow = 2)


#### double ####
# variables <- rep(0.0, times = 629)
# variables.priming <- rep(0.0, times = 629)
# variables[1:17] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
# variables.priming[1:17] <- scan(file = paste(path.parameters, "var-receptors-priming.txt", sep = ""))
# variables[44:52] <- scan(file = paste(path.parameters, "var-extrinsic.txt", sep = ""))
# variables.priming[44:52] <- variables[44:52]
# variables[27:43] <- varscale*(variables[1:17]^2)
# variables.priming[27:43] <- varscale*(variables.priming[1:17]^2)
# 
# 
# model.simulation.double <- do.call(run_model,
#                                    list(parameters = par.def,
#                                         variables = variables,
#                                         variables.priming = variables.priming,
#                                         tmesh = tmesh,
#                                         tmesh.list = tmesh.list,
#                                         stimulation.list = stimulation.list.all,
#                                         background = background))
# 
# result.double <- sapply(fun.likelihood.list,
#                         function(fun.likelihood){
#                           sum( likelihood(
#                             fun.likelihood = fun.likelihood,
#                             data.model = model.simulation.double$data.model %>% filter(stimulation %in% stimulation.list),
#                             data.exp.grouped = data.exp.grouped))
#                         })
# 
# 
# path.receptors <- paste(path.optimisation, "receptors", sep = "/")
# dir.create(path.receptors, recursive = TRUE, showWarnings = FALSE)
# 
# save_results(path.opt = path.receptors,
#              variables = variables,
#              variables.priming = variables.priming,
#              data.model.opt = model.simulation.double$data.model,
#              optimisation.opt = result.double,
#              par.opt = par.def,
#              res.list = list(),
#              data.exp.grouped = data.exp.grouped.all,
#              grid.ncol = 4,
#              grid.nrow = 2)
#### plots ####

optimisation.table <- data.table(
  id = character(),
  mvn.mean = numeric(),
  mvn = numeric(),
  lmvn = numeric(),
  mvn.sd_const = numeric(),
  lmvn.data = numeric()
)

data.model.list <- list()

path.receptors <- paste(path.optimisation, "receptors", sep = "/")
optimisation.table <- rbind(optimisation.table, read_optimisation(path = path.receptors, id = "receptors"))
data.model.list[["receptors"]] <- read.table(
  file = paste(path.receptors, "data_model.csv", sep = "/"),
  sep = ",",
  header = TRUE)


plot_results(data = data.exp.grouped,
             path.analysis = path.optimisation.data,
             path.output = path.optimisation.data,
             data.model.best.list = data.model.list,
             optimisation.best =  c("single"),
             plot.title = paste("comparison", id, sep ="-"), 
             filename.data_model = "data_model.csv",
             grid.ncol = 4)


