
#### pre processing ####
tmesh.exp <- unique(data.exp$time)
data.exp$logintensity <- log(data.exp$intensity)

data.exp.unique <- distinct(data.exp, priming, stimulation, time)
data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time)

path.parameters <- "resources/input/"
par.def <- scan(file = paste(path.parameters, "par.txt", sep = ""))
varscale <- 0.15
variables <- rep(0.0, times = 629)
variables.priming <- rep(0.0, times = 629)
variables[1:17] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
variables.priming[1:17] <- scan(file = paste(path.parameters, "var-priming.txt", sep = ""))
variables[44:52] <- scan(file = paste(path.parameters, "var-extrinsic.txt", sep = ""))
variables.priming[44:52] <- variables[44:52]
variables[27:43] <- varscale*(variables[1:17]^2)
variables.priming[27:43] <- varscale*(variables.priming[1:17]^2)

tmesh <- seq(from = 0, to = 100, by = 5)
tmesh.list <- which(tmesh %in% tmesh.exp)

stimulation.list <- unique(data.exp$stimulation)[-1]
background <- mean(data.exp[data.exp$time == 5,]$intensity)
background.var <- var(data.exp[data.exp$time == 5,]$intensity)


# data.exp.grouped.expand <-  rbind(rbind(data.exp.grouped, data.exp.grouped %>% filter(time == 5)), data.exp.grouped %>% filter(time == 5))
# data.exp.grouped.expand %>% ungroup() %>% group_by(time, priming) %>% summarise(count = n())

#### display model ####
res.def <- run_model(parameters = par.def,
                 variables = variables,
                 variables.priming = variables.priming,
                 tmesh = tmesh,
                 tmesh.list = tmesh.list,
                 stimulation.list = stimulation.list,
                 background = background)

likelihood(data.model = res.def$data.model,
           data.exp.grouped = data.exp.grouped,
           fun.likelihood = fun.likelihood.mvn)

res <- optimisation(par = par.def,
             variables = variables,
             variables.priming = variables.priming,
             tmesh = tmesh,
             tmesh.list = tmesh.list,
             stimulation.list = stimulation.list,
             background = background,
             data.exp.grouped = data.exp.grouped.expand,
             return.model = TRUE,
             fun.likelihood = fun.likelihood.mvn)


#### Gen SA ####

par <- par.def
par.lower <- par*10^(-2)
par.upper <- par*10^(2)

genSA.res <- GenSA(par = par.def,
                   fn = optimisation,
                   lower = par.lower,
                   upper = par.upper,
                   control = list(maxit = 1000, max.call = 1000, verbose = TRUE), 
                   fun_run_model = run_model,
                   variables = variables,
                   variables.priming = variables.priming,
                   tmesh = tmesh,
                   tmesh.list = tmesh.list,
                   stimulation.list = stimulation.list,
                   background = background,
                   data.exp.grouped = data.exp.grouped,
                   fun.likelihood = fun.likelihood.mvn)  

par.gensa <- scan(file = "resources/output/gensa/2017-01-26/par.txt")

res.opt <- optimisation(par = par.gensa,
                    variables = variables,
                    variables.priming  = variables.priming,
                    fun_run_model = run_model,
                    tmesh = tmesh,
                    tmesh.list = tmesh.list,
                    stimulation.list = stimulation.list,
                    background = background,
                    data.exp.grouped = data.exp.grouped,
                    return.model = TRUE,
                    fun.likelihood = fun.likelihood.mvn)

path.gensa <- paste("resources/output/extrinsic/gensa")
dir.create(path.gensa)
g <- compare_models(
  data = data.exp.grouped,
  data.model.list = list(def = res.def$data.model, opt = res.opt$data.model),
  plot.title = "Model compares",
  filename = paste(path.gensa, "2017-01-26-sd", sep ="/"),
  plot.save = TRUE,
  plot.return = TRUE
)

gensa.res <- save_results(path.opt = paste(path.gensa , "2017-01-26-2", sep ="/"),
                         par.opt = genSA.res$par,
                         res.list  = list(def = res.def$data.model),
                                          data.exp.grouped = data.exp.grouped,
                                          fun_run_model = run_model,
                                          fun.likelihood = fun.likelihood.mvn.mean,
                                          variables = variables,
                                          tmesh = tmesh,
                                          tmesh.list = tmesh.list,
                                          stimulation.list = stimulation.list,
                                          background = background
)

#### ####
par <- scan('')
res.par <- optimisation(fun_run_model = run_model_mean,
                          par = par,
                    variables = variables,
                    tmesh = tmesh,
                    tmesh.list = tmesh.list,
                    stimulation.list = stimulation.list,
                    background = background,
                    data.exp.grouped = data.exp.grouped,
                    return.model = TRUE,
                    fun.likelihood = fun.likelihood.mvn)

g <- compare_models(
  data = data.exp.grouped.expand,
  data.model.list = list(def = res$data.model, opt = res.par$data.model),
  plot.title = "Model compares",
  filename = "resources/output/par-def",
  plot.save = TRUE,
  plot.return = TRUE
)