### 
### double receptors
### 
tmesh.exp <- unique(data.exp$time)
data.exp$logintensity <- log(data.exp$intensity)

data.exp.unique <- distinct(data.exp, priming, stimulation, time)
data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time)

path.parameters <- "resources/input/"
par.def <- scan(file = paste(path.parameters, "par.txt", sep = ""))
varscale <- 0.15
variables.double <- rep(0.0, times = 629)
variables.double[1:34] <- scan(file = paste(path.parameters, "var-receptors.txt", sep = ""))
variables.double[35:68] <- varscale*(variables.double[1:34]^2)

tmesh <- seq(from = 0, to = 100, by = 5)
tmesh.list <- which(tmesh %in% tmesh.exp)

stimulation.list <- unique(data.exp$stimulation)[-1]
background <- mean(data.exp[data.exp$time == 5,]$intensity)
background.var <- var(data.exp[data.exp$time == 5,]$intensity)

#### display model ####
res.double <- optimisation(fun_run_model = run_model,
                        par = par.def,
                        variables = variables.double,
                        tmesh = tmesh,
                        tmesh.list = tmesh.list,
                        stimulation.list = stimulation.list,
                        background = background,
                        data.exp.grouped = data.exp.grouped,
                        return.model = TRUE,
                        fun.likelihood = fun.likelihood.mvn.mean)

path.double <- paste(path.output, "double", sep = "")
save_results(path.opt = paste(path.double , "2017-01-26", sep ="/"),
             par.opt = par.def,
             res.list  = list(def = res.def$data.model, gensa = gensa.res$data.model),
             data.exp.grouped = data.exp.grouped,
             fun_run_model = run_model,
             fun.likelihood = fun.likelihood.mvn.mean,
             variables = variables.double,
             tmesh = tmesh,
             tmesh.list = tmesh.list,
             stimulation.list = stimulation.list,
             background = background
)

g <- compare_models(
  data = data.exp.grouped,
  data.model.list = list(def = res.double$data.model, opt = res.par$data.model),
  plot.title = "Model compares",
  filename = "resources/output/double-opt",
  plot.save = TRUE,
  plot.return = TRUE
)