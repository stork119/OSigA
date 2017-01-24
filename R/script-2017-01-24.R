
#### ####
tmesh.exp <- unique(data.exp$time)
data.exp$logintensity <- log(data.exp$intensity)

data.exp.unique <- distinct(data.exp, priming, stimulation, time)
data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time)

varscale <- 0.15
variables <- rep(0.0, times = 629)
variables[1:34] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
variables[35:68] <- varscale*(variables[1:34]^2)

tmesh <- seq(from = 0, to = 100, by = 5)
tmesh.list <- which(tmesh %in% tmesh.exp)

stimulation.list <- unique(data.exp$stimulation)[-1]
background <- mean(data.exp[data.exp$time == 5,]$intensity)
background.var <- var(data.exp[data.exp$time == 5,]$intensity)

par.def <- scan(file = paste(path.parameters, "par.txt", sep = ""))
#### display model ####
res <- run_model(parameters = par.def,
                 variables = variables,
                 tmesh = tmesh,
                 tmesh.list = tmesh.list,
                 stimulation.list = stimulation.list,
                 background = background)

likelihood(data.model = res$data.model,
           data.exp.grouped = data.exp.grouped)

optimisation(par = par,
             variables = variables,
             tmesh = tmesh,
             tmesh.list = tmesh.list,
             stimulation.list = stimulation.list,
             background = background,
             data.exp.grouped = data.exp.grouped)

if(!res$error){
  ggplot(res$data.model, 
         mapping = aes(x = time,
                       y = res$data.model$m.norm,
                       color = factor(priming),
                       group = interaction(priming, stimulation))) +
    geom_line()
}


#### Gen SA ####

par <- par.def
par.lower <- par*10^(-1)
par.upper <- par*10^(1)

genSA.res <- GenSA(par = par.def,
                   fn = optimisation,
                   lower = par.lower,
                   upper = par.upper,
                   control = list(maxit = 10, max.call = 10, verbose = TRUE), 
                   variables = variables,
                   tmesh = tmesh,
                   tmesh.list = tmesh.list,
                   stimulation.list = stimulation.list,
                   background = background,
                   data.exp.grouped = data.exp.grouped)  