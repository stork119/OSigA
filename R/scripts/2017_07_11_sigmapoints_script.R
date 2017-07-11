### ###
### 2017-07-11-sigmapoints-script
### ###

setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")
source("R/sum_data_initialise.R")
source("R/model/ut/model_ut.R")

gplot.list <- list()
optimisation.conditions.toload <-
  LoadOptimisationConditions(path.optimisation = path.list$optimisation,
                             path.optimisation.data = path.list$optimisation.data,
                             maxit.tmp = Inf)
rm(list = labels(optimisation.conditions.toload))
attach(optimisation.conditions.toload)


sigmapoints <- LoadSigmapointsConditions()


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
                            fun_modify_input = PrepareModelArguments.ut,
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