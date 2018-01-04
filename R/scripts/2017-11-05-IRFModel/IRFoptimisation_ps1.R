### ###
### IRF model optimisation
### ###

source("R/scripts/2017-11-05-IRFModel/IRFcomputation.R")

#### optimisation initialisation ####

#par <- 0*ranges.opt 
#par <- par.new[ranges.opt]
#par[1:length(par.new)] <- par.new
#data.model.colnames <- c("pstat", "irf")
data.model.colnames <- c("pstat")

ranges.ps1 <- GetParametersRanges.ps1()
par.ps1 <- 0*ranges.ps1$opt # lub ustawic swoje par.ps1


#### optimisation running ####
stopfitness <- -1000
fun.optimisation = cma_es
maxit <- 1000

stimulations.pSTAT <- irfmodel.data.list$pSTAT %>% dplyr::distinct(stimulation) 
stimulations.irf <- irfmodel.data.list$irf %>% dplyr::distinct(stimulation)
stimulations  <- (stimulations.pSTAT %>% dplyr::inner_join(stimulations.irf))$stimulation

# data.raw.list <-
#   list(irfmodel.data.list$pSTATsum %>% dplyr::filter(stimulation %in% stimulations),
#        irfmodel.data.list$irfsum %>% dplyr::filter(stimulation %in% stimulations))

data.raw.list <-
    list(irfmodel.data.list$pSTATsum %>% dplyr::filter(stimulation %in% stimulations))

optimisation.res.ps1 <- do.call(
  fun.optimisation,
  list(par = par.ps1,
       fn = optimise.fun,
       control = list(maxit = maxit,
                      stopfitness = stopfitness,
                      diag.sigma  = TRUE,
                      keep.best  = TRUE,
                      diag.eigen  = TRUE,
                      diag.pop    = TRUE,
                      diag.value  = TRUE),
       lower = ranges.ps1$min[ranges.ps1$opt],
       upper = ranges.ps1$max[ranges.ps1$opt],
       data.raw.list = data.raw.list,
       stimulations = stimulations,
       ranges.base = ranges.ps1$base,
       ranges.factor = ranges.ps1$factor,
       ranges.opt = ranges.ps1$opt,
       model_fun = model_fun_stm_params.ps1,
       data.model.colnames = data.model.colnames)
)

#### plotting ####
#par.ps1 <- par.new
optimisation.res.ps1$value
par.ps1 <- optimisation.res.ps1$par

params <- ranges.ps1$factor
params[ranges.ps1$opt] <- ranges.ps1$factor[ranges.ps1$opt]*ranges.ps1$base[ranges.ps1$opt]^par.ps1

data.model <- model_fun_stm_params.ps1(
  stimulations = unique(data.raw.sum$stimulation),
  params = params)

data <- rbind(data.raw.sum %>% dplyr::select(-c(irf, irf.sd)),
              data.model %>% dplyr::select(-c(pstat.model, pstat.model.sd)))
data <- data[!is.na(data$pstat), ]
data <- data.frame(data)
data$type <- factor(data$type)

#data <- data[data$stimulation != 0.01, ]

g.list <- list()
g.list[["pstat"]] <- ggplot(data = data, 
                            mapping = aes(
                              x = stimulation,
                              y = pstat, 
                              ymin = pstat - pstat.sd,
                              ymax = pstat + pstat.sd,
                              group = type, 
                              color = type
                            ))+
  geom_errorbar() +
  geom_point() +
  geom_line() +
  ggtitle("pstat") + 
  do.call(theme_jetka, args = plot.args)
g.list[["pstat"]]


g.list[["pstat-factor"]] <- ggplot(data = data, 
                            mapping = aes(
                              x = factor(stimulation),
                              y = pstat, 
                              ymin = pstat - pstat.sd,
                              ymax = pstat + pstat.sd,
                              group = type, 
                              color = type
                            ))+
  geom_errorbar() +
  geom_point() +
  geom_line() +
  ggtitle("pstat") + 
  do.call(theme_jetka, args = plot.args)
g.list[["pstat-factor"]]

saveResults(
  model.type = "pstat",
  irfmodel.path.list = irfmodel.path.list,
  optimisation.res = optimisation.res.ps1,
  likelihood = optimisation.res.ps1$value, 
  par = optimisation.res.ps1$par,  
  ranges = ranges.ps1,
  stimulations = stimulations,
  stopfitness = stopfitness,
  fun.optimisation = fun.optimisation,
  maxit = maxit,
  g.list = g.list
)


#### ####
# optimise.fun(
#   par = par.new,
#   data.raw.list = data.raw.list,
#   stimulations = stimulations,
#   ranges.base = ranges.base,
#   ranges.factor = ranges.factor,
#   ranges.opt = ranges.opt
# )
# #### load RDS  ####
# irfmodel.path.list$optimisation.id <- "2017-12-13-summarise"
# irfmodel.path.list$output.path <-
#   paste(irfmodel.path.list$output.dir,
#         irfmodel.path.list$optimisation.id, sep = "/")
# results <- readRDS(
#   file = paste(irfmodel.path.list$output.path, "IRFmodel.RDS", sep = "/"))
# par <- results$par#[ranges.opt]
# 
# ## [1]  0.093594553 -0.001356297 -1.073560733 -0.103776065 -1.733805262 -0.523462295  0.640583619  3.960325670
# ##  [9] -3.225818492 -0.357628990
