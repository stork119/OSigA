### ###
### IRFoptimisation_irf1
### ###

params.ps1 <- ranges.ps1$factor
params.ps1[ranges.ps1$opt] <- 
  ranges.ps1$factor[ranges.ps1$opt]*ranges.ps1$base[ranges.ps1$opt]^par.ps1
data.model.ps1 <- model_fun_stm_params.ps1(
  stimulations = stimulations,
  params = params.ps1
)

par.irf1 <- ranges.irf1$par
stopfitness <- -10000
fun.optimisation = cma_es
maxit <- 1000
#### optimisation running ####
stimulations.pSTAT <- irfmodel.data.list$pSTAT %>% dplyr::distinct(stimulation) 
stimulations.irf <- irfmodel.data.list$irf %>% dplyr::distinct(stimulation)
stimulations  <- (stimulations.pSTAT %>% dplyr::inner_join(stimulations.irf))$stimulation

data.model.colnames <- c("irf")
 
data.raw.list <-
  list(irfmodel.data.list$irfsum %>% dplyr::filter(stimulation %in% stimulations))

optimisation.res <- do.call(
  fun.optimisation,
  list(par = par.irf1,
       fn = optimise.fun,
       control = list(maxit = maxit,
                      stopfitness = stopfitness,
                      diag.sigma  = FALSE,
                      #keep.best  = TRUE,
                      diag.eigen  = FALSE,
                      diag.pop    = FALSE,
                      diag.value  = FALSE),
       lower = ranges.irf1$min[ranges.irf1$opt],
       upper = ranges.irf1$max[ranges.irf1$opt],
       data.raw.list = data.raw.list,
       stimulations = stimulations,
       ranges.base = ranges.irf1$base,
       ranges.factor = ranges.irf1$factor,
       ranges.opt = ranges.irf1$opt,
       model_fun = model_fun_mean.irf1,
       data.model.colnames = data.model.colnames,
       data.model.ps1 = data.model.ps1)
)
par.irf1 <- optimisation.res$par 
#### ####
irfmodel.path.list$optimisation.id <- "2017-12-28-pSTAT"

params.ps1 <- ranges.ps1$factor
params.ps1[ranges.ps1$opt] <- 
  ranges.ps1$factor[ranges.ps1$opt]*ranges.ps1$base[ranges.ps1$opt]^par.ps1
data.model.ps1 <- model_fun_stm_params.ps1(
  stimulations = stimulations,
  params = params.ps1
)

par.irf1 <- optimisation.res.irf1$par
# par.irf1[1] <- -1
# par.irf1[4] <- 4
# par.irf1[5] <- 0.8

#### ####
params.irf1 <- ranges.irf1$factor
params.irf1[ranges.irf1$opt] <- 
  ranges.irf1$factor[ranges.irf1$opt]*ranges.irf1$base[ranges.irf1$opt]^par.irf1


data.model.irf1 <- model_fun_mean.irf1(params = params.irf1,
                                       data.model.ps1 = data.model.ps1)


data <- rbind(data.raw.sum %>% dplyr::select(-c(pstat, pstat.sd)),
              data.model.irf1 %>% dplyr::select(-c(pstat, pstat.sd, pstat.model, pstat.model.sd)))
data <- data[!is.na(data$irf), ]
data <- data.frame(data)
data$type <- factor(data$type)

g.list <- list()
g.list[["irf"]] <- ggplot(data = data,
                          mapping = aes(
                            x = stimulation,
                            y = irf,
                            ymin = irf - irf.sd,
                            ymax = irf + irf.sd,
                            group = type,
                            color = type))+
  geom_errorbar() +
  geom_point() +
  geom_line() +
  ggtitle("irf") +
  do.call(theme_jetka, args = plot.args)
g.list[["irf"]]

irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")
dir.create(irfmodel.path.list$output.path, 
           recursive = TRUE)

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(irfmodel.path.list$output.path, "IRFmodel-irf.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = g.list, ncol = 1, nrow = 1))))
saveRDS(
  file = paste(irfmodel.path.list$output.path, "IRFmodel-irf.RDS", sep = "/"),
  object = list(
    optimisation = optimisation.res,
    likelihood = optimisation.res$value,
    par = par.new, 
    ranges = ranges,
    #data.raw.list = data.raw.list,
    stimulations = stimulations,
    plots = g.list,
    stopfitness = stopfitness,
    fun.optimisation = fun.optimisation,
    maxit = maxit
  ))   


#### stochastic model ####
nsimulations <- 1000
ranges.ps1 <- GetParametersRanges.ps1()
data.sample.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
                                     par = par.ps1, 
                                     stimulations = stimulations, 
                                     nsimulations = nsimulations)

k <- 4

params <- ranges.irf1$factor
params[ranges.irf1$opt] <- ranges.irf1$factor[ranges.irf1$opt]*ranges.irf1$base[ranges.irf1$opt]^par.irf1
params.list <- GetParametersList.irf1(params = params)
ranges.irf1.stochastic <- GetParametersRanges.irf1(
  sd.max = 4, 
  hn.max = -4,
  theta.max = c(-4, -4, -4, -4),
  scale.max = -4,
  hn.factor = params.list$hn,
  theta.factor = params.list$theta,
  sd.factor = params.list$sd,
  scale.factor = params.list$scale
  )

maxit <- 1
stopfitness <- 0 

optimisation.res.irf1.stochastic <- do.call(
  fun.optimisation,
  list(par = ranges.irf1.stochastic$par,
       fn = optimise.fun.stochastic.irf1,
       control = list(maxit = maxit,
                      stopfitness = stopfitness,
                      diag.sigma  = FALSE,
                      #keep.best  = TRUE,
                      diag.eigen  = FALSE,
                      diag.pop    = FALSE,
                      diag.value  = FALSE),
       lower = ranges.irf1.stochastic$min[ranges.irf1.stochastic$opt],
       upper = ranges.irf1.stochastic$max[ranges.irf1.stochastic$opt],
       data.sample = data.sample.ps1,
       stimulations = stimulations,
       ranges = ranges.irf1.stochastic,
       irfmodel.data.list = irfmodel.data.list,
       k = k,
       nsimulations = nsimulations)
)
