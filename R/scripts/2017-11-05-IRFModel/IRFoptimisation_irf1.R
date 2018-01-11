### ###
### IRFoptimisation_irf1
### ###

source("R/scripts/2017-11-05-IRFModel/IRFcomputation.R")

#### ps1 simulations ####
params.ps1 <- ranges.ps1$factor
params.ps1[ranges.ps1$opt] <- 
  ranges.ps1$factor[ranges.ps1$opt]*ranges.ps1$base[ranges.ps1$opt]^par.ps1
data.model.ps1 <- model_fun_stm_params.ps1(
  stimulations = stimulations,
  params = params.ps1
)

#### optimisation running ####
stopfitness <- -10000
fun.optimisation = cma_es
maxit <- 10000
stimulations.pSTAT <- irfmodel.data.list$pSTAT %>% dplyr::distinct(stimulation) 
stimulations.irf <- irfmodel.data.list$irf %>% dplyr::distinct(stimulation)
stimulations  <- sort((stimulations.pSTAT %>% dplyr::inner_join(stimulations.irf))$stimulation)

data.model.colnames <- c("irf")
data.raw.list <-
  list(irfmodel.data.list$irfsum %>% dplyr::filter(stimulation %in% stimulations))

ranges.irf1 <- GetParametersRanges.irf1(scale.max = -4, sd.max = -4)
par.irf1 <- ranges.irf1$par #
optimisation.res.irf1 <- do.call(
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
par.irf1 <- optimisation.res.irf1$par 

#### simulations type I - improper ####

ranges.irf1 <- GetParametersRanges.irf1(scale.max = -4, sd.max = -4, ranges.max = 2, ranges.min = 2)


optimisation.res.irf1.simulations <- do.call(
  fun.optimisation,
  list(par = par.irf1.simulations,
       fn = optimise.fun,
       control = list(maxit = maxit,
                      stopfitness = stopfitness,
                      diag.sigma  = TRUE,
                      #keep.best  = TRUE,
                      diag.eigen  = TRUE,
                      diag.pop    = TRUE,
                      diag.value  = TRUE),
       lower = ranges.irf1$min[ranges.irf1$opt],
       upper = ranges.irf1$max[ranges.irf1$opt],
       data.raw.list = data.raw.list,
       stimulations = stimulations,
       ranges.base = ranges.irf1$base,
       ranges.factor = ranges.irf1$factor,
       ranges.opt = ranges.irf1$opt,
       model_fun = model_fun_simulations.irf1,
       data.model.colnames = data.model.colnames,
       data.model.ps1 = data.model.ps1,
       ranges = ranges.irf1,
       data.sample.ps1 = data.sample.sdnonconst.ps1)
)
par.irf1.simulations <- optimisation.res.irf1.simulations$par

# par.irf1[1] <- -1
# par.irf1[4] <- 4
# par.irf1[5] <- 0.8

#### ####

## ? o co mi chodzilo
#ranges.irf1 <- GetParametersRanges.irf1(scale.max = -4, sd.max = -4)
# params <- ranges.irf1$factor
# params[ranges.irf1$opt] <- ranges.irf1$factor[ranges.irf1$opt]*ranges.irf1$base[ranges.irf1$opt]^par.irf1
# params.list <- GetParametersList.irf1(params = params)
# params.list$sd <- par.irf.stochastic
# ranges.irf1 <- GetParametersRanges.irf1(
#   sd.max = 4, 
#   hn.max = -4,
#   theta.max = c(-4, -4, -4, -4),
#   scale.max = -4,
#   hn.factor = params.list$hn,
#   theta.factor = params.list$theta,
#   sd.factor = params.list$sd,
#   scale.factor = params.list$scale
# )

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

g.list[["irf-factor"]] <- ggplot(data = data,
                          mapping = aes(
                            x = factor(stimulation),
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
g.list[["irf-factor"]]

saveResults(
  model.type = "irf", #"pstat"
  irfmodel.path.list = irfmodel.path.list,
  optimisation.res = optimisation.res.irf1,
  likelihood = optimisation.res.irf1$value, 
  par = optimisation.res.irf1$par,  
  ranges = ranges.irf1,
  stimulations = stimulations,
  stopfitness = stopfitness,
  fun.optimisation = fun.optimisation,
  maxit = maxit,
  g.list = g.list
)

#### simulations type II - KL ####
## zkaladam ze nie ma szumu we procesie transkrypcji i translacji - czyli mamy do czynienia z ukladem deterministycznym
ranges.ps1 <- GetParametersRanges.ps1()
nsimulations <- 1000
data.sample.sdnonconst.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
                                                par = par.ps1, 
                                                stimulations = stimulations, 
                                                nsimulations = nsimulations,
                                                scaled = FALSE,
                                                sd.const = FALSE, 
                                                data.exp = irfmodel.data.list$pSTATsum)

data.sample.sdnonconst.ps1 <- data.sample.sdnonconst.ps1 %>%
  dplyr::mutate(response = response.pstat)


optimise.fun.simulations.irf1 <- 
  function(par = par.irf1,
           ranges = ranges.irf1,
           data.sample.ps1,
           nsimulations = 1000,
           no_cores = 6,
           irfmodel.data.list,
           k = 10,
           ...
           ){
    
    data.sample.irf1 <- 
  simulateMeanModel.irf(
    ranges = ranges,
    par = par,
    data.sample = data.sample.ps1,
    nsimulations = nsimulations,
    no_cores = no_cores) %>%
  dplyr::mutate(response.irf1 =  response.irf1.mean)

  stimulation.list <- (data.sample.irf1 %>% dplyr::distinct(stimulation))$stimulation
    
  registerDoParallel(no_cores)
  kl.list <- foreach(stm = stimulation.list) %dopar% {
    kl.div <- KL.divergence(
      X = (irfmodel.data.list$irf %>% 
             dplyr::filter(stimulation == stm))$logresponse,
      Y = (data.sample.irf1 %>% 
             dplyr::filter(stimulation == stm))$response.irf1,
      k = k)[k]
    return(kl.div)
  }
  stopImplicitCluster()
  
  kl <- do.call(sum, kl.list)
}



maxit <- 1000
nsimulations <- 1000
stopfitness <- 0 
no_cores <-8 
optimisation.res.irf1.simulations <- do.call(
  fun.optimisation,
  list(par = par,
       fn = optimise.fun.simulations.irf1,
       control = list(maxit = maxit,
                      stopfitness = stopfitness,
                      diag.sigma  = FALSE,
                      #keep.best  = TRUE,
                      diag.eigen  = TRUE,
                      diag.pop    = TRUE,
                      diag.value  = TRUE),
       lower = ranges.irf1$min[ranges.irf1$opt],
       upper = ranges.irf1$max[ranges.irf1$opt],
       data.sample = data.sample.sdnonconst.ps1,
       stimulations = stimulations,
       ranges = ranges.irf1,
       irfmodel.data.list = irfmodel.data.list,
       k = k,
       nsimulations = nsimulations,
       no_cores = no_cores)
)

optimisation.res.irf1.simulations$value
ranges.irf1
par <- optimisation.res.irf1.simulations$par
par.irf1.simulations <- par
data.sample.irf1 <- 
  simulateMeanModel.irf(
    ranges = ranges.irf1,
    par = par,
    data.sample = data.sample.sdnonconst.ps1,
    nsimulations = nsimulations,
    no_cores = 6) %>%
  dplyr::mutate(response.irf1 =  response.irf1.mean)

stimulation.list <- (data.sample.irf1 %>% dplyr::distinct(stimulation))$stimulation

data <-
  irfmodel.data.list$irf %>% 
  dplyr::select(stimulation, logresponse) %>%
  dplyr::mutate(type = "data") %>%
  rbind((data.sample.irf1 %>% 
           dplyr::mutate(logresponse = response.irf1.mean) %>%
           dplyr::select(stimulation,
                         logresponse) %>%
           dplyr::mutate(type = "model")))


g.plot.list <- foreach(stm = stimulation.list) %do% {
  ggplot(data %>% dplyr::filter(stimulation == stm),
         mapping = aes(x = logresponse, group = type, color = type))  + 
    geom_density() +
    do.call(theme_jetka, args = plot.args) +
    ggtitle(stm)
}
marrangeGrob(grobs = g.plot.list, nrow = 2, ncol = 3)


#### stochastic model ####
###
# 1. Zakaldamy ze istnieje szum w procesie pS1 -> IRF
# 2. Bedziemy sprawdzac, jak dobry fit dostaniemy, dy bedziemy losowac 
# 3. Dla optymalizacji dodajemy wagi - do implementacji !!!
# # losujemy probke, 
# # nastepnie dzilimy na q czesci 
# # z kadej grupy wyznaczamy element sredni
# # dodajemy wage do kazdego elementu zgodnie z licznoscia 
# Metoda jest empiryczna - majac rozklad normalny mozemy wyznaczyc jego kwantyle analitycznie !!!
# Ewentualnie jak mamy srednia i odchylenie standardowe to wyznaczmy k kwantyli i kazdy bedzie mial taka sama wage
# Kwantyle beda latwiejsze 
# 4. Optymalizujemy wszystkie parametry poza odchyleniem standarodwym ( bo jego optymalizacja nie ma sensu)


#### ps1 simulations ####
ranges.ps1 <- GetParametersRanges.ps1()
params.ps1 <- ranges.ps1$factor
params.ps1[ranges.ps1$opt] <- 
  ranges.ps1$factor[ranges.ps1$opt]*ranges.ps1$base[ranges.ps1$opt]^par.ps1
data.model.sdnonconst.ps1 <- model_fun_stm_params.ps1(
  stimulations = stimulations,
  params = params.ps1,
  sd.const = FALSE,
  data.exp = irfmodel.data.list$pSTATsum
)
data.sample.test.ps1 <- data.model.sdnonconst.ps1 %>%
  dplyr::mutate(response = exp(pstat.model)) %>%
  dplyr::mutate(response.pstat = exp(pstat.model)) %>%
  dplyr::select(stimulation, response, response.pstat)

q <- 50
i <- 1
data.model.ps1 <- data.model.sdnonconst.ps1
data.sample.quantiles.ps1.list <- foreach(i = 1:nrow(data.model.ps1)) %do% {
  df <- data.frame(logresponse = 
             qnorm(p = seq(from = 0, to = 1, length.out = q), 
                   mean = data.model.ps1[i,]$pstat.model, 
                   sd = data.model.ps1[i,]$pstat.model.sd),
           stimulation = data.model.ps1$stimulation[i]) %>% 
    dplyr::filter(!is.infinite(logresponse)) %>%
    dplyr::mutate(response = exp(logresponse),
                  response.pstat = exp(logresponse)) %>%
    dplyr::select(stimulation, response, response.pstat)
  return(df)
}
data.sample.quantiles.ps1 <- do.call(what = rbind, args = data.sample.quantiles.ps1.list)
data.sample.ps1 <- data.sample.quantiles.ps1

# nsimulations <- 1000
# ranges.ps1 <- GetParametersRanges.ps1()
# data.sample.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
#                                      par = par.ps1, 
#                                      stimulations = stimulations, 
#                                      nsimulations = nsimulations)

sd.list <- c(0.0001, 0.001, 0.01, 0.1, 1, 10)
fun.optimisation = "cma_es"
maxit <- 1000
nsimulations <- 1000
stopfitness <- 0 
no_cores <- 8 
k <- 10
no_cores <- 6
par.irf1.stochastic <- par.irf1
registerDoParallel(no_cores)
optimisation.res.irf1.stochastic.list <- foreach(sd.factor = sd.list) %dopar% {

  ranges.irf1.stochastic <- GetParametersRanges.irf1(scale.max = -4,
                                        sd.max = -4, 
                                        sd.factor = sd.factor,
                                        par = par.irf1.simulations,
                                        ranges.max = 2, ranges.min = 2)
  #par.irf1.stochastic <- ranges.irf1.stochastic$par # optimisation.res.irf1.stochastic$par
  #ranges.irf1.stochastic$par <- optimisation.res.irf1.stochastic$par

  optimisation.res.irf1.stochastic <- do.call(
    fun.optimisation,
    list(par = par.irf1,
         fn = optimise.fun.stochastic.irf1,
         control = list(maxit = maxit,
                        stopfitness = stopfitness,
                        diag.sigma  = TRUE,
                        #keep.best  = TRUE,
                        diag.eigen  = TRUE,
                        diag.pop    = TRUE,
                        diag.value  = TRUE),
         lower = ranges.irf1.stochastic$min[ranges.irf1.stochastic$opt],
         upper = ranges.irf1.stochastic$max[ranges.irf1.stochastic$opt],
         data.sample = data.sample.quantiles.ps1,
         stimulations = stimulations,
         ranges = ranges.irf1.stochastic,
         irfmodel.data.list = irfmodel.data.list,
         k = k,
         nsimulations = nsimulations,
         no_cores = no_cores)
  )
  
  par.irf1.stochastic <- optimisation.res.irf1.stochastic$par 

  saveResults(
    model.type = paste("irf-stochastic", 
                       gsub(pattern = "[.]", 
                            replacement = "_", 
                            x = as.numeric(sd.factor))),
    irfmodel.path.list = irfmodel.path.list,
    optimisation.res = optimisation.res.irf1.stochastic,
    likelihood = optimisation.res.irf1.stochastic$value, 
    par = optimisation.res.irf1.stochastic$par,  
    ranges = ranges.irf1.stochastic,
    stimulations = stimulations,
    stopfitness = stopfitness,
    fun.optimisation = fun.optimisation,
    maxit = maxit
  )
  return(optimisation.res.irf1.stochastic)
}
stopImplicitCluster()

par.irf1.simulations
par.irf1.stochastic <- optimisation.res.irf1.stochastic.list[[3]]$par
optimisation.res.irf1.stochastic.list[[3]]$value