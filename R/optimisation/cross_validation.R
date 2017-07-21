### ###
### cross validation ###
### ###


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


sigmapoints <- LoadSigmapointsConditions(path.optimisation = path.list$optimisation)
data.model <- read.table(paste("resources/output//optimisation/2017-07-19-ut-model_p1///data/", results.id, "data_model.csv", sep = "/"), header = TRUE, sep = ",")

data_sample_size <- 250
cv_times <- 100

no_cores <- 2
registerDoParallel(no_cores)
likelihood.list <- foreach(i = 1:cv_times) %dopar% {

  data.cv.list <- list()
  data.cv.list$grouped <- get_equal_data(data = data.list$data.exp, 
                                         sample_size = data_sample_size)
  
  data.cv.list$summarise  <- 
    data.cv.list$grouped %>% 
    dplyr::group_by(priming,
                    stimulation,
                    time) %>%
    dplyr::summarise(m.norm = mean(intensity),
              sd.norm   = var(intensity)) %>%
    dplyr::mutate(mean.lmvn = lmvn.mean(m = m.norm, sd = sd.norm),
                  sd.lmvn = lmvn.sd(m = m.norm, sd = sd.norm))
  
  sum(likelihood(data.model = data.model,# %>% dplyr::filter(time %in% tmesh[tmesh.list]),
               data.exp.grouped = data.cv.list$grouped,
               data.exp.summarise =   data.cv.list$summarise,
               fun.likelihood = fun.likelihood))

  
}
stopImplicitCluster()

likelihood.list <- unlist(likelihood.list)