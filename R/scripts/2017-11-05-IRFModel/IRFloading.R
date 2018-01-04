### ###
### Loading optimised parameters 
### ###

source("R/scripts/2017-11-05-IRFModel/IRFcomputation.R")

#### load pars ####

irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")

irfmodel.results.ps1 <- readRDS(file = 
                                  paste(irfmodel.path.list$output.path, "IRFmodel-pstat.RDS", sep = "/"))
par.ps1 <- irfmodel.results.ps1[[1]]$par 

irfmodel.results.irf1 <- readRDS(file = 
                                   paste(irfmodel.path.list$output.path, "IRFmodel-irf.RDS", sep = "/"))
par.irf1 <- irfmodel.results.irf1$par
