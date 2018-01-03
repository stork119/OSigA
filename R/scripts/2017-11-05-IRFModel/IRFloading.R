### ###
### Loading optimised parameters 
### ###


#### load pars ####

irfmodel.path.list$optimisation.id <- "2017-12-28-pSTAT"

irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")

irfmodel.results.ps1 <- readRDS(file = 
                                  paste(irfmodel.path.list$output.path, "IRFmodel-pstat.RDS", sep = "/"))
par.ps1 <- irfmodel.results.ps1$par

irfmodel.results.irf1 <- readRDS(file = 
                                   paste(irfmodel.path.list$output.path, "IRFmodel-irf.RDS", sep = "/"))
par.irf1 <- irfmodel.results.irf1$par
