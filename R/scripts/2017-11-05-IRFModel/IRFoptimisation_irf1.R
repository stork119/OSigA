### ###
### IRFoptimisation_irf1
### ###

ranges.ps1 <- GetParametersRanges.ps1()
data.sample.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
                                     par = par.ps1, 
                                     stimulations = stimulations, 
                                     nsimulations = nsimulations)
k <- 4
ranges.irf1 <- GetParametersRanges.irf1()
par.irf1 <- ranges.irf1$par

