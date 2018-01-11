### ###
### plot simulations
### ###

#### simulations ####
stimulations <- sort(stimulations)
nsimulations <- 10000
g.list <- list()
g.list[["simulations"]] <- list()
g.list[["errobars"]] <- list()
g.list[["variances"]] <- list()
#### ps1 simulations sdconst ####
data.sample.scaled.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
                                     par = par.ps1, 
                                     stimulations = stimulations, 
                                     nsimulations = nsimulations,
                                     scaled = TRUE,
                                     sd.const = TRUE, 
                                     data.exp = irfmodel.data.list$pSTATsum)

data.sample.scaled.ps1 <- data.sample.scaled.ps1 %>%
  dplyr::mutate(response = response.pstat)

plots.args.list <- list(data.exp = irfmodel.data.list$pSTAT %>% 
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "data"),
                        data.sample = data.sample.scaled.ps1 %>% 
                          dplyr::mutate(response = response.pstat)  %>%
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "model"),
                        stimulations = stimulations,
                        plot.title = "ps1",
                        plot.grob.title = "ps1 simulation",
                        plot.grob.nrow = 2,
                        plot.grob.ncol = 3)
g.list[["simulations"]][["ps1"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["ps1"]] <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["ps1"]] <-  do.call(plotLogSD, plots.args.list)


#### ps1 simulations sdconst ####
data.sample.sdnonconst.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
                                            par = par.ps1, 
                                            stimulations = stimulations, 
                                            nsimulations = nsimulations,
                                            scaled = TRUE,
                                            sd.const = FALSE, 
                                            data.exp = irfmodel.data.list$pSTATsum)

data.sample.sdnonconst.ps1 <- data.sample.sdnonconst.ps1 %>%
  dplyr::mutate(response = response.pstat)

plots.args.list <- list(data.exp = irfmodel.data.list$pSTAT %>% 
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "data"),
                        data.sample = data.sample.sdnonconst.ps1 %>% 
                          dplyr::mutate(response = response.pstat)  %>%
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "model"),
                        stimulations = stimulations,
                        plot.title = "ps1 -sdnonconst",
                        plot.grob.title = "ps1 simulation",
                        plot.grob.nrow = 2,
                        plot.grob.ncol = 3)
g.list[["simulations"]][["ps1-sdnonconst"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["ps1-sdnonconst"]] <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["ps1-sdnonconst"]] <-  do.call(plotLogSD, plots.args.list)


#### irf1 model ####
data.sample.irf1 <- 
  simulateMeanModel.irf(
    ranges = ranges.irf1,
    par = par.irf1,
    data.sample = data.sample.scaled.ps1,
    nsimulations = nsimulations,
    no_cores = 6)

plots.args.list <- list(data.exp = irfmodel.data.list$irf %>% 
                              dplyr::select(stimulation, response) %>%
                              dplyr::mutate(type = "data"),
                            data.sample = data.sample.irf1 %>% 
                              dplyr::mutate(response = exp(response.irf1.mean)) %>%
                              dplyr::select(stimulation, response) %>%
                              dplyr::mutate(type = "model"),
                            stimulations = stimulations,
                            plot.title = "irf1 simulations",
                            plot.grob.title = "irf1 simulation",
                            plot.grob.nrow = 2,
                            plot.grob.ncol = 3)

g.list[["simulations"]][["irf1"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["irf1"]] <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["irf1"]] <-  do.call(plotLogSD, plots.args.list)

#### irf sd non const ####
data.sample.sdnonconst.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
                                            par = par.ps1, 
                                            stimulations = stimulations, 
                                            nsimulations = nsimulations,
                                            scaled = FALSE,
                                            sd.const = FALSE, 
                                            data.exp = irfmodel.data.list$pSTATsum)

data.sample.sdnonconst.ps1 <- data.sample.sdnonconst.ps1 %>%
  dplyr::mutate(response = response.pstat)

data.sample.sdnonconst.irf1 <- 
  simulateMeanModel.irf(
    ranges = ranges.irf1,
    par = par.irf1,
    data.sample = data.sample.sdnonconst.ps1,
    nsimulations = nsimulations,
    no_cores = 6)

plots.args.list <- list(data.exp = irfmodel.data.list$irf %>% 
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "data"),
                        data.sample = data.sample.sdnonconst.irf1 %>% 
                          dplyr::mutate(response = exp(response.irf1.mean)) %>%
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "model"),
                        stimulations = stimulations,
                        plot.title = "irf1 simulations -sdnonconst",
                        plot.grob.title = "irf1 simulations",
                        plot.grob.nrow = 2,
                        plot.grob.ncol = 3)

g.list[["simulations"]][["irf1-sdnoncaled"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["irf1-sdnoncaled"]] <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["irf1-sdnoncaled"]] <-  do.call(plotLogSD, plots.args.list)

#### irf1 kl model ####
data.sample.irf1 <- 
  simulateMeanModel.irf(
    ranges = ranges.irf1,
    par = par.irf1.simulations,
    data.sample = data.sample.scaled.ps1,
    nsimulations = nsimulations,
    no_cores = 6)

plots.args.list <- list(data.exp = irfmodel.data.list$irf %>% 
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "data"),
                        data.sample = data.sample.irf1 %>% 
                          dplyr::mutate(response = exp(response.irf1.mean)) %>%
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "model"),
                        stimulations = stimulations,
                        plot.title = "KL irf1 simulation",
                        plot.grob.title = "KL irf1 simulation",
                        plot.grob.nrow = 2,
                        plot.grob.ncol = 3)

g.list[["simulations"]][["irf1-kl"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["irf1-kl"]] <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["irf1-kl"]] <-  do.call(plotLogSD, plots.args.list)

#### irf kl sd non const ####
data.sample.sdnonconst.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
                                                par = par.ps1, 
                                                stimulations = stimulations, 
                                                nsimulations = nsimulations,
                                                scaled = FALSE,
                                                sd.const = FALSE, 
                                                data.exp = irfmodel.data.list$pSTATsum)

data.sample.sdnonconst.ps1 <- data.sample.sdnonconst.ps1 %>%
  dplyr::mutate(response = response.pstat)

data.sample.sdnonconst.irf1 <- 
  simulateMeanModel.irf(
    ranges = ranges.irf1,
    par = par.irf1.simulations,
    data.sample = data.sample.sdnonconst.ps1,
    nsimulations = nsimulations,
    no_cores = 6)

plots.args.list <- list(data.exp = irfmodel.data.list$irf %>% 
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "data"),
                        data.sample = data.sample.sdnonconst.irf1 %>% 
                          dplyr::mutate(response = exp(response.irf1.mean)) %>%
                          dplyr::select(stimulation, response) %>%
                          dplyr::mutate(type = "model"),
                        stimulations = stimulations,
                        plot.title = "KL irf1 simulations -sdnonconst",
                        plot.grob.title = "KL irf1 simulations -sdnonconst",
                        plot.grob.nrow = 2,
                        plot.grob.ncol = 3)

g.list[["simulations"]][["irf1-kl-sdnoncaled"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["irf1-kl-sdnoncaled"]] <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["irf1-kl-sdnoncaled"]] <-  do.call(plotLogSD, plots.args.list)



#### irf1 model stochastic ####
i <- 4
par.irf1.stochastic <- optimisation.res.irf1.stochastic.list[[i]]$par
ranges.irf1.stochastic <- GetParametersRanges.irf1(scale.max = -4,
                                                   sd.max = -4, 
                                                   sd.factor = sd.list[[i]],
                                                   par = par.irf1.simulations,
                                                   ranges.max = 2, ranges.min = 2)
data.sample.irf1.stochastic <- 
  simulateModel.irf(
    ranges = ranges.irf1.stochastic,
    par = par.irf1.stochastic,
    data.sample = data.sample.ps1,
    nsimulations = nsimulations,
    no_cores = 16) %>% 
  dplyr::mutate(response.irf = exp(response.irf))

#KL.divergence((data.exp %>% dplyr::filter(stimulation == 1))$response, (data.sample %>% dplyr::filter(stimulation == 1))$response, k = 10)

plots.args.list <- list(data.exp = irfmodel.data.list$irf %>% 
                               dplyr::select(stimulation, response) %>%
                               dplyr::mutate(type = "data"),
                             data.sample = data.sample.irf1.stochastic %>% 
                               dplyr::mutate(response =  response.irf) %>%
                               dplyr::select(stimulation, response) %>%
                               dplyr::mutate(type = "model"),
                             stimulations = stimulations,
                             plot.title = "IRF1 STOCHASTIC",
                             plot.grob.title = "irf1 simulation",
                             plot.grob.nrow = 2,
                             plot.grob.ncol = 3)

g.list[["simulations"]][["irf1.stochastic"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["irf1.stochastic"]]   <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["irf1.stochastic"]] <-  do.call(plotLogSD, plots.args.list)
#### irf1 compare ####
plots.args.list <- list(data.exp = irfmodel.data.list$irf %>% 
                                          dplyr::select(stimulation, response) %>%
                                          dplyr::mutate(type = "data"),
                                        data.sample = data.sample.irf1.stochastic %>% 
                                          dplyr::mutate(response =  response.irf) %>%
                                          dplyr::select(stimulation, response) %>%
                                          dplyr::mutate(type = "model") %>%
                                          rbind((data.sample.irf1 %>% 
                                                  dplyr::mutate(response = exp(response.irf1.mean)) %>%
                                                  dplyr::select(stimulation, response) %>%
                                                  dplyr::mutate(type = "model-mean"))),
                                        stimulations = stimulations,
                                        plot.title = "IRF1 COMPARE",
                                        plot.grob.title = "irf1 compare simulation",
                                        plot.grob.nrow = 2,
                                        plot.grob.ncol = 3)

g.list[["simulations"]][["irf1.compare"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["irf1.compare"]]   <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["irf1.compare"]] <-  do.call(plotLogSD, plots.args.list)

#### ####
model.type <- "simulations"
irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")
dir.create(irfmodel.path.list$output.path, 
           recursive = TRUE)
do.call(what = pdf,
        args = append(
          plot.args.ggsave,
          list(file = paste(irfmodel.path.list$output.path, 
                            paste("IRFmodel-", 
                                  model.type,
                                  ".pdf", 
                                  sep = ""),
                            sep = "/"))))

l <- sapply(g.list[["simulations"]],
            function(x){
              print(x)
              return()
              })
dev.off()



do.call(what = pdf,
        args = append(
          plot.args.ggsave,
          list(file = paste(irfmodel.path.list$output.path, 
                            paste("IRFmodel-", 
                                  model.type,
                                  "-errorbars",
                                  ".pdf", 
                                  sep = ""),
                            sep = "/"))))

l <- sapply(g.list$errobars,
            function(x){
              print(x)
              return()
            })
dev.off()




do.call(what = pdf,
        args = append(
          plot.args.ggsave,
          list(file = paste(irfmodel.path.list$output.path, 
                            paste("IRFmodel-", 
                                  model.type,
                                  "-variances",
                                  ".pdf", 
                                  sep = ""),
                            sep = "/"))))

l <- sapply(g.list$variances,
            function(x){
              print(x)
              return()
            })
dev.off()
