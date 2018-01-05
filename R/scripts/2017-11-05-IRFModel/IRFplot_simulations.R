### ###
### plot simulations
### ###

#### plotSimulationsFun ####
plotSimulationsFun <- function(data.exp,
                               data.sample,
                               stimulations,
                               plot.title = "",
                               plot.grob.title = "",
                               plot.grob.nrow = 1,
                               plot.grob.ncol = 1,
                               ...){
  
  g.list <- foreach( stm =  stimulations ) %do% {
    
    ggplot(data = rbind(data.exp,
                        data.sample) %>%
             dplyr::filter_(paste("stimulation ==", stm)), 
           mapping = aes_string(x = "log(response)", 
                                group = "type",
                                fill = "type",
                                color = "type")) +
      do.call(theme_jetka, args = plot.args) +
      geom_density(alpha = 0.5) +
      ggtitle(paste(plot.title, stm))
  }
  
  g.plot <- marrangeGrob(grobs = g.list,
                         nrow = plot.grob.nrow, 
                         ncol = plot.grob.ncol,
                         top  = plot.grob.title)
  return(g.plot)
}

#### plot errorbars ####
plotSimulationsErrorbarsFun <- function(data.exp,
                               data.sample,
                               plot.title = "",
                               ...){
  
  
   g.plot <- ggplot(data = 
                      rbind(data.exp, data.sample) %>% 
                      dplyr::group_by_("type", "stimulation") %>%
                      dplyr::summarise_(logresponse = "mean(log(response))",
                                        logresponse.sd = "sd(log(response))"), 
           mapping = aes_string(x = "factor(stimulation)", 
                                y = "logresponse", 
                                ymin = "logresponse - logresponse.sd",
                                ymax = "logresponse + logresponse.sd",
                                group = "interaction(stimulation, type)",
                                color = "type")) +
      do.call(theme_jetka, args = plot.args) +
      geom_errorbar() +
      geom_point() +
      ggtitle(paste(plot.title))
  return(g.plot)
}


stimulations <- sort(stimulations)
#### plotLogSD ####
plotLogSD <- function(
  data.exp,
  data.sample,
  plot.title = "",
  ...
){
  data <- 
    rbind(data.exp, data.sample) %>% 
    dplyr::group_by_("type", "stimulation") %>%
    dplyr::summarise_(logresponse = "mean(log(response))",
                      logresponse.sd = "sd(log(response))")
  g.plot <- ggplot(data    = data, 
         mapping = aes(x = factor(stimulation),
                       y = logresponse.sd, 
                       group = factor(type),
                       fill = factor(type))) +
    geom_bar(stat = "identity", position = "dodge") +
    do.call(theme_jetka, args = plot.args) +
    ggtitle(paste(plot.title)) 
  return(g.plot)
}

#### simulations ####
stimulations <- sort(stimulations)
nsimulations <- 10000
g.list <- list()
g.list[["simulations"]] <- list()
g.list[["errobars"]] <- list()
g.list[["variances"]] <- list()
#### ps1 simulations ####
data.sample.scaled.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
                                     par = par.ps1, 
                                     stimulations = stimulations, 
                                     nsimulations = nsimulations,
                                     scaled = TRUE,
                                     sd.const = FALSE, 
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
                        plot.title = "",
                        plot.grob.title = "ps1 simulation",
                        plot.grob.nrow = 2,
                        plot.grob.ncol = 3)
g.list[["simulations"]][["ps1"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["ps1"]] <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["ps1"]] <-  do.call(plotLogSD, plots.args.list)


#### irf1 model ####
data.sample.irf1 <- 
  simulateMeanModel.irf(
    ranges = ranges.irf1.stochastic,
    par = 0,
    data.sample = data.sample.ps1,
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
                            plot.title = "",
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
                        plot.title = "",
                        plot.grob.title = "irf1 simulation",
                        plot.grob.nrow = 2,
                        plot.grob.ncol = 3)

g.list[["simulations"]][["irf1-sdnoncaled"]] <- do.call(plotSimulationsFun, plots.args.list)
g.list[["errobars"]][["irf1-sdnoncaled"]] <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
g.list[["variances"]][["irf1-sdnoncaled"]] <-  do.call(plotLogSD, plots.args.list)


#### irf1 model ####
data.sample.irf1.stochastic <- 
  simulateModel.irf(
    ranges = ranges.irf1.stochastic,
    par = par.irf1.stochastic,
    data.sample = data.sample.ps1,
    nsimulations = nsimulations,
    no_cores = 16)

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
