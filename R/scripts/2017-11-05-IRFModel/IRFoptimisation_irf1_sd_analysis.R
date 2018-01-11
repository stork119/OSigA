### ### ### ### ###
### sd analysis ###
### ### ### ### ###

#### paths ####
irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")
dir.create(irfmodel.path.list$output.path)
#### prepare ps1 ####
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
#### optimisation parameters ####
# nsimulations <- 1000
# ranges.ps1 <- GetParametersRanges.ps1()
# data.sample.ps1 <- simulateModel.ps1(ranges = ranges.ps1,
#                                      par = par.ps1, 
#                                      stimulations = stimulations, 
#                                      nsimulations = nsimulations)
sd.mean <- mean(data.raw.sum$irf.sd[!is.na(data.raw.sum$irf.sd)])
sd.list <- sd.mean*(2^seq(from = -8, to = 8, by = 0.5))

fun.optimisation = "cma_es"
maxit <- 1000
nsimulations <- 1000
stopfitness <- 0 
no_cores <- 8 
k <- 10
no_cores <- 10
par.irf1.stochastic <- par.irf1

#### optimisation arguments ####
ranges.max <- 8
ranges.irf1.stochastic <- GetParametersRanges.irf1(scale.max = -4,
                                                   sd.max = -4, 
                                                   ranges.max = 2, 
                                                   ranges.min = 2)

lhs.res <- randomLHS(10, length(ranges.irf1.stochastic$par))
lhs.res <- ranges.max*(2*rbind(matrix(0.5 + ranges.irf1.stochastic$par, nrow = 1), lhs.res) - 1)

write.table(x = lhs.res,
            file = paste(irfmodel.path.list$output.path, "parameters_list.csv", sep = "/"),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

write.table(x = sd.list,
            file = paste(irfmodel.path.list$output.path, "sd_list.csv", sep = "/"),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

df.exp.ids <- expand.grid(sd.id  = 1:length(sd.list), 
                          par.id = 1:nrow(lhs.res) )

write.table(x = df.exp.ids,
            file = paste(irfmodel.path.list$output.path, "df_ids.csv", sep = "/"),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

#### computations  ####
registerDoParallel(no_cores)
par.irf1.stochastic.list <- foreach(id = 1:nrow(df.exp.ids)) %dopar% {
  
  output.path <- paste(irfmodel.path.list$output.path, 
                       id, 
                       sep = "/")
  dir.create(output.path, recursive = TRUE)
  
  par.lhs <- lhs.res[df.exp.ids[id,]$par.id,]
  sd.factor <- sd.list[df.exp.ids[id,]$sd.id]
  
  ranges.irf1.stochastic <- GetParametersRanges.irf1(scale.max = -1,
                                                     sd.max = -1, 
                                                     sd.factor = sd.factor,
                                                     par = par.lhs,
                                                     ranges.max = ranges.max, 
                                                     ranges.min = ranges.max)
  par.irf1.stochastic <- ranges.irf1.stochastic$par # optimisation.res.irf1.stochastic$par
  #ranges.irf1.stochastic$par <- optimisation.res.irf1.stochastic$par
  
  optimisation.res.irf1.stochastic <- do.call(
    fun.optimisation,
    list(par = par.irf1.stochastic,
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
         data.sample = data.sample.quantiles.ps1,
         stimulations = stimulations,
         ranges = ranges.irf1.stochastic,
         irfmodel.data.list = irfmodel.data.list,
         k = k,
         nsimulations = nsimulations,
         no_cores = no_cores)
  )
  
  par.irf1.stochastic <- optimisation.res.irf1.stochastic$par 

  data.sample.irf1.stochastic <- 
    simulateModel.irf(
      ranges = ranges.irf1.stochastic,
      par = par.irf1.stochastic,
      data.sample = data.sample.sdnonconst.ps1,
      nsimulations = nsimulations,
      no_cores = 1) %>% 
    dplyr::mutate(response.irf = exp(response.irf))
  
  plots.args.list <- list(data.exp = irfmodel.data.list$irf %>% 
                            dplyr::select(stimulation, response) %>%
                            dplyr::mutate(type = "data"),
                          data.sample = data.sample.irf1.stochastic %>% 
                            dplyr::mutate(response =  response.irf) %>%
                            dplyr::select(stimulation, response) %>%
                            dplyr::mutate(type = "model"),
                          stimulations = stimulations,
                          plot.title = paste("IRF1 stochastic", id, sd.factor),
                          plot.grob.title = paste("IRF1 stochastic", id, sd.factor),
                          plot.grob.nrow = 2,
                          plot.grob.ncol = 3)
  g.list <- list()
  g.list[["simulations"]] <- do.call(plotSimulationsFun, plots.args.list)
  g.list[["errobars"]]   <-  do.call(plotSimulationsErrorbarsFun, plots.args.list)
  g.list[["variances"]] <-  do.call(plotLogSD, plots.args.list)
  
  model.type <- paste("irf-stochastic", id, sep = "_")
  
  saveResults(
    output.path = output.path,
    model.type =  model.type,
    id = id,
    irfmodel.path.list = irfmodel.path.list,
    optimisation.res = optimisation.res.irf1.stochastic,
    likelihood = optimisation.res.irf1.stochastic$value, 
    par.init = ranges.irf1.stochastic$par,  
    par = optimisation.res.irf1.stochastic$par, 
    par.lhs = par.lhs,
    sd.factor = sd.factor,
    ranges = ranges.irf1.stochastic,
    stimulations = stimulations,
    stopfitness = stopfitness,
    fun.optimisation = fun.optimisation,
    maxit = maxit
  )
  
  do.call(what = pdf,
          args = append(
            plot.args.ggsave,
            list(file = paste(output.path, 
                              paste("IRFmodel-", 
                                    model.type,
                                    ".pdf", 
                                    sep = ""),
                              sep = "/"))))
  
  l <- sapply(g.list,
              function(x){
                print(x)
                return()
              })
  dev.off()
  
  return(par.irf1.stochastic)
}
stopImplicitCluster()
