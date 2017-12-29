### ###
### IRF model optimisation
### ###


#### optimisation initialisation ####

#par <- 0*ranges.opt 
#par <- par.new[ranges.opt]
#par[1:length(par.new)] <- par.new
#data.model.colnames <- c("pstat", "irf")
data.model.colnames <- c("pstat")

ranges <- GetParametersRanges.ps1()
par <- 0*ranges$opt

stopfitness <- -1000
fun.optimisation = cma_es
maxit <- 1000
#### optimisation running ####
stimulations.pSTAT <- irfmodel.data.list$pSTAT %>% dplyr::distinct(stimulation) 
stimulations.irf <- irfmodel.data.list$irf %>% dplyr::distinct(stimulation)
stimulations  <- (stimulations.pSTAT %>% dplyr::inner_join(stimulations.irf))$stimulation

# data.raw.list <-
#   list(irfmodel.data.list$pSTATsum %>% dplyr::filter(stimulation %in% stimulations),
#        irfmodel.data.list$irfsum %>% dplyr::filter(stimulation %in% stimulations))

data.raw.list <-
    list(irfmodel.data.list$pSTATsum %>% dplyr::filter(stimulation %in% stimulations))

optimisation.res <- do.call(
  fun.optimisation,
  list(par = par,
       fn = optimise.fun,
       control = list(maxit = maxit,
                      stopfitness = stopfitness,
                      diag.sigma  = FALSE,
                      #keep.best  = TRUE,
                      diag.eigen  = FALSE,
                      diag.pop    = FALSE,
                      diag.value  = FALSE),
       lower = ranges$min[ranges$opt],
       upper = ranges$max[ranges$opt],
       data.raw.list = data.raw.list,
       stimulations = stimulations,
       ranges.base = ranges$base,
       ranges.factor = ranges$factor,
       ranges.opt = ranges$opt,
       model_fun = model_fun_stm_params.ps1,
       data.model.colnames = data.model.colnames)
)

#### plotting ####
irfmodel.path.list$optimisation.id <- "2017-12-28-pSTAT"
#par.ps1 <- par.new
optimisation.res$value
par.ps1 <- optimisation.res$par

ranges <- GetParametersRanges.ps1()
params <- ranges$factor
params[ranges$opt] <- ranges$factor[ranges$opt]*ranges$base[ranges$opt]^par.ps1

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

# g.list[["irf"]] <- ggplot(data = data, 
#                           mapping = aes(
#                             x = stimulation,
#                             y = irf, 
#                             ymin = irf - irf.sd,
#                             ymax = irf + irf.sd,
#                             group = type, 
#                             color = type))+
#   geom_errorbar() +
#   geom_point() +
#   geom_line() +
#   ggtitle("irf") + 
#   do.call(theme_jetka, args = plot.args)
# g.list[["irf"]]

irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")
dir.create(irfmodel.path.list$output.path, 
           recursive = TRUE)

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(irfmodel.path.list$output.path, "IRFmodel-pstat.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = g.list, ncol = 1, nrow = 1))))
saveRDS(
  file = paste(irfmodel.path.list$output.path, "IRFmodel-pstat.RDS", sep = "/"),
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
