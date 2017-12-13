### ###
### IRF1 model 
### ###
#### libraries ####
source("R/optimisation/initialise_optimisation.R")
source("R/scripts/2017-11-05-IRFModel/IRFlibrary.R")
#### organize data IRF ####
irfmodel.data.list <- list()
  
irfmodel.data.list$irf <- 
  poster.data.list$`Gamma-IRF-Nuclei` %>%
  dplyr::filter(time == 720)
irfmodel.data.list$irf <- 
  do.call(
    rbind,
    list(
      irfmodel.data.list$irf,
      poster.data.list$`Gamma-IRF-Nuclei` %>%
        dplyr::filter(time == 0) %>%
        dplyr::mutate(stimulation = 0)
    )
  ) %>%
  dplyr::mutate(response = Intensity_MeanIntensity_Alexa) %>%
  dplyr::mutate(logresponse = log(response))
irfmodel.data.list$irf <- get_equal_data(irfmodel.data.list$irf)

irfmodel.data.list$irfsum <- 
  irfmodel.data.list$irf %>% 
  dplyr::group_by(stimulation) %>% 
  dplyr::summarise(
    logresponse.mean = mean(logresponse),
    logresponse.sd = sd(logresponse)
    ) %>%
  dplyr::mutate(
    logresponse = logresponse.mean
  )

#### organize data pSTAT ####
irfmodel.data.list$pSTAT <- 
  poster.data.list$`Gamma-pStat-Nuclei` %>%
  dplyr::filter(time == 30) %>%
  dplyr::mutate(response = Intensity_MeanIntensity_Alexa) %>%
  dplyr::mutate(logresponse = log(response))

irfmodel.data.list$pSTATsum <- 
  irfmodel.data.list$pSTAT %>% 
  dplyr::group_by(stimulation) %>% 
  dplyr::summarise(
    logresponse.mean = mean(logresponse),
    logresponse.sd = sd(logresponse)
  ) %>%
  dplyr::mutate(
    logresponse = logresponse.mean
  )
irfmodel.data.list$pSTAT <- get_equal_data(irfmodel.data.list$pSTAT)
#### organize data all ####
data.raw.sum <-
  irfmodel.data.list$pSTATsum %>% 
  dplyr::full_join(irfmodel.data.list$irfsum, by = c("stimulation")) %>%
  dplyr::mutate(pstat = logresponse.x, 
                pstat.sd = logresponse.sd.x, 
                irf   = logresponse.y,
                irf.sd = logresponse.sd.y, 
                type  = "data") %>%
  dplyr::select(stimulation, pstat, pstat.sd, irf, irf.sd, type) %>%
  data.frame()
#### model definition ####
model_fun <- function(
  ymax = c(1, 1),
  ymin = c(0, 0),
  theta = c(1, 1, 1, 1),
  sd = c(1, 1),
  stm = 0,
  hn = 1,
  scale = ymax - ymin,
  bck = ymin,
  ...){
  
  hn[2] <- ifelse(length(hn) == 1, hn[1], hn[2])
  
  x <- c(0,0)
  x[1] <- (theta[2] * (stm^hn[1])/(stm^hn[1] + theta[1]^hn[1])) +  theta[3]
  x[2] <- (x[1]^hn[2])/(x[1]^hn[2] + theta[4]^hn[2])
  
  y <- c(0,0)
  y[1] <- scale[1]*x[1] + bck[1]
  y[2] <- scale[2]*x[2] + bck[2]
  
  return(list(x = x, y = y))
}
#### model computation ####
model_fun_stm <- function(
  stimulations,
  sd,
  ...
){
  y.list <- list()
  foreach( stm = stimulations ) %do% {
    result <-  model_fun(
      stm  = stm,
      sd   = sd,
      ...)
    y.list[[as.character(stm)]] <-
    data.frame( 
        pstat = result$y[1],
        pstat.sd = sd[1],
        pstat.model = result$x[1],
        irf   = result$y[2],
        irf.sd = sd[2],
        irf.model   = result$x[2],
        stimulation = stm,
        type = "model")
  }
  data.model <- do.call(rbind, y.list)
}
#### optimisation function ####
optimise.fun <- function(par,
                     stimulations,
                     data.raw.list,
                     ranges.factor,
                     ranges.base,
                     ranges.opt,
                     ...
                     ){
  params <- ranges.factor
  params[ranges.opt] <- ranges.factor[ranges.opt]*ranges.base[ranges.opt]^par
  sd <- params[c(3,4)]
  
  data.model <- model_fun_stm(
    stimulations = stimulations,
    hn = params[c(1, 2)],
    sd = sd,
    theta = params[c(5, 6, 7, 8)], 
    scale = params[c(9, 10)],
    bck = params[c(11, 12)]
  )
  likelihood.list <- foreach( data.i = 1:length(data.raw.list) ) %do% {
      data  <- data.raw.list[[data.i]]
      normalise <- (data %>%
                      dplyr::mutate(normalise = (logresponse^2)/(logresponse.sd^2)) %>%
                      dplyr::summarise(normalise = mean(normalise)))$normalise
      likelihood <-
        (data %>% 
        left_join(
          by =  "stimulation",
          ((data.model %>% 
              dplyr::mutate_(
                "logmodel" = data.model.colnames[data.i]))[, 
                                                           c("logmodel", 
                                                             "stimulation")])) %>%
        #dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(sd[data.i]^2)  + log(sd[data.i]^2))) #%>%
          dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(logresponse.sd^2)) %>%
        #  dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(logresponse^2)) %>%
        dplyr::summarise(likelihood = sum(likelihood)))$likelihood
      return(likelihood/normalise)
      #return(likelihood)         
  }
  likelihood <- sum(likelihood.list[[1]],likelihood.list[[2]])
  # print(likelihood)
  # print("\n")
  # print(par)
  return(as.numeric(likelihood))
}

#### optimisation initialisation ####
# hn_pstat hn_irf sd_pstat sd_irf theta1 theta2 theta3 theta4 scale_pstat scale_irf bck_stat bck_irf
ranges.min <- 2*c(-4, -4, 4, 4, -4, -4, -4, -4, -4, -4, -4, -4)
ranges.max <- 2*c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4)
ranges.base <- c(2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2)
ranges.factor <- c(1, 1, 1, 1, 1, 1, 1, 1, ymax - ymin, ymin)
ranges.opt <- which(ranges.min != ranges.max)
#par <- 0*ranges.opt 
par <- par.new
data.model.colnames <- c("pstat", "irf")

stopfitness <- 0
fun.optimisation = cma_es
maxit <- 1000
#### optimisation running ####
ymin <- c(min(irfmodel.data.list$irfsum$logresponse), 
          min(irfmodel.data.list$pSTATsum$logresponse))
ymax <- c(max(irfmodel.data.list$irfsum$logresponse), 
          max(irfmodel.data.list$pSTATsum$logresponse))
data.raw.list <-
  list(irfmodel.data.list$pSTATsum %>% dplyr::filter(stimulation %in% stimulations),
       irfmodel.data.list$irfsum %>% dplyr::filter(stimulation %in% stimulations))
stimulations.pSTAT <- irfmodel.data.list$pSTAT %>% dplyr::distinct(stimulation) 
stimulations.irf <- irfmodel.data.list$irf %>% dplyr::distinct(stimulation)
stimulations  <- (stimulations.pSTAT %>% dplyr::inner_join(stimulations.irf))$stimulation

optimisation.res <- do.call(
  fun.optimisation,
  list(par = par,
     fn = optimise.fun,
     control = list(maxit = maxit,
                    stopfitness = stopfitness,
                    diag.sigma = TRUE,
                    keep.best = TRUE,
                    diag.eigen = TRUE,
                    diag.pop = TRUE,
                    diag.value = TRUE),
     lower = ranges.min[ranges.opt],
     upper = ranges.max[ranges.opt],
     data.raw.list = data.raw.list,
     stimulations = stimulations,
     ranges.base = ranges.base,
     ranges.factor = ranges.factor,
     ranges.opt = ranges.opt)
  )
    
#### ####drtymoteuszzychchaoskiprparkowaniesejmmorawiecki
irfmodel.path.list$optimisation.id <- "2017-12-14-summarise-id"
optimisation.res$value
par.new <- optimisation.res$par
#par.new <- par
# par.new[c(12)] <- c(-100)
# ranges.factor[c(3,4)] <- c(4,4)
# par.new[c(3,4)] <- c(0,0)
#par.new <- par
#par.new[c(9,10)] <- par.new[c(9,10)] + c(1.3,0.2)
params <- ranges.factor
params[ranges.opt] <- ranges.factor[ranges.opt]*ranges.base[ranges.opt]^par.new
sd <- params[c(3,4)]

data.model <- model_fun_stm(
  stimulations = stimulations,
  hn = params[c(1, 2)],
  sd = sd,
  theta = params[c(5, 6, 7, 8)], 
  scale = params[c(9, 10)],
#  bck = c(0,0)
  bck = params[c(11, 12)] 
)

data <- rbind(data.raw.sum,
              data.model %>% dplyr::select(-c(pstat.model, irf.model)))
data <- data[!is.na(data$pstat), ]
data <- data[!is.na(data$irf), ]
data <- data.frame(data)
data$type <- factor(data$type)

#data <- data[data$stimulation != 0.01, ]

g.list <- list()
g.list[["pstat"]] <- ggplot(data = data, 
       mapping = aes(
         x = stimulation,
         y = pstat, 
         #ymin = pstat - pstat.sd,
         #ymax = pstat + pstat.sd,
         group = type, 
         color = type
         ))+
  #geom_errorbar() +
  geom_point() +
  geom_line() +
  ggtitle("pstat") + 
  do.call(theme_jetka, args = plot.args)
g.list[["pstat"]]

g.list[["irf"]] <- ggplot(data = data, 
       mapping = aes(
         x = stimulation,
         y = irf, 
         #ymin = irf - irf.sd,
         #ymax = irf + irf.sd,
         group = type, 
         color = type))+
  #geom_errorbar() +
  geom_point() +
  geom_line() +
  ggtitle("irf") + 
  do.call(theme_jetka, args = plot.args)
g.list[["irf"]]

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(irfmodel.path.list$output.path, "IRFmodel.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = g.list, ncol = 1, nrow = 1))))

irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")
dir.create(irfmodel.path.list$output.path, 
           recursive = TRUE)
saveRDS(
  file = paste(irfmodel.path.list$output.path, "IRFmodel.RDS", sep = "/"),
  object = list(
    optimisation = optimisation.res,
    likelihood = optimisation.res$value,
    par = par.new, 
    ranges.min = ranges.min,
    ranges.min = ranges.min,
    ranges.opt = ranges.opt,
    #data.raw.list = data.raw.list,
    stimulations = stimulations,
    ranges.base = ranges.base,
    ranges.factor = ranges.factor,
    plots = g.list,
    stopfitness = stopfitness,
    fun.optimisation = fun.optimisation,
    maxit = maxit
  ))   

#### ####
optimise.fun(
  par = par.new,
  data.raw.list = data.raw.list,
  stimulations = stimulations,
  ranges.base = ranges.base,
  ranges.factor = ranges.factor,
  ranges.opt = ranges.opt
)
#### load RDS  ####
irfmodel.path.list$optimisation.id <- "2017-12-13-summarise"
irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")
results <- readRDS(
  file = paste(irfmodel.path.list$output.path, "IRFmodel.RDS", sep = "/"))
  par <- results$par#[ranges.opt]
  
## [1]  0.093594553 -0.001356297 -1.073560733 -0.103776065 -1.733805262 -0.523462295  0.640583619  3.960325670
##  [9] -3.225818492 -0.357628990
  