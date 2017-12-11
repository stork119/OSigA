### ###
### IRF1 model 
### ###
#### libraries ####
source("R/optimisation/initialise_optimisation.R")
#### read input ####
poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"

poster.path.list$output.dir <- "resources/output/poster/"

irfmodel.path.list <- list()
irfmodel.path.list$output.dir <- "resources/output/IRFmodel"
dir.create(irfmodel.path.list$output.dir)

irfmodel.path.list$rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_joined.RDS"
poster.data.list <- readRDS(file = irfmodel.path.list$rds.path)

#### organize data ####
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

irfmodel.data.list$irfsum <- 
  irfmodel.data.list$irf %>% 
  dplyr::group_by(stimulation) %>% 
  dplyr::summarise(
    logresponse = mean(logresponse)
    )

#### ####
irfmodel.data.list$pSTAT <- 
  poster.data.list$`Gamma-pStat-Nuclei` %>%
  dplyr::filter(time == 30) %>%
  dplyr::mutate(response = Intensity_MeanIntensity_Alexa) %>%
  dplyr::mutate(logresponse = log(response))

irfmodel.data.list$pSTATsum <- 
  irfmodel.data.list$pSTAT %>% 
  dplyr::group_by(stimulation) %>% 
  dplyr::summarise(
    logresponse = mean(logresponse)
  )

#### ####
model_fun <- function(
  y = c(0,0),
  ymax = y,
  ymin = y,
  theta = c(1, 1,1),
  sd = c(1, 1),
  stm = 0,
  hn = 1,
  ...){
  hn[2] <- ifelse(length(hn) == 1, hn[1], hn[2])
  y[1] <- (ymax[1] - ymin[1])* ((stm^hn[1])/(stm^hn[1] + theta[1]^hn[1])) + ymin[1]
  y[2] <- (ymax[2] - ymin[2])*theta[3] * ((y[1]^hn[2])/(y[1]^hn[2] + theta[2]^hn[2])) + ymin[2]
  return(y)
}
#### ####
model_fun_stm <- function(
  stimulations,
  ...
){
  y.list <- list()
  foreach( stm = stimulations ) %do% {
    y <-  model_fun(
      stm  = stm,
      ...)
    y.list[[as.character(stm)]] <-
      data.frame( 
        pstat = y[1],
        irf   = y[2],
        stimulation = stm,
        type = "model")
  }
  data.model <- do.call(rbind, y.list)
}
#### ####
ymin <- c(min(irfmodel.data.list$irfsum$logresponse), 
          min(irfmodel.data.list$pSTATsum$logresponse))
ymax <- c(max(irfmodel.data.list$irfsum$logresponse), 
          max(irfmodel.data.list$pSTATsum$logresponse))

data.model <- model_fun_stm(
  stimulations = irfmodel.data.list$pSTATsum$stimulation,
  ymin = ymin,
  ymax = ymax,
  theta = c(1, 9.5, 5),
  hn = c(1,5))

data.raw <-
  irfmodel.data.list$pSTATsum %>% 
  dplyr::full_join(irfmodel.data.list$irfsum, by = c("stimulation")) %>%
  dplyr::mutate(pstat = logresponse.x, 
                irf   = logresponse.y,
                type  = "data") %>%
  dplyr::select(stimulation, pstat, irf, type)

data <- rbind(data.raw, data.model)
data <- data[!is.na(data$pstat), ]
data <- data[!is.na(data$irf), ]
data <- data.frame(data)
data$type <- factor(data$type)

data <- data[data$stimulation != 0.01, ]

ggplot(data = data, 
       mapping = aes(
         x = stimulation,
           y = pstat, 
           group = type, 
           color = type))+
  geom_point() +
  geom_line()


ggplot(data = data, 
       mapping = aes(
         x = stimulation,
         y = irf, 
         group = type, 
         color = type))+
  geom_point() +
  geom_line()