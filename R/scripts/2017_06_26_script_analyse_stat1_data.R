### ### ###
### 2017_05_27
### script_analyse_data
### ### ###
source("R/libraries.R")
source("R/libraries/library_data_manipulation.R")
source("R/graphics/libraries.R")
source("R/libraries/statistics/likelihood_lmvn.R")
source("R/libraries/statistics/lmvn.R")

#### initialise ###
path.stat1.list <- list()
path.stat1.list$data.input <- "resources/input/2015-04-23-KA00/ffc/data_quantify/"
path.stat1.list$output <- "resources/output/stat1/"
gplot.list$stat1 <- list()
#### ####

dir.create(path.stat1.list$output, showWarnings = FALSE, recursive = TRUE)

list.files(path.stat1.list$data.input)
data.list$stat1.cells <- 
  read.table(
    file = paste(
      path.stat1.list$data.input,
      "remove-Cells.csv",
      sep = "/"), 
    sep = ",", 
    header = TRUE) %>% 
  data.table()



#data.list$stat1.cells %>% dplyr::filter(protein.1.1 == 1) %>% dplyr::distinct(time.1.1)
gplot.list$stat1$boxplot <- ggplot(data = data.list$stat1.cells %>% dplyr::filter(protein.1.1 == 1),
       mapping = aes(x = factor(time.1.1), y = Intensity_MeanIntensity_Alexa, 
                     group = position.name, color = factor(prestimulation.1.1))) + 
  geom_boxplot() +
  facet_grid(~stimulation.1.1) +
  do.call(theme_jetka, plot.args)  +
  ggtitle("Boxplots compare")
  

data_log <- data.list$stat1.cells %>% dplyr::filter(protein.1.1 == 1, stimulation.1.1 <= 5) %>% dplyr::mutate(prestimulation.1.1 = factor(prestimulation.1.1), Intensity_MeanIntensity_Alexa_log = log(Intensity_MeanIntensity_Alexa))
gplot.list$stat1$density.log <- plot_density(data = data_log,
             group = "prestimulation.1.1", 
             compare_to_all = FALSE,
             x = "Intensity_MeanIntensity_Alexa_log") + ggtitle("logdensity")

gplot.list$stat1$density <- plot_density(data = data_log,
                                         group = "prestimulation.1.1", 
                                         compare_to_all = FALSE,
                                         x = "Intensity_MeanIntensity_Alexa") + ggtitle("density")


data_log.nonpriming.q95 <- quantile((data_log %>% dplyr::filter(prestimulation.1.1 == 0))$Intensity_MeanIntensity_Alexa_log, probs = 0.95)
data_log %>% filter(prestimulation.1.1 == 1000, Intensity_MeanIntensity_Alexa_log > data_log.nonpriming.q95)
g_dots <- plot_density(data = data_log %>% filter(prestimulation.1.1 == 1000, Intensity_MeanIntensity_Alexa_log > data_log.nonpriming.q95),
             group = "prestimulation.1.1", 
             compare_to_all = FALSE,
             x = "Intensity_MeanIntensity_Alexa")



do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.stat1.list$output, "compare_prm_noprm.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = gplot.list$stat1, ncol = 1, nrow =1 ))))


#### how many cells have similar stat level ####
data_log <- data.list$stat1.cells %>% dplyr::filter(protein.1.1 == 1, stimulation.1.1 <= 5) %>% dplyr::mutate(prestimulation.1.1 = factor(prestimulation.1.1), Intensity_MeanIntensity_Alexa_log = log(Intensity_MeanIntensity_Alexa))
data_log.nonpriming.q95 <- quantile((data_log %>% dplyr::filter(prestimulation.1.1 == 0))$Intensity_MeanIntensity_Alexa_log, probs = 0.95)

data_log.priming.q05 <- quantile((data_log %>% dplyr::filter(prestimulation.1.1 == 1000))$Intensity_MeanIntensity_Alexa_log, probs = 0.05)


df <- data.frame(priming_percentage = sum((data_log %>% dplyr::filter(prestimulation.1.1 == 1000))$Intensity_MeanIntensity_Alexa_log < data_log.nonpriming.q95)/nrow(data_log %>% dplyr::filter(prestimulation.1.1 == 1000)),
                 nonpriming_percentage = sum((data_log %>% dplyr::filter(prestimulation.1.1 == 0))$Intensity_MeanIntensity_Alexa_log > data_log.priming.q05)/nrow(data_log %>% dplyr::filter(prestimulation.1.1 == 0)))

write.table(x = df, file = paste(path.stat1.list$output, "cells_with_common_stat1.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = TRUE)


#### compare densities of experimental data  ####


data.compare.intersection <- data.list$data.exp %>% 
  data.table() %>% left_join( 
    y = (data.list$data.exp %>% 
           data.table() %>%
           dplyr::group_by(priming, stimulation, time) %>%
           dplyr::summarise(q05 = quantile(intensity, probs = 0.05), q95 = quantile(intensity, probs = 0.95))),
    by = c("stimulation", "time")) %>% 
  dplyr::filter(priming.x != priming.y) %>% 
  data.table() %>% 
  dplyr::mutate(q05_test = q05 < intensity, q95_test = q95 > intensity) %>%
  dplyr::group_by(priming.x, stimulation, time) %>% 
  dplyr::summarise(q05_test_sum = mean(q05_test), q95_test_sum = mean(q95_test) ) %>% 
  dplyr::mutate(priming = priming.x) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(priming, stimulation, time) %>% 
  dplyr::select(-priming.x)
              
write.table(file = paste(path.data.output, "data_common.csv", sep = ""),x = test, sep = ",", row.names = FALSE, col.names = TRUE)

data.compare.intersection.melt <-  reshape2::melt(data = data.compare.intersection, id.vars = c("priming", "stimulation", "time")) %>% 
  dplyr::filter((variable == "q05_test_sum" & priming == 0) | (variable == "q95_test_sum" & priming == 1000))

gplot.intersection.data <-  ggplot(data = data.compare.intersection.melt, mapping = aes(x = factor(time), y = value, group = variable, color = variable )) +
  geom_point() +
  facet_grid(priming~stimulation) +
  do.call(theme_jetka, plot.args) 
  
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "data_intersection.pdf", sep = ""),
                           plot = gplot.intersection.data)))

    
#### ####
time.list <- (data.list$data.exp %>% data.table() %>% distinct(time))$time 

gplot.list$density_compare <- list()
gplot.list$density_compare_log <- list()
for(t in time.list){
  gplot.list$density_compare[[paste(t, "nonpriming",   sep = "")]]  <- 
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(priming == 0, time == t) %>%
                   dplyr::mutate(stimulation = factor(stimulation)),
                 group = "stimulation", 
                 compare_to_all = FALSE,
                 x = "intensity") + 
    xlim(c(0,1000)) +  
    ggtitle(paste(t, "nonpriming",   sep = " "))
  
  gplot.list$density_compare[[paste(t, "priming",   sep = "")]]  <- 
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(priming == 1000, time == t) %>%
                   dplyr::mutate(stimulation = factor(stimulation)),
                 group = "stimulation", 
                 compare_to_all = FALSE,
                 x = "intensity") + 
    xlim(c(0,1000)) +  
    ggtitle(paste(t, "priming",   sep = " "))
  gplot.list$density_compare_log[[paste(t, "nonpriming",   sep = "")]]  <- 
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(priming == 0, time == t) %>%
                   dplyr::mutate(stimulation = factor(stimulation)),
                 group = "stimulation", 
                 compare_to_all = FALSE,
                 x = "logintensity") + 
    xlim(c(3.5,7.5)) +  
    ggtitle(paste(t, "nonpriming",   sep = " "))
  gplot.list$density_compare_log[[paste(t, "priming",   sep = "")]]  <- 
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(priming == 1000, time == t) %>%
                   dplyr::mutate(stimulation = factor(stimulation)),
                 group = "stimulation", 
                 compare_to_all = FALSE,
                 x = "logintensity") + 
    xlim(c(3.5,7.5)) +  
    ggtitle(paste(t, "priming",   sep = " "))
}


do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "density_compare_log.pdf", sep = ""),
                           plot = marrangeGrob( grobs = gplot.list$density_compare_log, ncol = 2, nrow = 1))))

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "density_compare.pdf", sep = ""),
                           plot = marrangeGrob( grobs = gplot.list$density_compare, ncol = 2, nrow = 1))))

## ##

gplot.list$density_compare_stm <- list()
gplot.list$density_compare_stm_log <- list()
time.list <- (data.list$data.exp %>% data.table() %>% distinct(time))$time 
stm.list  <- (data.list$data.exp %>% data.table() %>% distinct(stimulation))$stimulation
cond.list <- expand.grid(t = time.list, stm = stm.list)
for(i in 1:nrow(cond.list)){
  
  t <- cond.list[i,]$t
  stm <- cond.list[i,]$stm
  
  gplot.list$density_compare_stm[[i]]  <-
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(time == t, stimulation == stm) %>% 
                   dplyr::mutate(priming = factor(priming)),
                 group = "priming", 
                 compare_to_all = FALSE,
                 x = "intensity") + 
    xlim(c(0,1000)) +  
    ggtitle(paste(stm, t,   sep = " "))
  
  
  gplot.list$density_compare_stm_log[[i]]  <-
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(time == t, stimulation == stm) %>% 
                   dplyr::mutate(priming = factor(priming)),
                 group = "priming", 
                 compare_to_all = FALSE,
                 x = "intensity") + 
    xlim(c(0,1000)) +  
    ggtitle(paste(stm, t,   sep = " "))
  
}

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "density_compare_stm.pdf", sep = ""),
                           plot = marrangeGrob( grobs = gplot.list$density_compare_stm, ncol = 3, nrow = 2))))


do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "density_compare_stm_log.pdf", sep = ""),
                           plot = marrangeGrob( grobs = gplot.list$density_compare_stm_log, ncol = 3, nrow = 2))))

## density_compare_stm ##
stm.list <- (data.list$data.exp %>% data.table() %>% distinct(stimulation))$stimulation

gplot.list$density_compare_stm <- list()
gplot.list$density_compare_stm_log <- list()
for(stm in stm.list){
  gplot.list$density_compare_stm[[paste(stm, "nonpriming",   sep = "")]]  <- 
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(priming == 0, stimulation == stm)
                 %>% dplyr::mutate(time = factor(time)),
                 group = "time", 
                 compare_to_all = FALSE,
                 x = "intensity") + 
    xlim(c(0,1000)) +  
    ggtitle(paste(stm, "nonpriming",   sep = " "))
  
  gplot.list$density_compare_stm[[paste(stm, "priming",   sep = "")]]  <- 
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(priming == 1000, stimulation == stm)%>% dplyr::mutate(time = factor(time))
                 ,
                 group = "time", 
                 compare_to_all = FALSE,
                 x = "intensity") + 
    xlim(c(0,1000)) +  
    ggtitle(paste(stm, "priming",   sep = " "))
  gplot.list$density_compare_stm_log[[paste(stm, "nonpriming",   sep = "")]]  <- 
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(priming == 0, stimulation == stm) %>% dplyr::mutate(time = factor(time))
                 ,
                 group = "time", 
                 compare_to_all = FALSE,
                 x = "logintensity") + 
    xlim(c(3.5,7.5)) +  
    ggtitle(paste(stm, "nonpriming",   sep = " "))
  gplot.list$density_compare_stm_log[[paste(stm, "priming",   sep = "")]]  <- 
    plot_density(data = data.list$data.exp %>% 
                   data.table() %>%
                   dplyr::filter(priming == 1000, stimulation == stm) %>% dplyr::mutate(time = factor(time))
                 ,
                 group = "time", 
                 compare_to_all = FALSE,
                 x = "logintensity") + 
    xlim(c(3.5,7.5)) +  
    ggtitle(paste(stm, "priming",   sep = " "))
}


do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "density_compare_stm_log.pdf", sep = ""),
                           plot = marrangeGrob( grobs = gplot.list$density_compare_stm_log, ncol = 2, nrow = 1))))

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "density_compare_stm.pdf", sep = ""),
                           plot = marrangeGrob( grobs = gplot.list$density_compare_stm, ncol = 2, nrow = 1))))
