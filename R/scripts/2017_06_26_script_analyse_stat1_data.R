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


#### KZ 45 / 46 ####
source("R/data/normalize_data.R")
#### remopve data KZ45 ####
remove_data.KZ45 <- function(data){
  return(data)
}


remove_data.KZ45 <- function(data){
  return(dplyr::filter(data, !(id %in% c("G07", paste(LETTERS, 12, sep = "")))) %>% dplyr::filter(time.1.1 != 0))
}

fun.remove_data <- remove_data.KZ45
#### AnalyseConstaining ####
AnalyseConstaining <- function(path.input.costaining = "resources/input/2017-03-30-KZ45/raw/",
         filename = "ShrinkedNuclei.csv",
         path.output.costaining = paste(path.input.costaining, strsplit(filename, ".csv")[[1]], "/", sep = "/"),
         y = "Intensity_MeanIntensity_Alexa488",
         fun.remove_data =  function(data){ return(data)}){
  
  dir.create(path = path.output.costaining, showWarnings = FALSE, recursive = TRUE)

  data.list.costaining <- list()

  data.constaining <- 
    read.table(paste(path.input.costaining, filename , sep =  "/"), header = TRUE, sep = ",") %>% 
    normalize_data() %>% 
    fun.remove_data()
  
  gplot.list <- list()
  gplot.list$boxplot <- plot_boxplot_group(data = data.constaining,#x_factor = FALSE,
                                                  save_plot = FALSE, 
                                                  x = "time.1.1",
                                                  y = y, 
                                                  boxplot_group = "id",
                                                  facet_grid_group_y = "stimulation.1.1", 
                                                  facet_grid_group_x = "priming.1.1")# + xlim(c(-15,45))


  gplot.list$boxplot_labels <- 
    gplot.list$boxplot + 
    #ggplot() +
    geom_text(data = data.constaining %>% 
                dplyr::select(id, priming.1.1, time.1.1, stimulation.1.1, Intensity_MeanIntensity_Alexa488) %>% 
                dplyr::group_by(id, priming.1.1, time.1.1, stimulation.1.1) %>% 
                dplyr::summarise(Intensity_MeanIntensity_Alexa488 = mean(Intensity_MeanIntensity_Alexa488)) %>% 
                dplyr::mutate(id_num = (as.numeric(id)-6)/12),  
              aes(x = factor(time.1.1), group = id, label = id, y = 3*Intensity_MeanIntensity_Alexa488)) 
  
    
  do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.output.costaining, "boxplot.pdf", sep = "/"),
                           plot = gplot.list$boxplot)))


  stm.list <- data.constaining %>% dplyr::distinct(stimulation.1.1, time.1.1) %>% dplyr::arrange(stimulation.1.1, time.1.1)
  df.common.list <- list()
  gplot.list$density <- list()
  gplot.list$density_log <- list()
  for(stm.i in 1:nrow(stm.list)){
    stm <- stm.list[stm.i,]$stimulation.1.1
    t <- stm.list[stm.i,]$time.1.1
    #stm <- stm.list[1]
    gplot.list$density[[as.character(stm.i)]] <-
      plot_density(data = data.constaining %>%
                     dplyr::filter(stimulation.1.1 == stm, time.1.1 == t) %>%
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)),
                   x = y,
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste("stimulation", stm, "time", t))
    
    
    #stm <- stm.list[1]
    gplot.list$density_log[[as.character(stm.i)]] <-
      plot_density(data = data.constaining %>%
                     dplyr::filter(stimulation.1.1 == stm, time.1.1 == t) %>%
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)) %>% 
                     dplyr::mutate_(logint = paste("log(", y, ")")),
                   x = "logint",
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste("stimulation", stm, "time", t)) + xlab(paste("log(", y, ")"))
    
    do.call(what = ggsave,
            args = append(plot.args.ggsave,
                          list(filename = paste(path.output.costaining, "density_log.pdf", sep = "/"),
                               plot = marrangeGrob(grobs = gplot.list$density_log, ncol = 1, nrow = 1 ))))
    
    data.nonpriming.q95 <- quantile((data.constaining %>%
                                        dplyr::filter(stimulation.1.1 == stm, time.1.1 == t)  %>% 
                                        dplyr::filter(priming.1.1 == 0))[,y],
                                    probs = 0.95, na.rm = TRUE)
    data.priming.q05    <- quantile(( data.constaining %>%
                                        dplyr::filter(stimulation.1.1 == stm, time.1.1 == t)  %>% 
                                        dplyr::filter(priming.1.1 == 1000))[,y],
                                    probs = 0.05, na.rm = TRUE)
    
    df.common.list[[as.character(stm.i)]] <- data.frame(
      stm  = stm,
      time = t,
      type = y,
      filename = filename,
      priming_percentage = sum(( data.constaining %>%
                                   dplyr::filter(stimulation.1.1 == stm, time.1.1 == t) %>% 
                                   dplyr::filter(priming.1.1 == 1000))[,y] < data_log.nonpriming.q95)/
        nrow(data.constaining %>%
               dplyr::filter(stimulation.1.1 == stm, time.1.1 == t) %>% 
               dplyr::filter(priming.1.1 == 1000)),
      nonpriming_percentage = sum(( data.constaining %>%
                                      dplyr::filter(stimulation.1.1 == stm, time.1.1 == t) %>% 
                                      dplyr::filter(priming.1.1 == 0))[,y]> data_log.priming.q05)/
        nrow(data.constaining %>%
               dplyr::filter(stimulation.1.1 == stm, time.1.1 == t) %>% 
               dplyr::filter(priming.1.1 == 0)))
  }
  # do wyrzucenia wyrzucić kolumnę 7 i 12
  do.call(what = ggsave,
          args = append(plot.args.ggsave,
                        list(filename = paste(path.output.costaining, "density.pdf", sep = "/"),
                             plot = marrangeGrob(grobs = gplot.list$density, ncol = 1, nrow = 1 ))))

  do.call(what = ggsave,
          args = append(plot.args.ggsave,
                        list(filename = paste(path.output.costaining, "density_log.pdf", sep = "/"),
                             plot = marrangeGrob(grobs = gplot.list$density_log, ncol = 1, nrow = 1 ))))
  
  # do wyrzucenia wyrzucić kolumnę 7 i 12
  
  df.common <- do.call(what = rbind, args = df.common.list)
  write.table(x = df.common, 
              file = paste(path.output.costaining, "cells_common.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = TRUE)
  
  return(list(data.constaining = data.constaining,
              gplot.list = gplot.list, 
              df.common = df.common))
  
}
####  raw  ####
constainig.list.raw <- list()
constainig.list.raw[["nuclei"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/raw/",
                     filename = "ShrinkedNuclei.csv",
                     y = "Intensity_MeanIntensity_Alexa488")

constainig.list.raw[["cells.555"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/raw/",
                     filename = "Cells555.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

constainig.list.raw[["cells.masked.555"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/raw/",
                     filename = "CellsMasked555.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

constainig.list.raw[["cells.filtered.555"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/raw/",
                     filename = "CellsFiltered555.csv",
                     y = "Intensity_MeanIntensity_Alexa555")


constainig.list.raw[["cells.488"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/raw/",
                     filename = "Cells488.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

constainig.list.raw[["cells.masked.488"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/raw/",
                     filename = "CellsMasked488.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

constainig.list.raw[["cells.filtered.488"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/raw/",
                     filename = "CellsFiltered488.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

####  ffc ####
constainig.list.ffc <- list()
constainig.list.ffc[["nuclei"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/ffc/",
                     filename = "ShrinkedNuclei.csv",
                     y = "Intensity_MeanIntensity_Alexa488")

constainig.list.ffc[["cells.555"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/ffc/",
                     filename = "Cells555.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

constainig.list.ffc[["cells.masked.555"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/ffc/",
                     filename = "CellsMasked555.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

constainig.list.ffc[["cells.filtered.555"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/ffc/",
                     filename = "CellsFiltered555.csv",
                     y = "Intensity_MeanIntensity_Alexa555")


constainig.list.ffc[["cells.488"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/ffc/",
                     filename = "Cells488.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

constainig.list.ffc[["cells.masked.488"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/ffc/",
                     filename = "CellsMasked488.csv",
                     y = "Intensity_MeanIntensity_Alexa555")

constainig.list.ffc[["cells.filtered.488"]] <- 
  AnalyseConstaining(path.input.costaining = "resources/input/2017-03-30-KZ45/ffc/",
                     filename = "CellsFiltered488.csv",
                     y = "Intensity_MeanIntensity_Alexa555")


#### KZ58 ####
#### KZ 58/ KZ 51 preapre data ####
type <- "raw"
path.input.costaining <- paste("resources/input/2017-07-13-KZ62", type, "/", sep = "/")
for(filename in list.files(path.input.costaining, pattern = ".csv")){
  tryCatch({
    dt <- read.table(paste(path.input.costaining, filename, sep = "/"), header = TRUE, sep = ",") %>% data.table()
    dt <- dt %>% dplyr::mutate(stimulation.1.1 = ifelse(stimulation.1.1 == 1001, 1000, stimulation.1.1)) %>%
      dplyr::mutate(priming.1.1 = ifelse(time.1.1 == 0, 0, stimulation.1.1)) %>%
      dplyr::mutate(stimulation.2.1 = stimulation.1.1, 
                    time.2.1 = time.1.1 ) %>%
      dplyr::mutate(stimulation.1.1 = 0,
                    time.1.1 = 0)
    
    write.table(x = dt, 
                file = paste(path.input.costaining, filename, sep = "/"), 
                col.names = TRUE, 
                row.names = FALSE, 
                sep = ",")},
    error = function(e){
      print(e) 
      print(filename)})
}
#### CompareIFNBpriming####
CompareIFNBpriming <- function(data.constaining,
                       y,
                       path.output.costaining
                       ){
  dir.create(path = path.output.costaining, showWarnings = FALSE, recursive = TRUE)
  gplot.list <- list()
  gplot.list$boxplot <- plot_boxplot_group(data = data.constaining,#x_factor = FALSE,
                                           save_plot = FALSE, 
                                           x = "time.2.1",
                                           y = y, 
                                           boxplot_group = "time.2.1",
                                           facet_grid_group_y = "stimulation.2.1") + ylim(c(0,1000))
  
  
  
  do.call(what = ggsave,
          args = append(plot.args.ggsave,
                        list(filename = paste(path.output.costaining, "boxplot.pdf", sep = "/"),
                             plot = gplot.list$boxplot)))
  
  
  stm.list <- data.constaining %>% dplyr::distinct(stimulation.2.1) %>% dplyr::arrange(stimulation.2.1) %>% data.table()
  df.common.list <- list()
  gplot.list$density <- list()
  gplot.list$density_log <- list()
  for(stm.i in 1:nrow(stm.list)){
    stm <- stm.list[stm.i,]$stimulation.2.1
    #stm <- stm.list[1]
    gplot.list$density[[as.character(stm.i)]] <-
      plot_density(data = data.constaining %>%
                     dplyr::filter(stimulation.2.1 == stm) %>%
                     dplyr::mutate(time.2.1 = factor(time.2.1)),
                   x = y,
                   group = "time.2.1",
                   color = "time.2.1",
                   title = paste("stimulation", stm))
    
    
    #stm <- stm.list[1]
    gplot.list$density_log[[as.character(stm.i)]] <-
      plot_density(data = data.constaining %>%
                     dplyr::filter(stimulation.2.1 == stm) %>%
                     dplyr::mutate(time.2.1 = factor(time.2.1)) %>% 
                     dplyr::mutate_(logint = paste("log(", y, ")")),
                   x = "logint",
                   group = "time.2.1",
                   color = "time.2.1",
                   title = paste("stimulation", stm)) + xlab(paste("log(", y, ")"))
    
  }
  # do wyrzucenia wyrzucić kolumnę 7 i 12
  do.call(what = ggsave,
          args = append(plot.args.ggsave,
                        list(filename = paste(path.output.costaining, "density.pdf", sep = "/"),
                             plot = marrangeGrob(grobs = gplot.list$density, ncol = 1, nrow = 1 ))))
  
  do.call(what = ggsave,
          args = append(plot.args.ggsave,
                        list(filename = paste(path.output.costaining, "density_log.pdf", sep = "/"),
                             plot = marrangeGrob(grobs = gplot.list$density_log, ncol = 1, nrow = 1 ))))
  
  
  
  time.list <- data.constaining %>% dplyr::distinct(time.2.1) %>% dplyr::arrange(time.2.1) %>% data.table()
  gplot.list$density_time <- list()
  gplot.list$density_time_log <- list()
  for(t.i in 1:nrow(time.list)){
    t <- time.list[t.i,]$time.2.1
    #stm <- stm.list[1]
    gplot.list$density_time[[as.character(t.i)]] <-
      plot_density(data = data.constaining %>%
                     dplyr::filter(time.2.1 == t) %>%
                     dplyr::mutate(stimulation.2.1 = factor(stimulation.2.1)),
                   x = y,
                   group = "stimulation.2.1",
                   color = "stimulation.2.1",
                   title = paste("time", t))
    
    
    #stm <- stm.list[1]
    gplot.list$density_time_log[[as.character(t.i)]] <-
      plot_density(data = data.constaining %>%
                     dplyr::filter(time.2.1 == t) %>%
                     dplyr::mutate(stimulation.2.1 = factor(stimulation.2.1)) %>%
                     dplyr::mutate_(logint = paste("log(", y, ")")),
                   x = "logint",
                   group = "stimulation.2.1",
                   color = "stimulation.2.1",
                   title = paste("time", t)) + xlab(paste("log(", y, ")"))
    
  }
  
  do.call(what = ggsave,
          args = append(plot.args.ggsave,
                        list(filename = paste(path.output.costaining, "density_time.pdf", sep = "/"),
                             plot = marrangeGrob(grobs = gplot.list$density_time, ncol = 1, nrow = 1 ))))
  
  do.call(what = ggsave,
          args = append(plot.args.ggsave,
                        list(filename = paste(path.output.costaining, "density_time_log.pdf", sep = "/"),
                             plot = marrangeGrob(grobs = gplot.list$density_time_log, ncol = 1, nrow = 1 ))))
}

#### KZ 58 cells ####
type <- "ffc"
type.list <- c("ffc", "raw")
for(type in type.list){
path.input.costaining <- paste("resources/input/2017-07-13-KZ62/", type, "/", sep = "/")

filename.list <- list.files(path.input.costaining, pattern = ".csv", recursive = FALSE)
for(filename in filename.list[c(1,2,5,7)]){
  path.output.costaining <- paste(path.input.costaining, strsplit(filename, ".csv")[[1]], "/", sep = "/")
  y <- "Intensity_MeanIntensity_Alexa555"
  fun.remove_data  <-  function(data){ return(data)}
  dir.create(path = path.output.costaining, showWarnings = FALSE, recursive = TRUE)
  data.constaining <- 
    read.table(paste(path.input.costaining, filename , sep =  "/"), header = TRUE, sep = ",") %>% 
    normalize_data() %>% 
    fun.remove_data()
  CompareIFNBpriming(data.constaining = data.constaining, y = y, path.output.costaining = path.output.costaining)
     # do wyrzucenia wyrzucić kolumnę 7 i 12
    
}}
#### KZ 58 cytoplasm  ####
# 
# fun.remove_data <- function(data){
#   return(data)
# }
# filename <- "CytoplasmFiltered555.csv"
# type <- "ffc"
# path.input.costaining <- paste("resources/input/2017-06-29-KZ58", type, "/", sep = "/")
# data.cells <- 
#   read.table(paste(path.input.costaining, "Cells555.csv" , sep =  "/"), header = TRUE, sep = ",") %>% 
#   normalize_data() %>% 
#   fun.remove_data() %>% 
#   data.table()
# data.cells.filtered <- 
#   read.table(paste(path.input.costaining, "CellsFiltered555.csv" , sep =  "/"), header = TRUE, sep = ",") %>% 
#   normalize_data() %>% 
#   fun.remove_data() %>% 
#   data.table()
# data.nuclei <- 
#   read.table(paste(path.input.costaining, "Nuclei.csv" , sep =  "/"), header = TRUE, sep = ",") %>% 
#   normalize_data() %>% 
#   fun.remove_data() %>% 
#   data.table()
# 
# 
# 
# data.cells_nuclei <-  
#   data.cells %>% 
#   dplyr::ungroup() %>% 
#   dplyr::select_("well.name","ObjectNumber", "Parent_Nuclei") %>% 
#   left_join(data.nuclei,
#             by = c("well.name" = "well.name",
#                    "Parent_Nuclei" = "ObjectNumber")) %>%
#   data.table()
# 
# data.cytoplasm <- 
#   data.cells.filtered %>%
#   left_join((data.cells_nuclei), 
#             by = c("well.name" = "well.name",
#                    "Parent_Cells555" = "ObjectNumber")) %>% 
#   data.table()
# y <- "Intensity_MeanIntensity_Alexa555"
# 
# data.cytoplasm <- 
#   data.cytoplasm %>% 
#   dplyr::mutate_("IntensityMeanCytoplasm" = 
#                    paste("(Intensity_IntegratedIntensity_Alexa555.x - Intensity_IntegratedIntensity_Alexa555.y)",
#                          "(AreaShape_Area.x - AreaShape_Area.y)",
#                          sep = "/"),
#                  "IntensityIntegratedCytoplasm" = 
#                    paste("(Intensity_IntegratedIntensity_Alexa555.x - Intensity_IntegratedIntensity_Alexa555.y)", 
#                          sep = "/")) %>% 
#   dplyr::mutate_("IntensityRatioNucleiCytoplasm" = 
#                    paste("(Intensity_IntegratedIntensity_Alexa555.y)",
#                          "(IntensityIntegratedCytoplasm)",
#                          sep = "/"))  
# 
# data.cytoplasm <- data.cytoplasm %>% mutate(stimulation.1.1 = stimulation.1.1.x, time.1.1 = time.1.1.x)
# 
# write.table(x = data.cytoplasm, 
#             file = paste(path.input.costaining, filename, sep = "/"),
#             sep = ",",
#             row.names = FALSE, 
#             col.names = TRUE)
# 
# y = "IntensityMeanCytoplasm"
# CompareIFNBpriming(data.constaining = data.cytoplasm,
#            y = y,
#            path.output.costaining = paste(path.input.costaining, strsplit(filename, ".csv")[[1]], y, "/", sep = "/"))
# 
# y = "IntensityRatioNucleiCytoplasm"
# CompareIFNBpriming(data.constaining = data.cytoplasm,
#            y = y,
#            path.output.costaining = paste(path.input.costaining, strsplit(filename, ".csv")[[1]], y, "/", sep = "/"))
# 


#### 2017-07-06 - analysis ####
path.input <- "resources/input/"
experiment.list <-c("2017-06-29-KZ58", "2017-04-06-KZ46", "2017-03-30-KZ45", "2017-05-26-KZ51", "2017-07-13-KZ62")
gplot.list <- list()

i <- 5

gplot.list[[as.character(i)]]  <-list()
experiment <- experiment.list[i]
path.input.costaining <- paste(path.input, experiment, sep = "/")
type.list <- list.dirs(path.input.costaining, recursive = FALSE, full.names = FALSE)
type.i <- 1
for(type.i in 1:length(type.list)){
  gplot.list[[as.character(i)]][[type]]  <- list()
  y <- "Intensity_MeanIntensity_Alexa555"
  
  type <- type.list[type.i]
  path.input.costaining.type <- paste( path.input.costaining, type, sep = "/")
  
  data.cells <- 
    read.table(paste(path.input.costaining.type, "Cells555.csv" , sep =  "/"), header = TRUE, sep = ",") %>% 
    normalize_data() %>% 
    fun.remove_data() %>% 
    data.table()
  data.cells.filtered <- 
    read.table(paste(path.input.costaining.type, "CellsFiltered555.csv" , sep =  "/"), header = TRUE, sep = ",") %>% 
    normalize_data() %>% 
    fun.remove_data() %>% 
    data.table()
  data.nuclei <- 
    read.table(paste(path.input.costaining.type, "Nuclei.csv" , sep =  "/"), header = TRUE, sep = ",") %>% 
    normalize_data() %>% 
    fun.remove_data() %>% 
    data.table()
  
  if(!("stimulation.2.1" %in% colnames(data.cells))){
    data.cells <-   data.cells %>% dplyr::mutate(stimulation.2.1 = priming.1.1, time.2.1 = 24)
    data.cells.filtered <-   data.cells.filtered %>% dplyr::mutate(stimulation.2.1 = priming.1.1, time.2.1 = 24)
    data.nuclei <-   data.nuclei %>% dplyr::mutate(stimulation.2.1 = priming.1.1, time.2.1 = 24)
  }
  
    data.cells.summarise <- 
      data.cells.filtered %>%
    #  dplyr::group_by(priming.1.1, stimulation.1.1, time.1.1) %>%
      dplyr::group_by(priming.1.1, stimulation.1.1, time.1.1, stimulation.2.1, time.2.1) %>%
      dplyr::mutate_(intensity = y) %>%
      dplyr::summarise(intensity_mean = mean(intensity), 
                      intensity_median = median(intensity))
  
    priming <- data.cells.summarise %>% 
      dplyr::filter(priming.1.1 == 0) %>% 
      dplyr::ungroup() %>% 
      dplyr::summarise(intensity_mean = mean(intensity_mean),
                       intensity_median = mean(intensity_median))
    priming.mean <- priming$intensity_mean
    priming.median <- priming$intensity_median
    
    data.cells.summarise <- 
      data.cells.summarise  %>%
      dplyr::mutate(intensity_mean_ratio = intensity_mean/priming.mean,
                    intensity_median_ratio = intensity_median/priming.median)
    
    
    write.table(x = data.cells.summarise, file = paste(path.input.costaining.type,
                                                       "CellsFiltered555",
                                                       "data_ratio.csv", sep = "/"),
                sep = ",", 
                row.names = FALSE, 
                col.names = TRUE)
  
    write.table(x = data.cells.summarise, file = paste(path.input.costaining.type,
                                                       "data_ratio.csv", sep = "/"),
                sep = ",", 
                row.names = FALSE, 
                col.names = TRUE)
    
    ### cytoplasm
    
    data.cells_nuclei <-  
      data.cells %>% 
      dplyr::ungroup() %>% 
      dplyr::select_("well.name","ObjectNumber", "Parent_Nuclei") %>% 
      left_join(data.nuclei,
                by = c("well.name" = "well.name",
                       "Parent_Nuclei" = "ObjectNumber")) %>%
      data.table()
     
    
    data.cytoplasm <- 
      data.cells.filtered %>%
      left_join((data.cells_nuclei), 
                by = c("well.name" = "well.name",
                       "Parent_Cells555" = "ObjectNumber")) %>% 
      data.table()
    
    y <- "Intensity_MeanIntensity_Alexa555"
    
    data.cytoplasm <- 
      data.cytoplasm %>% 
      dplyr::mutate_("IntensityMeanCytoplasm" = 
                       paste("(Intensity_IntegratedIntensity_Alexa555.x - Intensity_IntegratedIntensity_Alexa555.y)",
                             "(AreaShape_Area.x - AreaShape_Area.y)",
                             sep = "/"),
                     "IntensityIntegratedCytoplasm" = 
                       paste("(Intensity_IntegratedIntensity_Alexa555.x - Intensity_IntegratedIntensity_Alexa555.y)", 
                             sep = "/")) %>% 
      dplyr::mutate_("IntensityRatioNucleiCytoplasm" = 
                       paste("(Intensity_IntegratedIntensity_Alexa555.y)",
                             "(IntensityIntegratedCytoplasm)",
                             sep = "/"))  %>%
      dplyr::mutate_("IntensityPercentageNucleiCytoplasm" = 
                       paste("100*(Intensity_IntegratedIntensity_Alexa555.y)",
                             "(IntensityIntegratedCytoplasm + Intensity_IntegratedIntensity_Alexa555.y)",
                             sep = "/"))  
    
    data.cytoplasm <- data.cytoplasm %>% 
      dplyr::mutate(priming.1.1 = priming.1.1.x, 
              stimulation.1.1 = stimulation.1.1.x, 
             stimulation.2.1 = stimulation.2.1.x,
             time.1.1 = time.1.1.x,
             time.2.1 = time.2.1.x) %>% 
      data.table()
    
    
    
    y <- "Intensity_MeanIntensity_Alexa555"
    gplot.list[[as.character(i)]][[type]][["cells"]] <-
      plot_density(data = data.cells.filtered %>%
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)),
                   x = y,
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "cells")) +
      xlim(c(0,1500))
    
    gplot.list[[as.character(i)]][[type]][["cells_log"]] <-
      plot_density(data = data.cells.filtered %>% 
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)) %>% 
                     dplyr::mutate_(logint = paste("log(", y, ")")),
                   x = "logint",
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "cells", "log")) + 
      xlim(c(2,8))
    
    y <- "Intensity_MeanIntensity_Alexa555"
    gplot.list[[as.character(i)]][[type]][["nuclei"]] <-
      plot_density(data = data.nuclei %>%
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)),
                   x = y,
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "nuclei")) +
      xlim(c(0,1500))
    
    gplot.list[[as.character(i)]][[type]][["nuclei_log"]] <-
      plot_density(data = data.nuclei %>% 
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)) %>% 
                     dplyr::mutate_(logint = paste("log(", y, ")")),
                   x = "logint",
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "nuclei", "log")) + 
      xlim(c(2,8))
    
  
    y <- "IntensityMeanCytoplasm"
    gplot.list[[as.character(i)]][[type]][["cytoplasm"]] <-
      plot_density(data = data.cytoplasm %>%
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)),
                   x = y,
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "cytoplasm")) +
      xlim(c(0,1500))
    
    gplot.list[[as.character(i)]][[type]][["cytoplasm_log"]] <-
      plot_density(data = data.cytoplasm %>% 
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)) %>% 
                     dplyr::mutate_(logint = paste("log(", y, ")")),
                   x = "logint",
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "cytoplasm", "log")) + 
      xlim(c(2,8))
      
    
    y <- "IntensityRatioNucleiCytoplasm"
    gplot.list[[as.character(i)]][[type]][["ratio"]] <-
      plot_density(data = data.cytoplasm %>%
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)),
                   x = y,
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "ratio nuclei/cytoplasm")) +
      xlim(c(0,20))
    
    gplot.list[[as.character(i)]][[type]][["ratio_log"]] <-
      plot_density(data = data.cytoplasm %>% 
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)) %>% 
                     dplyr::mutate_(logint = paste("log(", y, ")")),
                   x = "logint",
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "ratio nuclei/cytoplasm", "log")) + 
      xlim(c(-2,7))
    
    
    y <- "IntensityPercentageNucleiCytoplasm"
    gplot.list[[as.character(i)]][[type]][["percentage"]] <-
      plot_density(data = data.cytoplasm %>%
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)),
                   x = y,
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "Percentage in nuclei. nuclei/(nuclei+cytoplasm)")) 
    
    gplot.list[[as.character(i)]][[type]][["percentage_log"]] <-
      plot_density(data = data.cytoplasm %>% 
                     dplyr::mutate(priming.1.1 = factor(priming.1.1)) %>% 
                     dplyr::mutate_(logint = paste("log(", y, ")")),
                   x = "logint",
                   group = "priming.1.1",
                   color = "priming.1.1",
                   title = paste(experiment, "Percentage in nuclei. nuclei/(nuclei+cytoplasm)", "log")) 
    
    do.call(what = ggsave,
            args = append(plot.args.ggsave,
                          list(filename = paste(path.input.costaining.type, "density_comparison.pdf", sep = "/"),
                               plot = marrangeGrob(grobs = gplot.list[[as.character(i)]][[type]],
                                                   ncol = 1,
                                                   nrow = 1 ))))
}