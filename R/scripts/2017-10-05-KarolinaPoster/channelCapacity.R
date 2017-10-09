###
### 
###

# install.packages("openssl")
# install.packages("curl")
# install.packages("git2r")
# install.packages("libcurl", dependencies = TRUE)
# install.packages("devtools")
# library(devtools)
source("R/optimisation/initialise_optimisation.R")
#install.packages("e1071")
library(e1071)
library(CapacityLogReg)

# rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc.RDS"
# poster.data.list <- readRDS(file = rds.path)
# poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"
poster.path.list$output.dir <- "resources/output/poster/joined_zeroadded/"

#### plots by wells ####
gplot.list <- list()
r <-foreach(poster.label = labels(poster.data.list)) %do% {
  #poster.label <- labels(poster.data.list)[3]
  print(poster.label)
  data <- poster.data.list[[poster.label]] %>% data.frame()
  
  if("Intensity_MeanIntensity_Alexa" %in% colnames(data)){
    col_response <- "Intensity_MeanIntensity_Alexa"
  } else if("Intensity_MeanIntensity_Alexa488" %in% colnames(data)){
    col_response <- "Intensity_MeanIntensity_Alexa488"
  } else {
    col_response <- "Intensity_MeanIntensity_Alexa555"
  }
  
  path <- paste(poster.path.list$output.dir, poster.label, sep = "/")
  dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
  
  gplot.list[[poster.label]]  <- 
    plot_boxplot_group(
      data = data, 
      x = "well.name", 
      y = col_response, 
      save_plot = FALSE,
      ylim_max_const = TRUE,
      plot_title = poster.label,
      ylim_max = quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]])
  return()
}

ggsave(filename = paste(poster.path.list$output.dir, paste("experiments", ".pdf", sep = ""), sep = "/"), 
       plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
       width = plot.args$width,
       height =plot.args$height,
       useDingbats = plot.args$useDingbats)


#### channel capacity compuation ####
#no_cores <- 6
#registerDoParallel(no_cores)
foreach(poster.label = labels(poster.data.list)) %do% {
  #poster.label <- labels(poster.data.list)[2]
  tryCatch({
    print(poster.label)
    data <- poster.data.list[[poster.label]] %>% data.frame()
    
    if("time" %in% colnames(data)){
      col_time <- "time"
    } else if ("time.2.1" %in% colnames(data)){
      col_time <- "time.2.1"
    } else {
      col_time <- "time.1.1"
    }
    
    if("stimulation" %in% colnames(data)){
      col_stimulation <- "stimulation"
    } else if ("stimulation.2.1" %in% colnames(data)){
      col_stimulation <- "stimulation.2.1"
    } else {
      col_stimulation <- "stimulation.1.1"
    }
    
    
    if("Intensity_MeanIntensity_Alexa" %in% colnames(data)){
      col_response <- "Intensity_MeanIntensity_Alexa"
    } else if("Intensity_MeanIntensity_Alexa488" %in% colnames(data)){
      col_response <- "Intensity_MeanIntensity_Alexa488"
    } else {
      col_response <- "Intensity_MeanIntensity_Alexa555"
    }
    
    path <- paste(poster.path.list$output.dir, poster.label, sep = "/")
    dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
    gplot.list <- list()
    gplot.list$ids  <- 
      plot_boxplot_group(
        data = data, 
        x = "well.name", 
        y = col_response, 
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.label,
        ylim_max = quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]])
    
    gplot.list$x_time_mean  <- 
      plot_boxplot_group(
        data = data, 
        x = col_time, 
        y = col_response, 
        facet_grid_group_y = col_stimulation,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.label,
        ylim_max = max(quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]], 1500))
  
    gplot.list$x_time  <- 
      plot_boxplot_group(
        data = data, 
        x = col_time, 
        y = col_response, 
        boxplot_group = "well.name",
        facet_grid_group_y = col_stimulation,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.label,
        ylim_max = max(quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]], 1500))
    
    gplot.list$x_stimulation_mean  <- 
      plot_boxplot_group(
        data = data, 
        x =col_stimulation, 
        y = col_response, 
        facet_grid_group_y =  col_time,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.label,
        ylim_max = max(quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]], 1500))
      
    gplot.list$x_stimulation  <- 
      plot_boxplot_group(
        data = data, 
        x = col_stimulation, 
        y = col_response, 
        boxplot_group = "well.name",
        facet_grid_group_y =  col_time,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.label,
        ylim_max = max(quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]], 1500))
    
    ggsave(filename = paste(path, paste(poster.label, "-experiment", ".pdf", sep = ""), sep = "/"), 
           plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
           width = plot.args$width,
           height =plot.args$height,
           useDingbats = plot.args$useDingbats)
    cc.list <- list()
    no_cores <- 6
    registerDoParallel(no_cores)
    cc.list <- foreach(t = (data %>% dplyr::distinct_(col_time))[,col_time]) %dopar% {
      print(paste(poster.label, t))
    #for(t in (data %>% dplyr::distinct_(col_time))[,col_time]){
      output_path <- paste(path, "channel_capacity", t, "/", sep = "/")
      dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
      cc.output <- capacity_logreg_main(data %>% dplyr::filter_(paste(col_time, "==", t)),
                                                         graphs = TRUE,plot_width = plot.args$width,
                                                         plot_height = plot.args$height,
                                                          
                                      signal = col_stimulation,
                                      response = col_response,
                                      #                                forumla_string =
                                      cc_maxit = 50,
                                      lr_maxit = 2000,
                                      output_path = output_path)
    #  cc.list[[as.character(t)]] <- 
      return(list(time = t, cc = cc.output$cc))
    }
    stopImplicitCluster()
    cc.df <- do.call(rbind, cc.list)  
    write.table(x = cc.df, 
                file = paste(path, "channel_capacity.csv", sep = "/"), 
                sep = ",",
                row.names = FALSE, 
                col.names = TRUE)
    print(cc.df)
  }, error = function(e){print(e)})
}
#stopImplicitCluster()

#### channel capacity plot ####
gplot.list <- list()
r <- foreach(poster.label = labels(poster.data.list)) %do% {
  #poster.label <- labels(poster.data.list)[3]
  print(poster.label)
  data <- poster.data.list[[poster.label]] %>% data.frame()
  
  if("time" %in% colnames(data)){
    col_time <- "time"
  } else if ("time.2.1" %in% colnames(data)){
    col_time <- "time.2.1"
  } else {
    col_time <- "time.1.1"
  }
  
  if("stimulation" %in% colnames(data)){
    col_stimulation <- "stimulation"
  } else if ("stimulation.2.1" %in% colnames(data)){
    col_stimulation <- "stimulation.2.1"
  } else {
    col_stimulation <- "stimulation.1.1"
  }
  
  if("Intensity_MeanIntensity_Alexa" %in% colnames(data)){
    col_response <- "Intensity_MeanIntensity_Alexa"
  } else if("Intensity_MeanIntensity_Alexa488" %in% colnames(data)){
    col_response <- "Intensity_MeanIntensity_Alexa488"
  } else {
    col_response <- "Intensity_MeanIntensity_Alexa555"
  }
  path <- paste(poster.path.list$output.dir, poster.label, sep = "/")
  cc.df <- read.table(paste(path,"channel_capacity.csv", sep =  "/"), header = TRUE, sep = ",")
  stimulation.list <- unique(data[,col_stimulation])
  cc.df$stimulation <- floor(length(stimulation.list)/2)
  cc.df[,col_time] <- cc.df$time
  
  gplot.list[[poster.label]]  <- 
    plot_boxplot_group(
      data = data, 
      x = col_stimulation, 
      y = col_response, 
      facet_grid_group_y =  col_time,
      save_plot = FALSE,
      ylim_max_const = TRUE,
      plot_title = poster.label,
      ylim_max = max(quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]], 1500)) +
    geom_text(data = cc.df, mapping = aes(label = round(cc,2), y = 1200), size = 20, color = "red")

  return()
}
ggsave(filename = paste(poster.path.list$output.dir, paste("channel_capacity", ".pdf", sep = ""), sep = "/"), 
       plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
       width = plot.args$width,
       height =plot.args$height,
       useDingbats = plot.args$useDingbats)
