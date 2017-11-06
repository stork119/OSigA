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

rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_filtered.RDS"
poster.data.list <- readRDS(file = rds.path)
poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"
poster.path.list$output.dir <- "resources/output/poster/filtered/"
dir.create(path = poster.path.list$output.dir, 
           showWarnings = FALSE,
           recursive = TRUE)
#### plots by wells ####
gplot.list <- list()
r <-foreach(poster.label = labels(poster.data.list)[c(10)]) %do% {
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

sample.num <- 100
foreach(poster.label = labels(poster.data.list)[c(10)]) %do% {
  #poster.label <- labels(poster.data.list)[1]
  tryCatch({
    print(poster.label)
    data <- poster.data.list[[poster.label]] %>% 
      data.frame() 
    
    col_well <- "well.name"
    
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
  
    
    poster.label <- paste(poster.label, sep = "")
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
      data.time <- data %>% dplyr::filter_(paste(col_time, "==", t))
      cc.list.samples <- 
        foreach(sample.i = 1:sample.num) %do% {
          wells.df <- data.time %>% dplyr::group_by_(col_stimulation) %>% dplyr::distinct(well.name)
          wells.list <- c()
          for(stm in (wells.df %>% dplyr::distinct_(col_stimulation) %>% data.frame())[,col_stimulation]){
            wells <- as.character((wells.df %>% dplyr::filter_(paste(col_stimulation, '==', stm)) %>% data.frame())[, col_well])
            wells.list <- append(x = wells.list, values = sample(x = wells, size = 1))
          }
          
          output_path <- paste(path, "channel_capacity", t, sample.i, "/", sep = "/")
          dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
          write.table(x = data.frame(wells = wells.list), 
                      file = paste(output_path, "wells.csv", sep = "/"),
                      row.names = FALSE,
                      col.names = FALSE,
                      sep = ",")
          cc.output <- capacity_logreg_main(data.time %>% 
                                              dplyr::filter_(
                                                paste(
                                                  col_well, 
                                                  '%in%', 
                                                  'c(', 
                                                  paste("'", wells.list,"'", collapse = ",", sep = ""),
                                                  ')')),
                                            graphs = FALSE,
                                            plot_width = plot.args$width,
                                            plot_height = plot.args$height,
                                            signal = col_stimulation,
                                            response = col_response,
                                          #                                forumla_string =
                                          # cc_maxit = 50,
                                          # lr_maxit = 2000,
                                          output_path = output_path)
          file.remove(paste(output_path, "output.rds", sep = "/"))
          write.table(x = data.frame(cc = cc.output$cc, time = t, sample = sample.i), 
                      row.names = FALSE, 
                      col.names = TRUE, 
                      sep = ",", 
                      file = paste(output_path, "channel_capacity.csv", sep = "/"))
          return(list(time = t, cc = cc.output$cc, sample = sample.i))
      }
      #  cc.list[[as.character(t)]] <- 
      sample.i <- 0
      output_path <- paste(path, "channel_capacity", t, sample.i, "/", sep = "/")
      dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
      cc.output <- capacity_logreg_main(data.time,
                                        graphs = TRUE,
                                        plot_width = plot.args$width,
                                        plot_height = plot.args$height,
                                        signal = col_stimulation,
                                        response = col_response,
                                        output_path = output_path)
      write.table(x = data.frame(cc = cc.output$cc, time = t, sample = sample.i), 
                  row.names = FALSE, 
                  col.names = TRUE, 
                  sep = ",", 
                  file = paste(output_path, "channel_capacity.csv", sep = "/"))
      cc.list.samples[[sample.num + 1]] <- list(time = t, cc = cc.output$cc, sample = sample.i)
      cc.df <- do.call(rbind, cc.list.samples) 
      output_path <- paste(path, "channel_capacity", t, "/", sep = "/")
      write.table(x = cc.df, 
                  row.names = FALSE, 
                  col.names = TRUE, 
                  sep = ",", 
                  file = paste(output_path, "channel_capacity.csv", sep = "/"))
      return(cc.df)
      ### prepare output 
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

plot_fun <- "geom_boxplot"
plot_fun <- "geom_violin"

poster.labels.df <- data.frame(label = sort(labels(poster.data.list)))
poster.labels.df[,col_well] <- 0
poster.labels.df$title[1] = "pSTAT in nuclei; stm: gamma"
poster.labels.df$position[1] = 1
poster.labels.df$title[2] = "pSTAT in nuclei; stm: beta + gamma"
poster.labels.df$position[2] = 3
poster.labels.df$title[3] = "pSTAT in nuclei; stm: gamma + gamma"
poster.labels.df$position[3] = 4
poster.labels.df$title[4] = "pSTAT in nuclei; stm: gamma"
poster.labels.df$position[4] = 1
poster.labels.df$title[5] = "pSTAT in nuclei; stm: beta + gamma"
poster.labels.df$position[5] = 3
poster.labels.df$title[6] = "pSTAT in nuclei; stm: gamma + gamma"
poster.labels.df$position[6] = 4
poster.labels.df$title[7] = "IRF1 in nuclei; stm: gamma"
poster.labels.df$position[7] = 8
poster.labels.df$title[8] = "IRF1 in nuclei; stm: beta + gamma"
poster.labels.df$position[8] = 9
poster.labels.df$title[10] = "IRF1 in nuclei; stm: gamma"
poster.labels.df$position[10] = 8
poster.labels.df$title[11] = "IRF1 in nuclei; stm: beta + gamma"
poster.labels.df$position[11] = 9
poster.labels.df$title[9] = "pSTAT in nuclei; stm: beta"
poster.labels.df$position[9] = 2
poster.labels.df$title[12] = "STAT in cells; stm: beta"
poster.labels.df$position[12] = 6
poster.labels.df$title[13] = "STAT in cytoplasm; stm: beta"
poster.labels.df$position[13] = 7
poster.labels.df$title[14] = "STAT in nuclei; stm: beta"
poster.labels.df$position[14] = 5
poster.labels.df$title[15] = "STAT in cells; stm: beta"
poster.labels.df$position[15] = 6
poster.labels.df$title[16] = "STAT in cytoplasm; stm: beta"
poster.labels.df$position[16] = 7
poster.labels.df$title[17] = "STAT in nuclei; stm: beta"
poster.labels.df$position[17] = 5
poster.labels.df <- poster.labels.df %>% arrange(position, label)

# write.table(x = poster.labels.df,
#             file = 
#               paste(
#                 poster.path.list$output.dir, 
#                 "experiments_description.csv", 
#                 sep = "/"), 
#             sep = ",", 
#             row.names = FALSE, 
#             col.names = TRUE) 
              

gplot.list <- list()
gplot.list.wells <- list()
r <- foreach(poster.label.i = 1:nrow(poster.labels.df)) %do% {
  #poster.label <- labels(poster.data.list)[3]
  poster.label <- as.character(poster.labels.df[poster.label.i,]$label)
  poster.title <- as.character(poster.labels.df[poster.label.i,]$title)
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
  cc.df <- read.table(paste(path,"channel_capacity.csv", sep =  "/"), header = TRUE, sep = ",") %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(cc.mean = mean(cc), cc.sd = sqrt(var(cc)))
  stimulation.list <- unique(data[,col_stimulation])
  cc.df[,col_stimulation] <- floor(length(stimulation.list)/2)
  cc.df[,col_time] <- cc.df$time
  cc.df[,col_well] <- 0
  
  cc.df.0 <-  read.table(paste(path,"channel_capacity.csv", sep =  "/"), header = TRUE, sep = ",") %>%
    dplyr::filter(sample == 0)
  cc.df.0[,col_stimulation] <- floor(length(stimulation.list)/2)
  cc.df.0[,col_time] <- cc.df.0$time
  cc.df.0[,col_well] <- 0
  
  d <- data %>% 
    dplyr::group_by_(
      col_stimulation,
      col_time,
      col_well) %>% 
    dplyr::summarise(n = n()) %>%
    data.frame()
  d[,col_stimulation] <- factor(d[,col_stimulation])
  
  gplot.list.wells[[poster.label]]  <- 
    plot_boxplot_group(
      data = data, 
      x = col_stimulation, 
      y = col_response, 
      boxplot_group = col_well,
      facet_grid_group_y =  col_time,
      save_plot = FALSE,
      ylim_max_const = TRUE,
      plot_title = paste(poster.title, poster.label),
      plot_fun = plot_fun,
      ylim_max = max(quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]], 1500)) +
    geom_text(data = cc.df, 
              mapping = aes(
                label = paste(round(cc.mean,2), "+",  round(cc.sd,3), sep = ""),
                y = 1400),
              size = 10,
              color = "red") +
    geom_text(data = cc.df.0,
              mapping = aes(
                label = paste(round(cc,2), sep = ""),
                y = 1300),
              size = 10,
              color = "blue") +
    geom_text(data =  d,
              mapping = aes_string(
                label = "n",
                y = 1000,
                x = col_stimulation,
                group = col_well),
              inherit.aes = TRUE,
              position = position_dodge(width = 1),
              angle = 90,
              size = 5,
              color = "blue")
  

  gplot.list[[poster.label]]  <- 
    plot_boxplot_group(
      data = data, 
      x = col_stimulation, 
      y = col_response, 
      facet_grid_group_y =  col_time,
      save_plot = FALSE,
      ylim_max_const = TRUE,
      plot_title = paste(poster.title, poster.label),
      plot_fun = plot_fun,
      ylim_max = max(quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]], 1500)) +
    geom_text(data = cc.df, 
              mapping = aes(
                label = paste(round(cc.mean,2), "+",  round(cc.sd,3), sep = ""),
                y = 1200),
              size = 10,
              color = "red") +
    geom_text(data = cc.df.0,
              mapping = aes(
                label = paste(round(cc,2), sep = ""),
                y = 1000),
              size = 10,
              color = "blue")
  
  return()
}
ggsave(filename = paste(poster.path.list$output.dir, paste("channel_capacity_wells_", plot_fun, ".pdf", sep = ""), sep = "/"), 
       plot = marrangeGrob(grobs = gplot.list.wells, ncol = 1, nrow =1),
       width = plot.args$width,
       height =plot.args$height,
       useDingbats = plot.args$useDingbats)

ggsave(filename = paste(poster.path.list$output.dir, paste("channel_capacity_", plot_fun, ".pdf", sep = ""), sep = "/"), 
       plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
       width = plot.args$width,
       height =plot.args$height,
       useDingbats = plot.args$useDingbats)
