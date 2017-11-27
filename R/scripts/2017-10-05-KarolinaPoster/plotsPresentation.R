###
# data plots 
###


source("R/optimisation/initialise_optimisation.R")
#install.packages("e1071")
library(e1071)
library(CapacityLogReg)

rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_joined.RDS"
poster.data.list <- readRDS(file = rds.path)
poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"
poster.path.list$output.dir <- "resources/output/poster/joined/"



#foreach(poster.label = labels(poster.data.list)[c(10)]) %do% {
 plot.fill <- data.frame(fill = c("#ffffff", "#a4d09e", "#a4d09e", "#7abd6f", "#7abd6f", "#50af43", "#50af43", "#1fa637", "#177d34", "#175929"),
                        stm  =  c(0, 1.8, 3.6, 12, 9,   18, 36, 90,  180, 900),
                        stmnew  =  c(0, 2, 4,  10, 10,   20, 40,  100,  200, 1000))  
  
  #### Beta-Stat-Nuclei ####
  poster.label <- labels(poster.data.list)[1]
  ylimmax <- 1000 #1000
  stimulation.list <- list(c(10,1000))  
  plot.fill <- data.frame(fill = c("#ffffff", "#a4d09e", "#a4d09e", "#7abd6f", "#7abd6f", "#50af43", "#50af43", "#1fa637", "#177d34", "#175929"),
                          stm  =  c(0, 1.8, 3.6, 12, 9,   18, 36, 90,  180, 900),
                          stmnew  =  c(0, 2, 4,  10, 10,   20, 40,  100,  200, 1000))  
  
  #### Beta-Stat-Cytoplasm ####
  poster.label <- labels(poster.data.list)[2]
  ylimmax <- 600 #1000
  stimulation.list <- list(c(1000))  
  plot.fill <- data.frame(fill = c("#ffffff", "#a4d09e", "#a4d09e", "#7abd6f", "#7abd6f", "#50af43", "#50af43", "#1fa637", "#177d34", "#175929"),
                          stm  =  c(0, 1.8, 3.6, 12, 9,   18, 36, 90,  180, 900),
                          stmnew  =  c(0, 2, 4,  10, 10,   20, 40,  100,  200, 1000))  
  
  #### Beta-pStat-Nuclei ####  
  poster.label <- labels(poster.data.list)[4]
  ylimmax <- 300 #1000
  stimulation.list <- NULL #list(c(10,1000))  
  plot.fill <- data.frame(fill = c("#ffffff", "#a4d09e", "#a4d09e", "#7abd6f", "#7abd6f", "#50af43", "#50af43", "#1fa637", "#177d34", "#175929"),
                          stm  =  c(0, 1.8, 3.6, 12, 9,   18, 36, 90,  180, 900),
                          stmnew  =  c(0, 2, 4,  10, 10,   20, 40,  100,  200, 1000))  
  
  #### Gamma-pStat-Nuclei ####  
  poster.label <- labels(poster.data.list)[5]
  ylimmax <- 300 #1000
  stimulation.list <- NULL# list(c(0.1,1))  
  plot.fill <- data.frame(fill = c("#ffffff", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6", "#325aa6"),
                          stm  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 10),
                          stmnew  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 5))  
  
  #### Gamma-Gamma-pStat-Nuclei ####  
  poster.label <- labels(poster.data.list)[6]
  ylimmax <- 1000 #1000
  stimulation.list <- NULL #list(c(10,1000))  
  plot.fill <- data.frame(fill = c("#ffffff", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6", "#325aa6"),
                          stm  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 10),
                          stmnew  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 5))  
  
  #### Beta-Gamma-pStat-Nuclei ####  
  poster.label <- labels(poster.data.list)[7]
  ylimmax <- 1000 #1000
  stimulation.list <- NULL #list(c(10,1000))  
  plot.fill <- data.frame(fill = c("#ffffff", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6", "#325aa6"),
                          stm  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 10),
                          stmnew  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 5))  
  
    #### Gamma-IRF-Nuclei ####
  poster.label <- labels(poster.data.list)[8]
  ylimmax <- 600 #1000
  plot.fill <- data.frame(fill = c("#ffffff", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6"),
                          stm  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 10),
                          stmnew  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5))  
  stimulation.list <- NULL 
  
  #### Beta-Gamma-IRF-Nuclei ####
  poster.label <- labels(poster.data.list)[9]
  ylimmax <- 600 #1000
  plot.fill <- data.frame(fill = c("#ffffff", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6"),
                          stm  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 10),
                          stmnew  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5))  
  stimulation.list <- NULL 
  
  
  #### ####
  plot_title <- poster.labels.df[poster.labels.df$label  == poster.label,]$title
  tryCatch({
    print(poster.label)
    data <- poster.data.list[[poster.label]] %>% 
      data.frame() 
    data <- data %>% dplyr::filter(stimulation != 0)
    # data <- data %>% dplyr::filter_(paste(col_time, "!=", 42))  %>% dplyr::filter_(paste(col_time, "!=", 36))
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
    }
    
    if("Intensity_MeanIntensity_Alexa" %in% colnames(data)){
      col_response <- "Intensity_MeanIntensity_Alexa"
    } else if("Intensity_MeanIntensity_Alexa488" %in% colnames(data)){
      col_response <- "Intensity_MeanIntensity_Alexa488"
    } else {
      col_response <- "Intensity_MeanIntensity_Alexa555"
    }
    
    
    plot.fill[,col_stimulation] <- plot.fill[,"stm"] 
    
    data.0 <-  data %>% 
      dplyr::filter_(paste(col_time, "==", 0)) 
    if(nrow(data.0) > 0){
      data.0[,col_stimulation] <- 0
    }
    data <- data %>% rbind(data.0)
    gplot.list <- list()
    colors.list.list <- list()
    for( t in (data %>%  dplyr::distinct_(col_time) %>% 
               data.frame())[,col_time])
    {
      data.time <- data %>% 
        dplyr::filter_(paste(col_time, "==", t)) %>%
        rbind(data.0) %>%
        dplyr::arrange_(col_stimulation) %>%
        dplyr::left_join(plot.fill)
      data.time[,col_stimulation] <- data.time[,"stmnew"]
      #dplyr::mutate(stimulation.1.1 = stmnew)
      #dplyr::mutate_(paste(col_stimulation, "=", "stmnew"))
      colors.list <- as.character(
        (plot.fill %>% 
           dplyr::filter(stmnew %in%
                           as.character(unique(data.time[, col_stimulation]))))[,"fill"])
      data.time[, col_stimulation] <- factor(data.time[, col_stimulation])
      g <- ggplot(data = data.time, 
                  mapping = aes_string(
                    x = col_stimulation,
                    y = col_response,
                    group = col_stimulation,
                    fill = col_stimulation)) +
        geom_boxplot()  + 
        ylim(c(0,ylimmax)) + 
        theme_jetka() +
        ggtitle(
          paste("time:", t, ";", plot_title, ";", poster.label)) + 
        xlab("stimulation [U/ml]") +
        ylab("Fluorescent Intensity (A.U.)")
      if(length(colors.list) > 1){
        g <- g + scale_fill_manual(values = colors.list)
      }
      gplot.list[[as.character(t)]] <- g
      colors.list.list[[as.character(t)]] <- colors.list
    }
    dir.create(path =  paste(poster.path.list$output.dir, 
                             poster.label,
                             sep = "/"),
               showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = paste(poster.path.list$output.dir, 
                            poster.label,
                            paste(poster.label,
                                  "_Y_fluorescent_X_stimulation.pdf",
                                  sep = ""),
                            sep = "/"), 
           plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
           width = 8,
           height = 8,
           useDingbats = plot.args$useDingbats)
    
    data <- data %>% dplyr::left_join(plot.fill)
    data[,col_stimulation] <- data[,"stmnew"]
    if(is.null(stimulation.list)){
      stimulation.list <- sort(unique(data[,col_stimulation]))
    }
    gplot.list <- list()
    for(stimulation.i in 1:length(stimulation.list)){
      stimulations <- stimulation.list[[stimulation.i]]
      colors.list <- 
        unique(
          as.character(
        (plot.fill %>% 
           dplyr::filter(stmnew %in% stimulations))[,"fill"]))
      data[, col_time] <- factor(data[,col_time])
      g <- ggplot(data = 
                    data %>% 
                    dplyr::filter_(paste(col_stimulation, 
                                         "%in%", 
                                         "c(", 
                                         paste(stimulations, collapse = ","),
                                         ")")),
                  mapping = 
                    aes_string(
                      x = col_time,
                      y = col_response,
                      group = paste("interaction(", col_time, ",", col_stimulation,")"),
                      fill = paste("factor(", col_stimulation, ")"))) +
        #geom_boxplot()  + 
        ylim(c(0,ylimmax)) + 
        theme_jetka() +
        ggtitle(
          paste("stimulation:",
                paste(stimulations, collapse = ","),
                ";", 
                plot_title, ";",
                poster.label)) + 
        xlab("time") +
        ylab("Fluorescent Intensity (A.U.)")
      if(length(colors.list) == 1){
        g <- g +  geom_boxplot(fill = colors.list) 
      }  else {
        g <- g + geom_boxplot() + scale_fill_manual(values = colors.list)
     }
     gplot.list[[stimulation.i]] <- g
    }
    ggsave(filename = paste(poster.path.list$output.dir, 
                            poster.label,
                            paste(poster.label, 
                                  "_Y_fluorescent_X_time_", 
                                  ".pdf", sep = ""), sep = "/"), 
           plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
           width = 8,
           height = 8,
           useDingbats = plot.args$useDingbats)
  })

    
    
    
  
  
  
  
  
  #### ###