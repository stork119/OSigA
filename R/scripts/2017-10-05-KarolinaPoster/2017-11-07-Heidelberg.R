### ###
### Heidelberg script
### ###


rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_joined.RDS"
poster.data.list <- readRDS(file = rds.path)
poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"
poster.path.list$output.dir <- "resources/output/poster/joined/"

#### Gamma-pStat-Nuclei ####  
poster.label.list <- NULL
poster.label <- labels(poster.data.list)[5]
ylimmax <- 1000 #1000
stimulation.list <- NULL #list(c(10,1000))  
plot.fill <- data.frame(fill = c("#e6e9ed", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6"),
                        stm  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5),
                        stmnew  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5))  

#### Beta-Stat-Cytoplasm ####
poster.label.list <- NULL
poster.label <- labels(poster.data.list)[2]
ylimmax <- 600 #1000
stimulation.list <- list(c(1000))  
plot.fill <- data.frame(fill = c("#ffffff", "#a4d09e", "#a4d09e", "#7abd6f", "#7abd6f", "#50af43", "#50af43", "#1fa637", "#177d34", "#175929"),
                        stm  =  c(0, 1.8, 3.6, 12, 9,   18, 36, 90,  180, 900),
                        stmnew  =  c(0, 2, 4,  10, 10,   20, 40,  100,  200, 1000))  

#### Gamma-pStat-Nuclei & Beta-Gamma-pStat-Nuclei ####  
poster.label.list <- c(labels(poster.data.list)[5], labels(poster.data.list)[7])
ylimmax <- 600 #1000
stimulation.list <- NULL #list(c(10,1000))  
plot.fill <- data.frame(fill = c("#e6e9ed", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6"),
                        stm  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5),
                        stmnew  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5))  



#### histogram ####

if(is.null(poster.label.list)){
  poster.label.list <- c(poster.label)
}

data <- do.call(rbind, poster.data.list[poster.label.list])
poster.label <- paste(poster.label.list, collapse = "_")
# data <- data %>% dplyr::filter(stimulation != 0)
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
data <- data %>% 
  dplyr::left_join(plot.fill)
data[,col_stimulation] <- data[,"stmnew"]
# t <- 30
# data.time <- data %>% 
#   dplyr::filter_(paste(col_time, "==", t))
# 
# data.time <- data %>% 
#   dplyr::filter_(paste(col_time, "%in%", "c(",0,",",24,")")) %>%
#   dplyr::filter_(paste(col_stimulation, "%in%", "c(",0,",",1000,")"))
# 

data <- data %>% dplyr::filter_(paste(col_stimulation, "==", 1))

data.time <- data.time %>% dplyr::filter_(paste(col_stimulation, "!=", 0.01))

colors.list <- as.character(
  (plot.fill %>% 
     dplyr::filter(stmnew %in%
                     as.character(unique(data.time[, col_stimulation]))))[,"fill"])

g <- ggplot(data = data.time, 
            mapping = aes_string(
              #x = paste("log(", col_response, ")"),
              x = paste("log10(", col_response, ")"),
              group = paste("factor(", col_time, ")"),
              color = paste("factor(", col_stimulation, ")"))) +
  geom_density(size = 1.5)  + 
  theme_jetka() +
  ggtitle(
    paste("time:", t, ";", plot_title, ";", poster.label)) + 
  #xlab("stimulation [U/ml]") +
  #xlim(c(3,6.5)) +
  ylab("density") +
  xlab("log(fluorescent intensity [a.u.])") +
  scale_color_manual(values = colors.list)
g
ggsave(filename = paste(poster.path.list$output.dir, 
                        poster.label,
                        paste(poster.label,
                              "_density_plot.pdf",
                              sep = ""),
                        sep = "/"), 
       plot = g,
       width = 12,
       height = 6,
       useDingbats = plot.args$useDingbats)

#### ####
poster.label.nopriming <- labels(poster.data.list)[5]
poster.label.priming   <- labels(poster.data.list)[7]
data.stm <- rbind(poster.data.list[[poster.label.nopriming]],  
                  poster.data.list[[poster.label.priming]])
data.stm <- data.stm %>% dplyr::filter_(paste(col_stimulation, "!=", 0.01))
data.stm <- data.stm %>% dplyr::mutate(lwd = ifelse(priming == 1000, 0.5, 0.45))

colors.list <- as.character(
  (plot.fill %>% 
     dplyr::filter(stmnew %in%
                     as.character(unique(data.time[, col_stimulation]))))[,"fill"])

g <- ggplot(
  data = data.stm,
  mapping = aes_string(
    x = paste("factor(", col_stimulation, ")"),
    y = col_response,
    group = paste("interaction(", "priming",",",col_stimulation, ")"),
    #size  = "lwd", 
    fill  = paste("factor(", col_stimulation, ")"))) +
  geom_boxplot() +
  ylim(c(0,600)) +
  theme_jetka() +
  ggtitle(
    paste("time:", t, ";", plot_title, ";", poster.label)) + 
  xlab("stimulation [U/ml]") +
  ylab("Fluorescent Intensity (A.U.)") +
  scale_fill_manual(values = colors.list) 
g
ggsave(filename = paste(poster.path.list$output.dir, 
                        poster.label,
                        paste(poster.label,
                              "_sensing_boxplots_comparison.pdf",
                              sep = ""),
                        sep = "/"), 
       plot = g,
       width = 8,
       height = 6,
       useDingbats = plot.args$useDingbats)


#### ####
colors.list <- as.character(
  (plot.fill %>% 
     dplyr::filter(stmnew %in%
                     as.character(unique(data[, col_stimulation]))))[,"fill"])

g <- ggplot(
  data = data,
  mapping = aes_string(
    x = paste("factor(", col_time, ")"),
    y = col_response,
    group = paste("interaction(", col_time, ",", "priming",",",col_stimulation, ")"),
    #size  = "lwd", 
    fill  = paste("factor(", col_stimulation, ")"))) +
  geom_boxplot() +
  ylim(c(0,600)) +
  theme_jetka() +
  ggtitle(
    paste("time:", t, ";", plot_title, ";", poster.label)) + 
  xlab("time [min]") +
  ylab("Fluorescent Intensity (A.U.)") +
  scale_fill_manual(values = colors.list) 
g
dir.create(paste(poster.path.list$output.dir, 
                 poster.label, 
                 sep = "/"))
ggsave(filename = paste(poster.path.list$output.dir, 
                        poster.label,
                        paste(poster.label,
                              "_sensing_boxplots_comparison.pdf",
                              sep = ""),
                        sep = "/"), 
       plot = g,
       width = 8,
       height = 6,
       useDingbats = plot.args$useDingbats)
