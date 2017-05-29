### ### ###
### 2017_05_27
### script_analyse_data
### ### ###
source("R/libraries.R")
source("R/libraries/library_data_manipulation.R")
source("R/graphics/boxplot.R")

#### initialise ###
path.data.output <- "resources/output/data/"
path.data.input  <- "resources/input/data"

#### ####
dir.create(path = path.data.output, 
           recursive = TRUE, 
           showWarnings = TRUE)
plot.args <- list(theme.title_size = 9,
                  theme.margins =  rep(x = 0, times = 4))

#### ####
data.list <- read_data(path = path.data.input)

data.list$data.exp %>% dplyr::distinct(file)
g.list <- list()
#### boxplots all ####
g.list[["boxplots"]] <- list()
for(prm in (data.list$data.exp %>% dplyr::distinct_("priming"))$priming){
  g.list[["boxplots"]][[as.character(prm)]] <- 
    plot_boxplot_group(data = 
                     data.list$data.exp %>% 
                     dplyr::filter_("priming" != prm),
                   save_plot = FALSE,
                   x = "time",
                   y = "intensity",
                   boxplot_group = "position",
                   facet_grid_group_y = "stimulation",
                   ylim_max = 1500,
                   ylab = "pSTAT1 in nucleus [u]",
                   plot_title = paste("priming", prm))
}



#### histograms ####
data.groups <- 
  expand.grid(
    priming = (data.list$data.exp %>% 
                 dplyr::distinct_("priming"))$priming,
    stimulation = (data.list$data.exp %>% 
                     dplyr::distinct_("stimulation"))$stimulation,
    time = (data.list$data.exp %>% 
                     dplyr::distinct_("time"))$time) %>% 
  dplyr::arrange(priming, stimulation, time)

g.list[["histograms"]] <- list()
g.list[["normal"]] <- list()
g.list[["log"]] <- list()


for(i in 1:nrow(data.groups)){
  data <- data.list$data.exp %>% 
    dplyr::filter(priming == data.groups[i,"priming"],
                  stimulation == data.groups[i,"stimulation"],
                  time == data.groups[i,"time"]
                   )
  
  title <- paste("prm", data.groups[i,]$priming, 
                 "stm", data.groups[i,]$stimulation, 
                 "time", data.groups[i,]$time,
                 sep = " ") 
  
  x <- "intensity"
  group <- "position"
  xlim_max <- quantile(data$intensity, probs = 0.99)
  
  g.list[["normal"]][[as.character(i)]] <- 
    do.call(what = plot_density, 
            args = append(plot.args,
                          list(data = data,
                               x = "intensity", 
                               group = group, 
                               title = title, 
                               xlim_max = xlim_max, 
                               compare_to_all = TRUE)
                          ))
  
  g.list[["log"]][[as.character(i)]] <- 
    do.call(what = plot_density, 
            args = append(plot.args,
                          list(
                            data = data,
                            x = "log(intensity)", 
                            group = group, 
                            title = title, 
                            xlim_max = NULL, 
                            compare_to_all = TRUE
                          )))
      
}

#### save plots ####
ggsave(plot = marrangeGrob(grobs = g.list[["log"]],
                             nrow = 2,
                             ncol = 3,
                             top = NULL),
       filename = paste(path.data.output, 
                      "density-log.pdf", 
                      sep = "/" ),
         width = 24, 
         height = 12, 
         useDingbats = FALSE)

ggsave(plot = marrangeGrob(grobs = g.list[["normal"]],
                           nrow = 2,
                           ncol = 3,
                           top = NULL),
       filename = paste(path.data.output, 
                        "density-normal.pdf", 
                        sep = "/" ),
       width = 24, 
       height = 12, 
       useDingbats = FALSE)

ggsave(plot = marrangeGrob(grobs = g.list[["boxplots"]],
                           nrow = 1,
                           ncol = 1,
                           top = NULL),
       filename = paste(path.data.output, 
                        "boxplots.pdf", 
                        sep = "/" ),
       width = 24, 
       height = 12, 
       useDingbats = FALSE)
#### #### 

data.list$data.exp.norm <- get_equal_data(data = data.list$data.exp,
                                sample_size = 1000)

data.list$data.exp.summarise <- 
  data.list$data.exp.norm %>% 
  dplyr::group_by(priming,
                  stimulation,
                  time) %>%
  summarise(mean = mean(intensity),
            sd   = var(intensity))

data.list$data.exp.summarise$likelihood <- 
  (data.list$data.exp.norm %>%
     left_join(data.list$data.exp.summarise,
               by = c("priming",
                      "stimulation",
                      "time")) %>%
     dplyr::mutate(likelihood =
                     do.call(fun.likelihood.lmvn.data,
                             args = list(m = mean,
                                         sd = sd,
                                         model.m = intensity
                             ))) %>%
     dplyr::group_by(priming,
                     stimulation,
                     time) %>%
     summarise(likelihood = sum(likelihood)))$likelihood
    

g.list[["data_likelihood"]] <- 
  ggplot(data = data.list$data.exp.summarise, 
       mapping = aes(x = time, y = likelihood)) +
    geom_point() +
    facet_grid(priming ~ stimulation) +
    do.call(theme_jetka, args = plot.args)
ggsave(plot = g.list[["data_likelihood"]],
       filename = paste(path.data.output, 
                        "data_likelihood.pdf", 
                        sep = "/" ),
       width = 24, 
       height = 12, 
       useDingbats = FALSE)