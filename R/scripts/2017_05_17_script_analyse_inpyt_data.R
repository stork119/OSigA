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
path.data.output <- "resources/output/data/"
path.data.input  <- "resources/input/data"

#### ####
dir.create(path = path.data.output, 
           recursive = TRUE, 
           showWarnings = TRUE)
plot.args <- list(theme.title_size = 9,
                  theme.margins =  rep(x = 0, times = 4),
                  width = 24, 
                  height = 12, 
                  useDingbats = FALSE)

#### ####
data.list <- read_data(path = path.data.input)

data.list$data.exp %>% dplyr::distinct(file)
g.list <- list()
#### boxplots all ####
g.list[["boxplots"]] <- list()
for(prm in (data.list$data.exp %>% dplyr::distinct_("priming"))$priming){
  g.list[["boxplots"]][[as.character(prm)]] <- 
    do.call(what = plot_boxplot_group,
            args = append(plot.args,
                          list(
                            data = 
                              data.list$data.exp %>% 
                              dplyr::filter_("priming" != prm),
                            save_plot = FALSE,
                            x = "time",
                            y = "intensity",
                            boxplot_group = "position",
                            facet_grid_group_y = "stimulation",
                            ylim_max = 1500,
                            ylab = "pSTAT1 in nucleus [u]",
                            plot_title = paste("priming", prm)
                          )))
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
data.list$data.exp.summarise<-
  data.list$data.exp.summarise %>%
  dplyr::mutate(lmvn.mean = lmvn.mean(m = mean, sd = sd),
                lmvn.sd = lmvn.sd(m = mean, sd = sd))

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

#### plot mean and variance ####

PlotMeanVariance <- function(
  y,
  yerror){
  g.list.tmp <- list()
  for(prm in 
      (data.list$data.exp.summarise %>%
       dplyr::ungroup() %>%
       dplyr::distinct_("priming"))$priming){
    
    
    g.list.tmp[[as.character(prm)]] <-
      do.call(what = plot_points,
              args = append(plot.args,
                            list(
                              data = data.list$data.exp.summarise %>% 
                                dplyr::filter(priming == prm) %>%
                                dplyr::mutate_(var_mutate = yerror) %>%
                                dplyr::mutate(sqrt_yerror = sqrt(var_mutate)),
                              x = "time",
                              y = y,
                              yerror = "sqrt_yerror",
                              facet_grid_group_y = "stimulation",
                              title = prm
                            )))
    
    
  }
  return(g.list.tmp)
}

g.list[["mean_variance"]] <- PlotMeanVariance(y = "mean", yerror = "sd")
ggsave(plot = marrangeGrob(grobs = g.list[["mean_variance"]], nrow = 1, ncol = 1),
       filename = paste(path.data.output, 
                        "mean_variance.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)

g.list[["lmvn_mean_variance"]] <- PlotMeanVariance(y = "lmvn.mean", yerror = "lmvn.sd")
ggsave(plot = marrangeGrob(grobs = g.list[["lmvn_mean_variance"]], nrow = 1, ncol = 1),
       filename = paste(path.data.output, 
                        "lmvn_mean_variance.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)

#### LMVN ####
y <- "lmvn.sd"
g.list[["lmvn_variance"]] <- do.call(what = plot_points,
          args = append(plot.args,
                        list(
                          data = data.list$data.exp.summarise %>%
                            dplyr::mutate_(var_mutate = y) %>%
                            dplyr::mutate(sqrt_y = sqrt(var_mutate)),
                          x = "time",
                          y = "sqrt_y",
                          facet_grid_group_y = "stimulation",
                          facet_grid_group_x = "priming",
                          title = "variance"
                        )))
ggsave(plot = g.list[["lmvn_variance"]],
       filename = paste(path.data.output, 
                        "lmvn_variance.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)


#### ####
ggplot(data = data.list$data.exp.summarise %>% 
         dplyr::mutate(likelihood = lmvn.mean^2/lmvn.mean), 
       mapping = aes(x = time, y = likelihood)) 

g.list[["data_likelihood_lmvn.mean_lmvn.sd"]] <- do.call(what = plot_points,
        args = append(plot.args,
                      list(
                        data.list$data.exp.summarise %>% 
                          dplyr::mutate(likelihood = lmvn.mean^2/lmvn.mean),
                        x = "time",
                        y = "likelihood",
                        facet_grid_group_y = "stimulation",
                        facet_grid_group_x = "priming",
                        title = "lmvn.mean^2/lmvn.sd"
                      )))
ggsave(plot = g.list[["data_likelihood_lmvn.mean_lmvn.sd"]],
       filename = paste(path.data.output, 
                        "data_likelihood_lmvn.mean_lmvn.sd.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)
