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

data.list$data.exp.norm <- get_equal_data(data = data.list$data.exp,
                                          sample_size = 1000)

data.list$data.exp.summarise <- 
  data.list$data.exp.norm %>% 
  dplyr::group_by(priming,
                  stimulation,
                  time) %>%
  summarise(m.norm = mean(intensity),
            sd.norm   = var(intensity))

data.list$data.exp.summarise<-
  data.list$data.exp.summarise %>%
  dplyr::mutate(mean.lmvn = lmvn.mean(m = m.norm, sd = sd.norm),
                sd.lmvn = lmvn.sd(m = m.norm, sd = sd.norm))

data.list$data.exp.summarise$likelihood <- 
  (data.list$data.exp.norm %>%
     left_join(data.list$data.exp.summarise,
               by = c("priming",
                      "stimulation",
                      "time")) %>%
     dplyr::mutate(likelihood =
                     do.call(ComputeLikelihood.lmvn,
                             args = list(m = m.norm,
                                         sd = sd.norm,
                                          X  = intensity
                             ))) %>%
     dplyr::group_by(priming,
                     stimulation,
                     time) %>%
     summarise(likelihood = sum(likelihood)))$likelihood

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
#### likelihood plots #### 

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
#### mean variance ####
g.list[["mean_variance"]] <- PlotMeanVariance(y = "m.norm", yerror = "sd.norm")
ggsave(plot = marrangeGrob(grobs = g.list[["mean_variance"]], nrow = 1, ncol = 1),
       filename = paste(path.data.output, 
                        "mean_variance.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)

g.list[["lmvn_mean_variance"]] <- PlotMeanVariance(y = "mean.lmvn", yerror = "sd.lmvn")
ggsave(plot = marrangeGrob(grobs = g.list[["lmvn_mean_variance"]], nrow = 1, ncol = 1),
       filename = paste(path.data.output, 
                        "lmvn_mean_variance.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)

#### LMVN ####
y <- "sd.lmvn"
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
# ggplot(data = data.list$data.exp.summarise %>% 
#          dplyr::mutate(likelihood = mean.lmvn^2/mean.lmvn), 
#        mapping = aes(x = time, y = likelihood)) 

g.list[["data_likelihood_mean.lmvn_sd.lmvn"]] <- do.call(what = plot_points,
        args = append(plot.args,
                      list(
                        data.list$data.exp.summarise %>% 
                          dplyr::mutate(likelihood = mean.lmvn^2/mean.lmvn),
                        x = "time",
                        y = "likelihood",
                        facet_grid_group_y = "stimulation",
                        facet_grid_group_x = "priming",
                        title = "mean.lmvn^2/sd.lmvn"
                      )))
ggsave(plot = g.list[["data_likelihood_mean.lmvn_sd.lmvn"]],
       filename = paste(path.data.output, 
                        "data_likelihood_mean.lmvn_sd.lmvn.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)

#### normalisation testing ####
attach(LoadOptimisationConditions(path.optimisation,
                                       path.optimisation.data,
                                       maxit.tmp = Inf))

data.exp.summarise.optimisation$likelihood <-
  likelihood(data.model = data.exp.summarise.optimisation,
          data.exp.grouped = data.exp.grouped.optimisation,
          data.exp.summarise = data.exp.summarise.optimisation,
           fun.likelihood = fun.likelihood.list.sd_data
)


data.model.list[[1]]$likelihood_2 <- likelihood(data.model = data.model.list[[1]],
           data.exp.grouped = data.exp.grouped.optimisation,
           data.exp.summarise = data.exp.summarise.optimisation,
           fun.likelihood = fun.likelihood.list.sd_data
)




sum(data.exp.summarise.optimisation$likelihood)
sum(data.model.list[[1]]$likelihood)
gplot.list.likelihood <- list()
gplot.list.likelihood[["likelihood_model_vs_data"]] <- ggplot(data = rbind(data.frame(data.exp.summarise.optimisation %>%
         dplyr::mutate(type = "data")),
         data.model.list[[1]][,c(colnames(data.exp.summarise.optimisation), "type")]),
       mapping = aes(x = factor(time), y = likelihood, colour = type)) +
  geom_point() +
  #geom_line() +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  #geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  #geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data"),
  #              mapping = aes(ymin = mean.lmvn - sqrt(sd.lmvn),
  #                            ymax = mean.lmvn + sqrt(sd.lmvn)),
  #              color = "black") +
  ggtitle(paste("Compare best model with data means", collapse = " "))

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "likelihood_lmvn_compare_data_with_best_model.pdf", sep = ""),
                           plot = gplot.list.likelihood[["likelihood_model_vs_data"]])))


# data.model.list[[1]]$likelihood_2 <-
#   likelihood(data.model = data.model.list[[1]],
#              data.exp.grouped = data.exp.grouped.optimisation,
#              data.exp.summarise = data.exp.summarise.optimisation,
#              fun.likelihood = fun.likelihood.list.sd_data
#   )
# 
# sum(data.model.list[[1]]$likelihood)
# # data.exp.grouped.optimisation %>%
# #   dplyr::group_by(priming, stimulation, time) %>%
# #   dplyr::filter(priming == 1000, stimulation == 1, time == 30) %>%
# #   dplyr::summarise(m.norm = mean(intensity),
# #                    mean.lmvn = mean(logintensity),
# #                    sd.norm = var(intensity),
# #                    sd.lmvn = var(logintensity))
# 
# model.tmp <- data.exp.grouped.optimisation %>%
#   dplyr::filter(priming == 1000, stimulation == 1, time == 30) %>%
#   mutate(data = (logintensity - 6.09351)^2,
#          model = ((logintensity - 5.986634)^2)/0.1249113) %>%
#   summarise(sum_data = sum(data)/0.1249113, sum_model = sum(model))
# 
# 
# fun.likelihood.list.sd_data(logintensity = (data.list$data.exp.norm %>%
#                                               dplyr::filter(priming == 1000, stimulation == 1, time == 30))$logintensity, data.model.tmp = tmp.model, data.exp.summarise = data.exp.summarise.optimisation)
# 
# tmp.data <- data.exp.summarise.optimisation %>% dplyr::filter(priming == 1000, stimulation == 1, time == 30)
# tmp.model <- data.model.list[[1]] %>% dplyr::filter(priming == 1000, stimulation == 1, time == 30)
# mean.lmvn(m.norm = tmp.model$m.norm, sd.norm = tmp.data$sd.norm)
# sd.lmvn(m.norm = tmp.model$m.norm, sd.norm = tmp.data$sd.norm)

gplot.list.likelihood[["density"]] <- list()
for(i in 1:nrow(data.exp.summarise.optimisation)){
  prm <- data.exp.summarise.optimisation[i,]$priming
  t <- data.exp.summarise.optimisation[i,]$time
  stm <- data.exp.summarise.optimisation[i,]$stimulation
  data <- data.exp.grouped.optimisation %>%
    data.table() %>%
    dplyr::filter(priming == prm &
                  stimulation == stm &
                  time == t
    )
  
  
  title <- paste("prm", data.exp.summarise.optimisation[i,]$priming,
                 "stm", data.exp.summarise.optimisation[i,]$stimulation,
                 "time", data.exp.summarise.optimisation[i,]$time,
                 sep = " ")

  x <- "intensity"
  group <- "position"
#  xlim_max <- quantile(data$intensity, probs = 0.99)

  data.rand <- data.frame(logintensity =
                            rnorm(n = 500,
                                  mean = data.exp.summarise.optimisation[i,]$mean.lmvn,
                                  sd = data.exp.summarise.optimisation[i,]$sd.lmvn))
  data.rand <- cbind(data.rand, data.frame(data.exp.summarise.optimisation[i,]))
  data.rand <- data.rand %>% dplyr::mutate(intensity = exp(logintensity), position = "rand") %>% data.table()
  data.rand$cell <- 1:nrow(data.rand)
  data.rand$file <- ""
  data <- rbind(data.frame(data), data.frame(data.rand)[,colnames(data)]) %>% data.table()
  # g.list[["normal"]][[as.character(i)]] <-
  #   do.call(what = plot_density,
  #           args = append(plot.args,
  #                         list(data = data,
  #                              x = "intensity",
  #                              group = group,
  #                              title = title,
  #                              xlim_max = xlim_max,
  #                              compare_to_all = TRUE)
  #           ))

  gplot.list.likelihood[["density"]][[as.character(i)]] <-
    do.call(what = plot_density,
            args = append(plot.args,
                          list(
                            data = data,
                            x = "logintensity",
                            group = group,
                            title = title,
                            xlim_max = NULL,
                            compare_to_all = FALSE
                          )))

}
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "likelihood_lmvn_compare_data_with_rand.pdf", sep = ""),
                           plot = marrangeGrob(grobs = gplot.list.likelihood[["density"]], ncol = 1, nrow = 1))))

do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.data.output, "likelihood_lmvn_compare_data_with_rand_twopages.pdf", sep = ""),
                           plot = marrangeGrob(grobs = gplot.list.likelihood[["density"]], ncol = 4, nrow = 2))))
