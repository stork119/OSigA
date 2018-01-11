### ###
### IRFplot_library
### ###

#### plotSimulationsFun ####
plotSimulationsFun <- function(data.exp,
                               data.sample,
                               stimulations,
                               plot.title = "",
                               plot.grob.title = "",
                               plot.grob.nrow = 1,
                               plot.grob.ncol = 1,
                               ...){
  
  g.list <- foreach( stm =  stimulations ) %do% {
    
    ggplot(data = rbind(data.exp,
                        data.sample) %>%
             dplyr::filter_(paste("stimulation ==", stm)), 
           mapping = aes_string(x = "log(response)", 
                                group = "type",
                                fill = "type",
                                color = "type")) +
      do.call(theme_jetka, args = plot.args) +
      geom_density(alpha = 0.5) +
      ggtitle(paste(plot.title, stm))
  }
  
  g.plot <- marrangeGrob(grobs = g.list,
                         nrow = plot.grob.nrow, 
                         ncol = plot.grob.ncol,
                         top  = plot.grob.title)
  return(g.plot)
}

#### plot errorbars ####
plotSimulationsErrorbarsFun <- function(data.exp,
                                        data.sample,
                                        plot.title = "",
                                        ...){
  
  
  g.plot <- ggplot(data = 
                     rbind(data.exp, data.sample) %>% 
                     dplyr::group_by_("type", "stimulation") %>%
                     dplyr::summarise_(logresponse = "mean(log(response))",
                                       logresponse.sd = "sd(log(response))"), 
                   mapping = aes_string(x = "factor(stimulation)", 
                                        y = "logresponse", 
                                        ymin = "logresponse - logresponse.sd",
                                        ymax = "logresponse + logresponse.sd",
                                        group = "interaction(stimulation, type)",
                                        color = "type")) +
    do.call(theme_jetka, args = plot.args) +
    geom_errorbar() +
    geom_point() +
    ggtitle(paste(plot.title))
  return(g.plot)
}


stimulations <- sort(stimulations)
#### plotLogSD  ####
plotLogSD <- function(
  data.exp,
  data.sample,
  plot.title = "",
  ...
){
  data <- 
    rbind(data.exp, data.sample) %>% 
    dplyr::group_by_("type", "stimulation") %>%
    dplyr::summarise_(logresponse = "mean(log(response))",
                      logresponse.sd = "sd(log(response))")
  g.plot <- ggplot(data    = data, 
                   mapping = aes(x = factor(stimulation),
                                 y = logresponse.sd, 
                                 group = factor(type),
                                 fill = factor(type))) +
    geom_bar(stat = "identity", position = "dodge") +
    do.call(theme_jetka, args = plot.args) +
    ggtitle(paste(plot.title)) 
  return(g.plot)
}

