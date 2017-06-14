### ###
### model plot
### ###

#### plot trajectories ####
plot_trajectories <- function(path,
                              data.trajectory,
                              plot.args = list(),
                              plot.args.ggsave = list(),
                              save = TRUE,
                              filename = "trajectories.pdf"
){
  var.list <- list()
  for(var in (data.trajectory %>% 
              dplyr::distinct(var))$var){
    var.list[[paste("y",as.character(var), sep = "")]] <- var
  }
  var.list$p16p17 <- c(16,17)
  var.list$stat_cyto <- c(1,2,3)
  
  gplot.trajectory.list <- list()
  plot.args.tmp <- plot.args
  plot.args.tmp$theme.title_size <- 18
  for(var.i in 1:length(var.list)){
    gplot.trajectory.list[[labels(var.list)[var.i]]] <- 
      ggplot(data.trajectory %>%  
               dplyr::filter(var %in% var.list[[var.i]]) %>%
               dplyr::group_by(priming, stimulation, time) %>%
               dplyr::summarise(m = sum(m)) %>%
               dplyr::mutate(type = "type"),
             mapping = aes(x = factor(time), y = m, group = type)) +
      geom_point() +
      geom_line() +
      facet_grid(priming ~ stimulation) + 
      do.call(what = theme_jetka, args = append(plot.args, list(theme.text_size = 18))) +
      ggtitle(labels(var.list)[var.i]) +
      expand_limits(y= 0)
  }
  if(save){
    do.call(what = ggsave,
            args = append(plot.args.ggsave,
                          list(filename = paste(path, filename, sep = "/"),
                               plot = marrangeGrob(grobs = gplot.trajectory.list, nrow = 1, ncol = 1))))
  }
  return(gplot.trajectory.list)
}
#### compare model with experimental data ####