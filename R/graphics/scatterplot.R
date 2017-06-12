### ###
### scatterplot
### ###

scatterplot_lieklihood <- function(
  data,
  colnames.list = c("p1", "p6", "stat1"),
  path.list,
  filename = "likelihood_scatterplot.pdf",
  best_percentage = 1.01,
  text_size = 32,
  width = 12,
  height = 9
){
  grob.list <- list()
  
  grobs.grid <- expand.grid(x = colnames.list, y = colnames.list)
  
  data.best <- data %>% dplyr::filter(likelihood < ((data %>% dplyr::summarise(minlik = min(likelihood)))$minlik*best_percentage))
  for(grobs.grid.i in 1:nrow(grobs.grid)){
    grob <- grobs.grid[grobs.grid.i,]
    if(grob$x == grob$y){
      grob.list[[grobs.grid.i]] <- 
        ggplot(
          data = data.frame(text = grob$x, x = 1, y = 1),
          aes(x = x, y = y, label = text)) +
        geom_text(size = text_size) +
        do.call(what = theme_jetka, args = plot.args) 
    } else {
      x_character  <- paste("log(", as.character(grob$x), ")")
      y_character  <- paste("log(", as.character(grob$y), ")")
      grob.list[[grobs.grid.i]] <- 
        ggplot(data = data,
               mapping = aes_string(x = x_character,
                                    y = y_character,
                                    color = "log(likelihood)")) + 
        geom_point(size = 2 ) + 
        geom_point(data = data.best,
                   aes_string(
                     x = x_character,
                     y = y_character),
                   color = "red")  +
        do.call(what = theme_jetka, args = plot.args) 
    }
  }
  scatterplot <- marrangeGrob(grobs = grob.list, nrow = length(colnames.list), ncol = length(colnames.list))
  plot.args.ggsave.tmp <- plot.args.ggsave
  plot.args.ggsave.tmp$width <- width
  plot.args.ggsave.tmp$height <- height
  do.call(what = ggsave,
          args = append(plot.args.ggsave.tmp,
                        list(filename = paste(path.list$optimisation.analysis, filename, sep = "/"),
                             plot = scatterplot)))
  return(scatterplot)
  
}