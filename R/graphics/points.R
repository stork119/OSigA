### ###
### error plot
### ###


plot_points <- function(data,
                        ...,
                        x,
                        y,
                        yerror = "",
                        facet_grid_group_x = "",
                        facet_grid_group_y = "",
                        title = "",
                        x.as_factor = TRUE,
                        ylim_min = 0,
                        ylim_max = NULL
){
  if(x.as_factor){
    data[,x] <- factor(x = data[,x][[1]])
  }
  
  gplot <- ggplot(data = data,
                  mapping = aes_string(x = x,
                                       y = y)) +
    geom_point() +
    theme_jetka(...) +
    ggtitle(title)
  if(facet_grid_group_y != ""){
    gplot <- gplot + 
      facet_grid(paste(facet_grid_group_x, "~", facet_grid_group_y, sep = " "),
                 scale ="free",
                 space = "free") 
  }
  
  if(yerror != ""){
    gplot <- gplot +
      geom_errorbar(mapping = aes_string(x = x,
                                         ymin = paste(y, "-", yerror), 
                                         ymax = paste(y, "+", yerror)))
  }
  
  if(is.null(ylim_max)){
    ylim_max <- max(data[,y])
    if(yerror != ""){
      ylim_max <- ylim_max + max(data[,yerror])
    }
    ylim_max <- 1.1*ylim_max
  }
  
  gplot <- gplot + ylim(c(ylim_min, ylim_max))
  
  return(gplot)
}