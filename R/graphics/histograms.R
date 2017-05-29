### ###
### histograms 
### ###
source("R/graphics/theme_jetka.R")

#### plot_boxplot_group ####
plot_density <- function(data,
                         ...,
                         x,
                         group,
                         color = group,
                         title = "",
                         xlim_max = NULL,
                         compare_to_all = FALSE
){
  gplot <- 
    ggplot(data, aes_string(x = x,
                            group = group,
                            color = color)) +
    geom_density() +
    theme_jetka(...)  + 
    ggtitle(title)
  
  if(!is.null(xlim_max)){
    gplot <- gplot + xlim(c(0,xlim_max))
  }
  
  if(compare_to_all){
    gplot <- gplot + geom_density(data = 
                                    data %>% 
                                    dplyr::mutate_(
                                      .dots = 
                                        setNames(
                                          list(lazyeval::interp(~ a, a = 'all')),
                                          group)),
                                  mapping = aes_string(x = x))
  }
  return(gplot)
}