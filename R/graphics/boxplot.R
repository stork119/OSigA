### ###
### boxplot 
### ###
source("R/graphics/theme_jetka.R")

#### plot_boxplot_group ####
plot_boxplot_group <- function(data,
                               ...,
                               output_path = NULL,
                               filename = NULL,
                               x = "time",
                               y = "intensity",
                               boxplot_group = x,
                               facet_grid_group_y = "",
                               facet_grid_group_x = "",
                               ylab = y,
                               xlab = x,
                               ylim_min = 0,
                               ylim_max = 2000,
                               plot_width = 24,
                               plot_height = 8,
                               plot_title = "",
                               xlab_angle = 90,
                               xlab_hjust = 0,
                               legend_position = "bottom",
                               plot_fun = "geom_boxplot",
                               normalize_data   = FALSE, 
                               normalize_factor = 65535,
                               ylim_max_const = TRUE,
                               x_factor = TRUE,
                               save_plot = TRUE){
  CheckColumnExistence <- function(data, columns.list = list()){
    columns_existance <- (unlist(columns.list) %in% colnames(data))
    if_exists <- sum(!columns_existance) == 0
    if(!if_exists){
      print(unlist(columns.list)[-which(columns_existance)])
    }
    return(if_exists)
  }
  ylim_min <- as.integer(ylim_min)
  ylim_max  <- as.integer(ylim_max)
  plot_width <- as.integer(plot_width)
  plot_height <- as.integer(plot_height)
  xlab_angle <- as.integer(xlab_angle)
  xlab_hjust <- as.integer(xlab_hjust)
  normalize_data <- as.integer(normalize_data)
  normalize_factor <- as.integer(normalize_factor)
  ylim_max_const <- as.integer(ylim_max_const)
  x_factor  <- as.integer(x_factor)
  
  if(!CheckColumnExistence(data = data, list(x,y,boxplot_group))){
    return()
  }
  
  if(normalize_data){
    data[,y] <- normalize_factor*data[,y]
  }
  if(!ylim_max_const){
    ylim_max <- 1.2*max(data[,y])
  }
  if(x_factor){
    data[,x] <- factor(data[,x])
  }
  
  #  data$bg <- factor(sapply(1:nrow(data),function(i){paste(as.character(data[i,boxplot_group]), sep = " ", collapse = " ")}))
  
  gplot <- ggplot(data = data, 
                  aes_string(x = x,
                             y = y,
                             group = boxplot_group)
  ) + 
    do.call(plot_fun, args = list()) +
    ylim(ylim_min, ylim_max) +
    xlab(xlab) +
    ylab(ylab) + 
    ggtitle(plot_title) +
    theme_jetka(...)
  if(facet_grid_group_x != "" || facet_grid_group_y != ""){
    gplot <- 
      gplot + 
      facet_grid(paste(facet_grid_group_x, "~", facet_grid_group_y, sep = " "),
                 scale ="free",
                 space = "free")
  }
  try({
    if(save_plot){
      output_path <- normalizePath(output_path, "/")
      dir.create(path = output_path, recursive = TRUE, showWarnings = FALSE)
      ggsave(filename = paste(output_path, "/", filename, ".pdf", sep = ""),
             plot = gplot,
             width = plot_width,
             height = plot_height,
             useDingbats = FALSE)
    }
  })
  return(gplot)
}