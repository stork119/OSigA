### ###
### plotsScripts
### ###
source("R/scripts/2017-10-05-KarolinaPoster/scriptsLibrary.R")
#### prepare data####
poster.path.list <- getPathsList(normalisation = "ffc_filtered")
poster.data.list <- readRDS(file = poster.path.list$rds.path)

#### plots by wells ####
gplot.list <- list()
r <-foreach(poster.label = labels(poster.data.list)) %do% {
  #poster.label <- labels(poster.data.list)[1]
  tryCatch({
    print(poster.label)
    data <- poster.data.list[[poster.label]]$data %>% data.frame()
    column.names <- getColumnNames(data)
    
    gplot.list[[poster.label]]  <-
      plot_boxplot_group(
        data = data,
        x = column.names$well,
        y = column.names$response,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.label,
        ylim_max = quantile(x = data[,column.names$response], na.rm = TRUE, probs = 0.95)[[1]]
      )
  }, error = function(e){print(e)})
  return()
}

ggsave(filename = paste(poster.path.list$output.dir, paste("experiments", ".pdf", sep = ""), sep = "/"),
       plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
       width = plot.args$width,
       height =plot.args$height,
       useDingbats = plot.args$useDingbats)

#### plotting combinations ####
no_cores <- 6
registerDoParallel(no_cores)
foreach(poster.label = labels(poster.data.list)[19]) %dopar% {
  #poster.label <- labels(poster.data.list)[1]
  tryCatch({
    print(poster.label)
    data <- poster.data.list[[poster.label]]$data %>% 
      data.frame() 
    
    column.names <- getColumnNames(data)
    
    poster.title <- paste(poster.data.list[[poster.label]]$id,
                          poster.data.list[[poster.label]]$title, sep = "")
    path <- paste(poster.path.list$output.dir, poster.data.list[[poster.label]]$id, sep = "/")
    dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
    gplot.list <- list()
    gplot.list$ids  <- 
      plot_boxplot_group(
        data = data, 
        x = column.names$well, 
        y = column.names$response, 
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.title,
        ylim_max = quantile(x = data[,column.names$response], na.rm = TRUE, probs = 0.95)[[1]])
    
    gplot.list$x_time_mean  <- 
      plot_boxplot_group(
        data = data, 
        x = column.names$time, 
        y = column.names$response, 
        facet_grid_group_y = column.names$stimulation,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.title,
        ylim_max = max(quantile(x = data[,column.names$response], na.rm = TRUE, probs = 0.95)[[1]], 1500))
    
    gplot.list$x_time  <- 
      plot_boxplot_group(
        data = data, 
        x = column.names$time, 
        y = column.names$response, 
        boxplot_group = column.names$well,
        facet_grid_group_y = column.names$stimulation,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.title,
        ylim_max = max(quantile(x = data[,column.names$response], na.rm = TRUE, probs = 0.95)[[1]], 1500))
    
    gplot.list$x_stimulation_mean  <- 
      plot_boxplot_group(
        data = data, 
        x = column.names$stimulation, 
        y = column.names$response, 
        facet_grid_group_y =  column.names$time,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.title,
        ylim_max = max(quantile(x = data[,column.names$response], na.rm = TRUE, probs = 0.95)[[1]], 1500))
    
    gplot.list$x_stimulation  <- 
      plot_boxplot_group(
        data = data, 
        x = column.names$stimulation, 
        y = column.names$response, 
        boxplot_group = column.names$well,
        facet_grid_group_y =  column.names$time,
        save_plot = FALSE,
        ylim_max_const = TRUE,
        plot_title = poster.title,
        ylim_max = max(quantile(x = data[,column.names$response], na.rm = TRUE, probs = 0.95)[[1]], 1500))
    
    ggsave(filename = paste(path, 
                            paste(poster.data.list[[poster.label]]$id,
                                  "-experiment",
                                  ".pdf", sep = ""), sep = "/"), 
           plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
           width = plot.args$width,
           height =plot.args$height,
           useDingbats = plot.args$useDingbats)
  }, error = function(e){print(e)})
}
stopImplicitCluster()

#### plots divided by time/stimulation ####
plots.arguments <- list()
plots.arguments$geom_fun <- geom_boxplot
plots.arguments$geom_fun <- geom_violin
plots.arguments$group_fun <- function(column.names){
  return(column.names$stimulation)
}
response_fun <- function(x){log(x)}
plots.arguments$ylimmax <- 10
plots.arguments$title <- ""
plots.arguments$factors <- list(time = FALSE)
plots.arguments$fill <-
  data.frame(fill = c("#ffffff", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6", "#213f72"),
                        stm  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 10),
                        stmnew  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 10))  
plots.arguments$stimulation.list <- NULL 

### args 1
plots.arguments$group_fun <- function(column.names){
  return(column.names$well)
}
title.group_by_well <- "_wells"


plots.arguments$group_fun <- function(column.names){
  return(column.names$stimulation)
}
title.group_by_well <- ""

### args 2
plots.arguments$plot.type <- paste(title.group_by_well, "_violin", sep = "")
plots.arguments$geom_fun <- geom_violin

plots.arguments$plot.type <- paste(title.group_by_well, "_boxplot", sep = "")
plots.arguments$geom_fun <- geom_boxplot

### args 3
plots.arguments$height <- 6
plots.arguments$width  <- 10
plots.arguments$plot_fun.list <- list(plotByTime)

plots.arguments$height <- 6
plots.arguments$width  <- 12
plots.arguments$plot_fun.list <- list(plotByStm)

combinePlot(poster.data.list = poster.data.list, 
            plots.arguments = plots.arguments, 
            plot.args = plot.args,
            poster.path.list = poster.path.list)

#### STAT1 ####
plots.arguments$stimulation.list <- NULL
plots.arguments$fill <- data.frame(fill = c("#ffffff", "#a4d09e", "#a4d09e", "#7abd6f",  "#50af43", "#1fa637", "#177d34", "#175929"),
                        stm  =  c(0, 1.8, 9, 18, 90, 180,  900, 1800),
                        stmnew  =  c(0, 2, 10, 20, 100, 200, 1000,  2000))  

#plots.arguments$ylimmax <- 10
plots.arguments$factors <- list(time = FALSE)

### args 1
plots.arguments$group_fun <- function(column.names){
  return(column.names$well)
}
title.group_by_well <- "_wells"


plots.arguments$group_fun <- function(column.names){
  return(column.names$stimulation)
}
title.group_by_well <- ""

### args 2
plots.arguments$plot.type <- paste(title.group_by_well, "_violin", sep = "")
plots.arguments$geom_fun <- geom_violin
response_fun <- function(x){log(x)}
plots.arguments$ylimmax <- 10 #1000


plots.arguments$plot.type <- paste(title.group_by_well, "_boxplot", sep = "")
plots.arguments$geom_fun <- geom_boxplot
response_fun <- function(x){x}
plots.arguments$ylimmax <- 2000 #1000

### args 3
plots.arguments$height <- 6
plots.arguments$width  <- 10
plots.arguments$plot_fun.list <- list(plotByTime)

plots.arguments$height <- 6
plots.arguments$width  <- 12
plots.arguments$plot_fun.list <- list(plotByStm)


###

###
combinePlot(poster.data.list = poster.data.list, 
            poster.data.list.ids = c(19),
            plots.arguments = plots.arguments, 
            plot.args = plot.args,
            poster.path.list = poster.path.list,
            response_fun = response_fun)
