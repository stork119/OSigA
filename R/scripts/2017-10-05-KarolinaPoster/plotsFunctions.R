### ###
### plots functions
### ###

#### plotByTime ####
plotByTime <- function(
  data,
  column.names,
  plots.arguments,
  poster.label,
  data.0,
  plot.args,
  poster.path.list,
  response_fun,
  ...
){
  
  data <- data %>% 
    dplyr::mutate_(.dots = 
                     setNames(object = 
                                list(paste("response_fun(", column.names$response, ")")), 
                              nm = column.names$response))
  
  data.0 <- data.0 %>% 
    dplyr::mutate_(.dots = 
                     setNames(object = 
                                list(paste("response_fun(", column.names$response, ")")), 
                              nm = column.names$response))
  
  gplot.list <- list()  
  colors.list <- list()
  for( t in (data %>%  dplyr::distinct_(column.names$time) %>% 
             data.frame())[,column.names$time])
  {    tryCatch({
      print(t)
      # t <- 30
      data.time <- data %>% 
        dplyr::filter_(paste(column.names$time, "==", t)) 
      
      if(nrow(data.time %>% dplyr::filter_(paste(column.names$stimulation, "==", 0))) == 0){
        data.time <- rbind(
          data.time,
          (data.0 %>% 
             dplyr::mutate_(.dots = 
                              setNames(object = 
                                         list(t), 
                                       nm = column.names$time)) %>%
             dplyr::mutate_(.dots = 
                              setNames(object = 
                                         list(0), 
                                       nm = column.names$stimulation) )
          ))
      }
      data.time <- data.time %>% 
        dplyr::arrange_(column.names$stimulation) %>%
        dplyr::left_join(plots.arguments$fill)
      
      data.time[,column.names$stimulation] <- data.time[,"stmnew"]
      
      colors.list[[as.character(t)]] <- as.character(
        (plots.arguments$fill %>% 
           dplyr::filter(stmnew %in%
                           as.character(unique(data.time[, column.names$stimulation]))))[,"fill"])
      
      data.time[, column.names$stimulation] <-
        factor(data.time[, column.names$stimulation])
      
      
      g<-
        ggplot(data = data.time, 
                  mapping = aes_string(
                    x = column.names$stimulation,
                    y = column.names$response,
                    group =  paste(
                      "interaction(",
                      paste(do.call(plots.arguments$group_fun, list(column.names)),
                            collapse = ", "),
                      ")"),
                    fill = column.names$stimulation)) +
        ylim(c(0,plots.arguments$ylimmax)) + 
        # do.call(what = plots.arguments$geom_fun,
        #                                              args = list()) +
        theme_jetka() +
        ggtitle(
          paste("time:", t, ";", plots.arguments$title)) + 
        xlab("stimulation") +
        ylab("Fluorescent Intensity (A.U.)")
      if(length(colors.list[[as.character(t)]]) > 1){
        x <- colors.list[[as.character(t)]]
        force(x)
        g <- g+
          do.call(what = plots.arguments$geom_fun,
                               args = list()) +
          scale_fill_manual(values = x)
      } else {
        x <- colors.list[[as.character(t)]]
        force(x)
         g <- g +
           do.call(what = plots.arguments$geom_fun,
                   args = list(fill = x))
      }
      gplot.list[[as.character(t)]] <- ggplot_build(g)$plot
    }, error = function(e){print(e)})
    
  }
  dir.create(path =  paste(poster.path.list$output.dir, 
                           poster.label,
                           sep = "/"),
             showWarnings = FALSE,
             recursive = TRUE)
  
  ggsave(filename = paste(poster.path.list$output.dir,
                          poster.label,
                          paste(poster.label,
                                plots.arguments$plot.type,
                                "_Y_fluorescent_X_stimulation.pdf",
                                sep = ""),
                          sep = "/"),
         plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
         width = plots.arguments$width,
         height = plots.arguments$height,
         useDingbats = plot.args$useDingbats)
  return()#list(plot = gplot.list, colors = colors.list))
}

#### plotByStm ####
plotByStm <- function(
  data,
  column.names,
  plots.arguments,
  poster.label,
  plot.args,
  poster.path.list,
  response_fun,
  ...
){
  
  data <- data %>% 
    dplyr::mutate_(.dots = 
                     setNames(object = 
                                list(paste("response_fun(", column.names$response, ")")), 
                              nm = column.names$response))
  
  data <- data %>% dplyr::left_join(plots.arguments$fill)
  data[,column.names$stimulation] <- data[,"stmnew"]
  if(is.null(plots.arguments$stimulation.list)){
    plots.arguments$stimulation.list <- sort(unique(data[,column.names$stimulation]))
  }
  
  gplot.list <- list()
  for(stimulation.i in 1:length(plots.arguments$stimulation.list)){
    stimulations <- plots.arguments$stimulation.list[[stimulation.i]]
    colors.list <- 
      unique(
        as.character(
          (plots.arguments$fill %>% 
             dplyr::filter(stmnew %in% stimulations))[,"fill"]))
    if(!is.null(plots.arguments$factors$time) & plots.arguments$factors$time){
      data[, column.names$time] <- factor(data[,column.names$time])
    }
    g <- ggplot(data = 
                  data %>% 
                  dplyr::filter_(paste(column.names$stimulation, 
                                       "%in%", 
                                       "c(", 
                                       paste(stimulations, collapse = ","),
                                       ")")),
                mapping = 
                  aes_string(
                    x = column.names$time,
                    y = column.names$response,
                    group = paste(
                      "interaction(", 
                      column.names$time, ",", 
                      paste(
                        do.call(plots.arguments$group_fun, 
                                list(column.names)),
                        collapse = ","),
                      ")"),
                    fill = paste("factor(", column.names$stimulation, ")"))) +
      ylim(c(0,plots.arguments$ylimmax)) + 
      theme_jetka() +
      ggtitle(
        paste("stimulation:",
              paste(stimulations, collapse = ","),
              ";", 
              plots.arguments$title)) + 
      xlab("time") +
      ylab("Fluorescent Intensity (A.U.)")
    if(length(colors.list) == 1){
      g <- g + 
        do.call(plots.arguments$geom_fun, list(fill = colors.list))
    }  else {
      g <- g +
        do.call(plots.arguments$geom_fun, list()) + 
        scale_fill_manual(values = colors.list)
    }
    gplot.list[[stimulation.i]] <- (ggplot_build(g))$plot
  }
  ggsave(filename = paste(poster.path.list$output.dir, 
                          poster.label,
                          paste(poster.label, 
                                plots.arguments$plot.type,
                                "_Y_fluorescent_X_time_", 
                                ".pdf", sep = ""), sep = "/"), 
         plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow =1),
         width = plots.arguments$width,
         height = plots.arguments$height,
         useDingbats = plot.args$useDingbats)
}


#### combinePlot ####
combinePlot <- 
  function(poster.data.list, 
           plots.arguments, 
           plot.args,
           poster.path.list,
           response_fun = function(x){x},
           poster.data.list.ids = NULL,
           no_cores = 6){
    if(is.null(poster.data.list.ids)){
      poster.data.list.ids <- 1:length(poster.data.list)
    }
    registerDoParallel(no_cores)
    foreach(poster.label = labels(poster.data.list)[poster.data.list.ids]) %dopar% {
      tryCatch({
        print(poster.label)
        
        data <- poster.data.list[[poster.label]]$data %>% 
          data.frame() 
        
        column.names <- getColumnNames(data)
        
        plots.arguments$title <- paste(#poster.data.list[[poster.label]]$title,
          poster.data.list[[poster.label]]$id)
        
        plots.arguments$fill[,column.names$stimulation] <- plots.arguments$fill[,"stm"] 
        
        data.0 <-  data %>% 
          dplyr::filter_(paste(column.names$time, "==", 0)) 
        
        if(nrow(data.0) > 0){
          data.0[,column.names$stimulation] <- 0
        }
        
        data <- data %>% rbind(data.0)
        
        # data <- data %>% 
        #   dplyr::mutate_(.dots = 
        #                    setNames(object = 
        #                               list(paste("response_fun(", column.names$response, ")")), 
        #                             nm = column.names$response))
        # 
        # data.0 <- data.0 %>% 
        #   dplyr::mutate_(.dots = 
        #                    setNames(object = 
        #                               list(paste("response_fun(", column.names$response, ")")), 
        #                             nm = column.names$response))
        
        ### iterate: time X: stimulation
        l <- lapply(plots.arguments$plot_fun.list,
                    do.call, args = list(data = data,
                                         column.names = column.names,
                                         plots.arguments = plots.arguments,
                                         poster.label = poster.label,
                                         data.0 = data.0,
                                         poster.path.list = poster.path.list,
                                         plot.args = plot.args,
                                         response_fun = response_fun))
      })
    }
    stopImplicitCluster()
  }
