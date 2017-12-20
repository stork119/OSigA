### ###
### plots Manuscript
### ###
title.manuscript <- "_manuscript"
#### gamma pSTAT KA11
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
plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_violin", sep = "")
plots.arguments$geom_fun <- geom_violin

plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_boxplot", sep = "")
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

data.manuscript <- list()
data.0.manuscript <- list()
#### gamma pSTAT ####
poster.label  <- "2016-01-28-KA11-C"
data <- poster.data.list[[poster.label]]$data %>% 
      data.frame() 

column.names <- getColumnNames(data)

labels.0 <- c("C01", "C02")
    
data.0 <-  data %>% 
  dplyr::filter_(
    paste(column.names$time, 
          "==", 0, "|",
          column.names$stimulation, 
          "==", 0)) %>%
  dplyr::filter_(
    paste(column.names$well,
          "%in%",
          "c(", 
          paste(
            sapply(labels.0, function(lab){paste("'",lab,"'", sep = "")}),
            collapse = ","),
          ")", 
          sep = " "))

if(nrow(data.0) > 0){
  data.0[,column.names$stimulation] <- 0
  data.0[,column.names$time] <- 0
}

  data <- data %>% 
    dplyr::filter_(
      paste(column.names$time, 
            "!=", 0, "&",
            column.names$stimulation, 
            "!=", 0)) %>%
    rbind(data.0)

  plots.arguments$title <- paste(#poster.data.list[[poster.label]]$title,
    poster.data.list[[poster.label]]$id)
  
  plots.arguments$fill[,column.names$stimulation] <- plots.arguments$fill[,"stm"] 
  
  data.manuscript[[poster.label]] <- data
  data.0.manuscript[[poster.label]] <- data.0
### iterate: time X: stimulation
  plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_violin", sep = "")
  plots.arguments$geom_fun <- geom_violin
  plots.arguments$height <- 6
  plots.arguments$width  <- 8
  plots.arguments$ylimmax <- 10
  
  g.list <-  do.call(
    what = plotByTime,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args)
    )
#### beta + gamma pSTAT ####
  poster.label  <- "2016-01-28-KA11-IFNB"
  data <- poster.data.list[[poster.label]]$data %>% 
    data.frame() 
  
  column.names <- getColumnNames(data)
  
  labels.0 <- c("A07", "A08")
  
  data.0 <-  data %>% 
    dplyr::filter_(
      paste(column.names$time, 
            "==", 0, "|",
            column.names$stimulation, 
            "==", 0)) %>%
    dplyr::filter_(
      paste(column.names$well,
            "%in%",
            "c(", 
            paste(
              sapply(labels.0, function(lab){paste("'",lab,"'", sep = "")}),
              collapse = ","),
            ")", 
            sep = " "))
  
  if(nrow(data.0) > 0){
    data.0[,column.names$stimulation] <- 0
    data.0[,column.names$time] <- 0
  }
  
  data <- data %>% 
    dplyr::filter_(
      paste(column.names$time, 
            "!=", 0, "&",
            column.names$stimulation, 
            "!=", 0)) %>%
    rbind(data.0)
  
  plots.arguments$title <- paste(#poster.data.list[[poster.label]]$title,
    poster.data.list[[poster.label]]$id)
  
  plots.arguments$fill[,column.names$stimulation] <- plots.arguments$fill[,"stm"] 
  
  
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

  data.manuscript[[poster.label]] <- data
  data.0.manuscript[[poster.label]] <- data.0
  ### iterate: time X: stimulation
  plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_violin", sep = "")
  plots.arguments$geom_fun <- geom_violin
  plots.arguments$height <- 6
  plots.arguments$width  <- 8
  plots.arguments$ylimmax <- 10
  
  g.list <-  do.call(
    what = plotByTime,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args)
  )
  
  plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_violin", sep = "")
  plots.arguments$geom_fun <- geom_violin
  plots.arguments$height <- 6
  plots.arguments$width  <- 12
  plots.arguments$ylimmax <- 10
  
  g.list <-  do.call(
    what = plotByStm,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args)
  )
  
#### beta + gamma vs gamma pSTAT ####
  data <- rbind.smart(data.manuscript$`2016-01-28-KA11-C`, data.manuscript$`2016-01-28-KA11-IFNB`)
  data.0 <- rbind.smart(data.0.manuscript$`2016-01-28-KA11-C`, data.0.manuscript$`2016-01-28-KA11-IFNB`)
  
  column.names <- getColumnNames(data)
  
  
  ### iterate: time X: stimulation
  plots.arguments$plot.type <- paste("_compare", title.manuscript, title.group_by_well, "_violin", sep = "")
  plots.arguments$geom_fun <- geom_violin
  plots.arguments$height <- 6
  plots.arguments$width  <- 8
  plots.arguments$ylimmax <- 10
  plots.arguments$group_fun <- function(column.names){
    return(c(column.names$stimulation, column.names$priming))
  }  
  g.list <-  do.call(
    what = plotByTime,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args)
  )
  
  plots.arguments$plot.type <- paste("_compare", title.manuscript, title.group_by_well, "_violin", sep = "")
  plots.arguments$geom_fun <- geom_violin
  plots.arguments$height <- 6
  plots.arguments$width  <- 12
  plots.arguments$ylimmax <- 10
  
  g.list <-  do.call(
    what = plotByStm,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args)
  )
  
  plots.arguments$plot.type <- paste("_compare", title.manuscript, title.group_by_well, "_boxplot", sep = "")
  plots.arguments$geom_fun <- geom_boxplot
  plots.arguments$height <- 6
  plots.arguments$width  <- 8
  plots.arguments$ylimmax <- 10
  plots.arguments$group_fun <- function(column.names){
    return(c(column.names$stimulation, column.names$priming))
  }  
  g.list <-  do.call(
    what = plotByTime,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args)
  )
  
  plots.arguments$plot.type <- paste("_compare", title.manuscript, title.group_by_well, "_boxplot", sep = "")
  plots.arguments$geom_fun <- geom_boxplot
  plots.arguments$height <- 6
  plots.arguments$width  <- 12
  plots.arguments$ylimmax <- 10
  
  g.list <-  do.call(
    what = plotByStm,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args)
  )
  
#### beta stat1 ####

  
  
  poster.label  <- "2017-10-13-KZ83"
  data <- poster.data.list[[poster.label]]$data %>% 
    data.frame() 
  
  column.names <- getColumnNames(data)
  
  labels.0 <- c("A10", "A11", "A12")
  
  data.0 <-  data %>% 
    dplyr::filter_(
      paste(column.names$time, 
            "==", 0, "|",
            column.names$stimulation, 
            "==", 0)) %>%
    dplyr::filter_(
      paste(column.names$well,
            "%in%",
            "c(", 
            paste(
              sapply(labels.0, function(lab){paste("'",lab,"'", sep = "")}),
              collapse = ","),
            ")", 
            sep = " "))
  
  if(nrow(data.0) > 0){
    data.0[,column.names$stimulation] <- 0
    data.0[,column.names$time] <- 0
  }
  
  data <- data %>% 
    dplyr::filter_(
      paste(column.names$time, 
            "==", 24, "&",
            column.names$stimulation, 
            "==", 900)) %>%
    rbind(data.0)
  
  plots.arguments$title <- paste("STAT1 copy number")
  
  
  
  ### iterate: time X: stimulation
  
  plots.arguments$stimulation.list <- NULL
  plots.arguments$fill <- data.frame(fill = c("#ffffff", "#a4d09e", "#a4d09e", "#7abd6f",  "#50af43", "#1fa637", "#177d34", "#175929"),
                                     stm  =  c(0, 1.8, 9, 18, 90, 180,  900, 1800),
                                     stmnew  =  c(0, 2, 10, 20, 100, 200, 1000,  2000))  
  plots.arguments$fill[,column.names$stimulation] <- plots.arguments$fill[,"stm"] 
  
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
  plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_violin", sep = "")
  plots.arguments$geom_fun <- geom_violin
  response_fun <- function(x){log(x)}
  plots.arguments$ylimmax <- 10 #1000
  
  
  plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_boxplot", sep = "")
  plots.arguments$geom_fun <- geom_boxplot
  response_fun <- function(x){x}
  plots.arguments$ylimmax <- 1500 #1000
  
  ### args 3
  plots.arguments$height <- 6
  plots.arguments$width  <- 10
  plots.arguments$plot_fun.list <- list(plotByTime)
  
  g.list <-  do.call(
    what = plotByTime,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args,
                response_fun = response_fun)
  )
  
  plots.arguments$height <- 6
  plots.arguments$width  <- 12
  plots.arguments$plot_fun.list <- list(plotByStm)
  
  g.list <-  do.call(
    what = plotByStm,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args,
                response_fun = response_fun)
  )
  
  
  ### KZ73 KZ78
  poster.label  <- "2017-08-25-KZ73-Nuclei"
  data.1 <- poster.data.list[[poster.label]]$data %>% 
    data.frame() 
  
  poster.label  <- "2017-09-15-KZ78-Nuclei"
  data.2 <- poster.data.list[[poster.label]]$data %>% 
    data.frame() 
  
  poster.label <- "KZ73-KZ78"
  
  data <- rbind.smart(data.1, data.2)
  
  column.names <- getColumnNames(data)
  
  # labels.0 <- c("A10", "A11", "A12")

  data.0 <-  data %>%
    dplyr::filter_(
      paste(column.names$time,
            "==", 0, "|",
            column.names$stimulation,
            "==", 0))# %>%
    # dplyr::filter_(
    #   paste(column.names$well,
    #         "%in%",
    #         "c(",
    #         paste(
    #           sapply(labels.0, function(lab){paste("'",lab,"'", sep = "")}),
    #           collapse = ","),
    #         ")",
    #         sep = " "))

  if(nrow(data.0) > 0){
    data.0[,column.names$stimulation] <- 0
    data.0[,column.names$time] <- 0
  }
  
  data <- data %>% 
    dplyr::filter()
  
  data <- data %>%
    dplyr::filter_(
      paste(column.names$time,
            "==", 24, "&",
            column.names$stimulation,
            "==", 900)) %>%
    rbind(data.0)

   plots.arguments$title <- paste("STAT1 copy number")
  
  
  
  ### iterate: time X: stimulation
  
  plots.arguments$stimulation.list <- NULL
  plots.arguments$fill <- data.frame(fill = c("#ffffff", "#a4d09e", "#a4d09e", "#7abd6f",  "#50af43", "#1fa637", "#177d34", "#175929"),
                                     stm  =  c(0, 1.8, 9, 18, 90, 180,  900, 1800),
                                     stmnew  =  c(0, 2, 10, 20, 100, 200, 1000,  2000))  
  plots.arguments$fill[,column.names$stimulation] <- plots.arguments$fill[,"stm"] 
  
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
  plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_violin", sep = "")
  plots.arguments$geom_fun <- geom_violin
  response_fun <- function(x){log(x)}
  plots.arguments$ylimmax <- 10 #1000
  
  
  plots.arguments$plot.type <- paste(title.manuscript, title.group_by_well, "_boxplot", sep = "")
  plots.arguments$geom_fun <- geom_boxplot
  response_fun <- function(x){x}
  plots.arguments$ylimmax <- 2000 #1000
  
  ### args 3
  plots.arguments$height <- 4
  plots.arguments$width  <- 6
  plots.arguments$plot_fun.list <- list(plotByTime)
  
  g.list <-  do.call(
    what = plotByTime,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args,
                response_fun = response_fun)
  )
  
  plots.arguments$height <- 6
  plots.arguments$width  <- 12
  plots.arguments$plot_fun.list <- list(plotByStm)
  
  g.list <-  do.call(
    what = plotByStm,
    args = list(data = data,
                column.names = column.names,
                plots.arguments = plots.arguments,
                poster.label = poster.label,
                data.0 = data.0,
                poster.path.list = poster.path.list,
                plot.args = plot.args,
                response_fun = response_fun)
  )
  
  data.sum <- data %>% dplyr::group_by(time.1.1) %>% dplyr::summarise(mean(Intensity_MeanIntensity_Alexa555))