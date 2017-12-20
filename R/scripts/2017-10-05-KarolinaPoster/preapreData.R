### ### 
### script - preparing data 
### ###

source("R/optimisation/initialise_optimisation.R")

poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"
poster.path.list$input.rds.dir <- "resources/input/poster/2017-12-18/rds/"
dir.create(poster.path.list$input.rds.dir, recursive = TRUE)

poster.path.list$output.dir <- "resources/output/poster/2017-12-18/"
dir.create(poster.path.list$output.dir)

#### reading data ####
poster.data.list <- list()
poster.path.list$experiments <- 
  list(
    list(id = "2015-12-13-KA07-C",    title = "pSTAT in nuclei; stm: gamma"),
    list(id = "2015-12-13-KA07-IFNB", title = "pSTAT in nuclei; prm: beta; stm: gamma"),
    list(id = "2015-12-13-KA07-IFNG", title = "pSTAT in nuclei; prm: gamma; stm: gamma"),
    list(id = "2015-12-29-KA08-C",    title = "pSTAT in nuclei; stm: gamma"),
    list(id = "2015-12-29-KA08-IFNB", title = "pSTAT in nuclei; prm: beta; stm: gamma"),
    list(id = "2015-12-29-KA08-IFNG", title = "pSTAT in nuclei; prm: gamma; stm: gamma"),
    list(id = "2016-01-21-KA09-C",    title = "pSTAT in nuclei; stm: gamma"),
    list(id = "2016-01-21-KA09-IFNB", title = "pSTAT in nuclei; prm: beta; stm: gamma"),
    # list(id = "2016-01-21-KA09-IFNG", title = "pSTAT in nuclei; prm: gamma; stm: gamma"),
    list(id = "2016-01-26-KA10-C",    title = "pSTAT in nuclei; stm: gamma"),
    list(id = "2016-01-26-KA10-IFNB", title = "pSTAT in nuclei; prm: beta; stm: gamma"),
    list(id = "2016-01-26-KA10-IFNG", title = "pSTAT in nuclei; prm: gamma; stm: gamma"),
    list(id = "2016-01-28-KA11-C",    title = "pSTAT in nuclei; stm: gamma"),
    list(id = "2016-01-28-KA11-IFNB", title = "pSTAT in nuclei; prm: beta; stm: gamma"),
    list(id = "2016-01-28-KA11-IFNG", title = "pSTAT in nuclei; prm: gamma; stm: gamma"),
    list(id = "2017-07-18-KS25", path = "2017-07-18-KS25/488", title = "IRF1 in nuclei; stm: gamma"),
    list(id = "2017-07-18-KS26", path = "2017-07-18-KS26/488", title = "IRF1 in nuclei; prm: beta; stm: gamma"),
    list(id = "2017-07-25-KS27", path = "2017-07-25-KS27/488", title = "IRF1 in nuclei; stm: gamma"),
    list(id = "2017-07-25-KS28", path = "2017-07-25-KS28/488", title = "IRF1 in nuclei; prm: beta; stm: gamma"),
    list(id = "2017-10-13-KZ83", path = "2017-10-13-KZ83/555", title = "STAT1 in nuclei; stm: beta"),
    # list(id = "2017-07-19-KZ63", path = "2017-07-19-KZ63/488"),
    list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", output_id = "2017-08-25-KZ73-Nuclei"),
    list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", output_id = "2017-09-15-KZ78-Nuclei")
    # list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", file = "Cells555.csv", output_id = "2017-08-25-KZ73-Cells"),
    # list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", file = "Cells555.csv", output_id = "2017-09-15-KZ78-Cells"),
    # list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", file = "Cytoplasm555.csv", output_id = "2017-08-25-KZ73-Cytoplasm"),
    # list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", file = "Cytoplasm555.csv", output_id = "2017-09-15-KZ78-Cytoplasm")
    )


g <- foreach(experiment = poster.path.list$experiments) %do% {
  tryCatch({
    if(is.null(experiment[["file"]])){
      experiment$file <- "ShrinkedNuclei.csv"
    }
    if(is.null(experiment[["path"]])){
      experiment$path <- experiment$id
    }
    path <- paste(poster.path.list$input.dir,
                  experiment$path,
                  "ffc/data_quantify/",
                  experiment$file,
                  sep = "/")
   # print(path)
    if(is.null(experiment[["output_id"]])){
      experiment$output_id <- experiment$id
    }
    poster.data.list[[experiment$output_id]] <- list()
    poster.data.list[[experiment$output_id]]$data <-
      read.table(file = path,
                 header = TRUE,
                 sep = ",") %>%
      normalize_data() %>%
      data.table()
    poster.data.list[[experiment$output_id]]$id <- experiment$id
    poster.data.list[[experiment$output_id]]$title <- experiment$title
  }, error = function(e){ print(e) })
  return()
}


#### daving RDS ####
saveRDS(object = poster.data.list, 
        file = paste(poster.path.list$input.rds.dir, "data_ffc.RDS", sep = "/"))


#### filter data ####
rds.path <- paste(poster.path.list$input.rds.dir, "data_ffc.RDS", sep = "/")
poster.data.list <- readRDS(file = rds.path)

poster.filtered.list <- list(
  list(id = "2016-01-26-KA10-C", well.names = c('A01', 'B01', 'H12', 'D10')),
  list(id = "2016-01-26-KA10-IFNB" , well.names = c('A01', 'A02', 'A03', 'A04', 'D09')),
  list(id = "2016-01-26-KA10-IFNG", well.names = c('A01', 'D02', 'H12')),
  list(id = '2016-01-28-KA11-C', well.names = c(paste('A0', 1:9, sep = ''), 'A10', 'A11', 'A12', 'B01', 'B11', 'B12', 'D09')),
  list(id = '2016-01-28-KA11-IFNB', well.names = c('A01', 'C01','D01', 'D02', 'D10', 'E01', 'F01', 'B01', 'G01', 'H01')),
  list(id = '2016-01-28-KA11-IFNG', well.names = c(paste('A0', 1:3, sep = ''), 'B01', 'C01', 'D01', 'E01', 'E10', 'F01')),
  list(id = '2017-07-18-KS26', well.names = c('D01')),
  list(id = '2017-07-25-KS27', well.names = c('D02','F01', 'F02', 'F10', 'H04')),
  list(id = '2017-07-25-KS28', well.names = c('D02', 'D03', 'H01', 'H02')),
  list(id = '2017-10-13-KZ83', well.names = c('A13', 'A14', 'A15', 'A16', 'A17', 'A18', 'A19',
                                              'A20', 'A21', 'A22', 'A23', 'B21', 'C15', 'C18',
                                              'C22', 'D18', 'D22', 'E06', 'E17', 'E19', 'E22',
                                              'G13', 'G19', 'H17')),
  # list(id = '2017-08-25-KZ73-Nuclei', well.names = c('A11', paste(LETTERS[1:8], '01', sep = ""))),
  # list(id = '2017-08-25-KZ73-Cytoplasm', well.names = c(paste('A0', 1:4, sep = ''), 'A11', 'A12', paste('B0', 1:3, sep = ''), 'H01', 'H12',"well.name")),
  # list(id = '2017-08-25-KZ73-Cells', well.names = c(paste('A0', 1:4, sep = ''), 'A11', 'A12', paste('B0', 1:3, sep = ''), 'H01', 'H12', "well.name"),
   list(id = '2017-08-25-KZ73-Nuclei', 
        well.names = c(
          as.vector(matrix(sapply(LETTERS[c(1,3:8)], 
            function(l) {
              paste(l, c(paste('0', 3:9, sep = ""),
                         "10", "11", "12"), sep = "")}), ncol = 1)), "A01", "A02")),
  # list(id = '2017-08-25-KZ73-Cells', well.names = c(as.vector(matrix(sapply(LETTERS[c(1,3:8)], function(l) {paste(l, c(paste('0', 3:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), "A01", "A02")),
  # list(id = '2017-08-25-KZ73-Cytoplasm', well.names = c(as.vector(matrix(sapply(LETTERS[c(1,3:8)], function(l) {paste(l, c(paste('0', 3:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), "A01", "A02")),
   list(id = '2017-09-15-KZ78-Nuclei',
        well.names = 
          c(as.vector(matrix(sapply(LETTERS[c(2)], function(l) {
            paste(l, c(paste('0', 1:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)),
          as.vector(matrix(sapply(LETTERS[1:8], function(l) {paste(l, c("01", "02"), sep = "")}), ncol = 1))))
  # list(id = '2017-09-15-KZ78-Cells', well.names = c(as.vector(matrix(sapply(LETTERS[c(2)], function(l) {paste(l, c(paste('0', 1:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), 
  #                                                   as.vector(matrix(sapply(LETTERS[1:8], function(l) {paste(l, c("01", "02"), sep = "")}), ncol = 1)))),
  # list(id = '2017-09-15-KZ78-Cytoplasm', well.names = c(as.vector(matrix(sapply(LETTERS[c(2)], function(l) {paste(l, c(paste('0', 1:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), 
  #                                                       as.vector(matrix(sapply(LETTERS[1:8], function(l) {paste(l, c("01", "02"), sep = "")}), ncol = 1))))
  
  )


### filters to text file
filtered.file <- paste(poster.path.list$input.dir, "ffc_filtered.txt", sep = "/")
file.remove(filtered.file)
lapply(poster.filtered.list, 
       function(filtered.list, ...){
         write(x = unlist(filtered.list), ...)},
       file = filtered.file,
       append =TRUE,
       ncolumns = 1000
       )

###
poster.data.list.2 <- poster.data.list
for(poster.filter in poster.filtered.list){
  #poster.filter <- poster.filtered.list[[1]]
  poster.data.list.2[[poster.filter$id]]$data <-
    poster.data.list[[poster.filter$id]]$data[
      -which(poster.data.list[[poster.filter$id]]$data$well.name %in% poster.filter$well.names),]
}
poster.data.list <- poster.data.list.2
saveRDS(object = poster.data.list, 
        file = paste(poster.path.list$input.rds.dir, "data_ffc_filtered.RDS", sep = "/"))

#### connect data ####
rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_filtered.RDS"
poster.data.list <- readRDS(file = rds.path)

poster.data.list.new <- list()
# BetaStat1
part.list <- list("Nuclei", "Cytoplasm", "Cells")
for(part in part.list){
#   exp.id <- "2017-09-15-KZ78"
#   data.KZ78 <- poster.data.list[[paste(exp.id, part, sep = "-")]] %>% 
#     dplyr::mutate(Intensity_MeanIntensity_Alexa = Intensity_MeanIntensity_Alexa555) %>%
#     dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
#     dplyr::mutate(experiment = exp.id) %>%
#     dplyr::mutate(time = time.1.1) %>%
#     dplyr::mutate(stimulation = stimulation.1.1) %>%
#     dplyr::select(experiment, well.name, time, stimulation, Intensity_MeanIntensity_Alexa) %>%
#     data.table()
#     
#   exp.id <- "2017-08-25-KZ73"
#   data.KZ73 <- poster.data.list[[paste(exp.id, part, sep = "-")]] %>% 
#     dplyr::mutate(Intensity_MeanIntensity_Alexa = Intensity_MeanIntensity_Alexa555) %>%
#     dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
#     dplyr::mutate(experiment = exp.id) %>%
#     dplyr::mutate(time = time.1.1) %>%
#     dplyr::mutate(stimulation = stimulation.1.1) %>%
#     dplyr::select(experiment, well.name, time, stimulation, Intensity_MeanIntensity_Alexa)%>%
#     data.table()
#   
#   poster.data.list.new[[paste("Beta-Stat", part, sep = "-")]] <- 
#   rbind(data.KZ73, data.KZ78)
# }
# 
# # Beta-pStat
# exp.id <- "2017-07-19-KZ63"
# poster.data.list.new[["Beta-pStat-Nuclei"]] <- 
#   poster.data.list[[exp.id]] %>%
#   dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
#   dplyr::mutate(experiment = exp.id)  %>%
#   dplyr::mutate(priming = 0)  %>%
#   dplyr::mutate(time = time.1.1) %>%
#   dplyr::mutate(stimulation = stimulation.1.1) %>%
#   dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa)%>%
#   data.table()
  
# Beta-gamma

exp.id <- "2016-01-26-KA10-C"
poster.data.list.new[["Gamma-pStat-Nuclei"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(priming = 0)  %>%
  dplyr::mutate(time = time.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa)%>%
  data.table()

exp.id <- "2016-01-28-KA11-C"
poster.data.list.new[["Gamma-pStat-Nuclei"]] <- 
  poster.data.list.new[["Gamma-pStat-Nuclei"]] %>%
  rbind(
    poster.data.list[[exp.id]] %>%
    dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
    dplyr::mutate(experiment = exp.id)  %>%
    dplyr::mutate(priming = 0)  %>%
      dplyr::mutate(time = time.1.1) %>%
      dplyr::mutate(stimulation = stimulation.1.1) %>%
    dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa)%>%
    data.table()
  )

# Gamma-Stat
# exp.id <- "2016-01-26-KA10-C"
# poster.data.list.new[["Gamma-pStat"]] <- 
#   poster.data.list[[exp.id]] %>%
#   dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
#   dplyr::mutate(experiment = exp.id)  %>%
#   dplyr::select(experiment, well.name, time.1.1, stimulation.1.1, Intensity_MeanIntensity_Alexa)%>%
#   data.table()

# Beta-Gamma-pStat
exp.id <- "2016-01-26-KA10-IFNG"
poster.data.list.new[["Gamma-Gamma-pStat-Nuclei"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(time = time.1.1) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa)%>%
  data.table()

exp.id <- "2016-01-28-KA11-IFNG"
poster.data.list.new[["Gamma-Gamma-pStat-Nuclei"]] <- 
  poster.data.list.new[["Gamma-Gamma-pStat-Nuclei"]] %>%
  rbind(
    poster.data.list[[exp.id]] %>%
      dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
      dplyr::mutate(experiment = exp.id)  %>%
      dplyr::mutate(time = time.1.1) %>%
      dplyr::mutate(priming = priming.1.1) %>%
      dplyr::mutate(stimulation = stimulation.1.1) %>%
      dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa)%>%
      data.table()
  )

# Beta-Beta-pStat
exp.id <- "2016-01-26-KA10-IFNB"
poster.data.list.new[["Beta-Gamma-pStat-Nuclei"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(time = time.1.1) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
  data.table()

exp.id <- "2016-01-28-KA11-IFNB"
poster.data.list.new[["Beta-Gamma-pStat-Nuclei"]] <- 
  poster.data.list.new[["Beta-Gamma-pStat-Nuclei"]] %>%
  rbind(
    poster.data.list[[exp.id]] %>%
      dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
      dplyr::mutate(experiment = exp.id)  %>%
      dplyr::mutate(time = time.1.1) %>%
      dplyr::mutate(priming = priming.1.1) %>%
      dplyr::mutate(stimulation = stimulation.1.1) %>%
      dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
      data.table()
  )

# Beta-IRF
exp.id <- "2017-07-18-KS25"
poster.data.list.new[["Gamma-IRF-Nuclei"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(time = time.2.1) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
  data.table()

exp.id <- "2017-07-25-KS27"
poster.data.list.new[["Gamma-IRF-Nuclei"]] <- 
  poster.data.list.new[["Gamma-IRF-Nuclei"]] %>%
  rbind(
    poster.data.list[[exp.id]] %>%
      dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
      dplyr::mutate(experiment = exp.id)  %>%
      dplyr::mutate(time = time.2.1) %>%
      dplyr::mutate(priming = priming.1.1) %>%
      dplyr::mutate(stimulation = stimulation.1.1) %>%
      dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
      data.table()
  )
# Beta-Gamma-IRF
exp.id <- "2017-07-18-KS26"
poster.data.list.new[["Beta-Gamma-IRF-Nuclei"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(time = time.2.1) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
  data.table()

exp.id <- "2017-07-25-KS28"
poster.data.list.new[["Beta-Gamma-IRF-Nuclei"]] <- 
  poster.data.list.new[["Beta-Gamma-IRF-Nuclei"]] %>%
  rbind(
    poster.data.list[[exp.id]] %>%
      dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
      dplyr::mutate(experiment = exp.id)  %>%
      dplyr::mutate(time = time.2.1) %>%
      dplyr::mutate(priming = priming.1.1) %>%
      dplyr::mutate(stimulation = stimulation.1.1) %>%
      dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
      data.table()
  )

poster.data.list <- poster.data.list.new
saveRDS(object = poster.data.list, 
        file = paste(poster.path.list$input.dir, "data_ffc_joined.RDS", sep = "/"))

#### add zero ####
rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_joined.RDS"
poster.data.list <- readRDS(file = rds.path)
poster.data.list.new <- poster.data.list
for(poster.label in labels(poster.data.list)){
  #experiment.list <- unique(poster.data.list[[poster.label]]$experiment)
  #data.list.tmp <- list()
  #for( experiment_ in experiment.list){
    data.tmp <- poster.data.list[[poster.label]]# %>% dplyr::filter(experiment == experiment_)
    time.list <- unique(data.tmp$time)
    time.list <- time.list[which(time.list != 0)]
    data.zero <- data.tmp  %>% 
      dplyr::filter(time == 0 | stimulation == 0)
    for( t in time.list){
      poster.data.list.new[[poster.label]] <-
        poster.data.list.new[[poster.label]] %>%
        rbind(
          data.zero %>% 
            dplyr::mutate(stimulation = 0) %>%
            dplyr::mutate(time = t)
        )
    }
  }
#}

g <- plot_boxplot_group(
  data = poster.data.list.new$`Beta-Stat-Nuclei` %>% data.frame(),# %>% dplyr::filter(experiment =="2017-09-15-KZ78" ), 
  x = col_stimulation,#"well.name", 
  y = col_response, 
 # boxplot_group = col_time,#"well.name", 
  facet_grid_group_y = col_time,
  save_plot = FALSE
#  ylim_max_const = TRUE,
#  plot_title = poster.label
  #ylim_max = max(quantile(x = data[,col_response], na.rm = TRUE, probs = 0.95)[[1]], 1500)
  )


poster.data.list <- poster.data.list.new
saveRDS(object = poster.data.list.new, 
        file = paste(poster.path.list$input.dir, "data_ffc_zero_added.RDS", sep = "/"))


#### ####
rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_filtered.RDS"
poster.data.list <- readRDS(file = rds.path)
poster.data.labels <- labels(poster.data.list)
r <- foreach(poster.label.i = 1:length(poster.data.labels) %do% {
  
 poster.label <- poster.data.labels[poster.label.i]
 data <- poster.data.list[[poster.label]] %>% data.frame()
  
  if("time" %in% colnames(data)){
    col_time <- "time"
  } else if ("time.2.1" %in% colnames(data)){
    col_time <- "time.2.1"
  } else {
    col_time <- "time.1.1"
  }
  data[,"time"]  <- data[,col_time]
 
  if("stimulation" %in% colnames(data)){
    col_stimulation <- "stimulation"
  } else if ("stimulation.2.1" %in% colnames(data)){
    col_stimulation <- "stimulation.2.1"
  } else {
    col_stimulation <- "stimulation.1.1"
  }
  data[,"stimulation"]  <- data[,col_stimulation]
  
  if("Intensity_MeanIntensity_Alexa" %in% colnames(data)){
    col_response <- "Intensity_MeanIntensity_Alexa"
  } else if("Intensity_MeanIntensity_Alexa488" %in% colnames(data)){
    col_response <- "Intensity_MeanIntensity_Alexa488"
  } else {
    col_response <- "Intensity_MeanIntensity_Alexa555"
  }
  data[,"intensity"]  <- data[,col_response]
  
  poster.data.list.new[[poster.label]] <- data
}
  