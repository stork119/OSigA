### ### 
### script - preparing data 
### ###

source("R/optimisation/initialise_optimisation.R")

poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"

poster.path.list$output.dir <- "resources/output/poster/"
dir.create(poster.path.list$output.dir)

#### reading data ####
poster.data.list <- list()
poster.path.list$experiments <- 
  list(
    list(id = "2016-01-26-KA10-C"),
    list(id = "2016-01-26-KA10-IFNB"),
    list(id = "2016-01-26-KA10-IFNG"),
    list(id = "2016-01-28-KA11-C"),
    list(id = "2016-01-28-KA11-IFNB"),
    list(id = "2016-01-28-KA11-IFNG"),
    list(id = "2017-07-18-KS25", path = "2017-07-18-KS25/488"),
    list(id = "2017-07-18-KS26", path = "2017-07-18-KS26/488"),
    list(id = "2017-07-25-KS27", path = "2017-07-25-KS27/488"),
    list(id = "2017-07-25-KS28", path = "2017-07-25-KS28/488"),
    list(id = "2017-07-19-KZ63", path = "2017-07-19-KZ63/488"),
    list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", output_id = "2017-08-25-KZ73-Nuclei"),
    list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", output_id = "2017-09-15-KZ78-Nuclei"),
    list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", file = "Cells555.csv", output_id = "2017-08-25-KZ73-Cells"),
    list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", file = "Cells555.csv", output_id = "2017-09-15-KZ78-Cells"),
    list(id = "2017-08-25-KZ73", path = "2017-08-25-KZ73/555", file = "Cytoplasm555.csv", output_id = "2017-08-25-KZ73-Cytoplasm"),
    list(id = "2017-09-15-KZ78", path = "2017-09-15-KZ78/555", file = "Cytoplasm555.csv", output_id = "2017-09-15-KZ78-Cytoplasm")
    )


for(experiment in poster.path.list$experiments){
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
  
  if(is.null(experiment[["output_id"]])){
    experiment$output_id <- experiment$id
  }
  poster.data.list[[experiment$output_id]] <-
    read.table(file = path, 
               header = TRUE, 
               sep = ",") %>% 
    normalize_data() %>%
    data.table()  
}


#### daving RDS ####
saveRDS(object = poster.data.list, 
        file = paste(poster.path.list$input.dir, "data_ffc.RDS", sep = "/"))


#### filter data ####
rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc.RDS"
poster.data.list <- readRDS(file = rds.path)

poster.filtered.list <- list(
  list(id = "2016-01-26-KA10-C", well.names = c('A01', 'B01', 'H12')),
  list(id = "2016-01-26-KA10-IFNB" , well.names = c('A02', 'A03', 'A04')),
  list(id = "2016-01-26-KA10-IFNG", well.names = c('D02', 'H12')),
  list(id = '2016-01-28-KA11-C', well.names = c(paste('A0', 1:9, sep = ''), 'A10', 'A11', 'A12', 'B01', 'B11', 'B12')),
  list(id = '2016-01-28-KA11-IFNB', well.names = c('B01', 'G01', 'H01')),
  list(id = '2016-01-28-KA11-IFNG', well.names = c(paste('A0', 1:3, sep = ''), 'B01', 'C01')),
  list(id = '2017-07-18-KS26', well.names = c('D01')),
  list(id = '2017-07-25-KS27', well.names = c('D02','F01', 'F02', 'F10', 'H04')),
  list(id = '2017-07-25-KS28', well.names = c('D02', 'D03')),
  # list(id = '2017-08-25-KZ73-Nuclei', well.names = c('A11', paste(LETTERS[1:8], '01', sep = ""))),
  # list(id = '2017-08-25-KZ73-Cytoplasm', well.names = c(paste('A0', 1:4, sep = ''), 'A11', 'A12', paste('B0', 1:3, sep = ''), 'H01', 'H12',"well.name")),
  # list(id = '2017-08-25-KZ73-Cells', well.names = c(paste('A0', 1:4, sep = ''), 'A11', 'A12', paste('B0', 1:3, sep = ''), 'H01', 'H12', "well.name"),
  list(id = '2017-08-25-KZ73-Nuclei', well.names = c(as.vector(matrix(sapply(LETTERS[c(1,3:8)], function(l) {paste(l, c(paste('0', 3:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), "A01", "A02")),
  list(id = '2017-08-25-KZ73-Cells', well.names = c(as.vector(matrix(sapply(LETTERS[c(1,3:8)], function(l) {paste(l, c(paste('0', 3:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), "A01", "A02")),
  list(id = '2017-08-25-KZ73-Cytoplasm', well.names = c(as.vector(matrix(sapply(LETTERS[c(1,3:8)], function(l) {paste(l, c(paste('0', 3:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), "A01", "A02")),
  list(id = '2017-09-15-KZ78-Nuclei', well.names = c(as.vector(matrix(sapply(LETTERS[c(2)], function(l) {paste(l, c(paste('0', 1:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), 
                                                     as.vector(matrix(sapply(LETTERS[1:8], function(l) {paste(l, c("01", "02"), sep = "")}), ncol = 1)))),
  list(id = '2017-09-15-KZ78-Cells', well.names = c(as.vector(matrix(sapply(LETTERS[c(2)], function(l) {paste(l, c(paste('0', 1:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), 
                                                    as.vector(matrix(sapply(LETTERS[1:8], function(l) {paste(l, c("01", "02"), sep = "")}), ncol = 1)))),
  list(id = '2017-09-15-KZ78-Cytoplasm', well.names = c(as.vector(matrix(sapply(LETTERS[c(2)], function(l) {paste(l, c(paste('0', 1:9, sep = ""), "10", "11", "12"), sep = "")}), ncol = 1)), 
                                                        as.vector(matrix(sapply(LETTERS[1:8], function(l) {paste(l, c("01", "02"), sep = "")}), ncol = 1))))
  
  )

poster.data.list.2 <- poster.data.list
for(poster.filter in poster.filtered.list){
  #poster.filter <- poster.filtered.list[[1]]
  poster.data.list.2[[poster.filter$id]] <-
    poster.data.list[[poster.filter$id]][
      -which(poster.data.list[[poster.filter$id]]$well.name %in% poster.filter$well.names),]
}
poster.data.list <- poster.data.list.2
saveRDS(object = poster.data.list, 
        file = paste(poster.path.list$input.dir, "data_ffc_filtered.RDS", sep = "/"))


#### connect data ####
rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_filtered.RDS"
poster.data.list <- readRDS(file = rds.path)

poster.data.list.new <- list()
# BetaStat1
part.list <- list("Nuclei", "Cytoplasm", "Cells")
for(part in part.list){
  exp.id <- "2017-09-15-KZ78"
  data.KZ78 <- poster.data.list[[paste(exp.id, part, sep = "-")]] %>% 
    dplyr::mutate(Intensity_MeanIntensity_Alexa = Intensity_MeanIntensity_Alexa555) %>%
    dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
    dplyr::mutate(experiment = exp.id) %>%
    dplyr::mutate(time = time.1.1) %>%
    dplyr::mutate(stimulation = stimulation.1.1) %>%
    dplyr::select(experiment, well.name, time, stimulation, Intensity_MeanIntensity_Alexa) %>%
    data.table()
    
  exp.id <- "2017-08-25-KZ73"
  data.KZ73 <- poster.data.list[[paste(exp.id, part, sep = "-")]] %>% 
    dplyr::mutate(Intensity_MeanIntensity_Alexa = Intensity_MeanIntensity_Alexa555) %>%
    dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
    dplyr::mutate(experiment = exp.id) %>%
    dplyr::mutate(time = time.1.1) %>%
    dplyr::mutate(stimulation = stimulation.1.1) %>%
    dplyr::select(experiment, well.name, time, stimulation, Intensity_MeanIntensity_Alexa)%>%
    data.table()
  
  poster.data.list.new[[paste("Beta-Stat", part, sep = "-")]] <- 
  rbind(data.KZ73, data.KZ78)
}

# Beta-pStat
exp.id <- "2017-07-19-KZ63"
poster.data.list.new[["Beta-pStat"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(priming = 0)  %>%
  dplyr::mutate(time = time.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa)%>%
  data.table()
  
# Beta-gamma

exp.id <- "2016-01-26-KA10-C"
poster.data.list.new[["Gamma-pStat"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(priming = 0)  %>%
  dplyr::mutate(time = time.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa)%>%
  data.table()

exp.id <- "2016-01-28-KA11-C"
poster.data.list.new[["Gamma-pStat"]] <- 
  poster.data.list.new[["Gamma-pStat"]] %>%
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
poster.data.list.new[["Gamma-Gamma-pStat"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(time = time.1.1) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa)%>%
  data.table()

exp.id <- "2016-01-28-KA11-IFNG"
poster.data.list.new[["Gamma-Gamma-pStat"]] <- 
  poster.data.list.new[["Gamma-Gamma-pStat"]] %>%
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
poster.data.list.new[["Beta-Gamma-pStat"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(time = time.1.1) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
  data.table()

exp.id <- "2016-01-28-KA11-IFNB"
poster.data.list.new[["Beta-Gamma-pStat"]] <- 
  poster.data.list.new[["Beta-Gamma-pStat"]] %>%
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
poster.data.list.new[["Beta-IRF"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(time = time.2.1) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
  data.table()

# Beta-Gamma-IRF
exp.id <- "2017-07-18-KS26"
poster.data.list.new[["Beta-Gamma-IRF"]] <- 
  poster.data.list[[exp.id]] %>%
  dplyr::mutate(well.name = paste(exp.id, well.name, sep = "-")) %>%
  dplyr::mutate(experiment = exp.id)  %>%
  dplyr::mutate(time = time.2.1) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::select(experiment, well.name, time, priming, stimulation, Intensity_MeanIntensity_Alexa) %>%
  data.table()

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
