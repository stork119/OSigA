### ###
### library_data_manipulation
### ###

#### normalise data ####
### function returns data with equal number of elements 
### in all subrgroups (priming, stimulation, time)
get_equal_data  <- function(
  data,
  sample_size = 1000,
  ...){
  
  d.distinct <- data %>%
    dplyr::distinct(priming,
             stimulation,
             time) 
  d.i <- d.distinct[1,]
  data.equal <- data %>% 
    filter(stimulation == d.i$stimulation,
           priming == d.i$priming,
           time == d.i$time) %>% 
    sample_n(sample_size, replace = TRUE)
  
  for(i in 2:nrow(d.distinct)){
    d.i <- d.distinct[i,]
    data.equal <- data.equal %>%
      rbind(data %>% 
              filter(stimulation == d.i$stimulation,
                     priming == d.i$priming,
                     time == d.i$time) %>% 
              sample_n(sample_size, replace = TRUE))
  }
  return(data.equal)
}

#### reading data ####
read_data <- function(
  path = "resources/input/data",
  files.list = list("KA10_C.csv",
                    "KA10_IFNB.csv"),
  file.header = TRUE,
  file.sep = "\t",
  columns.list = list(position    = "PositionName",
                      cell        = "ObjectNumber",
                      priming     = "group.1.1",
                      stimulation = "group.2.1",
                      time        = "compare.1.1",
                      intensity   = "ShrinkedNuclei.Intensity"
                      ),
  ...
){
  data.exp <- data.table(position = character(),
                         cell  = numeric(),
                         priming = numeric(),
                         stimulation = numeric(),
                         time = numeric(),
                         intensity = numeric(),
                         file = character())
  for(filename in files.list){
    path.file <- paste(path, filename, sep = "/")
    data.exp.part <- read.table(
      file = path.file,
      header = file.header,
      sep = file.sep)
    data.exp <- rbind(data.exp,
                      data.table(position    = data.exp.part[,columns.list$position],
                                 cell        = data.exp.part[,columns.list$cell],
                                 priming     = data.exp.part[,columns.list$priming],
                                 stimulation = data.exp.part[,columns.list$stimulation],
                                 time        = data.exp.part[,columns.list$time],
                                 intensity   = data.exp.part[,columns.list$intensity],
                                 file        = strsplit(filename, split = "[.]")[[1]][1]))
    
  }
  data.exp$logintensity <- log(data.exp$intensity)
  data.exp.unique <- distinct(data.exp, priming, stimulation, time)
  data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time)
  
  return( list(data.exp = data.exp#,
               #data.exp.grouped = data.exp.grouped,
               #data.exp.unique = data.exp.unique
               )) 
}

#### ####
prepare_data <- function(...){
  
}