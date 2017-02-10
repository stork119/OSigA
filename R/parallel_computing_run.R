###
### parallel_computing_run
###
setwd("~/Documents/modelling/")
source("R/parallel_computing.R")

path.optimisation <- paste(path.output, "cmaes/normalized/2017-02-04-5/", sep = "/")

data.exp.grouped <-  read.table(
  file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
  sep = ",",
  header = TRUE)


run_parallel_computations(path.optimisation = path.optimisation,
                          data.exp.grouped = data.exp.grouped,
                          no_cores = 18,
                          #maxit.tmp   =  10,
                          fun.optimisation = pureCMAES,
                          optimisation.res.par = "xmin")


#### ####
# d <- data.exp.grouped %>% filter(stimulation %in% stimulation.list)
# d.distinct <- d %>% distinct(priming, stimulation, time) 
# 
# d.i <- d.distinct[1,]
# data.exp.grouped.equal <- data.exp.grouped %>% 
#   filter(stimulation == d.i$stimulation, priming == d.i$priming, time == d.i$time) %>% 
#   sample_n(1000, replace = TRUE)
# 
# for(i in 2:nrow(d.distinct)){
#   d.i <- d.distinct[i,]
#   data.exp.grouped.equal <- data.exp.grouped.equal %>%
#     rbind(data.exp.grouped %>% 
#             filter(stimulation == d.i$stimulation, priming == d.i$priming, time == d.i$time) %>% 
#             sample_n(1000, replace = TRUE))
# }
# 
# write.table(x = data.exp.grouped.equal,
#             file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)
#### ####
