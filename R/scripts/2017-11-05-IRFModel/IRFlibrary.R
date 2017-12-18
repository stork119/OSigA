###
### IRF library 
###

library(e1071)
library(CapacityLogReg)
#### read input ####
poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/"

poster.path.list$output.dir <- "resources/output/poster/"

irfmodel.path.list <- list()
irfmodel.path.list$output.dir <- "resources/output/IRFmodel"
dir.create(irfmodel.path.list$output.dir)

irfmodel.path.list$rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_joined.RDS"
poster.data.list <- readRDS(file = irfmodel.path.list$rds.path)
