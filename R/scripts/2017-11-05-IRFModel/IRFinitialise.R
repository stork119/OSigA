###
### IRF library 
###
# setwd("../../../")
# source("R/optimisation/initialise_optimisation.R")
#source("R/scripts/2017-11-05-IRFModel/IRFlibrary.R")
source("R/scripts/2017-11-05-IRFModel/IRFlibraries.R")
source("R/graphics/libraries.R")
source("R/libraries/library_data_manipulation.R")

#### read input ####
poster.path.list <- list()
poster.path.list$input.dir <- "resources/input/poster/2017-12-18/rds/"

poster.path.list$output.dir <- "resources/output/poster/"

irfmodel.path.list <- list()
irfmodel.path.list$output.dir <- "resources/output/IRFmodel"
#dir.create(irfmodel.path.list$output.dir, recursive = TRUE)

irfmodel.path.list$rds.path <- "/home/knt/Documents/modelling/resources/input/poster/data_ffc_filtered.RDS"
poster.data.list <- readRDS(file = irfmodel.path.list$rds.path)

source("R/scripts/2017-11-05-IRFModel/IRFlibraries.R")
source("R/scripts/2017-11-05-IRFModel/IRFdata.R")
source("R/scripts/2017-11-05-IRFModel/IRFmodel_irf1.R")
source("R/scripts/2017-11-05-IRFModel/IRFmodel_ps1.R")
source("R/scripts/2017-11-05-IRFModel/IRFoptimise_fun.R")
source("R/scripts/2017-11-05-IRFModel/IRFsave.R")

setwd(dir = "knt/Documents/modelling/")
