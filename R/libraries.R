### ###
### libraries 
### ###

setwd("~/Documents/modelling/")

#### graphics ####
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)

#### dev ####
library(doParallel)
library(foreach)
library(Rcpp)
#install.packages("dtplyr")
#library(dtplyr)
library(dplyr)
#library(plyr)
library(data.table)

#### math ####
library(lhs)

#### optimisation ####
#install.packages("GenSA")
library(GenSA)
library(adagio)

#### SOURCES ####
source("R/jakstat_data.R")
source("R/jakstat_model.R")
source("R/lmvn_library.R")
source("R/likelihood_library.R")
source("R/jakstat_estimation.R")
source("R/theme_jetka.R")
source("R/jakstat_plot.R")
source("R/jakstat_data_library.R")