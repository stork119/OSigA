### ###
### libraries 
### ###


### ###
### libraries install
### ###
packages.list <- c(
  
  #### graphics ####
  "ggplot2",
  "corrplot",
  "Hmisc",
  "reshape2",
  "grid",
  "gridExtra",
  "ggthemes",
  
  #### optimisation algorithm ####
  "adagio",
  "cmaes",
  "lhs",
  
  #### ####
  "stringr",
  
  #### mathsstats ####
  # "MASS",
  # "GPfit",
  # "kernlab",
  # "FNN",
  # "Hotelling",
  # "MCMCpack", ### wishart distribution
  # "invgamma",
  # "akima",
  # "psych",
  # "mvtnorm",
  
  #### dev ####
  "doParallel",
  "parallel",
  "foreach",
  "Rcpp",
  
  #### data ####
  "data.table",
  "Matrix",
  
  #### data manipulation ####
  "dplyr",
  "lazyeval"
  
  #### channel capacity ####
  # "caret",
  # "glmnet",
  # "nnet",
  # "e1071"
)

packages.list.new <- packages.list[!(packages.list %in%
                                       installed.packages()[,"Package"])]
sapply(packages.list.new, remove.packages)
if(length(packages.list.new)){
  install.packages(packages.list.new)
}
sapply(packages.list, require, character.only = TRUE)
