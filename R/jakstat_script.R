library(Rcpp)
#install.packages("dtplyr")
library(dtplyr)
#install.packages("ggthemes")
library(ggplot2)
#install.packages("GenSA")
library(GenSA)

####TODO: Read from Makefile ####
  
  #### copmpiling static libraries ####
  setwd("~/Documents/modelling/")
  source("R/jakstat_data.R")
  source("R/jakstat_estimation.R")

#### ####
  source_cpp_filename <- "C++/test/jakstat/jakstatRcpp.cc"
  
  PKG_CXXFLAGS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -mtune=core2 -std=c++11"
  PKG_LIBS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib  -lsundials_cvodes -lsundials_nvecserial -Wl,-rpath,/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -pthread"
  
  Sys.setenv("PKG_CXXFLAGS" = PKG_CXXFLAGS)
  Sys.setenv("PKG_LIBS" = PKG_LIBS)
  
  sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)
  
  path.parameters <- "resources/input/"
    
  par <- scan(file = paste(path.parameters, "par.txt", sep = ""))
  #### ####
likelihood <- function(par, variables, tmesh, tmesh.exp.i.list){
  
    data.model <- data.table(time = numeric(),
                             m = numeric(),
                             sd = numeric(),
                             priming = numeric(), 
                             stimulation = numeric())
    print(par)
    for(stm in unique(data.exp$stimulation)){
      res <- rmain(parameters = par, 
                   variables = variables, 
                   stm = stm, 
                   tmesh = tmesh, 
                   time_interval = 100, 
                   time_computation = 1000*60*5)
      if(res$success){
        for(tmesh.exp.i in tmesh.exp.i.list){
          data.model <- rbind(data.model,
                            data.table(time = c(tmesh[tmesh.exp.i], tmesh[tmesh.exp.i]),
                                       m = c(res$output[[tmesh.exp.i]][14],
                                             res$output[[tmesh.exp.i]][31]),
                                       sd =  c(res$output[[tmesh.exp.i]][48],
                                             res$output[[tmesh.exp.i]][65]),
                                       priming = c(0, 1000), 
                                       stimulation = c(stm, stm))
                            )
        }
      } else {
        return(Inf)
      }
    }
    
    data.model <- normalization(data.model, background = data.exp.background)
    data.model <- lmvn(data.model)
    
    result <- sum(sapply(13:nrow(data.model),
    function(data.model.i){
      data.model.tmp <- data.model[data.model.i,]
      return((data.exp.grouped %>%
           filter(priming == data.model.tmp$priming,
                  time == data.model.tmp$time,
                  stimulation == data.model.tmp$stimulation) %>%
           mutate(likelihood = fun.likelihood(logintensity, intensity, data.model.tmp = data.model.tmp)) %>%
           summarise(likelihood.sum = sum(likelihood)))$likelihood.sum)
    }
    ))
    print(result)
    
    return(result)
}

#### ####        
varscale <- 0.15
variables <- rep(0.0, times = 629)
variables[1:34] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
variables[35:68] <- varscale*(variables[1:34]^2)
tmesh <- seq(from = 0, to = 100, by = 5)
tmesh.exp.i.list <- which(tmesh %in% tmesh.exp)

#### optimisation ####

par.lower <- par*10^(-1)
par.upper <- par*10^(1)

#likelihood(par, variables = variables, tmesh = tmesh, tmesh.exp.i.list = tmesh.exp.i.list )    

genSA.res <- GenSA(par = par,
                   fn = likelihood,
                   lower = par.lower,
                   upper = par.upper,
                   control = list(maxit = 10, max.call = 10, verbose = TRUE), 
                   variables = variables, 
                   tmesh = tmesh,
                   tmesh.exp.i.list = tmesh.exp.i.list)    



  
