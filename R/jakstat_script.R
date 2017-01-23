library(Rcpp)
library(data.table)
#install.packages("ggthemes")
library(ggplot2)
#install.packages("GenSA")
library(GenSA)

####TODO: Read from Makefile ####
  
  #### copmpiling static libraries
  setwd("~/Documents/modelling/")

  source_cpp_filename <- "C++/test/jakstat/jakstatRcpp.cc"
  
  PKG_CXXFLAGS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -mtune=core2 -std=c++11"
  PKG_LIBS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib  -lsundials_cvodes -lsundials_nvecserial -Wl,-rpath,/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -pthread"
  
  Sys.setenv("PKG_CXXFLAGS" = PKG_CXXFLAGS)
  Sys.setenv("PKG_LIBS" = PKG_LIBS)
  
  sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)
  
  path.parameters <- "resources/input/"
    
  par <- scan(file = paste(path.parameters, "par.txt", sep = ""))
  varscale <- 0.15
  variables <- rep(0.0, times = 629)
  variables[1:34] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
  variables[35:68] <- varscale*(variables[1:34]^2)
  tmesh <- seq(from = 0, to = 100, by = 5)
  
  tmesh.exp.i.list <- which(tmesh %in% tmesh.exp)
  data.model <- data.frame()
    
  data.model <- data.frame(t = numeric(),
                           mean = numeric(),
                           sd = numeric(),
                           cond = character(),
                           stm = numeric())
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
                          data.frame(t = c(tmesh[tmesh.exp.i], tmesh[tmesh.exp.i]),
                                     m = c(res$output[[tmesh.exp.i]][14],
                                           res$output[[tmesh.exp.i]][31]),
                                     sd =  c(res$output[[tmesh.exp.i]][48],
                                           res$output[[tmesh.exp.i]][65]),
                                     cond = c("C", "IFNB"), 
                                     stm = c(stm, stm))
                          )
      }
    } else {
      data.model <- rbind(data.model,
                          data.frame(t = c(tmesh[tmesh.exp.i], tmesh[tmesh.exp.i]),
                                     m  = 0,
                                     sd = 0,
                                     cond = c("C", "IFNB"), 
                                     stm = c(stm, stm))
      )
    }
  }
  
  data.model <- normalization(data.model, background = data.exp.background)
  data.model <- 
    
     
  ?GenSA
  
  
  
  data <- data.table(y = numeric(), t = numeric(), val = numeric(), cond = numeric())
  for(i in 1:length(res)){
    data.tmp <- data.table(y = 1:34,
                           t = rep(tmesh[i], 34),
                           val = res[[i]][1:34],
                           cond = rep(cond.tmp, 34))
    data <- rbind(data, data.tmp)
  }
  
  ggplot(data.frame(data[data$y == 13,]), aes(x = t, y= val)) + geom_line() 