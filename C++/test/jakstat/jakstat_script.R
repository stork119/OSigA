library(Rcpp)
library(data.table)
#install.packages("ggthemes")
library(ggplot2)


####TODO: Read from Makefile ####
  
  #### copmpiling static libraries
  source_cpp_filename <- "jakstatRcpp.cc"
  
  PKG_CXXFLAGS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -mtune=core2 -std=c++11"
  PKG_LIBS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib  -lsundials_cvodes -lsundials_nvecserial -Wl,-rpath,/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -pthread"
  
  Sys.setenv("PKG_CXXFLAGS" = PKG_CXXFLAGS)
  Sys.setenv("PKG_LIBS" = PKG_LIBS)
  
  sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)
  
  par <- scan(file = "par.txt")
  variables <- rep(0.0, times = 170)
  variables[1:17] <- scan(file = "var.txt")
  cond.tmp <- 1
  tmesh <- seq(from = 0, to = 100, by = 5)
  res <- rmain(parameters = par, variables = variables)
  data <- data.table(y = numeric(), t = numeric(), val = numeric(), cond = numeric())
  for(i in 1:length(res)){
    data.tmp <- data.table(y = 1:34,
                           t = rep(tmesh[i], 34),
                           val = res[[i]][1:34],
                           cond = rep(cond.tmp, 34))
    data <- rbind(data, data.tmp)
  }
  
  ggplot(data.frame(data[data$y == 13,]), aes(x = t, y= val)) + geom_line() 