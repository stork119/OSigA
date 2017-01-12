library(Rcpp)

#setwd("D:/KN/ClustSense/ClustDesign/R")


wd <- dirname(parent.frame(2)$ofile)

####TODO: Read from Makefile ####
  
  #### copmpiling static libraries
  source_cpp_filename <- "jakstatRcpp.cc"
  PKG_CXXFLAGS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib"
  PKG_LIBS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib  -lsundials_cvodes -lsundials_nvecserial"
  
  Sys.setenv("PKG_CXXFLAGS" = PKG_CXXFLAGS)
  Sys.setenv("PKG_LIBS" = PKG_LIBS)
  
  sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)
  
  PKG_CXXFLAGS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib"
  PKG_LIBS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -lsundials_cvodes -lsundials_nvecserial -Wl,-rpath,/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib "
  
  Sys.setenv("PKG_CXXFLAGS" = PKG_CXXFLAGS)
  Sys.setenv("PKG_LIBS" = PKG_LIBS)
  
  sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)
  
  par <- scan(file = "par.txt")
  var <- scan(file = "var.txt")
  
  rmain(par)