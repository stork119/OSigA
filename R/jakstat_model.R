### ###
### jakstat_model ###
### ###


PKG_CXXFLAGS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -mtune=core2 -std=c++11"
PKG_LIBS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib  -lsundials_cvodes -lsundials_nvecserial -Wl,-rpath,/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -pthread"

Sys.setenv("PKG_CXXFLAGS" = PKG_CXXFLAGS)
Sys.setenv("PKG_LIBS" = PKG_LIBS)

source_cpp_filename <- "C++/test/jakstat/jakstatrcpp_means.cc"
sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)

source_cpp_filename <- "C++/test/jakstat/jakstatrcpp_extrinsic.cc"
sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)
