### ###
### jakstat_model ###
### ###


if( Sys.info()[["sysname"]] == "Linux"){
  if(Sys.info()[["nodename"]] == "piotrk-pc"){
    PKG_CXXFLAGS="-I/usr/share/sundials/sundials-2.7.0/instdir/include -I/usr/share/sundials/sundials-2.7.0/instdir/include/cvode -I/usr/share/sundials/sundials-2.7.0/instdir/include/cvodes -L/usr/share/sundials/sundials-2.7.0/instdir/lib -mtune=core2 -std=c++11      -Wl,-rpath,/usr/share/sundials/sundials-2.7.0/instdir/lib"
    PKG_LIBS="-I/usr/share/sundials/sundials-2.7.0/instdir/include -L/usr/share/sundials/sundials-2.7.0/instdir/lib  -lsundials_cvodes -lsundials_nvecserial -Wl,-rpath,/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -pthread      -Wl,-rpath,/usr/share/sundials/sundials-2.7.0/instdir/lib"
  } else {
    PKG_CXXFLAGS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -mtune=core2 -std=c++11"
    PKG_LIBS="-I/home/knt/Programs/Sundials/cvodes-2.9.0-inst/include -L/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib  -lsundials_cvodes -lsundials_nvecserial -Wl,-rpath,/home/knt/Programs/Sundials/cvodes-2.9.0-inst/lib -pthread"
  }
} else if(Sys.info()[["sysname"]] == "Windows" ){
  print("NOT WORKING !!!")
  break()
} else {
  PKG_CXXFLAGS="-I/usr/local/include -I/usr/local/include/cvode -I/usr/local/include/cvodes -L/usr/local/lib  -mtune=core2 -std=c++11"
  PKG_LIBS="-I/usr/local/include -I/usr/local/include/cvode -I/usr/local/include/cvodes -L/usr/local/lib  -lsundials_cvodes -lsundials_nvecserial -Wl,-rpath,/usr/local/lib -pthread"
}

Sys.setenv("PKG_CXXFLAGS" = PKG_CXXFLAGS)
Sys.setenv("PKG_LIBS" = PKG_LIBS)

source_cpp_filename <- "C++/test/jakstat/jakstatrcpp_means.cc"
sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)

source_cpp_filename <- "C++/test/jakstat/jakstatrcpp_extrinsic.cc"
sourceCpp(file = source_cpp_filename, rebuild = TRUE, showOutput = TRUE)
