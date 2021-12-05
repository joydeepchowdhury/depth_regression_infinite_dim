require(Rcpp)
Rcpp.package.skeleton("projmed", cpp_files = "projmed.cpp")
install.packages("projmed", repos=NULL, type="source")