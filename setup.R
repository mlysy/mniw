#--- compile the package ---------------------------------------------------

require(Rcpp)
require(RcppEigen)
require(devtools)
pkg.name <- "mniw"

compileAttributes()
document()
install()

test()

test_path()

build()

#--- initial compile -------------------------------------------------------

# default package skeleton
require(Rcpp)
# Rcpp.package.skeleton(name = "test", example_code = FALSE)
