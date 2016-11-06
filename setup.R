#--- compile the package --------------------------------------------------------

# default package skeleton
require(Rcpp)
# Rcpp.package.skeleton(name = "test", example_code = FALSE)

require(RcppEigen)
require(devtools)
pkg.name <- "mniw"

compileAttributes()
document()
install()

build(pkg.name)

test(pkg.name)
