#--- compile the package ---------------------------------------------------

require(Rcpp)
require(RcppEigen)
require(devtools)
pkg.path <- "C:/Users/Jerome/Documents/R/mniw/mniw"
#pkg.name <- "mniw"

compileAttributes(pkgdir = pkg.path)
document(pkg = pkg.path)
install(pkg = pkg.path)

test()

test_path()

build()

#--- initial compile -------------------------------------------------------

# default package skeleton
require(Rcpp)
# Rcpp.package.skeleton(name = "test", example_code = FALSE)
