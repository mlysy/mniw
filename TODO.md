- [x]  *Test hierarchical model*.  Make sure the **R** wrapper does what it's supposed to and that the underlying **C++** code and testing **Stan** code give the same output. Put all `.stan` and `.R` test files into `tests/dontrun`.

- [ ] *Documentation*: same template as **`LMN`**, except we may not need a vignette since the package is somewhat self-explanatory.  Perhaps a short vignette illustrating how much longer it is to sample from an MNIW distribution using other packages (e.g., **`MCMCpack`**), and perhaps the Gibbs sampler for the hierarchical model as well.

- [x]  *Convert **C++** code as Header-only library*.  To do this, put everything into header files *except* files ending in `Export.cpp`.  Those are the interface functions between **C++** and **R**.  Don't edit `RcppExports.cpp` since that's created by the **`Rcpp`** package, which itself facilitates the interface process.

- [x]  Standardize and coordinate argument names in R and C++, including unit tests.

- [ ]  Create package vignette illustrating package functionality.  Might also be worth comparing to [**mvnfast**](https://mfasiolo.github.io/mvnfast/articles/mvnfast.html).

- [x] Convenience function `.get_PQ` names dimensions.  To fix this change its output to unamed arguments, but then must make sure they are never referred to by name.

- [ ] Clean up really old commented code.

- [ ] Remove default values for distribution parameters.  This probably does more harm than good, e.g., if users define but forget to pass in a particular input...

- [ ] Add example for `MatrixT`.

- [ ] Add multivariate-t distribution for good measure.

- [ ] Finish documentation and example for Gibbs sampler.

- [ ] Add vectorized matrix inversion to `crossprodV`.

- [ ] Add `Omega` precision input to `rMNIW`?  How about `rmniw`?
