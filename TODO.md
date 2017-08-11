1.  *Documentation*: same template as **`LMN`**, except we may not need a vignette since the package is somewhat self-explanatory.  Perhaps a short vignette illustrating how much longer it is to sample from an MNIW distribution using other packages (e.g., **`MCMCpack`**), and perhaps the Gibbs sampler for the hierarchical model as well.

2.  *Convert **C++** code as Header-only library*.  To do this, put everything into header files *except* files ending in `Export.cpp`.  Those are the interface functions between **C++** and **R**.  Don't edit `RcppExports.cpp` since that's created by the **`Rcpp`** package, which itself facilitates the interface process.

3.  *Links of interest*. 
    + Hadley Wickham's free e-book on [package authoring](http://r-pkgs.had.co.nz/).
    + The [**`devtools`**](https://cran.r-project.org/web/packages/devtools/index.html) package which goes hand-in-hand with the book.
    + The [**`Rcpp`**](https://cran.r-project.org/web/packages/Rcpp/index.html) package for interfacing **C++**/**R**.  In particular, if you change the exported **C++** functions (which are preceded by `//[[Rcpp::export]]` in the `.cpp` source code), then you should run `Rcpp::compileAttributes()` with the working directory set to the package folder in order to update the `src/RcppExports.cpp` and `R/RcppExports.R` files.
