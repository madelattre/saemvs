
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saemvs

<!-- badges: start -->

<!-- badges: end -->

The goal of saemvs is to …

## Installation

This package contains C++ code using **Rcpp** and **RcppArmadillo**,
which requires a working C++ and Fortran toolchain. Please follow the
instructions for your operating system:

### macOS

- R ≥ 4.5 is required.
- gfortran ≥ 14.x must be installed for compiling C++ code.  
  Download the official gfortran 14.x installer here:  
  <https://mac.r-project.org/tools/>  
- Verify your installation:

``` bash
gfortran --version
```

### Windows

`Rtools` must be installed for compiling packages with C++ and Fortran
code. Make sure to use the version of `Rtools` corresponding to your R
version: <https://cran.r-project.org/bin/windows/Rtools/>

### Linux

Ensure gfortran and a C++ compiler are installed via your package
manager. Example for Debian/Ubuntu:

``` bash
sudo apt update
sudo apt install gfortran build-essential
```

### Installing the package

After ensuring your system meets the requirements:

``` r
# Install dependencies first
install.packages(c("Rcpp", "RcppArmadillo"), dependencies = TRUE)

# Then install the package (from a local .tar.gz or source)
remotes::install_local("path/to/your_package.tar.gz", dependencies = TRUE)

## Example

This is a basic example which shows you how to solve a common problem:


``` r
library(saemvs)
## basic example code
```
