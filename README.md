# PhiMRF
R package for **P**oisson **Hi**erarchical **M**arkov **R**andom **F**ield model for analysis of spatial dependency on count data.

## Package Features
See vignetts


## Getting Started
### Prerequisites
#### 1. [Lapacke](https://www.netlib.org/lapack/lapacke.html)

Lapacke is a C interface to LAPACK routines.
To install on linux if you don't already have lapack and openblas:
```
sudo apt install liblapack3 
sudo apt install liblapack-dev 
sudo apt install libopenblas-base 
sudo apt install libopenblas-dev 
sudo apt install liblapacke-dev 
sudo apt install liblapack-dev
``` 
See [here](https://icl.cs.utk.edu/lapack-for-windows/lapack/#lapacke) for instructions for Windows.

#### 2. [Rmath Standalone library](https://cran.r-project.org/doc/manuals/R-exts.html#Standalone-Mathlib)

Rmath enables the use of common R functions such as `rnorm()` in C.

Detailed installation instructions can be found [here](https://colinfay.me/r-installation-administration/the-standalone-rmath-library.html) for linux, mac and Windows.

Here is how to build the Rmath library:

   1. Find your R home directory by running `R.home()` within R.

   2. Go to your R home directory (`$R_HOME`) and see if the the directory `$R_HOME/src/nmath/standalone/` exits 

    *If yes, run `make` in that directory.

    *If not, [install R from source](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Top), and then go to the steps above.


#### 3. Check if you have the required libraries

```
ldconfig -p | grep liblapacke
ldconfig -p | grep libRmath
ldconfig -p | grep liblapack
ldconfig -p | grep libblas
```
### Install



