# PhiMRF
R package for **P**oisson **Hi**erarchical **M**arkov **R**andom **F**ield model for analysis of spatial dependency on count data.

## Package Features

#### **Spatial dependency for count data**. 

Spatial neighborhoods are represented by networks, and thus not limited to regular lattices.

#### **Fast and parallel computation**. 

The MCMC process is written in C and takes advantage of BLAS and LAPACK routines. OpenMP enables parallel computation.

#### **Large matrices friendly**. 

The package has special functions that handles large matrices that would otherwise crash your R. 

## Model Description
See Vignettes.

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

Here is a brief outline:

   1. Find your R home directory by running `R.home()` within R.

   2. Go to your R home directory (`$R_HOME`) and see if the the directory `$R_HOME/src/nmath/standalone/` exists 

    If yes, run `make` in that directory.

    If not, [install R from source](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Top), and then go to the steps above.
    
Detailed installation instructions can be found [here](https://colinfay.me/r-installation-administration/the-standalone-rmath-library.html) for Linux, MacOS and Windows.



#### 3. Check if you have the required libraries
In a command-line console, type:
```
ldconfig -p | grep liblapacke
ldconfig -p | grep libRmath
ldconfig -p | grep liblapack
ldconfig -p | grep libblas
```

### Installation

Install the package in R using the [`devtools` package](https://cran.r-project.org/web/packages/devtools/index.html).

```
library(devtools)
devtools::install_github("https://github.com/ashleyzhou972/PhiMRF")
```



