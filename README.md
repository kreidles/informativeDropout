The informativeDropout Package
=========================

The informativeDropout package for R (>3.2.0) implements several models to account
for dropout in longitudinal studies.  

The package provides companion code for the manuscripts:

Moore, C.M., Carlson, N.E., MaWhinney, S. A Bayesian Natural Cubic B Spline Varying Coefficient Method
for Non-Ignorable Dropout, In review.

Moore, C.M., Carlson, N.E., MaWhinney, S. A Dirichlet Process Mixture Model for
Non-Ignorable Dropout, In review.

### Installation instructions

The results in the above manuscript were produced using R version 3.6.2. To reproduce the results,
perform the following steps:

* Install R version 3.6.x or higher by following the instructions at http://www.r-project.org
* From the R environment, install and load the "devtools" package
```R
install.packages("devtools")
library(devtools)
```
* Install the "informativeDropout" package directly from Github.com
```R
install_github(repo="kreidles/informativeDropout", ref="develop")
```
* Load the library
```R
library(informativeDropout)
```

### Instructions for replicating results in "A Dirichlet Process Mixture Model for Non-Ignorable Dropout"

Install the informativeDropout package as described above.

Run the DP model demo program

```R
demo(ba_dp_paper)
```
The demo will run for 200,000 iterations may take up to 24 hours to run, depending on processor speed. 

### Instructions for replicating results in "A Bayesian Natural Cubic B Spline Varying Coefficient Method for Non-Ignorable Dropout"

Install the informativeDropout package as described above.

Run the DP model demo program

```R
demo(BNSV_paper)
```




