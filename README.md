## Usage instructions

### dependencies
System packages (Fedora): `R`
R packages: `Rcpp dplyr RColorBrewer ggplot2 tidyr gridExtra`

For the moment, copy Q1_binary.rds into /simulations (or add it to the repo if it's allowed to be here).

fits_xyz runs the simulation, arguments are:

`n p SNR num_bi num_bij perc_viol regression_alpha`

e.g. `./fits_xyz.R 1000 100 10 0 20 0 0.9`

Results will be saved to /fits_proper, to be read by the graph generating code.
No warning is given on incorrect usage, and unwanted results will be saved to fits_proper, so be sure to remove these if they don't look right.
This currently rebuilds core.cpp on every run. To avoid this you will probably need to install the local copy of xyz, and replace `sourceCPP('xyz/core.cpp')` with `library(xyz)`

To generate graphs, run `./PrecRecF1.R $regression_alpha`. Previous simulations will be read from results_proper, and graphs will be saved in /PrecRecF1. To allow experimentation with different values of regression_alpha, only results with the chosen value will be read.
