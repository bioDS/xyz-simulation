### dependencies
System packages (Fedora): `R`
R packages: `Rcpp dplyr RColorBrewer ggplot2 tidyr gridExtra stringr`

For the moment, copy Q1_binary.rds into /simulations (or add it to the repo if it's allowed to be here).

### Usage
`./generate_data.sh $multiplier $count`

Will generate a $count sets of data to be graphed, with sizes (n,p, etc.) scaled by $multiplier. Graph scripts assume the multiplier is either 10 or 1. Output is stored in ./simulated_data.

`./generate_lethal_data.sh $multiplier $count`

Will do the same, but will generate data with lethal pairs. For the moment this is only used in one graph. Output is stored in ./simulated_data.

`./fit_all_1204.sh $threads`

Will check whether each generated dataset has already been (fitted/checked?) by xyz, and, if it hasn't, do so now. This will be done on `$threads` separate threads (2 by default).

`./generate_graphs.sh [large/small/store]`

Will (re)create every graph. 'large' will use n=10,000, p=1,000. 'small' will use n=1,000, p=100. 'store' will move every graph currently present into a subdirectory "pdfs_`date +%s`" (e.g. pdfs_1548386299). To find out when such a directory was made (in a more readable fashion), use `date --date=@1548386299` (adjusting the date as required).

N.B. attempts using different values of L (number of projections, for l_diff) need to be run manually (using `./fits_xyz_only $simulated_data/filename $L write`).

