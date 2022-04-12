# snpR 1.2.1

## UI
* Added a verbose option to many functions to avoid putting quite so much noise onto the console.

## Features
### Major
### Minor
* Renamed calc_SFS and calc_pairwise_LD_heatmap to lowercase for consistency.
* Got rid of sfs arguments in plot_sfs and calc_direcitonality. They now just take provided sfs objects as x, consistent with other functions.
* Reworked calc_pairwise_fst bootstrapping to be more memory efficient.
* Added a chr_order argument to plot_manhattan() to allow for manual resorting of chromosomes (since factors are coerced away in snpRdata objects).

## Documentation
* Added a `NEWS.md` file to track changes to the package.
* Cleaned up the association testing vignette.


## Bug fixes
* Fixed a bug in run_random_forest() where formulas would be incorrectly specified the first time a line of code was run due to a weird environment scope issue.
* Fixed a bug in get.snpR.stats() when requesting fst values from both a facet with and without fst calculated would throw an error during fst matrix construction. Implemented a test.
* Added ggtree citation to plot_tree.
* Fixed a bug where running no formula but generating importance estimates would throw an error in run_random_forest()