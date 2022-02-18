# snpR 1.2.1

## UI
* Added a verbose option to many functions to avoid putting quite so much noise onto the console.

## Features

## Documentation
* Added a `NEWS.md` file to track changes to the package.
* Cleaned up the association testing vignette.


## Bug fixes
* Fixed a bug in run_random_forest() where formulas would be incorrectly specified the first time a line of code was run due to a wierd environment scope issue.