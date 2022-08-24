## Test environments

* macOS-latest (release) -- via github actions
* windows-latest (release) -- via github actions
* ubuntu-latest (oldrel-1) -- via github actions
* ubuntu-latest (release) -- via github actions

## R CMD check results

0 errors | 0 warnings | 2 notes
```
N  checking CRAN incoming feasibility
   Maintainer: 'William Hemstrom <hemstrow@gmail.com>'
   
   New submission
   
   Possibly misspelled words in DESCRIPTION:
     Nucleotide (3:56, 12:45)
     Polymorphism (4:9, 12:56)
```

* This is a new release.
* Nucleotide and polymorphism are both correctly spelled genetics terminology essential to describing the package.

```
N  checking for detritus in the temp directory
   Found the following files/directories:
     'lastMiKTeXException'
```
* See [R-hub issue #503](https://github.com/r-hub/rhub/issues/503). This is probably due to bug/crash in MiKTeX and can be ignored.

## Downstream Dependencies

There are currently no downstream dependencies for this package.
