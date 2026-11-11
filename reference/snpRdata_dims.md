# Get the dimensions of a snpRdata object

Get the dimensions of a snpRdata object

## Usage

``` r
# S4 method for class 'snpRdata'
nsnps(x)

# S4 method for class 'snpRdata'
nrow(x)

# S4 method for class 'snpRdata'
nsamps(x)

# S4 method for class 'snpRdata'
ncol(x)

# S4 method for class 'snpRdata'
dim(x)
```

## Arguments

- x:

  snpRdata object

## Functions

- `nsnps(snpRdata)`: get the number of SNPs

- `nrow(snpRdata)`: identical to nsnps

- `nsamps(snpRdata)`: get the number of samples

- `ncol(snpRdata)`: identical to nsamps

- `dim(snpRdata)`: get the number of SNPs and samples
