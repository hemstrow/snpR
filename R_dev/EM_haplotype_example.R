# read in the example file:
ex <- readRDS("R_dev/haplotype_frequency_estimation_example.RDS")

# this contains counts of all of the "phenotypes", or the visible genotypes of two different haplotypes
ex <- ex[,-grep("NN", colnames(ex))]
x <- ex[1,]
phenos <- which(x != 0)
x <- x[phenos] # the cleaned version is the data for one row, only the cells that actually have data.

# haptable, containing the tabulated unambigious haplotypes
hap.ex <- readRDS("R_dev/hapmat_example.RDS")
haptable <- hap.ex[1,] # these are the counts from the unambigious loci
haps <- c("AC", "AG", "GC", "GG")
names(haptable) <- haps


# so here's an example funciton. Right now, it only works for the first row of ex!
single_haplotype_estimation <- function(x, haptable, sigma = 0.0001){

  # find the double het. Should be able to use an approach like this when this gets extended to work with everything.
  # cj values for each possible genotype:
  s1 <- substr(names(x), 1, 1)
  s2 <- substr(names(x), 2, 2)
  s3 <- substr(names(x), 3, 3)
  s4 <- substr(names(x), 4, 4)
  het.1 <- s1 != s2
  het.2 <- s3 != s4



  # First, make a guess at the starting haplotype frequencies. We'll do this by taking the unambigious haplotype frequencies,
  # then making a guess at the haplotype composition in the double heterozygote assuming that all possible haplotypes are equally likely
  doub.het <- which(het.1 + het.2 == 2) # identify double heterozygotes

  nhap.counts <- haptable # grab the haplotypes
  ehap.counts <- nhap.counts + .5*x[doub.het] # assuming that both haplopairs are equaly likely in the double het
  shap.freqs <- ehap.counts/sum(ehap.counts) # get the starting haplotype frequencies



  # now that we have our starting conditions, we will do the EM loops.
  # 1) First, we find out how many of each haplotype
  # we expect to get from our double heterozygotes given the initial haplotype frequencies we guessed above.
  # 2) Then, we use those expected frequencies to update our estimates of the haplotype frequencies.
  # 3) repeat 1 and 2 until the difference between the haplotype frequencies between loop iterations is less than sigma.


  # we'll use a while loop, which will run as long as the conditions are met. Note that this can freeze your computer if
  # the conditions are NEVER met! Try just hitting the stop sign, if that doesn't work you'll need to restart Rstudio.


  diff <- sigma + 1 # initialize the diff. Doesn't matter what it is as long as it's larger than sigma.
  while(diff > sigma){

    # 1)
    # expectation, which is that we are drawing haplotypes (aka alleles) from a pool of options. Follows HWE, essentially,
    # but the "alleles" are actually haplotypes
    op1.e <- (2*shap.freqs["AC"]*shap.freqs["GG"])/
      ((2*shap.freqs["AC"]*shap.freqs["GG"])+(2*shap.freqs["AG"]*shap.freqs["GC"])) # percentage of AC/GG haplo pairs
    op2.e <- 1 - op1.e

    # maximization: given the expected haplotype frequencies, how many of each haplotype should we have? get new frequencies
    n1hap.freqs <- haptable # grab the known haplotype frequencies form the unambigious phenotypes again.
    n1hap.freqs[c("AC", "GG")] <- n1hap.freqs[c("AC", "GG")] + (x[4]*op1.e*.5) # we basically add the expected number of haplotypes for the double heterozygotes
    n1hap.freqs[c("AG", "GC")] <- n1hap.freqs[c("AG", "GC")] + (x[4]*op2.e*.5)
    n1hap.freqs <- n1hap.freqs/sum(n1hap.freqs)

    # calculate the diff and update
    diff <- sum(abs(n1hap.freqs - shap.freqs))

    shap.freqs <- n1hap.freqs
  }

  # return the output
  return(shap.freqs)
}

# run
out <- single_haplotype_estimation(x, haptable)



# The goal now is to extend this so that the function runs when provided with the full ex and hap.ex files.
# Note that each row will stop at a different point, since the difference in haplotype frequencies between
# iterations will reach sigma at different numbers of iterations. What we can do is either keep going until
# everything hits sigma, or just stop each row when they cross sigma and only iterate the remainder. The latter
# should be faster.

# Paper source for EM:
# Maximum-likelihood estimation of molecular haplotype frequencies in a diploid population. Excoffier, L., and Slatkin, M. (1995)
# The math in this paper IS NOT HELPFUL in general. Here's maybe a better source:
# https://homes.cs.washington.edu/~suinlee/genome541/lecture3-genetics-Lee.pdf

# Note: the haplotype table is provided by the tabulate_haplotypes function in the calc_pairwise_ld function in the package.
# Go to the R/stat_functions.R file, then ctrl + F "tabulate_haplotypes" to find it. That function calls "GtoH", a function that is defined just above
# it in the script!
