armitage_example
# task: conduct an armitage test for every SNP (row) in the example dataset.
# get example data
readRDS("R_dev/armitage_example.RDS")
head(readRDS)


# here's a function I wrote to do this for a single SNP, taking data in a 2x3 matrix:
# see Armitage (1955), Tests for Linear Trends in Proportions and Frequencies, Biometrics for equations

# table with cases on top, controls on bottom, then genotypes pp, pq, and qq in columns 1, 2, 3
m <- matrix(c(19, 29, 24, 497, 560, 269), 2, byrow = T) # example table, you probably won't need to remake this.
colnames(m) <- c("pp", "pq", "qq")
rownames(m) <- c("case", "control")
wv <- c(0,1,2) # these are weights. These should be 0 for all "pp", 1 for all "pq", and 2 for all "qq" genotypes
wv

# the function
calc_single_amitage <- function(m, w){

  # define variables, names relate to Armitage 1955
  n <- m[1,] # cases
  N <- colSums(m) # Number of individuals with each genotype
  bT <- sum(m) # total number of individuals sequenced at this loci
  t <- sum(m[1,]) # total number of "case" individuals

  # get the sums that go into the equations, see the paper
  s1 <- sum(n*wv)
  s2 <- sum(N*wv)
  s3 <- sum(N*wv^2)

  # equation 5
  b <- (bT*s1 - t*s2)/(bT*s3 - (s2^2))

  # equation 6
  Vb <- (t*(bT - t))/(bT*(bT*s3 - s2^2))

  # equation 7
  chi <- (b^2)/Vb

  # use the pchisq function with 1 degree of freedom to get a p-value and return.
  return(pchisq(chi, 1, lower.tail = F))
}

# the output, returns a p-value
calc_single_amitage(m, wv)




# you'll need to expand this function to work on the whole cast_gs object!
# I haven't done this yet, but my tips are:
#    1) split the data into two sets, one with cases and one with controls
#    2) figure out weights for each column (pp = 0, pq = 1, qq = 2)
#    3) do the calculation. Remember that any column with a 0 in it won't effect the resulting p-values at all,
#       so there's no reason to remove them or worry about them for each row!
head(armitage_example)
a<- armitage_example
w<- wv
calc_armitage<- function(a, w){#where a is the matrix you want to run the test on, and w is the weights
  #separate controls
  control<- matrix(nrow = nrow(a), ncol=10)
  control<- a[, c("AA_control", "AC_control", "AG_control", "AT_control", "CC_control", "CG_control"
              , "CT_control", "GG_control", "GT_control", "TT_control")]
  #separate cases
  case<- matrix(nrow = nrow(a), ncol=10)
  case<- a[, c("AA_case", "AC_case", "AG_case", "AT_case", "CC_case", "CG_case"
                  , "CT_case", "GG_case", "GT_case", "TT_case")]
  #identify ps qs and pqs
  p_cont<- matrixStats::rowMaxs(control)
  p_case<- matrixStats::rowMaxs(case)
  q_cont<- matrixStats::rowMins(control)
  q_case<- matrixStats::rowMins(case)
  pq_cont<- matrixStats::rowMedians(control)
  pq_case<- matrixStats::rowMedians(case)

  #assign weight to p q pq?
  p_cont<- w[, 1]

  }
calc_armitage(a)
head(calc_armitage(a))
m
p_cont
