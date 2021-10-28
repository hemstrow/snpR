# dat <- stickSNPs
# 
# #=================make the function=================
# # here's an example of the structure of the function
# # abba_baba(dat, facet, levels = c(p1, p2, p3, pO)) # data, then the sample metadata column that defines populations/species, then the identity of the p1, p2, p3, pO populations
# # abba_baba(dat, "pop", levels = c("ASP", "CLF", "SMR", "OPL")) # might look like this!
# 
# 
# # start by making a table of the ancestral and derived allele counts in each population
# 
# dat <- add.facets.snpR.data(dat, "pop")
# temp <- cbind(dat@facet.meta, dat@geno.tables$as)
# # temp$anc_count <- 0
# # temp$der_count <- 0
# # # if you can fill these two columns, we should be good to calculate the abba/baba stats!
# # 
# # 
# # # this might help, this code does something really similar, but for major and minor alleles
# # get.ac <- function(x, maj, min, mis.al){
# #   # initialize:
# #   if(is.null(nrow(x))){
# #     temp.x <- matrix(x, ncol = length(x))
# #     colnames(temp.x) <- names(x)
# #     x <- temp.x
# #   }
# #   
# #   out <- data.frame(n_total = numeric(length(maj)),
# #                     n_alleles = numeric(length(maj)),
# #                     ni1 = numeric(length(maj)),
# #                     ni2 = numeric(length(maj)))
# #   
# #   
# #   # get the column from as matching the target allele.
# #   maj.col.match <- match(maj, colnames(x))
# #   out$ni1 <- t(x)[maj.col.match + seq(0, length(x) - ncol(x), by = ncol(x))]
# #   
# #   # ni1 is the rowsums minus this
# #   out$ni2 <- rowSums(x) - out$ni1
# #   out$n_total <- rowSums(x)
# #   out$n_alleles <- rowSums(ifelse(out[,3:4] != 0, 1, 0))
# #   out[is.na(out)] <- 0
# #   
# #   return(out)
# # }
# # 
# # 
# # #finding anc count for each allele
# # table(temp[,8])
# # 
# # #what we're trying to get: anc and der count for each level and each site
# # 
# # #test
# # #setup
# # anc<-temp[,8]
# # anc_countA <-length(which(anc=="A"))
# # anc_countT <-length(which(anc=="T"))
# # anc_countC <-length(which(anc=="C"))
# # anc_countG <-length(which(anc=="G"))
# # 
# # #functions
# # get.anc <-function(anc){
# #   if(anc=="A") return(anc_countA)
# #   if(anc=="T")return(anc_countT)
# #   if(anc=="C")return(anc_countC)
# #   if(anc=="G")return(anc_countT)
# #   
# # }
# # 
# # get.der <-function(anc){
# #   if(anc=="T") return(anc_countA)
# #   if(anc=="A")return(anc_countT)
# #   if(anc=="G")return(anc_countC)
# #   if(anc=="C")return(anc_countT)
# #   
# # }
# # 
# # #apply function to data, create vectors
# # anctotals <-sapply(anc, get.anc)
# # dertotals <-sapply(anc, get.der)
# # #add vector to table
# # temp$anc_count <- anctotals
# # temp$der_count <- dertotals
# # 
# # #values needed:
# # anc and der count at each locus. ex. row 15, anc is 28
# # locus is position on the table
# # 
# # #11/4/2020 looping each nucleotide
# # 
# A<-temp[,9]
# C<-temp[,10]
# G<-temp[,11]
# t<-temp[,12]
# # 
# # 
# # get.anc <-function(anc){
# #   if(anc=="A") return(A)
# #   if(anc=="T")return(t)
# #   if(anc=="C")return(C)
# #   if(anc=="G")return(G)
# # }
# # #tapply good to use
# # 
# # get.der <-function(anc){
# #   if(anc=="T") return(A)
# #   if(anc=="A")return(t)
# #   if(anc=="G")return(C)
# #   if(anc=="C")return(G)
# #   
# # }
# # 
# # #apply function to data, create vectors
# # anctotals <-sapply(anc, get.anc)
# # dertotals <-sapply(anc, get.der)
# # #add vector to table
# # temp$anc_count <- anctotals
# # temp$der_count <- dertotals
# # 
# # 
# # #subset data to pop for calculation
# # test <- sample(c("A", "C", "G", "T"), 1000, T)
# # test <- data.frame(anc = test,
# #                    A = sample(100, 1000, T),
# #                    T = sample(100, 1000, T),
# #                    C = sample(100, 1000, T),
# #                    G = sample(100, 1000, T))
# # 
# temp$anc_count <- 0

temp <- readRDS("R_dev/abba_baba_test_data.RDS")
#colnames(temp)[12:13] <- c("der","ref")
#saveRDS(temp,"R_dev/abba_baba_test_data.RDS")
#temp$ref_count <-NULL
#temp$der_count <- NULL




# temp[which(temp$anc == "A"),]$anc_count <- temp[which(temp$anc == "A"),]$A



#using temp bc i don't know how to put sample data into a table w/columns like temp
temp$der_count <-0
temp[which(temp$der == "A"),]$der_count <- temp[which(temp$der == "A"),]$A
temp[which(temp$der == "C"),]$der_count <- temp[which(temp$der == "C"),]$C
temp[which(temp$der == "T"),]$der_count <- temp[which(temp$der == "T"),]$T
temp[which(temp$der == "G"),]$der_count <- temp[which(temp$der == "G"),]$G

temp$ref_count <-0
temp[which(temp$ref == "A"),]$ref_count <- temp[which(temp$ref == "A"),]$A
temp[which(temp$ref == "C"),]$ref_count <- temp[which(temp$ref == "C"),]$C
temp[which(temp$ref == "T"),]$ref_count <- temp[which(temp$ref == "T"),]$T
temp[which(temp$ref == "G"),]$ref_count <- temp[which(temp$ref == "G"),]$G

totals <- temp[,8:11]
rowSums(temp[,8:11], na.rm=T)

#temp$der_count <- rowSums(temp[,8:11]) - temp$anc_count
#subfacet is the diff species/populations that are being compared

#user input: define p1, p2, p3, and O


# input=make a combo c() with the necessary pops, output=abba baba calculaltions
temp$der_freq <- (temp$der_count)/(temp$der_count + temp$ref_count)

p1 <- "ASP"
p2 <-"CLF"
p3 <-"OPL"
pO <-"PAL"

A <- temp$der_freq[temp$subfacet == p1]
B <- temp$der_freq[temp$subfacet == p2]
C <- temp$der_freq[temp$subfacet == p3]
D <- temp$der_freq[temp$subfacet == pO]

#abba
abba <- (1-A)*B*C*(1-D)
#baba
baba <- A*(1-B)*C*(1-D)

D_stat <- (sum(abba) - sum(baba)) / (sum(abba) + sum(baba))

devtools::load_all(".")

#create a function that has 2 arguments: the dataset and a vector of the population names
#in the order p1, p2,p3, and pO. It returns the D-stat

popnames <- c("ASP", "CLF", "OPL", "PAL")
Dstat_calculate <- function (data, popnames ) {
  data <- data[which(data$subfacet %in% popnames),]
  data$der_count <-0
  data[which(data$der == "A"),]$der_count <- data[which(data$der == "A"),]$A
  data[which(data$der == "C"),]$der_count <- data[which(data$der == "C"),]$C
  data[which(data$der == "T"),]$der_count <- data[which(data$der == "T"),]$T
  data[which(data$der == "G"),]$der_count <- data[which(data$der == "G"),]$G
  
  data$ref_count <-0
  data[which(data$ref == "A"),]$ref_count <- data[which(data$ref == "A"),]$A
  data[which(data$ref == "C"),]$ref_count <- data[which(data$ref == "C"),]$C
  data[which(data$ref == "T"),]$ref_count <- data[which(data$ref == "T"),]$T
  data[which(data$ref == "G"),]$ref_count <- data[which(data$ref == "G"),]$G
  
  data$der_freq <- (data$der_count)/(data$der_count + data$ref_count)
  
  p1 <- popnames[1]
  p2 <- popnames[2]
  p3 <- popnames[3]
  pO <- popnames[4]
  A <- data$der_freq[data$subfacet == p1]
  B <- data$der_freq[data$subfacet == p2]
  C <- data$der_freq[data$subfacet == p3]
  D <- data$der_freq[data$subfacet == pO]
  
  abba <- (1-A)*B*C*(1-D)
  baba <- A*(1-B)*C*(1-D)
  
  D_stat1 <- (sum(abba) - sum(baba)) / (sum(abba) + sum(baba))
  return(D_stat1)
}
#for efficiency: subset the data such that the only calculations the function does are for the population being considered
# %in%
