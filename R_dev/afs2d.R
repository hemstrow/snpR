test<- readRDS("R_dev/2d_sfs_example_data.RDS")
pop1<- test$SMR.as
pop2<- test$OPL.as
do.2afs<- function(pop1, pop2, bins){
  #get p and q for each pop
  p.1<- matrixStats::rowMaxs(pop1)/rowSums(pop1)
  q.1<- 1-p.1
  p.2<- matrixStats::rowMaxs(pop2)/rowSums(pop2)
  q.2<- 1-p.2

  #make bins
  lbin<- seq(0, .5, length.out = bins + 1)

  #do for loop
  my.mat<- matrix(NA, bins, bins)
  for(i in 1:(length(lbin)-1)){
    out<- lbin[i] < q.1 & q.1 <= lbin[i+1]
    num<- sum(out, na.rm = T)

    for(j in 1:(length(lbin)-1)){
      out2<- lbin[j] < q.2 & q.2 <= lbin[j+1]
      num2<- sum(out2, na.rm = T)
      my.mat[i, j]<- num + num2

    }

  }

  #bin names
  binnames<- lbin[-length(lbin)]
  colnames(my.mat)<- binnames
  rownames(my.mat)<- binnames

  #melt data for graph
  my.mat<- cbind(frequency = binnames, my.mat)
  melt.mat<- reshape2::melt(my.mat, id.vars = "frequency")

  #graph data
  plot.dat<- melt.mat
  plot.dat$pop1<- melt.mat[, "Var1"]
  plot.dat$pop2<- melt.mat[,"Var2"]
  plot.dat$num_loci<- melt.mat[,"value"]
  library(ggplot2)
  plot<- ggplot(plot.dat, aes(x = pop1, y = pop2)) + geom_tile(aes(fill = num_loci)) + theme_bw() + scale_fill_gradient(low = "white", high = "blue")

  print(plot)
  return(my.mat)

}
