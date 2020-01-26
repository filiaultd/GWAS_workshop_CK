
results_file = "./data/subset_flowering_time_16_gwas.csv"
results <- read.csv(file=results_file, stringsAsFactors=FALSE)

results_file_nc = "./data/subset_flowering_time_16_gwas_nc.csv"
results_nc <- read.csv(file=results_file_nc, stringsAsFactors=FALSE)

results_file_maf1 = './data/subset_flowering_time_16_gwas_maf1.csv'
results_maf1<- read.csv(file=results_file_maf1, stringsAsFactors=FALSE)

# Make sure what we've loaded is as expected
dim(results)
head(results)

##################################################
### plotting -log10 p-value across the genome by chromosome
### takes results file from script 2
### uses the column headers, so don't change those!
###########################################
at.chr.lengths <- c(30427671,19698289,23459830,18585056,26975502) # the lengths of the 5 chromosomes in A. thaliana

manhattan.plot <- function(results, chr.lengths){
  # the plot will change colors to denote new chromsomes.  This section lays the groundwork for this
  chr.add <- c(0,cumsum(chr.lengths))[1:5]
  max.bp <- sum(chr.lengths)
  chr.colors <- rep(c("blue", "dodgerblue"),ceiling(length(chr.lengths/2)))
  chr.mids <- chr.add + (chr.lengths/2)
  results.s <- split(results, results$chr)
  
  # make the plot 
  # generate an empty plot with the right dimensions
  plot(results$pos,-log10(results$pvalue), xlim=c(0,max.bp), type="n", xlab="Chromosome", ylab="-log10 p-value", xaxt="n")
  # set up the x axis
  axis(1,at=chr.mids,labels=c(1:5))
  # does plotting by chromosome
  for(up.chr in 1:length(results.s)){
    print(up.chr)
    up.c <- results.s[[up.chr]]
    up.add <- chr.add[up.chr]
    up.c$pos.plot <- up.c$pos + up.add
    points(up.c$pos.plot, -log10(up.c$pvalue), col=chr.colors[up.chr])
  }
}

manhattan.plot(results=results, chr.lengths=at.chr.lengths)

manhattan.plot(results=results_nc, chr.lengths=at.chr.lengths)

manhattan.plot(results=results_maf1, chr.lengths=at.chr.lengths)

########### qq plots
library(lattice)

par(mfrow=c(1,2))

qqmath(~ -log10(results$pvalue), distribution=function(x){-log10(qunif(1-x))},xlab="expected -log10 pvalues", ylab="observed -log10 pvalues", main="QQplot of p-values",panel = function(x, ...) {
  panel.abline(a=0,b=1)
  panel.qqmath(x, ...)})

qqmath(~ -log10(results_nc$pvalue), distribution=function(x){-log10(qunif(1-x))},xlab="expected -log10 pvalues", ylab="observed -log10 pvalues", main="QQplot of p-values",panel = function(x, ...) {
  panel.abline(a=0,b=1)
  panel.qqmath(x, ...)})



exp.pvalues<-(rank(results_nc$pvalue, ties.method="first")+.5)/(length(results_nc$pvalue)+1)

#Make plot
plot(-log10(exp.pvalues), -log10(results_nc$pvalue), asp=1)
abline(0,1)
