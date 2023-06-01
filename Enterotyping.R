library(cluster)
library(clusterSim)
library(ade4)
library(tidyverse)

# used the following guide: https://enterotype.embl.de/enterotypes.html

#functions
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log(x/y))

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

# making genus abudance table from phyloseq


#filtering out samples that have multiple values (Keys are shared for 1376 rows)
row_names_df_to_remove<-c("g__uncultured_bacterium", "g__", "g__<empty>", "g__<empty>=*","g__uncultured", "g__uncultured_bacterium=*", "g__uncultured=*", "g__uncultured=*_1")


genus_abundance = subsetG %>% tax_glom(taxrank = "Genus") %>% 
                                transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
                                dplyr::select(Genus, Sample, Abundance) %>% 
                                filter(!Genus %in% row_names_df_to_remove) %>% spread(Sample, Abundance)


genus_abundance = as.data.frame(genus_abundance)

GA <- genus_abundance[,-1]
rownames(GA) <- genus_abundance[,1]

data.dist=dist.JSD(GA) 


nclusters = clusterSim::index.G1(t(GA), data.cluster, d = data.dist, centrotypes = "medoids")

nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=clusterSim::index.G1(t(GA),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index") # 3 is optimal number

data.cluster=pam.clustering(data.dist, k=3)

obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

data.denoized=noise.removal(GA, percent=0.01)


obs.pca=dudi.pca(data.frame(t(GA)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)

s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(4,2,3))


obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F)

s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))

