#### Load packages
library(mia) # Broad package, includes clustering functions
library(bluster) # Used for clustering
library(scater)
library(scran) # A wrapper for bluster and tse objects
library(NbClust) # To find out the optimal number of clusters
library(cobiclust)
library(dendextend) # For creating dendrograms with additional options
library(factoextra) # Visualize optomial number of clusters
library(patchwork)
library(pheatmap)
library(biclust)
library(ecodist)
library(sechm)
library(simpr)

# used the following guides: https://microbiome.github.io/OMA/clustering.html, https://microbiome.github.io/OMA/viz-chapter.html 


# Trying out different methods for finding optimal number of clusters:

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)

tse <- transformCounts(tse, method = "relabundance")

assay <- t(assay(tse, "relabundance"))

diss_jaccard <- vegdist(assay, method = "jaccard")

res_jaccard <- NbClust(
  diss = diss_jaccard, distance = NULL, method = "complete",
  index = "mcclain"
)
res_jaccard$Best.nc



# silhouette (ASW)
diss_jaccard <- as.matrix(diss_jaccard) 
fviz_nbclust(diss_jaccard, kmeans, method = "silhouette") # 2 seems optimal

# k-means
res_jaccard <- NbClust(data = diss_jaccard,
                       diss = diss_jaccard, distance = NULL, method = "kmeans",
                       index = "all")
res_jaccard$Best.nc

set.seed(1337)
km <- kmeans(diss_jaccard, 2, nstart = 25)
colData(tse)$clusters <- as.factor(km$cluster)
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "bray")
plotReducedDim(tse, "MDS", colour_by = "clusters") # optimal number is 2

# DMM (Laplace approximation) ASV
tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse_dmn <- mia::runDMN(tse, name = "DMN", k = 1:7) # calculate most likely number of clusters from 1 to 7
tse_dmn
getDMN(tse_dmn)
miaViz::plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace") # 2 again

# genus level

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse_genus <- agglomerateByRank(tse, rank = "Genus", agglomerateTree = TRUE)
tse_dmn <- mia::runDMN(tse_genus, name = "DMN", k = 1:7) # calculate most likely number of clusters from 1 to 7
tse_dmn
getDMN(tse_dmn)
miaViz::plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace") # Gives 3! as best fit for genus level data

# phylum level

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse_phylum <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree = TRUE)
tse_dmn <- mia::runDMN(tse_phylum, name = "DMN", k = 1:7) # calculate most likely number of clusters from 1 to 7
tse_dmn
getDMN(tse_dmn)
miaViz::plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace") # Gives 2 as best fit for phylum level data


# Hierarchal clustering

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse <- transformCounts(tse, method = "relabundance")
tse <- runMDS(tse,
              assay.type = "relabundance",
              FUN = vegan::vegdist,
              method = "bray"
)

assay <- t(assay(tse, "relabundance"))
diss_bray <- vegdist(assay, method = "bray")
hc_bray <- hclust(diss_bray, method = "complete")
plot(hc_bray)
hcd = as.dendrogram(hc_bray)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

colorCode <- c(Control=cbPalette[2], CRC = cbPalette[3])

grouping = cutree(hc_bray, k = 2) # most methods gave 2 clusters

labels_colors(hcd) <- colorCode[grouping][order.dendrogram(hcd)]
plot(hcd)

hclust.out <- clusterRows(assay, HclustParam(method = "complete"), full = TRUE)
colData(tse)$clusters <- hclust.out$clusters
dendro <- as.dendrogram(hclust.out$objects$hclust)
plot(dendro)


labels_colors(dendro) <- colorCode[grouping][order.dendrogram(dendro)]
plot(dendro)


pam.out <- clusterCells(tse,
                        assay.type = "relabundance",
                        BLUSPARAM = PamParam(centers = 2)
)

pam.out

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)

tse <- transformCounts(tse, method = "relabundance")

assay <- assay(tse, "relabundance")
assay <- t(assay)

diss_bray <- vegdist(assay, method = "bray")
hc_bray <- hclust(diss_bray, method = "complete")

res_bray <- NbClust(
  diss = diss_bray, distance = NULL, method = "kmeans",
  index = "silhouette"
)
res_bray$Best.nc
cutree(hc_bray, k = 3)

dendro

dend <- color_branches(dendro, k = 2)
labels(dend) <- NULL
plot(dend)

# Hierarchical clustering lijkt onzinnig, 2 clusters en gebruikt gewoon de root


# PCoA for ASV data, BC with DMM, euclidian ( make sure tse_dmn is on right taxonomic level)

dmn_group <- calculateDMNgroup(tse_dmn,
                               variable = "Age", assay.type = "counts",
                               k = 2, seed = .Machine$integer.max)

dmn_group <- calculateDMNgroup(tse_dmn,
                               variable = "Farm2", assay.type = "counts",
                               k = 2, seed = .Machine$integer.max)
dmn_group <- calculateDMNgroup(tse_dmn,
                               variable = "AB", assay.type = "counts",
                               k = 2, seed = .Machine$integer.max)
dmn_group
DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn)) # measure weights
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))) # sample-cluster assignment probablities
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn))) # taxa contribution
prob <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))
colnames(prob) <- c("comp1", "comp2")
vec <- colnames(prob)[max.col(prob, ties.method = "first")]
assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformCounts(tse, assay.type = "pseudo", method = "relabundance")
tse <- transformCounts(tse, "relabundance", method = "clr")
df <- calculateMDS(tse, assay.type = "clr", method = "euclidean")
euclidean_pcoa_df <- data.frame(
  pcoa1 = df[, 1],
  pcoa2 = df[, 2])
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = vec)

ggplot(
  data = euclidean_dmm_pcoa_df,
  aes(
    x = pcoa1, y = pcoa2,
    color = dmm_component
  )
) +
  geom_point() +
  labs(
    x = "Coordinate 1",
    y = "Coordinate 2",
    title = "PCoA with Aitchison distances")

# UMAP with different ks

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse <- transformCounts(tse, method = "rclr")
tse <- runUMAP(tse, name = "UMAP", assay.type = "rclr")
k <- c(2, 3, 5, 10)
ClustAndPlot <- function(x) {
  # Creating the graph and running the short random walks algorithm
  graph_clusters <- clusterRows(t(assays(tse)$rclr), NNGraphParam(k = x))
  
  # Results of the clustering as a color for each sample
  plotUMAP(tse, colour_by = I(graph_clusters)) +
    labs(title = paste0("k = ", x))
}
plots <- lapply(k, ClustAndPlot)
(plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])

# boxplots
ClustDiagPlot <- function(x) {
  # Getting the clustering results
  graph_clusters <- clusterRows(t(assays(tse)$rclr), NNGraphParam(k = x))
  
  # Computing the diagnostic info
  sil <- approxSilhouette(t(assays(tse)$rclr), graph_clusters)
  
  # Plotting as a boxlpot to observe cluster separation
  boxplot(split(sil$width, graph_clusters), main = paste0("k = ", x))
}
# Applying the function for different k values
res <- lapply(k, ClustDiagPlot)


# declutter R environment by removing objects that no longer serve a purpose
rm(ps_prim)