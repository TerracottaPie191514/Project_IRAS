library(mia)
library(bluster)
library(scater)
library(scran)
library(NbClust)
library(cobiclust)
library(dendextend)
library(factoextra)
library(miaViz)
library(patchwork)


install.packages("miaViz")

data("peerj13075", package = "mia")
tse <- peerj13075
rm(peerj13075)

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetG)

tse <- transformCounts(tse, method = "relabundance")

x <- t(assay(tse, "relabundance"))
hclust.out <- clusterRows(x, HclustParam())
colData(tse)$clusters <- hclust.out
hclust.out

tse <- runMDS(tse,
              assay.type = "relabundance",
              FUN = vegan::vegdist,
              method = "bray"
)

plotReducedDim(tse, "MDS", colour_by = "clusters")

hclust.out <- clusterRows(x, HclustParam(method = "complete"), full = TRUE)
colData(tse)$clusters <- hclust.out$clusters
dendro <- as.dendrogram(hclust.out$objects$hclust)
plot(dendro)


pam.out <- clusterCells(tse,
                        assay.type = "relabundance",
                        BLUSPARAM = PamParam(centers = 5)
)

assay <- assay(tse, "relabundance")
assay <- t(assay)
diss <- vegdist(assay, method = "bray")
hc <- hclust(diss, method = "complete")
dendro <- as.dendrogram(hc)
plot(dendro)

res <- NbClust(
  diss = diss, distance = NULL, method = "ward.D2",
  index = "silhouette"
)
res$Best.nc
cutree(hc, k = 2)
dendro

dend <- color_branches(dendro, k = 2)
labels(dend) <- NULL
plot(dend)

# Hierarchical clustering lijkt onzinnig, 2 clusters en gebruikt gewoon de root

diss <- as.matrix(diss)
fviz_nbclust(diss, kmeans, method = "silhouette")

set.seed(1337)
km <- kmeans(diss, 2, nstart = 25)
colData(tse)$clusters <- as.factor(km$cluster)
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "bray")
plotReducedDim(tse, "MDS", colour_by = "clusters")


tse = makeTreeSummarizedExperimentFromPhyloseq(subsetG)
tse <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree = TRUE)
tse_dmn <- mia::runDMN(tse, name = "DMN", k = 1:7)
tse_dmn
names(metadata(tse_dmn))
getDMN(tse_dmn)
miaViz::plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace")
dmn_group <- calculateDMNgroup(tse_dmn,
                               variable = "Age", assay.type = "counts",
                               k = 2, seed = .Machine$integer.max
)

dmn_group
DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))
prob <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))
colnames(prob) <- c("comp1", "comp2")
vec <- colnames(prob)[max.col(prob, ties.method = "first")]

assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformCounts(tse, assay.type = "pseudo", method = "relabundance")
tse <- transformCounts(tse, "relabundance", method = "clr")
df <- calculateMDS(tse, assay.type = "clr", method = "euclidean")
euclidean_pcoa_df <- data.frame(
  pcoa1 = df[, 1],
  pcoa2 = df[, 2]
)
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = vec
)
euclidean_dmm_plot <- ggplot(
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
    title = "PCoA with Aitchison distances"
  ) +
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_dmm_plot

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetG)
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

data("HintikkaXOData")
