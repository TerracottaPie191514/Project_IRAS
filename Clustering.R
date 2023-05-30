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
library(pheatmap)
library(ggplot2)
library(biclust)

install.packages("miaViz")


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
plotReducedDim(tse, "MDS", colour_by = "Farm2")


  
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

dmn_group <- calculateDMNgroup(tse_dmn,
                               variable = "Farm2", assay.type = "counts",
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
  pcoa2 = df[, 2],
  pcoa3 = df[, 3],
  pcoa4 = df[, 4]
)
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = vec
)

euclidean_dmm_pcoa_df

test = cbind(euclidean_dmm_pcoa_df,metadf)

euclidean_dmm_plot <- ggplot(
  data = test,
  aes(
    x = pcoa1, y = pcoa2,
    color = FeedType
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

tse.subsampled <- subsampleCounts(tse, 
                                  min_size = 60000, 
                                  name = "subsampled",
                                  replace = TRUE,
                                  seed = 1337)

mae = MultiAssayExperiment(c("tse" = tse, "sub_tse" = tse.subsampled))

mae[[1]] <- subsetByPrevalentTaxa(mae[[1]], rank = "Genus", prevalence = 0.2, detection = 0.001)
mae[[1]] <- transformCounts(mae[[1]], method = "relabundance")
mae[[1]] <- transformCounts(mae[[1]], "relabundance", method = "rclr")
clusters <- cobiclust(assay(mae[[1]], "counts"))
row_clusters <- clusters$classification$rowclass
col_clusters <- clusters$classification$colclass
rowData(mae[[1]])$clusters <- factor(row_clusters)
colData(mae[[1]])$clusters <- factor(col_clusters)
mae[[1]] <- mae[[1]][order(rowData(mae[[1]])$clusters), order(colData(mae[[1]])$clusters)]
clusters$classification

mae[[1]] <- transformCounts(mae[[1]],
                            assay.type = "rclr",
                            MARGIN = "features",
                            method = "z", name = "clr_z"
)
annotation_col <- data.frame(colData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_col) <- "col_clusters"

annotation_row <- data.frame(rowData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_row) <- "row_clusters"

pheatmap(assay(mae[[1]], "clr_z"),
         cluster_rows = F, cluster_cols = F,
         annotation_col = annotation_col,
         annotation_row = annotation_row
)
melt_assay <- meltAssay(mae[[1]], assay.type = "rclr", add_col_data = T, add_row_data = T)
p1 <- ggplot(melt_assay) +
  geom_boxplot(aes(x = clusters.x, y = rclr)) +
  labs(x = "Taxa clusters")

p2 <- ggplot(melt_assay) +
  geom_boxplot(aes(x = clusters.y, y = rclr)) +
  labs(x = "Sample clusters")

p1 + p2


mae[[1]] <- mae[[1]][, colnames(mae[[2]])]
rownames(mae[[1]]) <- make.unique(rownames(mae[[1]]))
corr <- getExperimentCrossCorrelation(mae, 1, 2,
                                      assay.type1 = "rclr",
                                      assay.type2 = "counts",
                                      mode = "matrix",
                                      cor_threshold = 0.2
)

set.seed(1337)
bc <- biclust(corr,
              method = BCPlaid(), fit.model = y ~ m,
              background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10,
              iter.startup = 10, iter.layer = 100, verbose = FALSE
)

.get_biclusters_from_biclust <- function(bc, assay) {
  # Get cluster information for columns and rows
  bc_columns <- t(bc@NumberxCol)
  bc_columns <- data.frame(bc_columns)
  bc_rows <- bc@RowxNumber
  bc_rows <- data.frame(bc_rows)
  
  # Get data into right format
  bc_columns <- .manipulate_bc_data(bc_columns, assay, "col")
  bc_rows <- .manipulate_bc_data(bc_rows, assay, "row")
  
  return(list(bc_columns = bc_columns, bc_rows = bc_rows))
}

.manipulate_bc_data <- function(bc_clusters, assay, row_col) {
  # Get right dimension
  dim <- ifelse(row_col == "col", ncol(assay), nrow(assay))
  # Get column/row names
  if (row_col == "col") {
    names <- colnames(assay)
  } else {
    names <- rownames(assay)
  }
  
  # If no clusters were found, create one. Otherwise create additional
  # cluster which
  # contain those samples that are not included in clusters that were found.
  if (nrow(bc_clusters) != dim) {
    bc_clusters <- data.frame(cluster = rep(TRUE, dim))
  } else {
    # Create additional cluster that includes those samples/features that
    # are not included in other clusters.
    vec <- ifelse(rowSums(bc_clusters) > 0, FALSE, TRUE)
    # If additional cluster contains samples, then add it
    if (any(vec)) {
      bc_clusters <- cbind(bc_clusters, vec)
    }
  }
  # Adjust row and column names
  rownames(bc_clusters) <- names
  colnames(bc_clusters) <- paste0("cluster_", 1:ncol(bc_clusters))
  return(bc_clusters)
}

bcs <- .get_biclusters_from_biclust(bc, corr)

bicluster_rows <- bcs$bc_rows
bicluster_columns <- bcs$bc_columns
head(bicluster_rows)

.sum_mean_median_var <- function(tse1, tse2, assay.type1, assay.type2, clusters1, clusters2) {
  list <- list()
  # Create a data frame that includes all the information
  for (i in 1:ncol(clusters1)) {
    # Subset data based on cluster
    tse_subset1 <- tse1[clusters1[, i], ]
    tse_subset2 <- tse2[clusters2[, i], ]
    # Get assay
    assay1 <- assay(tse_subset1, assay.type1)
    assay2 <- assay(tse_subset2, assay.type2)
    # Calculate sum, mean, median, and mean variance
    sum1 <- colSums2(assay1, na.rm = T)
    mean1 <- colMeans2(assay1, na.rm = T)
    median1 <- colMedians(assay1, na.rm = T)
    var1 <- colVars(assay1, na.rm = T)
    
    sum2 <- colSums2(assay2, na.rm = T)
    mean2 <- colMeans2(assay2, na.rm = T)
    median2 <- colMedians(assay2, na.rm = T)
    var2 <- colVars(assay2, na.rm = T)
    
    list[[i]] <- data.frame(
      sample = colnames(tse1), sum1, sum2, mean1, mean2,
      median1, median2, var1, var2
    )
  }
  
  return(list)
}

df <- .sum_mean_median_var(mae[[1]], mae[[2]], "rclr", "nmr", bicluster_rows, bicluster_columns)

pics <- list()
for (i in seq_along(df)) {
  pics[[i]] <- ggplot(df[[i]]) +
    geom_point(aes(x = median1, y = median2)) +
    labs(
      title = paste0("Cluster ", i),
      x = "Taxa (rclr median)",
      y = "Metabolites (abs. median)"
    )
  print(pics[[i]])
}

bicluster_columns <- data.frame(apply(bicluster_columns, 2, as.factor))
bicluster_rows <- data.frame(apply(bicluster_rows, 2, as.factor))

if (ncol(bicluster_rows) > ncol(bicluster_columns)) {
  cluster_names <- colnames(bicluster_rows)
} else {
  cluster_names <- colnames(bicluster_columns)
}
annotation_colors <- list()
for (name in cluster_names) {
  annotation_colors[[name]] <- c("TRUE" = "red", "FALSE" = "white")
}

pheatmap(corr,
         cluster_cols = F, cluster_rows = F,
         annotation_col = bicluster_columns,
         annotation_row = bicluster_rows,
         annotation_colors = annotation_colors
)
