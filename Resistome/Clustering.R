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
library(ecodist)
library(sechm)

# used the following guides: https://microbiome.github.io/OMA/clustering.html, https://microbiome.github.io/OMA/viz-chapter.html 

# Rps_new = subset_samples(Rps, Sample_Unique != "2_57")

# use this for color schemes in heatmaps?
library(viridis)

x <- y <- seq(-8*pi, 8*pi, len = 40)
r <- sqrt(outer(x^2, y^2, "+"))
filled.contour(cos(r^2)*exp(-r/(2*pi)), 
               axes=FALSE,
               color.palette=viridis,
               asp=1)



tse = makeTreeSummarizedExperimentFromPhyloseq(Rps)

tse <- transformCounts(tse, method = "relabundance")

x <- t(assay(tse, "relabundance"))
hclust.out <- clusterRows(x, HclustParam())
colData(tse)$clusters <- hclust.out
hclust.out$clusters

# Community type clustering

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

pam.out

tse = makeTreeSummarizedExperimentFromPhyloseq(Rps)

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
#cutree(hc, k = 3)
#dendro

#dend <- color_branches(dendro, k = 2)
#labels(dend) <- NULL
#plot(dend)

# Hierarchical clustering lijkt onzinnig, 2 clusters en gebruikt gewoon de root


# Trying out different methods for finding optimal number of clusters:

tse = makeTreeSummarizedExperimentFromPhyloseq(Rps)

tse <- transformCounts(tse, method = "relabundance")

assay <- assay(tse, "relabundance")
assay <- t(assay)

diss_jaccard <- vegdist(assay, method = "jaccard")

res_jaccard <- NbClust(
  diss = diss_jaccard, distance = NULL, method = "complete",
  index = "mcclain"
)
res_jaccard$Best.nc


res_jaccard <- NbClust(data = diss_jaccard,
  diss = diss_jaccard, distance = NULL, method = "kmeans",
  index = "all"
)
res_jaccard$Best.nc


diss_jaccard <- as.matrix(diss_jaccard)
fviz_nbclust(diss_jaccard, kmeans, method = "silhouette")

set.seed(1337)
km <- kmeans(diss_jaccard, 2, nstart = 25)
colData(tse)$clusters <- as.factor(km$cluster)
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "bray")
plotReducedDim(tse, "MDS", colour_by = "Farm2")



#testing
tse2 = makeTreeSummarizedExperimentFromPhyloseq(Rps)
tse2 <- agglomerateByRank(tse2, rank = "ARGCluster90", agglomerateTree = TRUE)



#
  
tse = makeTreeSummarizedExperimentFromPhyloseq(Rps)
#agglomerateByRank(tse) 
colnames(rowData(tse)) = c("Domain","Phylum","Class","Order") # Domain = AMR_class_primary, Phylum = AMR_class_secondary, Class = ARGCluster90, Order = ID_Clust_Refsequence
tse <- agglomerateByRank(tse, rank = "Class", agglomerateTree = TRUE) # this is not working
tse_dmn <- mia::runDMN(tse, name = "DMN", k = 1:7)
tse_dmn
names(metadata(tse_dmn))
getDMN(tse_dmn)
miaViz::plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace") # Gives 2 as best fit
dmn_group <- calculateDMNgroup(tse_dmn,
                               variable = "Age", assay.type = "counts",
                               k = 3, seed = .Machine$integer.max
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
  pcoa2 = df[, 2]
)
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = vec
)

euclidean_dmm_pcoa_df

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

tse = makeTreeSummarizedExperimentFromPhyloseq(Rps)
colnames(rowData(tse)) = c("Domain","Phylum","Class","Order") # Domain = AMR_class_primary, Phylum = AMR_class_secondary, Class = ARGCluster90, Order = ID_Clust_Refsequence

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
                                  min_size = 20000, 
                                  name = "subsampled",
                                  replace = TRUE,
                                  seed = 1337)

mae = MultiAssayExperiment(c("tse" = tse, "sub_tse" = tse.subsampled))


mae[[1]] <- subsetByPrevalentTaxa(mae[[1]], rank = "Domain", prevalence = 0.2, detection = 0.001)
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

#Jaccard
dis <- vegan::vegdist(t(assays(tse)$counts), method = "jaccard")
jaccard_pcoa <- ecodist::pco(dis)



tse_order <- agglomerateByRank(tse,
                               rank = "Order",
                               onRankOnly = TRUE)
tse_order <- transformCounts(tse_order,
                             assay.type = "counts",
                             method = "relabundance")
top_taxa <- getTopTaxa(tse_order,
                       top = 10,
                       assay.type = "relabundance")
order_renamed <- lapply(rowData(tse_order)$Order,
                        function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_order)$Order <- as.character(order_renamed)
miaViz::plotAbundance(tse_order,
              assay.type = "relabundance",
              rank = "Class",
              order_rank_by = "abund",
              order_sample_by = "c__Clostridia")


tse_order$Farm2 = as.factor(tse_order$Farm2)
tse_order$AB = as.factor(tse_order$AB)


plots <- miaViz::plotAbundance(tse_order,
                       assay.type = "relabundance",
                       rank = "Order",
                       order_rank_by = "abund",
#                       order_sample_by = "o__Clostridiales",
                        order_sample_by = "AB",
                       features = "AB")

plots[[1]] <- plots[[1]] +
  theme(legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))
plots[[2]] <- plots[[2]] +
  theme(legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.direction = "vertical")

legend <- wrap_plots(as_ggplot(get_legend(plots[[1]])), as_ggplot(get_legend(plots[[2]])), ncol = 1) 
plots[[1]] <- plots[[1]] + theme(legend.position = "none")
plots[[2]] <- plots[[2]] + theme(legend.position = "none", axis.title.x=element_blank()) 

plot <- wrap_plots(plots[[2]], plots[[1]], ncol = 1, heights = c(2, 10))
wrap_plots(plot, legend, nrow = 1, widths = c(2, 1))




tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

tse_phylum <- transformCounts(tse_phylum, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)
tse_phylum <- transformCounts(tse_phylum, assay.type = "clr",
                              MARGIN = "features", 
                              method = "z", name = "clr_z")


top_taxa <- getTopTaxa(tse_phylum, top = 20)
tse_phylum <- tse_phylum[top_taxa, ]

mat <- assay(tse_phylum, "clr_z")

pheatmap(mat)


taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

taxa_clusters <- cutree(tree = taxa_hclust, k = 3)

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters

rowData(tse_phylum)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_phylum))), ]

# Prints taxa and their clusters
rowData(tse_phylum)$clusters


sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

# to view the tree, run
sample_tree

# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 3))

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 

# Order data based on 
tse_phylum <- tse_phylum[ , rownames(sample_data)]

# Add sample type data
sample_data$sample_types <- colData(tse_phylum)$Farm2

sample_data


breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)

pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_data,
         breaks = breaks,
         color = colors)

sechm(tse_phylum, 
      features = rownames(tse_phylum), 
      assayName = "clr", 
      do.scale = TRUE, 
      top_annotation = c("AB"), 
      gaps_at = "AB",
      cluster_cols = TRUE, cluster_rows = TRUE)

sessionInfo()