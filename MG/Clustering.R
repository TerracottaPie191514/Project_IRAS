#### Load packages
library(mia) # Broad package, includes clustering functions
library(bluster) # Used for clustering
library(scater) # visualisation, reduced dimensions
library(scran) # A wrapper for bluster and tse objects
library(NbClust) # To find out the optimal number of clusters
library(dendextend) # For creating dendrograms with additional options
library(factoextra) # Visualize optomial number of clusters
library(cluster) # For clustering algorithms, specifically used for PAM.

# used the following guides: https://microbiome.github.io/OMA/clustering.html, https://microucph.github.io/amplicon_data_analysis/html/cluster.html, https://www.datacamp.com/tutorial/hierarchical-clustering-R, https://rpubs.com/TBrach/68544


# Trying out different distances, aggregation methods and indices for finding optimal number of clusters, on ASV level for jaccard:

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)

tse <- transformCounts(tse, method = "relabundance")

assay <- t(assay(tse, "relabundance"))

diss_jaccard <- vegdist(assay, method = "jaccard")

# different aggregation methods and indices will grant different amount of clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "complete", index = "mcclain")$Best.nc # 2 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "complete", index = "frey")$Best.nc # 2 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "complete", index = "cindex")$Best.nc # 15 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "complete", index = "silhouette")$Best.nc # 2 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "complete", index = "dunn")$Best.nc # 9 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "ward.D2", index = "silhouette")$Best.nc # 3 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "ward.D", index = "silhouette")$Best.nc # 14 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "single", index = "silhouette")$Best.nc # 2 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "average", index = "silhouette")$Best.nc # 2 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "mcquitty", index = "silhouette")$Best.nc # 12 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "median", index = "silhouette")$Best.nc # 2 clusters
NbClust(diss = diss_jaccard, distance = NULL, method = "centroid", index = "silhouette")$Best.nc # 2 clusters



# silhouette (ASW), different clustering methods
diss_jaccard <- as.matrix(diss_jaccard) 
fviz_nbclust(diss_jaccard, kmeans, method = "silhouette") # 2 seems optimal for k-means
fviz_nbclust(diss_jaccard, cluster::pam, method = "silhouette") # 2 seems optimal for PAM
fviz_nbclust(diss_jaccard, hcut, method = "silhouette") # 2 seems optimal for hcut

# k-means 
set.seed(1337)
km <- kmeans(diss_jaccard, 2, nstart = 25)
colData(tse)$clusters <- as.factor(km$cluster)
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "jaccard")
plotReducedDim(tse, "MDS", colour_by = "clusters")

# now, let's repeat this for BC

diss_bray <- vegdist(assay, method = "bray")

#diss_euk <- vegdist(assay, method = "euclidian")


NbClust(diss = diss_bray, distance = NULL, method = "complete", index = "mcclain")$Best.nc # two clusters
NbClust(diss = diss_bray, distance = NULL, method = "complete", index = "frey")$Best.nc # 3 clusters
NbClust(diss = diss_bray, distance = NULL, method = "complete", index = "cindex")$Best.nc # 15 clusters
NbClust(diss = diss_bray, distance = NULL, method = "complete", index = "silhouette")$Best.nc # two clusters
NbClust(diss = diss_bray, distance = NULL, method = "complete", index = "dunn")$Best.nc # 9 clusters
NbClust(diss = diss_bray, distance = NULL, method = "ward.D2", index = "silhouette")$Best.nc # 3 clusters
NbClust(diss = diss_bray, distance = NULL, method = "ward.D", index = "silhouette")$Best.nc # 3 clusters
NbClust(diss = diss_bray, distance = NULL, method = "single", index = "silhouette")$Best.nc # 2 clusters
NbClust(diss = diss_bray, distance = NULL, method = "average", index = "silhouette")$Best.nc # 2 clusters
NbClust(diss = diss_bray, distance = NULL, method = "mcquitty", index = "silhouette")$Best.nc # 6 clusters
NbClust(diss = diss_bray, distance = NULL, method = "median", index = "silhouette")$Best.nc # 2 clusters
NbClust(diss = diss_bray, distance = NULL, method = "centroid", index = "silhouette")$Best.nc # 2 clusters

# silhouette (ASW)
diss_bray <- as.matrix(diss_bray) 
fviz_nbclust(diss_bray, kmeans, method = "silhouette") # 2 seems optimal
fviz_nbclust(diss_bray, cluster::pam, method = "silhouette") # 2 seems optimal for PAM
fviz_nbclust(diss_bray, hcut, method = "silhouette") # 2 seems optimal for hcut

fviz_nbclust(diss_bray, kmeans, method = "gap_stat") # 1 seems optimal for k-means gap stat
fviz_nbclust(diss_bray, cluster::pam, method = "gap_stat") # 1 seems optimal for PAM gap stat
fviz_nbclust(diss_bray, hcut, method = "gap_stat") # 1 seems optimal for hcut gap stat

# k-means 
set.seed(1337)
km <- kmeans(diss_bray, 2, nstart = 25)
colData(tse)$clusters <- as.factor(km$cluster)
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "bray")
plotReducedDim(tse, "MDS", colour_by = "clusters") 


# DMM (Laplace approximation) - ASV level
tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse_dmn <- mia::runDMN(tse, name = "DMN", k = 1:7) # calculate most likely number of clusters from 1 to 7
tse_dmn
getDMN(tse_dmn)
miaViz::plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace") # 1 cluster seems optimal

# genus level

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse_genus <- agglomerateByRank(tse, rank = "Genus", agglomerateTree = TRUE)
tse_dmn <- mia::runDMN(tse_genus, name = "DMN", k = 1:7) # calculate most likely number of clusters from 1 to 7
tse_dmn
getDMN(tse_dmn)
miaViz::plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace") # Gives 1 as best fit for genus level data

# phylum level

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse_phylum <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree = TRUE)
tse_dmn <- mia::runDMN(tse_phylum, name = "DMN", k = 1:7) # calculate most likely number of clusters from 1 to 7
tse_dmn
getDMN(tse_dmn)
miaViz::plotDMNFit(tse_dmn, type = "laplace")
getBestDMNFit(tse_dmn, type = "laplace") # Gives 1 as best fit for phylum level data


# Hierarchal clustering BC asv

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse <- transformCounts(tse, method = "relabundance")
tse <- runMDS(tse,
              assay.type = "relabundance",
              FUN = vegan::vegdist,
              method = "bray"
)

hc_bray <- hclust(vegdist(t(assay(tse, "relabundance")), method = "bray"), method = "complete")
plot(hc_bray)
hcd = as.dendrogram(hc_bray)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

colorCode <- c(Control=cbPalette[2], CRC = cbPalette[3])

grouping = cutree(hc_bray, k = 2) # all methods gave 2 clusters, based on cuttree

labels_colors(hcd) <- colorCode[grouping][order.dendrogram(hcd)]
plot(hcd)

hclust.out <- clusterRows(assay, HclustParam(method = "complete"), full = TRUE) # cutting based on complete
colData(tse)$clusters <- hclust.out$clusters
dendro <- as.dendrogram(hclust.out$objects$hclust)
plot(dendro)


labels_colors(dendro) <- colorCode[grouping][order.dendrogram(dendro)]
plot(dendro)

col_val_map <- randomcoloR::distinctColorPalette(k) %>%
  as.list() %>% 
  setNames(paste0("clust_", seq(k)))

dend <- color_branches(dendro, k = 2, col = unlist(col_val_map))
labels(dend) <- NULL
plot(dend) # very similar to 16S scripts: based on two visualisations, only a few samples are clustered distinctly, 
# based on splitting at the root, which is not informative. this particular plot splits it more down the middle# PAM clustering
labels_colors(hcd) <- colorCode[pam.out][order.dendrogram(hcd)]
plot(hcd)


tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)

tse <- transformCounts(tse, method = "relabundance")

pam.out <- clusterCells(tse,
                        assay.type = "relabundance",
                        BLUSPARAM = PamParam(centers = 2)
)

pam.out


n_iterations <- 1000
previous_cluster_assignment <- NULL
cluster_assignments <- list()

# loop that runs PAM clusterings X times and stores the results in a list, additionally checks if any clusters have changed
for (i in 1:n_iterations) {
  result <- clusterCells(tse, assay.type = "relabundance", BLUSPARAM = PamParam(centers = 2))
  cluster_assignments[[i]] <- result
  
  # Check if cluster assignments have changed
  if (!is.null(previous_cluster_assignment)) {
    samples_changed <- which(result != previous_cluster_assignment)
    if (length(samples_changed) > 0) {
      cat(sprintf("In iteration %d, the following samples changed clusters: %s\n", i, paste(samples_changed, collapse = ", ")))
    }
  }
  previous_cluster_assignment <- result
}

# To see if all of the clusters are the same or not
if (all(sapply(cluster_assignments, identical, cluster_assignments[[1]]))) {
  cat("All cluster assignments are the same across iterations.\n")
} else {
  cat("Cluster assignments vary across iterations.\n")
}

# There are no differences in clusters when run 1000 times

# save to metadata and make original PCoA plot
subsetMG@sam_data$PAM_clust = pam.out
sample_data(subsetMG)$PAM_clust = as.factor(sample_data(subsetMG)$PAM_clust)
pcoa_bc = ordinate(subsetMG, "PCoA", "bray")

plot_pcoa_ordination(subsetMG, pcoa_bc, "PAM_clust", "PCoA Bray Curtis")
#plot_pcoa_ordination(subsetMG, pcoa_bc, "Cluster", "PCoA Bray Curtis")

# change shape to different variables, age
plot_ordination(subsetMG, pcoa_bc, color = "PAM_clust", shape = "Age") + 
  geom_point(size = 3)  + labs(title = "PCoA Bray curtis", color = "AMR_class_primary")

# change shape to different variables, farm
plot_ordination(subsetMG, pcoa_bc, color = "PAM_clust", shape = "Farm2") + 
  geom_point(size = 3)  + labs(title = "PCoA Bray curtis", color = "AMR_class_primary")



# Create PAM PCoA - from 2 to 10 clusters
phy_rel <- transform_sample_counts(subsetMG, function(x) log10(x+1/sum(x+1)))
UF <- UniFrac(phy_rel, weighted = TRUE)
n_clust <- 2:10
pam_list <- lapply(n_clust, function(x) pam(UF, k = x))

sil_width <- lapply(pam_list, function(x) mean(x$silinfo$widths[, "sil_width"]))
plot(n_clust, sil_width, type="l")
pcoa_data <- cmdscale(UF, eig = TRUE)
pcoa_df <- data.frame(PC1 = c(pcoa_data$points[,1]),
                      PC2 = c(pcoa_data$points[,2]),
                      Sample = rownames(pcoa_data$points))

# Add sample data
Samp <- data.frame(sample_data(subsetMG))
Samp$Sample <- sample_names(subsetMG)

pcoa_df <- merge(pcoa_df, Samp, by = "Sample")

# Add cluster information
clusters <- factor(pam_list[[which.max(sil_width)]]$clustering)
pcoa_df <- merge(pcoa_df, clusters, by.x = "Sample", by.y = "row.names")
colnames(pcoa_df)[ncol(pcoa_df)] <- "PAM"

# Variance explained
ve <- pcoa_data$eig/sum(pcoa_data$eig)

# Plot
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = PAM)) +
  theme_bw() +
  geom_point() +
  xlab(paste0("PCoA 1 (",round(ve[1]*100,1),"%)")) +
  ylab(paste0("PCoA 2 (",round(ve[2]*100,1),"%)"))

# facet by clusters and colour by farm
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Farm2)) +
  theme_bw() +
  geom_point() +
  xlab(paste0("PCoA 1 (",round(ve[1]*100,1),"%)")) +
  ylab(paste0("PCoA 2 (",round(ve[2]*100,1),"%)")) +
  facet_wrap(~PAM)

# facet by clusters and colour by AB
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = AB)) +
  theme_bw() +
  geom_point() +
  xlab(paste0("PCoA 1 (",round(ve[1]*100,1),"%)")) +
  ylab(paste0("PCoA 2 (",round(ve[2]*100,1),"%)")) +
  facet_wrap(~PAM)

# facet by clusters and colour by Age
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Age)) +
  theme_bw() +
  geom_point() +
  xlab(paste0("PCoA 1 (",round(ve[1]*100,1),"%)")) +
  ylab(paste0("PCoA 2 (",round(ve[2]*100,1),"%)")) +
  facet_wrap(~PAM)

# facet by clusters and colour by Agent
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Cox)) +
  theme_bw() +
  geom_point() +
  xlab(paste0("PCoA 1 (",round(ve[1]*100,1),"%)")) +
  ylab(paste0("PCoA 2 (",round(ve[2]*100,1),"%)")) +
  facet_wrap(~PAM)

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
colnames(prob) <- c("comp1")
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

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
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
rm(tse_phylum, tse_dmn, assay, diss_jaccard, diss_bray, k, dmn_group, euclidean_pcoa_df, euclidean_dmm_pcoa_df, plots, 
   res, ve, Samp, clusters, pam.out, phy_rel, UF, n_clust, pam_list, sil_width, pcoa_df, tse, tse_genus, hc_bray, km,
   hcd, dend, dendro, col_val_map, hclust.out, pcoa_data, prob, df, cbPalette, colorCode, grouping, vec)
