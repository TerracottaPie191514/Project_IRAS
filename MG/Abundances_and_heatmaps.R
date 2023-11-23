#### Load packages
library(ggpubr) # Publication quality figures, based on ggplot2.
library(RColorBrewer) # Color options.
library(microbiome) # For data analysis and visualisation.
library(microbiomeutilities) # Some utility tools for microbiome package.
library(mia) # microbiome analysis package, making tse objects.
library(sechm) # Used for plotting heatmaps.
library(ggtree) # For creating trees, hierarchical clustering for heatmaps
library(pheatmap) # Creating heatmaps.
library(viridis) # Creating colour pallettes.
library(patchwork) # Used to add plots together into the same plot.



# used the following guides:, https://microbiome.github.io/OMA/viz-chapter.html, https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/composition-plots.html#heatmaps, https://david-barnett.github.io/microViz/articles/web-only/compositions.html#sorting-the-barplot
# https://david-barnett.github.io/microViz/articles/web-only/compositions.html

# absolute abundances, since there are a lot of phyla (43), we will only include the top 5 phyla 

subsetMG %>% aggregate_top_taxa2("Phylum", top = 5) %>% plot_bar(fill="Phylum", title = "Absolute abundances per sample")

# for plotting abundances of specific stables
subsetMG %>% aggregate_top_taxa2("Phylum", top = 5) %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% plot_bar(fill="Phylum")


# visualisation on AB at Phylum level, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- subsetMG %>% aggregate_top_taxa2("Phylum", top = 11) %>% phyloseq::tax_glom("Phylum") 
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Phylum"]

psmelt(ps_prim) %>% #AB
  ggplot(aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

psmelt(ps_prim) %>% #Age
  ggplot(aes(x = Age, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

psmelt(ps_prim) %>% #Farm
  ggplot(aes(x = Farm2, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

psmelt(ps_prim) %>% #Stable
  ggplot(aes(x = Stables, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

psmelt(ps_prim) %>% #Agent
  ggplot(aes(x = Cox, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

# visualisation on AB at Genus level, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- subsetMG %>% aggregate_top_taxa2("Genus", top = 13) %>% phyloseq::tax_glom("Genus")
# some top hits will not work properly, there are in acutality only 11 genera being selected above
# this is because of problems within the taxonomy info, empty taxonomies etc will be found, with tax_fix() can replace these ranks with their higher orders
# the problem is that when replacing these taxonomies with their higher ranks you are not looking at the abundances of these genera but rather the higher rank
# therefore we opt to not include these taxonomies and rather skip these unknown taxonomies
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Genus"]

psmelt(ps_prim) %>% # AB
  ggplot(aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

psmelt(ps_prim) %>% #Age
  ggplot(aes(x = Age, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

psmelt(ps_prim) %>% #Farm
  ggplot(aes(x = Farm2, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

psmelt(ps_prim) %>% #Stable
  ggplot(aes(x = Stables, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

psmelt(ps_prim) %>% #Agent
  ggplot(aes(x = Cox, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

# Plots of relative abundances, fixing some genes that are clustered in the data twice, showing top 12 taxa and others are clustered

# relative abundance of Phyla shown by age and farm (deprecated, farm2 needs to be made relative again)

subsetMG %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "Phylum", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Farm2") +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age = ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Farm", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of genera by farm and age")



# Relative abundance for both stable and antibiotics used

subsetMG %>% tax_fix() %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "Phylum", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "FarmRoundStable") +
  facet_wrap(
    facets = vars(AB), labeller = as_labeller(~ paste("Antiobotics used: ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Stable", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of Phyla by stable and antibiotics used")


# Same plots but with Genus

subsetMG %>% aggregate_top_taxa2("Genus", top = 8) %>% phyloseq::tax_glom("Genus") %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "Genus", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Farm2") +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age = ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Farm", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of Genera by farm and age")


subsetMG %>% aggregate_top_taxa2("Genus", top = 10) %>% phyloseq::tax_glom("Genus") %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "Genus", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
    x = "FarmRoundStable",
    n_taxa = 12, other_name = "Other ARG", merge_other = F) +
  facet_wrap(
    facets = vars(AB), labeller = as_labeller(~ paste("Antiobotics used: ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Stable", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of Genera by stable and antibiotics used")


# relabundance with other category (deprecated - pretty ugly)

subsetMG.rel <- subsetMG  %>% aggregate_top_taxa2("Phylum", top = 5) %>% microbiome::transform("compositional")


plot_composition(subsetMG.rel, x.label = "Id") + theme(legend.position = "bottom") +
  scale_fill_brewer("Phylum", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))

# instead we will plot both 16S and MG data in the same figure

# get the samples in the same order
sample_names(subset16S) = sample_names(subsetMG)

dataset1 =  ps_filter(subset16S)
dataset2 =  ps_filter(subsetMG)

dataset1 %<>% ps_mutate(dataset = "16S")
dataset2 %<>% ps_mutate(dataset = "MG")

sample_names(dataset1) <- paste(sample_names(dataset1), "16S", sep="_")
sample_names(dataset2) <- paste(sample_names(dataset2), "MG", sep="_")

combined <- phyloseq::merge_phyloseq(
  dataset1 %>% tax_agg("Phylum") %>% ps_get(),
  dataset2 %>% tax_agg("Phylum") %>% ps_get()
)

combined %>%
  comp_barplot("Phylum", facet_by = "dataset", n_taxa = 12, palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
               other_name = "Other Phylum", merge_other = F, sample_order = "asis") +
  coord_flip() + ggtitle("Metataxonomic vs metagenomic relative abundances of Phyla")


# rel abundance on phylum level (old version without other phyla)

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)

tse <- transformCounts(tse, method = "relabundance")

tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)
tse_phylum <- transformCounts(tse_phylum,
                              assay.type = "counts",
                              method = "relabundance")
top_taxa <- getTopTaxa(tse_phylum,
                       top = 10,
                       assay.type = "relabundance")
phylum_renamed <- lapply(rowData(tse_phylum)$Phylum,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_phylum)$Phylum <- as.character(phylum_renamed)

# rel abundance figures, can order by specific taxa
miaViz::plotAbundance(tse_phylum,
                      assay.type = "relabundance",
                      rank = "Phylum",
                      order_rank_by = "abund")
                    
tse_phylum$Farm2 = as.factor(tse_phylum$Farm2)
tse_phylum$AB = as.factor(tse_phylum$AB)

# Add AB plot on top

plots <- miaViz::plotAbundance(tse_phylum,
                               assay.type = "relabundance",
                               rank = "Phylum",
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


# heatmaps on phylum level

tse = subsetMG %>% aggregate_top_taxa2("Phylum", top = 11) %>% phyloseq::tax_glom("Phylum") %>% makeTreeSummarizedExperimentFromPhyloseq()

tse <- transformCounts(tse, method = "relabundance")

tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

tse_phylum <- transformCounts(tse_phylum, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)
tse_phylum <- transformCounts(tse_phylum, assay.type = "clr",
                              MARGIN = "features", 
                              method = "z", name = "clr_z")


top_taxa <- getTopTaxa(tse_phylum, top = 10)
#tse_phylum <- tse_phylum[top_taxa, ]

# Phylum AB heatmap
tse_phylum@metadata$anno_colors$AB = c(yes = "darkred",no ="darkblue")

sechm(tse_phylum,
      features = rownames(tse_phylum),
      assayName = "clr",
      do.scale = TRUE,
      top_annotation = "AB", 
      gaps_at = "AB",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)

# heatmap with AB and stable
tse_phylum@metadata$anno_colors$AB = c(yes = "darkred",no ="darkblue")

sechm(tse_phylum,
      features = rownames(tse_phylum),
      assayName = "clr",
      do.scale = TRUE,
      top_annotation = c("AB"), 
      gaps_at = "Stables",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)

# AB and agent
sechm(tse_phylum,
      features = rownames(tse_phylum),
      assayName = "clr",
      do.scale = TRUE,
      top_annotation = c("AB"), 
      gaps_at = "Cox",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)
# AB and age, we do see interesting shifts here
sechm(tse_phylum,
      features = rownames(tse_phylum),
      assayName = "clr",
      do.scale = TRUE,
      top_annotation = c("AB"), 
      gaps_at = "Age",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)

# Phylum heatmap
mat <- assay(tse_phylum, "clr_z")

pheatmap(mat)

# Phylum heatmap hierarchal clustering with AB

# Clustering both samples and features hierarchically 

taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

taxa_tree # based on this three, we'll create two clusters

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

taxa_clusters <- cutree(tree = taxa_hclust, k = 3) # 3 clusters based on tree figure

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
sample_clusters <- factor(cutree(tree = sample_hclust, k = 2)) # 2 clusters based on methods in Clustering.R script

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 

# Order data based on 
tse_phylum <- tse_phylum[ , rownames(sample_data)]

# Add sample type data
sample_data$sample_types <- colData(tse_phylum)$AB

sample_data


breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
#colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1) replaced with viridis pallette

pheatmap(mat, annotation_row = taxa_clusters,
         annotation_col = sample_data,
         breaks = breaks,
         color = colorRampPalette(viridis(256))(length(breaks)-1))


# heatmaps on OTU level

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse <- transformCounts(tse, method = "relabundance")
tse <- transformCounts(tse, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)
tse <- transformCounts(tse, assay.type = "clr",
                       MARGIN = "features", 
                       method = "z", name = "clr_z")
top_taxa <- getTopTaxa(tse, top = 20)
tse <- tse[top_taxa, ]

# OTU heatmap AB
tse@metadata$anno_colors$AB = c(yes = "darkred",no ="darkblue")

sechm(tse, 
      features = rownames(tse), 
      assayName = "clr", 
      do.scale = TRUE, 
      top_annotation = c("AB"), 
      gaps_at = "AB",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)

# OTU heatmap
mat <- assay(tse, "clr_z")

pheatmap(mat)

# OTU heatmap hierarchal clustering with AB

# Clustering both samples and features hierarchically 

taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

taxa_tree # based on this three, we'll create two clusters

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

taxa_clusters <- cutree(tree = taxa_hclust, k = 2) # 2 clusters based on methods in Clustering.R script

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
sample_clusters <- factor(cutree(tree = sample_hclust, k = 2))

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 

# Order data based on 
tse_phylum <- tse_phylum[ , rownames(sample_data)]

# Add sample type data
sample_data$sample_types <- colData(tse_phylum)$AB

sample_data


breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )

pheatmap(mat, annotation_row = taxa_clusters,
         annotation_col = sample_data,
         breaks = breaks,
         color = colorRampPalette(viridis(256))(length(breaks)-1))


# declutter R environment by removing objects that no longer serve a purpose
rm(ps_prim, dataset1, dataset2, combined, subsetMG.rel, tse, tse_phylum, taxmat, taxic, taxa_tree, taxa_hclust, 
   taxa_clusters, sample_tree, sample_hclust, sample_data, plot, plots, phylum_renamed, mat, legend, breaks, dis, 
   sample_clusters, samples_ordered, taxa_ordered, top_taxa)
