#### Load packages
library(ggpubr) # Publication quality figures, based on ggplot2.
library(RColorBrewer) # Color options.
library(microbiomeutilities) # Some utility tools for microbiome package.
library(mia) # microbiome analysis package, making tse objects.
library(sechm) # Used for plotting heatmaps.
library(ggtree) # For creating trees, hierarchical clustering for heatmaps
library(pheatmap) # Creating heatmaps.
library(viridis) # Creating colour pallettes.
library(patchwork) # Used to add plots together into the same plot.


# used the following guides : https://david-barnett.github.io/microViz/articles/web-only/compositions.html, https://microbiome.github.io/OMA/viz-chapter.html

# absolute abundances - phylum
plot_bar(subset16S, fill="Phylum", title = "Absolute abundances per sample")

# for plotting abundances of specific stables
subset16S %>% ps_filter(Stables == c("Farm2R1S1")) %>% plot_bar(fill="Phylum")


# visualisation on AB at Phylum level, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- phyloseq::tax_glom(subset16S, "Phylum")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Phylum"]

psmelt(ps_prim) %>% # AB
  ggplot(data = ., aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free")

psmelt(ps_prim) %>% # Age
  ggplot(data = ., aes(x = Age, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free")

psmelt(ps_prim) %>% # Farm
  ggplot(data = ., aes(x = Farm2, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free")

psmelt(ps_prim) %>% # Stable
  ggplot(data = ., aes(x = Stables, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free")

psmelt(ps_prim) %>% # Agent
  ggplot(data = ., aes(x = Cox, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free")

# visualisation on AB at Genus level, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- subset16S %>% aggregate_top_taxa2("Genus", top = 5) %>% phyloseq::tax_glom("Genus")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Genus"]

psmelt(ps_prim) %>% # AB
  ggplot(data = ., aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Genus, scales = "free")

psmelt(ps_prim) %>% # Age
  ggplot(data = ., aes(x = Age, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Genus, scales = "free")

psmelt(ps_prim) %>% # Farm
  ggplot(data = ., aes(x = Farm2, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Genus, scales = "free")

psmelt(ps_prim) %>% # Stable
  ggplot(data = ., aes(x = Stables, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Genus, scales = "free")

psmelt(ps_prim) %>% # Agent
  ggplot(data = ., aes(x = Cox, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Genus, scales = "free")

# Check the amount of unique genera in samples which have and have not been treated with antibiotics
subset16S %>% ps_filter(AB == "no") %>% get_taxa_unique("Genus") # 93 different genera for non AB treated
subset16S %>% ps_filter(AB == "yes") %>% get_taxa_unique("Genus") # 74 different genera for AB treated
subset16S %>% get_taxa_unique("Genus") # 93 different genes in total, which are all present in non-treated

# Plots of relative abundances, fixing some genes that are clustered in the data twice, showing top 12 taxa and others are clustered

# relative abundance of Phyla shown by age and farm (deprecated, farm2 needs to be made relative again)

subset16S %>%
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
  ggtitle("Relative abundance of Phyla by farm and age")



# Relative abundance for both stable and antibiotics used

subset16S %>% 
  ps_arrange(Stables) %>%
  ps_mutate(
    Stables = factor(Stables, rev(unique(Stables)))
  ) %>%
  comp_barplot(
    tax_level = "Phylum", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Stables") +
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

subset16S %>% aggregate_top_taxa2("Genus", top = 8) %>% phyloseq::tax_glom("Genus") %>%
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


subset16S %>% aggregate_top_taxa2("Genus", top = 8) %>% phyloseq::tax_glom("Genus") %>% 
  ps_arrange(Stables) %>%
  ps_mutate(
    Stables = factor(Stables, rev(unique(Stables)))
  ) %>%
  comp_barplot(
    tax_level = "Genus", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
    x = "Stables",
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


# figures with other category (deprecated - ugly plot)

subset16S.rel <- microbiome::transform(subset16S, "compositional")

plot_composition(subset16S.rel, x.label = "Id") + theme(legend.position = "bottom") +
  scale_fill_brewer("Phylum", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))

# another relative abundance plot (very ugly)

ps1.com <- subset16S

# if you have dada2/deblur output and sequences as taxa names, then you can change them as follows
taxa_names(ps1.com) <- paste0("ASV_", rownames(tax_table(ps1.com)))

# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

# colourCount = length(unique(taxic$Family))  #define number of variable colors based on number of Family (change the level accordingly to phylum/class/order)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.


taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object


# now edit the unclassified taxa
tax_table(ps1.com)[tax_table(ps1.com)[, "Phylum"] == "", "Phylum"] <- "Unclassified phylum"

# it would be nice to have the Taxonomic names in italics.
# for that we set this

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))


## Now we need to plot at family level, we can do it as follows:

# first remove the phy_tree

ps1.com@phy_tree <- NULL

# Second merge at family level

ps1.com.fam <- microbiomeutilities::aggregate_top_taxa2(ps1.com, "Phylum", top = 10)

plot_composition(ps1.com.fam) + theme(legend.position = "bottom") +
  scale_fill_brewer("Family", palette = "Paired") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

# rel abundance on phylum level

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)

tse <- transformCounts(tse, method = "relabundance")

tse_phylum <- agglomerateByRank(tse,
                               rank = "Phylum",
                               onRankOnly = TRUE)
tse_phylum <- transformCounts(tse_order,
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
                      order_rank_by = "abund",
                      order_sample_by = "p__Firmicutes")


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

tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

tse_phylum <- transformCounts(tse_phylum, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)
tse_phylum <- transformCounts(tse_phylum, assay.type = "clr",
                              MARGIN = "features", 
                              method = "z", name = "clr_z")


#top_taxa <- getTopTaxa(tse_phylum, top = 20) there are few phyla in this data so no need to exclude some
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

taxa_clusters <- cutree(tree = taxa_hclust, k = 2) # 2 clusters based on tree figure

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


# heatmaps on ASV level

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse <- transformCounts(tse, method = "relabundance")
tse <- transformCounts(tse, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)
tse <- transformCounts(tse, assay.type = "clr",
                              MARGIN = "features", 
                              method = "z", name = "clr_z")
top_taxa <- getTopTaxa(tse, top = 20)
tse <- tse[top_taxa, ]

# ASV heatmap AB
tse@metadata$anno_colors$AB = c(yes = "darkred",no ="darkblue")

sechm(tse, 
      features = rownames(tse), 
      assayName = "clr", 
      do.scale = TRUE, 
      top_annotation = c("AB"), 
      gaps_at = "AB",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)

# ASV heatmap
mat <- assay(tse, "clr_z")

pheatmap(mat)

# ASV heatmap hierarchal clustering with AB

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
rm(ps_prim, tse, tse_phylum, taxmat, taxic, taxa_tree, taxa_hclust, taxa_clusters, subset16S.rel, sample_tree,
   sample_hclust, sample_data, ps1.com, ps1.com.fam, plots, plot, phylum_renamed, mat, legend, guide_italics,
   breaks, new.tax, sample_clusters, samples_ordered, taxa_ordered, top_taxa)
