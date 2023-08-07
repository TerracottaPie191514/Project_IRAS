#### Load packages
library(microbiome)
library(phyloseq)
library(microbiomeutilities)
library(RColorBrewer)
library(ggpubr)
library(DT)
library(data.table)
library(tidyverse)
library(pheatmap)
library(picante)
library(nlme)
library(scales)
library(readxl)


### loading a subset of metagenomic data into phyloseq format
subsetMG= readRDS("Phyloseq") # this reads a pre-existing phyloseq object containing OTU and tax tables

# reading in and combining metadata from 16S and metagenomic origins, adding missing underscores
firm_names = read_excel("./Metagenomic/FIRM_MetaNames.xlsx")
meta_data = read.csv("MetaData.csv", header = TRUE, sep = ",")
meta_data_MG = left_join(firm_names, meta_data, by="SampleID")
meta_data_MG$Raw_data_name = sub(" ", "_", meta_data_MG$Raw_data_name)

# subsetting the metadata and using the metagenomic names ([Ff]irm*) as rownames
meta_data_MG_subset <- meta_data_MG[ meta_data_MG$Raw_data_name %in% colnames(subsetMG@otu_table), ]
meta_data_MG_subset %<>% remove_rownames %>% column_to_rownames(var="Raw_data_name")

# creating tree and making phyloseq components, adding tree and sample data components to phyloseq
random_tree = rtree(ntaxa(subsetMG), rooted=TRUE, tip.label=taxa_names(subsetMG))
meta_data_MG_subset = sample_data(meta_data_MG_subset)
subsetMG = merge_phyloseq(subsetMG, meta_data_MG_subset, random_tree)
class(subsetMG)

# overview data
datatable(tax_table(subsetMG))
rank_names(subsetMG)
sort(get_taxa_unique(subsetMG, "AMR_class_primary"))
sort(sample_sums(subsetMG))
sample_variables(subsetMG)
taxa_names(subsetMG)
plot_bar(subsetMG, fill="AMR_class_primary")

subsetMG %>% ps_filter(AB == "no") %>% get_taxa_unique("ARGCluster90") 
subsetMG %>% ps_filter(AB == "yes") %>% get_taxa_unique("ARGCluster90") 


sort(get_taxa_unique(subsetMG, "ARGCluster90"))
library(microViz)

subsetMG %>% tax_fix(unknowns = c("lnu(B)")) %>% ps_filter(AB == "no") %>% 
  ps_calc_dominant(rank = "ARGCluster90") %>% comp_barplot(tax_level = "ARGCluster90", n_taxa = 12) + coord_flip()


ps_calc_dominant(
  subsetMG,
  "ARGCluster90",
  threshold = 0.3,
  n_max = 6,
  var = paste("dominant", rank, sep = "_"),
  none = "none",
  other = "other"
)

# factorizing variables as not to create problems with visualisation later down the line
sample_data(subsetMG)$Cluster = as.factor(sample_data(subsetMG)$Cluster)
sample_data(subsetMG)$FlockSize = as.factor(sample_data(subsetMG)$FlockSize)
sample_data(subsetMG)$AgeParentStock = as.factor(sample_data(subsetMG)$AgeParentStock)
sample_data(subsetMG)$Age = as.factor(sample_data(subsetMG)$Age)
sample_data(subsetMG)$LibraryNumber = as.factor(sample_data(subsetMG)$LibraryNumber)