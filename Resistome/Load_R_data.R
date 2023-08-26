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
library(microViz)


### loading a subset of metagenomic data into phyloseq format
Rps= readRDS("Phyloseq") # this reads a pre-existing phyloseq object containing OTU and tax tables, with [fF]irm_x_x names as sample_names
Rps_tpm = readRDS("Phyloseq_tpm") # 

#We rewrite the sample names to a format filtering out Firm and firm and the first underscore so that it lines up with the column of our meta data
sample_names(Rps) = sapply(regmatches(sample_names(Rps), regexpr("_", sample_names(Rps)), invert = TRUE), "[[", 2) 

# Because the names in both metadata sets do not completely overlap, we need to manually edit one of the samples whose name was not included in the FIRM metadata file

sample_names(Rps)[68] = "4_65"

# reading in and combining metadata from 16S and metagenomic origins, adding missing underscores
firm_names = read_excel("./Metagenomic/FIRM_MetaNames.xlsx")
firm_names = firm_names[,-2] # Remove wrongful Raw_data_name column, to avoid confusion

meta_data = read.csv("MetaData.csv", header = TRUE, sep = ",")
meta_data_R = dplyr::right_join(firm_names, meta_data, by="SampleID")

# using Sample_Unique as rownames so we can match the two sets in phyloseq
meta_data_R %<>% remove_rownames %>% column_to_rownames(var="Sample_Unique")

# creating tree and making phyloseq components, adding tree and sample data components to phyloseq
random_tree = rtree(ntaxa(Rps), rooted=TRUE, tip.label=taxa_names(Rps))
meta_data_R = sample_data(meta_data_R)
Rps = merge_phyloseq(Rps, meta_data_R, random_tree)
class(Rps)


# repeat for tpm

random_tree2 = rtree(ntaxa(Rps_tpm), rooted=TRUE, tip.label=taxa_names(Rps_tpm))
Rps_tpm = merge_phyloseq(Rps_tpm, meta_data_R, random_tree2)

#rounding the "counts" for phyloseq to function

otu_table(Rps) = otu_table(round(as((otu_table(Rps)), "matrix")), taxa_are_rows(Rps))
otu_table(Rps_tpm) = otu_table(round(as((otu_table(Rps_tpm)), "matrix")), taxa_are_rows(Rps_tpm))


# overview data
datatable(tax_table(Rps))
rank_names(Rps)
sort(get_taxa_unique(Rps, "AMR_class_primary"))
sort(sample_sums(Rps)) #min is 118 and max 93167, enormous difference
sample_variables(Rps)
taxa_names(Rps)


# absolute abundances
plot_bar(Rps, fill="AMR_class_primary")
plot_bar(Rps_tpm, fill="AMR_class_primary")

table(tax_table(Rps)[, "AMR_class_primary"])
table(tax_table(Rps)[, "ARGCluster90"])

# for plotting abundances of specific stables
Rps %>% ps_filter(FarmRoundStable == c("Farm2R1S2")) %>% plot_bar(fill="AMR_class_primary")


# relative abundances
ps_rel_abund = transform_sample_counts(Rps, function(x){x / sum(x)})
plot_bar(ps_rel_abund, fill = "AMR_class_primary") +
  geom_bar(aes(color = AMR_class_primary, fill = AMR_class_primary), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ AB, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


ps_prim <- phyloseq::tax_glom(Rps, "AMR_class_primary")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "AMR_class_primary"]

psmelt(ps_prim) %>%
  ggplot(data = ., aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

# Check the amount of uniqe ARGs in samples which have and have not been treated with antibiotics
Rps %>% ps_filter(AB == "no") %>% get_taxa_unique("ARGCluster90") 
Rps %>% ps_filter(AB == "yes") %>% get_taxa_unique("ARGCluster90") 

Rps_tpm %>% ps_filter(AB == "no") %>% get_taxa_unique("ARGCluster90") 
Rps_tpm %>% ps_filter(AB == "yes") %>% get_taxa_unique("ARGCluster90") 


# plots of relative abundances
Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)")) %>% ps_filter(AB == "no") %>% 
  ps_calc_dominant(rank = "ARGCluster90") %>% comp_barplot(tax_level = "ARGCluster90", n_taxa = 12) + coord_flip()


Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)")) %>% ps_filter(AB == "no") %>% 
  ps_calc_dominant(rank = "ARGCluster90") %>% comp_barplot(tax_level = "ARGCluster90", n_taxa = 12) + coord_flip()

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)")) %>% ps_filter(AB == "yes") %>% 
  ps_calc_dominant(rank = "ARGCluster90") %>% comp_barplot(tax_level = "AMR_class_primary", n_taxa = 12) + coord_flip()

Rps %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>% ps_calc_dominant(
  "ARGCluster90",
  threshold = 0.3,
  n_max = 6,
  var = paste("dominant", rank, sep = "_"),
  none = "none",
  other = "other"
)


# factorizing variables as not to create problems with visualisation later down the line
sample_data(Rps)$Cluster = as.factor(sample_data(Rps)$Cluster)
sample_data(Rps)$FlockSize = as.factor(sample_data(Rps)$FlockSize)
sample_data(Rps)$AgeParentStock = as.factor(sample_data(Rps)$AgeParentStock)
sample_data(Rps)$Age = as.factor(sample_data(Rps)$Age)
sample_data(Rps)$LibraryNumber = as.factor(sample_data(Rps)$LibraryNumber)
