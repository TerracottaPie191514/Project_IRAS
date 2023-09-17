#### Load packages
library(microbiome) # For data analysis and visualisation.
library(phyloseq) # Data analysis and visualisation, also the basis of data object.
library(microbiomeutilities) # Some utility tools for microbiome package.
library(RColorBrewer) # Color options.
library(ggpubr) # Publication quality figures, based on ggplot2.
library(DT) # Interactive tables in html and markdown.
library(data.table) # Giving overview of data.
library(tidyverse) # Data handling and much more.
library(readxl) # Reading in excel files.
library(ape) # Phylogenetic package, used for creating random trees and as dependency for other packages.
library(magrittr) # Data handling, specifically assignment pipes
library(microViz) # Both analysis and visualisation


### loading a subset of metagenomic data into phyloseq format
Rps= readRDS("Phyloseq") # this reads a pre-existing phyloseq object containing OTU and tax tables, with [fF]irm_x_x names as sample_names
Rps_tpm = readRDS("Phyloseq_tpm") # also read in TPM data instead of FPKM

#We rewrite the sample names to a format filtering out Firm and firm and the first underscore so that it lines up with the column of our meta data
sample_names(Rps) = sapply(regmatches(sample_names(Rps), regexpr("_", sample_names(Rps)), invert = TRUE), "[[", 2) 
sample_names(Rps_tpm) = sapply(regmatches(sample_names(Rps_tpm), regexpr("_", sample_names(Rps_tpm)), invert = TRUE), "[[", 2) 

# Because the names in both metadata sets do not completely overlap, we need to manually edit one of the samples whose name was not included in the FIRM metadata file
sample_names(Rps)[68] = "4_65"
sample_names(Rps_tpm)[68] = "4_65"


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

# we will scale the otu table by a factor 1000 to preserve data before the rounding step 

otu_table(Rps) = otu_table(Rps) * 1000
#otu_table(Rps_tpm) = otu_table(Rps_tpm) * 1000


#rounding the "counts" for phyloseq to function

otu_table(Rps) = otu_table(round(as((otu_table(Rps)), "matrix")), taxa_are_rows(Rps))
otu_table(Rps_tpm) = otu_table(round(as((otu_table(Rps_tpm)), "matrix")), taxa_are_rows(Rps_tpm))


# overview data
datatable(tax_table(Rps))
rank_names(Rps) # Shows classes and ARGs
sort(get_taxa_unique(Rps, "AMR_class_primary")) # Shows primary AMR classes
sort(sample_sums(Rps)) # Amount of unique "taxa" per sample, the min is 118 and max 93167, which is an enormous difference
summary(sample_sums(Rps)) # summary of the sampling depths
sample_variables(Rps) # metadata variables
taxa_names(Rps) # ARGs

# Stable "Farm2R1S1 has 11/12 lowest sampling depths of the dataset, and there are only 2 other samples with a lower depth than the 12th  sample
Rps %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% sample_sums() %>% sort()

# absolute abundances
plot_bar(Rps, fill="AMR_class_primary")
plot_bar(Rps_tpm, fill="AMR_class_primary")

# Amount of different AMR classes present.
sort(table(tax_table(Rps)[, "AMR_class_primary"]))
sort(table(tax_table(Rps)[, "ARGCluster90"]))

# for plotting abundances of specific stables
Rps %>% ps_filter(FarmRoundStable == c("Farm2R1S2")) %>% plot_bar(fill="AMR_class_primary")


# relative abundances of primary AMR, split by AB (very ugly)
ps_rel_abund = transform_sample_counts(Rps, function(x){x / sum(x)})
plot_bar(ps_rel_abund, fill = "AMR_class_primary") +
  geom_bar(aes(color = AMR_class_primary, fill = AMR_class_primary), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ AB, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# visualisation on primary AMR classes, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- phyloseq::tax_glom(Rps, "AMR_class_primary")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "AMR_class_primary"]

psmelt(ps_prim) %>%
  ggplot(data = ., aes(x = FarmRoundStable, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

# Check the amount of uniqe ARGs in samples which have and have not been treated with antibiotics
Rps %>% ps_filter(AB == "no") %>% get_taxa_unique("ARGCluster90") # 163 different genes for non AB treated
Rps %>% ps_filter(AB == "yes") %>% get_taxa_unique("ARGCluster90") # 81 different genes for AB treated
Rps %>% get_taxa_unique("ARGCluster90") # 188 different genes in total, meaning 25 do not overlap


Rps_tpm %>% ps_filter(AB == "no") %>% get_taxa_unique("ARGCluster90") # 168 different genes for non AB treated
Rps_tpm %>% ps_filter(AB == "yes") %>% get_taxa_unique("ARGCluster90") # 97 different genes for AB treated
Rps_tpm %>% get_taxa_unique("ARGCluster90") # 186 different genes in total, meaning 18 do not overlap




# Plots of relative abundances, fixing some genes that are clustered in the data twice, showing top 12 taxa and others are clustered

# TPM and FPKM in the same plot, comparing non AB
dataset1 =  ps_filter(Rps)
dataset2 =  ps_filter(Rps_tpm)

dataset1 %<>% ps_mutate(dataset = "FPKM")
dataset2 %<>% ps_mutate(dataset = "TPM")

sample_names(dataset1) <- paste(sample_names(dataset1), "FPKM", sep="_")
sample_names(dataset2) <- paste(sample_names(dataset2), "TPM", sep="_")

combined <- phyloseq::merge_phyloseq(
  dataset1 %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>% tax_agg("ARGCluster90") %>% ps_get(),
  dataset2 %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>% tax_agg("ARGCluster90") %>% ps_get()
)

combined %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>%
  comp_barplot("ARGCluster90", facet_by = "dataset", n_taxa = 12, palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
               other_name = "Other ARG", merge_other = F, sample_order = "asis") +
  coord_flip() + ggtitle("FPKM vs TPM relative abundance of ARGs")

# uitprobeersel

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>% ps_arrange(Age) %>% 
  comp_barplot(tax_level = "ARGCluster90", n_taxa = 12, sample_order = "asis", palette = colorRampPalette(brewer.pal(8,"Accent"))(13)) +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age", ., "days")),
    scales = "fixed") +
  coord_flip()

# relative abundance of AMR primary class shown by age and farm

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of primary antimicrobial class by farm and age (FPKM)")

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of primary antimicrobial class by farm and age (TPM)")

# uitprobeersel (lelijk)
Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Farm2") +
  facet_grid(
    cols = vars(Age), rows = vars(AB),
    labeller = labeller(.cols = as_labeller(~ paste("Age = ", .))),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Farm", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of ARGs by farm and age (FPKM)")

# Relative abundance for both stable and antibiotics used

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of primary antimicrobial class by stable and antibiotics used (FPKM)")

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of primary antimicrobial class by stable and antibiotics used (TPM)")


# Same plots but with ARG

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "ARGCluster90", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of ARG by farm and age (FPKM)")

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "ARGCluster90", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of ARG by farm and age (TPM)")


Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "ARGCluster90", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of ARG by stable and antibiotics used (FPKM)")


Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "ARGCluster90", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of ARG by stable and antibiotics used (TPM)")


# 2_57 is a big outlier, for some reason barely has tet(W) and other common ARGs, and therefore a lot more other abundant ARGs

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)")) %>% ps_filter(SampleID == "7.F1S2.21.09") %>% ps_calc_dominant(rank = "ARGCluster90") %>% 
  comp_barplot(tax_level = "ARGCluster90", n_taxa = 12, palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
               other_name = "Other ARG", merge_other = F) +
  coord_flip()



# factorizing variables as not to create problems with visualisation later down the line
sample_data(Rps)$Cluster = as.factor(sample_data(Rps)$Cluster)
sample_data(Rps)$FlockSize = as.factor(sample_data(Rps)$FlockSize)
sample_data(Rps)$AgeParentStock = as.factor(sample_data(Rps)$AgeParentStock)
sample_data(Rps)$Age = as.factor(sample_data(Rps)$Age)
sample_data(Rps)$LibraryNumber = as.factor(sample_data(Rps)$LibraryNumber)

# repeat for TPM
sample_data(Rps_tpm)$Cluster = as.factor(sample_data(Rps)$Cluster)
sample_data(Rps_tpm)$FlockSize = as.factor(sample_data(Rps)$FlockSize)
sample_data(Rps_tpm)$AgeParentStock = as.factor(sample_data(Rps)$AgeParentStock)
sample_data(Rps_tpm)$Age = as.factor(sample_data(Rps)$Age)
sample_data(Rps_tpm)$LibraryNumber = as.factor(sample_data(Rps)$LibraryNumber)

# declutter R environment by removing objects that no longer serve a purpose
rm(meta_data, firm_names, random_tree, random_tree2, ps_prim, ps_rel_abund, dataset1, dataset2, combined)
