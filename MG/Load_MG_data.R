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
subsetMG= import_biom("kraken2_output.biom") # this imports a .biom created by kraken2-biom  containing OTU and tax tables, with [fF]irm_x_x names as sample_names

#We rewrite the sample names to a format filtering out Firm and firm and the first underscore so that it lines up with the column of our meta data
sample_names(subsetMG) = sapply(regmatches(sample_names(subsetMG), regexpr("_", sample_names(subsetMG)), invert = TRUE), "[[", 2) 

# Because the names in both metadata sets do not completely overlap, we need to manually edit one of the samples whose name was not included in the FIRM metadata file
sample_names(subsetMG)[68] = "4_65"

# reading in and combining metadata from 16S and metagenomic origins, adding missing underscores
firm_names = read_excel("./Metagenomic/FIRM_MetaNames.xlsx")
firm_names = firm_names[,-2] # Remove wrongful Raw_data_name column, to avoid confusion

meta_data = read.csv("MetaData.csv", header = TRUE, sep = ",")
meta_data_MG = dplyr::right_join(firm_names, meta_data, by="SampleID")

# using Sample_Unique as rownames so we can match the two sets in phyloseq
meta_data_MG %<>% remove_rownames %>% column_to_rownames(var="Sample_Unique")

# creating tree and making phyloseq components, adding tree and sample data components to phyloseq
random_tree = rtree(ntaxa(subsetMG), rooted=TRUE, tip.label=taxa_names(subsetMG))
meta_data_MG = sample_data(meta_data_MG)
subsetMG = merge_phyloseq(subsetMG, meta_data_MG, random_tree)
class(subsetMG)

# set Rank names
colnames(tax_table(subsetMG)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(subsetMG) 

### overview data
datatable(tax_table(subsetMG))
subsetMG # 7058 taxa

# filter out non bacterial domains (no chloroplast, mitochondrial "taxa" present)
subsetMG <- subset_taxa(subsetMG, Domain!="k__Archaea")
subsetMG <- subset_taxa(subsetMG, Domain!="k__Viruses")
subsetMG <- subset_taxa(subsetMG, Domain!="k__Eukaryota")

subsetMG # 6355 taxa

# Stable "Farm2R1S1"  has the three lowest sampling depths of the dataset, the other nine samples are fairly averagehttp://127.0.0.1:30549/graphics/plot_zoom_png?width=1504&height=963
subsetMG %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% sample_sums() %>% sort()

# Amount of different taxa present.
sort(table(tax_table(subsetMG)[, "Phylum"]))
sort(table(tax_table(subsetMG)[, "Order"]))
sort(table(tax_table(subsetMG)[, "Family"]))

# absolute abundances, since there are a lot of phyla (43), we will only include phyla which are present at least 50 times 


#HIER NOG ALLES SAMEN VOEGEN WAT OTHERS IS

subsetMG %>% subset_taxa(Phylum == c("p__Pseudomonadota", "p__Actinomycetota", "p__Bacillota", "p__Bacteroidota",
                                     "p__Thermodesulfobacteriota", "p__Cyanobacteriota", "p__Campylobacterota", 
                                     "p__Mycoplasmatota", "p__Planctomycetota", "p__Spirochaetota")) %>% 
  plot_bar(fill="Phylum", title = "Absolute abundances per sample")



# for plotting abundances of specific stables
subsetMG %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% plot_bar(fill="Phylum")


# visualisation on Phyla, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- phyloseq::tax_glom(subsetMG, "Phylum")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Phylum"]

psmelt(ps_prim) %>%
  ggplot(aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

# Check the amount of unique Orders in samples which have and have not been treated with antibiotics
subsetMG %>% ps_filter(AB == "no") %>% get_taxa_unique("Order") # 200 different orders for non AB treated
subsetMG %>% ps_filter(AB == "yes") %>% get_taxa_unique("Order") # 159 different orders for AB treated
subsetMG %>% get_taxa_unique("Order") # 203 different order in total, so 3 orders are not found in non AB


# Plots of relative abundances, fixing some genes that are clustered in the data twice, showing top 12 taxa and others are clustered

# relative abundance of Phyla shown by age and farm

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

subsetMG %>% 
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
  ggtitle("Relative abundance of genera by stable and antibiotics used")


# Same plots but with Order

subsetMG %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "Order", bar_width = 0.7, sample_order = "asis", 
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


subsetMG %>% 
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "Order", bar_width = 0.7, sample_order = "asis", 
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


# relabundance with other category
subsetMG.rel <- microbiome::transform(subsetMG, "compositional")

plot_composition(subsetMG.rel, sample.sort = "FarmRoundStable", x.label = "Id") + theme(legend.position = "bottom") +
 scale_fill_brewer("Phylum", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun)

# factorizing variables as not to create problems with visualization later down the line
sample_data(subsetMG)$Cluster = as.factor(sample_data(subsetMG)$Cluster)
sample_data(subsetMG)$FlockSize = as.factor(sample_data(subsetMG)$FlockSize)
sample_data(subsetMG)$AgeParentStock = as.factor(sample_data(subsetMG)$AgeParentStock)
sample_data(subsetMG)$Age = as.factor(sample_data(subsetMG)$Age)
sample_data(subsetMG)$LibraryNumber = as.factor(sample_data(subsetMG)$LibraryNumber)

# declutter R environment by removing objects that no longer serve a purpose
rm(pseq, ps, ps_prim, subset, treefile)