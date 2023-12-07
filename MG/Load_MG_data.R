#### Load packages
library(phyloseq) # Data analysis and visualisation, also the basis of data object.
library(DT) # Interactive tables in html and markdown, overview of data.
library(tidyverse) # Data handling and much more.
library(readxl) # Reading in excel files.
library(ape) # Phylogenetic package, used for creating random trees and as dependency for other packages.
library(magrittr) # Data handling, specifically assignment pipes
library(microViz) # Both analysis and visualisation
library(plyr) # to apply functions, transform data

### Loading a subset of metagenomic data into phyloseq format
subsetMG= import_biom("kraken2_output.biom") # this imports a .biom created by kraken2-biom  containing OTU and tax tables, with [fF]irm_x_x names as sample_names

# We rewrite the sample names to a format filtering out Firm and firm and the first underscore so that it lines up with the column of our meta data
sample_names(subsetMG) = sapply(regmatches(sample_names(subsetMG), regexpr("_", sample_names(subsetMG)), invert = TRUE), "[[", 2) 

# Because the names in both metadata sets do not completely overlap, we need to manually edit one of the samples whose name was not included in the FIRM metadata file
sample_names(subsetMG)[68] = "4_65"

# reading in and combining metadata from 16S and metagenomic origins, adding missing underscores
firm_names = read_excel("./Metagenomic/FIRM_MetaNames.xlsx")
firm_names = firm_names[,-2] # Remove wrongful Raw_data_name column, to avoid confusion

meta_data = read.csv("MetaData.csv", header = TRUE, sep = ",")
meta_data_MG = dplyr::right_join(firm_names, meta_data, by="SampleID")

# using Sample_Unique as rownames so we can match the two sets in phyloseq
rownames(meta_data_MG) = meta_data_MG$Sample_Unique

# now we'll also add in microbial load
microbial_load = read.table("bacterial_load_kraken2.tab", sep = "\t", header = TRUE)
microbial_load$Sample_Unique = sapply(regmatches(microbial_load$Sample_Unique, regexpr("_",microbial_load$Sample_Unique), invert = TRUE), "[[", 2) 
microbial_load$Sample_Unique[68] = "4_65"
meta_data_MG = dplyr::right_join(meta_data_MG, microbial_load, by="Sample_Unique")

# creating tree and making phyloseq components, adding tree and sample data components to phyloseq
set.seed("878") # setting seed for reproducibility purposes
random_tree = rtree(ntaxa(subsetMG), rooted=TRUE, tip.label=taxa_names(subsetMG))
meta_data_MG = sample_data(meta_data_MG)
rownames(meta_data_MG) = meta_data_MG$Sample_Unique
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

# Amount of different taxa present.
sort(table(tax_table(subsetMG)[, "Phylum"]))
sort(table(tax_table(subsetMG)[, "Order"]))
sort(table(tax_table(subsetMG)[, "Family"]))
sort(table(tax_table(subsetMG)[, "Genus"]))


# Check the amount of unique Orders in samples which have and have not been treated with antibiotics
subsetMG %>% ps_filter(AB == "no") %>% get_taxa_unique("Order") # 200 different orders for non AB treated
subsetMG %>% ps_filter(AB == "yes") %>% get_taxa_unique("Order") # 159 different orders for AB treated
subsetMG %>% get_taxa_unique("Order") # 203 different order in total, so 3 orders are not found in non AB

# Check the amount of unique Species in samples which have and have not been treated with antibiotics
subsetMG %>% ps_filter(AB == "no") %>% get_taxa_unique("Species") # 4464 different orders for non AB treated
subsetMG %>% ps_filter(AB == "yes") %>% get_taxa_unique("Species") # 2347 different orders for AB treated
subsetMG %>% get_taxa_unique("Species") # 4706 different order in total, so 242 species are not found in non AB

# Check the amount of unique taxa in samples which have and have not been treated with antibiotics
subsetMG %>% ps_filter(AB == "no")  # 6014 different taxa for non AB treated
subsetMG %>% ps_filter(AB == "yes") # 3148 different taxa for AB treated
# 6355 different taxa in total, so 341 taxa are not found in non AB


# Stable "Farm2R1S1"  has the three lowest sampling depths of the dataset, the other nine samples are fairly average 
subsetMG %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% sample_sums() %>% sort()

# factorizing variables as not to create problems with visualization later down the line
sample_data(subsetMG)$Cluster = as.factor(sample_data(subsetMG)$Cluster)
sample_data(subsetMG)$FlockSize = as.factor(sample_data(subsetMG)$FlockSize)
sample_data(subsetMG)$AgeParentStock = as.factor(sample_data(subsetMG)$AgeParentStock)
sample_data(subsetMG)$Age = as.factor(sample_data(subsetMG)$Age)
sample_data(subsetMG)$LibraryNumber = as.factor(sample_data(subsetMG)$LibraryNumber)

# add stable column with shorter names
sample_data(subsetMG)$FarmRoundStable = as.factor(sample_data(subsetMG)$FarmRoundStable)
subsetMG@sam_data$Stables = revalue(sample_data(subsetMG)$FarmRoundStable, c("Farm1R1S1"="Stable1", "Farm1R1S2"="Stable2", "Farm2R1S1"="Stable3", "Farm2R1S2"="Stable4",
                                                                              "Farm2R2S1"="Stable5", "Farm2R2S2"="Stable6", "Farm3R1S1"="Stable7", "Farm3R1S2"="Stable8",
                                                                              "Farm4R1S1"="Stable9", "Farm4R1S2"="Stable10"))
# Shortening agent names
subsetMG@sam_data$Cox[subsetMG@sam_data$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
subsetMG@sam_data$Cox[subsetMG@sam_data$Cox == "narasin(monteban)"] = "Monteban"
subsetMG@sam_data$Cox[subsetMG@sam_data$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

# declutter R environment by removing objects that no longer serve a purpose
rm(meta_data, firm_names, meta_data_MG, microbial_load, random_tree) 
