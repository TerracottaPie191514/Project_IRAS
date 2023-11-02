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
library(plyr) # to apply functions, transform data


### loading a subset of metagenomic data into phyloseq format
Rps= readRDS("Phyloseq_k2") # this reads a pre-existing phyloseq object containing OTU and tax tables, with [fF]irm_x_x names as sample_names, kraken2 count data
Rps_mp= readRDS("Phyloseq") # reads in the data with count data corrected with metaphlan bacterial counts
Rps_tpm = readRDS("Phyloseq_tpm") # also read in TPM data instead of FPKM

#We rewrite the sample names to a format filtering out Firm and firm and the first underscore so that it lines up with the column of our meta data
sample_names(Rps) = sapply(regmatches(sample_names(Rps), regexpr("_", sample_names(Rps)), invert = TRUE), "[[", 2) 
sample_names(Rps_mp) = sapply(regmatches(sample_names(Rps_mp), regexpr("_", sample_names(Rps_mp)), invert = TRUE), "[[", 2) 
sample_names(Rps_tpm) = sapply(regmatches(sample_names(Rps_tpm), regexpr("_", sample_names(Rps_tpm)), invert = TRUE), "[[", 2) 

# Because the names in both metadata sets do not completely overlap, we need to manually edit one of the samples whose name was not included in the FIRM metadata file
sample_names(Rps)[68] = "4_65"
sample_names(Rps_mp)[68] = "4_65"
sample_names(Rps_tpm)[68] = "4_65"


# reading in and combining metadata from 16S and metagenomic origins, adding missing underscores
firm_names = read_excel("./Metagenomic/FIRM_MetaNames.xlsx")
firm_names = firm_names[,-2] # Remove wrongful Raw_data_name column, to avoid confusion

meta_data = read.csv("MetaData.csv", header = TRUE, sep = ",")
meta_data_R = dplyr::right_join(firm_names, meta_data, by="SampleID")

# using Sample_Unique as rownames so we can match the two sets in phyloseq
rownames(meta_data_R) = meta_data_R$Sample_Unique

# now we'll also add in microbial load
microbial_load = read.table("bacterial_load_kraken2.tab", sep = "\t", header = TRUE)
microbial_load$Sample_Unique = sapply(regmatches(microbial_load$Sample_Unique, regexpr("_",microbial_load$Sample_Unique), invert = TRUE), "[[", 2) 
microbial_load$Sample_Unique[68] = "4_65"
meta_data_R = dplyr::right_join(meta_data_R, microbial_load, by="Sample_Unique")

# creating tree and making phyloseq components, adding tree and sample data components to phyloseq
set.seed("877") # setting seed for reproducibility purposes
random_tree = rtree(ntaxa(Rps), rooted=TRUE, tip.label=taxa_names(Rps))
meta_data_R = sample_data(meta_data_R)
rownames(meta_data_R) = meta_data_R$Sample_Unique
Rps = merge_phyloseq(Rps, meta_data_R, random_tree)
class(Rps)

# repeat for mp

random_tree2 = rtree(ntaxa(Rps_mp), rooted=TRUE, tip.label=taxa_names(Rps_mp))
Rps_mp = merge_phyloseq(Rps_mp, meta_data_R, random_tree2)

# repeat for tpm

random_tree3 = rtree(ntaxa(Rps_tpm), rooted=TRUE, tip.label=taxa_names(Rps_tpm))
Rps_tpm = merge_phyloseq(Rps_tpm, meta_data_R, random_tree3)

# we will scale the otu table by a factor 1000 to preserve data before the rounding step 

otu_table(Rps) = otu_table(Rps) * 1000
otu_table(Rps_mp) = otu_table(Rps_mp) * 1000


#rounding the "counts" for phyloseq to function

otu_table(Rps) = otu_table(round(as((otu_table(Rps)), "matrix")), taxa_are_rows(Rps))
otu_table(Rps_mp) = otu_table(round(as((otu_table(Rps_mp)), "matrix")), taxa_are_rows(Rps_mp))
otu_table(Rps_tpm) = otu_table(round(as((otu_table(Rps_tpm)), "matrix")), taxa_are_rows(Rps_tpm))


# sample_data(Rps)$Sample_Unique = sample_names(Rps)


# overview data
datatable(tax_table(Rps))
rank_names(Rps) # Shows classes and ARGs
sort(get_taxa_unique(Rps, "AMR_class_primary")) # Shows primary AMR classes
sort(sample_sums(Rps)) # Amount of unique "taxa" per sample, the min is 1365913 and max 44483138, which is a big difference
summary(sample_sums(Rps)) # summary of the sampling depths
summary(sample_sums(Rps_mp)) # there are big differences between kraken2 and metaphlan counts data, with metaphlan having lower min and higher max and a much higer mean and median
sample_variables(Rps) # metadata variables
taxa_names(Rps) # ARGs

# Stable "Farm2R1S1"  has the five lowest sampling depths of the dataset, some other samples are also very low, but others not so much 
Rps %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% sample_sums() %>% sort()

# For metaphlan data, Stable "Farm2R1S1 has 11/12 lowest sampling depths of the dataset, and there are only 2 other samples with a lower depth than the 12th  sample
Rps_mp %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% sample_sums() %>% sort()


# factorizing variables as not to create problems with visualisation later down the line
sample_data(Rps)$Cluster = as.factor(sample_data(Rps)$Cluster)
sample_data(Rps)$FlockSize = as.factor(sample_data(Rps)$FlockSize)
sample_data(Rps)$AgeParentStock = as.factor(sample_data(Rps)$AgeParentStock)
sample_data(Rps)$Age = as.factor(sample_data(Rps)$Age)
sample_data(Rps)$LibraryNumber = as.factor(sample_data(Rps)$LibraryNumber)

# repeat for MP
sample_data(Rps_mp)$Cluster = as.factor(sample_data(Rps_mp)$Cluster)
sample_data(Rps_mp)$FlockSize = as.factor(sample_data(Rps_mp)$FlockSize)
sample_data(Rps_mp)$AgeParentStock = as.factor(sample_data(Rps_mp)$AgeParentStock)
sample_data(Rps_mp)$Age = as.factor(sample_data(Rps_mp)$Age)
sample_data(Rps_mp)$LibraryNumber = as.factor(sample_data(Rps_mp)$LibraryNumber)

# repeat for TPM
sample_data(Rps_tpm)$Cluster = as.factor(sample_data(Rps_tpm)$Cluster)
sample_data(Rps_tpm)$FlockSize = as.factor(sample_data(Rps_tpm)$FlockSize)
sample_data(Rps_tpm)$AgeParentStock = as.factor(sample_data(Rps_tpm)$AgeParentStock)
sample_data(Rps_tpm)$Age = as.factor(sample_data(Rps_tpm)$Age)
sample_data(Rps_tpm)$LibraryNumber = as.factor(sample_data(Rps_tpm)$LibraryNumber)

# add stable column with shorter names
sample_data(Rps)$FarmRoundStable = as.factor(sample_data(Rps)$FarmRoundStable)
Rps@sam_data$Stable = revalue(sample_data(Rps)$FarmRoundStable, c("Farm1R1S1"="Stable1", "Farm1R1S2"="Stable2", "Farm2R1S1"="Stable3", "Farm2R1S2"="Stable4",
                                                                              "Farm2R2S1"="Stable5", "Farm2R2S2"="Stable6", "Farm3R1S1"="Stable7", "Farm3R1S2"="Stable8",
                                                                              "Farm4R1S1"="Stable9", "Farm4R1S2"="Stable10"))
# Shortening agent names
Rps@sam_data$Cox[Rps@sam_data$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
Rps@sam_data$Cox[Rps@sam_data$Cox == "narasin(monteban)"] = "Monteban"
Rps@sam_data$Cox[Rps@sam_data$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

# repeat
sample_data(Rps_mp)$FarmRoundStable = as.factor(sample_data(Rps_mp)$FarmRoundStable)
Rps_mp@sam_data$Stable = revalue(sample_data(Rps_mp)$FarmRoundStable, c("Farm1R1S1"="Stable1", "Farm1R1S2"="Stable2", "Farm2R1S1"="Stable3", "Farm2R1S2"="Stable4",
                                                                              "Farm2R2S1"="Stable5", "Farm2R2S2"="Stable6", "Farm3R1S1"="Stable7", "Farm3R1S2"="Stable8",
                                                                              "Farm4R1S1"="Stable9", "Farm4R1S2"="Stable10"))
Rps_mp@sam_data$Cox[Rps_mp@sam_data$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
Rps_mp@sam_data$Cox[Rps_mp@sam_data$Cox == "narasin(monteban)"] = "Monteban"
Rps_mp@sam_data$Cox[Rps_mp@sam_data$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

sample_data(Rps_tpm)$FarmRoundStable = as.factor(sample_data(Rps_tpm)$FarmRoundStable)
Rps_tpm@sam_data$Stables = revalue(sample_data(Rps_tpm)$FarmRoundStable, c("Farm1R1S1"="Stable1", "Farm1R1S2"="Stable2", "Farm2R1S1"="Stable3", "Farm2R1S2"="Stable4",
                                                                              "Farm2R2S1"="Stable5", "Farm2R2S2"="Stable6", "Farm3R1S1"="Stable7", "Farm3R1S2"="Stable8",
                                                                              "Farm4R1S1"="Stable9", "Farm4R1S2"="Stable10"))
Rps_tpm@sam_data$Cox[Rps_tpm@sam_data$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
Rps_tpm@sam_data$Cox[Rps_tpm@sam_data$Cox == "narasin(monteban)"] = "Monteban"
Rps_tpm@sam_data$Cox[Rps_tpm@sam_data$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"
                       
# declutter R environment by removing objects that no longer serve a purpose
rm(meta_data, meta_data_R, firm_names, microbial_load, random_tree, random_tree2, random_tree3) 
