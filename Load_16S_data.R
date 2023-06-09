##  For Martijn Melissen - April 2023
### ImportData subset n=120 field study, n=6 per farm
#### n=120 METAGENOMIC and 16S

#### background:  https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html   #Introbackground
####              https://microbiome.github.io/OMA/    # chapter 10, clustering
#### Preprocces NG-tax:
####              https://www.frontiersin.org/articles/10.3389/fgene.2019.01366/full
####              https://f1000research.com/articles/5-1791

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

### create phyloseq object, https://joey711.github.io/phyloseq/
pseq <- read_phyloseq(otu.file= "ASV.biom1",
                      taxonomy.file = NULL,
                      metadata.file = "MetaData.csv",
                      type="biom", sep =";" )


treefile <- read_tree("all_asvTREE.tree")
ps <-merge_phyloseq(pseq, treefile)
ps

sort(sample_sums(ps))
#> ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6249 taxa and 180 samples ]
#sample_data() Sample Data:       [ 180 samples by 26 sample variables ]
#tax_table()   Taxonomy Table:    [ 6249 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 6249 tips and 6248 internal nodes ]

### 180 samples

### overview data
datatable(tax_table(ps))

### remove some contamination
subset <- subset_taxa(ps, Domain !="NA")
subset <- subset_taxa(subset,Family !="f__Mitochondria=*")
subset <- subset_taxa(subset,Family !="f__Mitochondria")
subset <- subset_taxa(subset, Order !="o__Chloroplast")
#do you know why? filter out plant and eukaryote data that had chloroplast and mitochondrium bacterial dna
subset <- subset_taxa(subset, Domain!="k__Archaea") # can be discussed about it

### remove taxa with zeros
subset <- prune_taxa(taxa_sums(subset) > 0, subset)

### subset phyloseq object n=120 metagenomics data
subsetG  <- subset_samples(subset,  Metagenomics == "yes" )    #n=120
subsetG <- prune_taxa(taxa_sums(subsetG) > 0, subsetG)
subsetG

datatable(tax_table(subsetG))

#subsetG
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1536 taxa and 120 samples ]
#sample_data() Sample Data:       [ 120 samples by 26 sample variables ]
#tax_table()   Taxonomy Table:    [ 1536 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1536 tips and 1535 internal nodes ]

rank_names(subsetG)
sort(get_taxa_unique(subsetG, "Genus"))
sample_sums(subsetG)
sample_variables(subsetG)


taxa_names(subsetG)

subset16S = subsetG
subset16S@tax_table = gsub("=\\*|~\\*|\\*|<empty>","",subsetG@tax_table)

datatable(tax_table(subset16S))



sample_variables(subset16S)

sample_data(subset16S)$Cluster = as.factor(sample_data(subset16S)$Cluster)
sample_data(subset16S)$FlockSize = as.factor(sample_data(subset16S)$FlockSize)
sample_data(subset16S)$AgeParentStock = as.factor(sample_data(subset16S)$AgeParentStock)
sample_data(subset16S)$Age = as.factor(sample_data(subset16S)$Age)
sample_data(subset16S)$LibraryNumber = as.factor(sample_data(subset16S)$LibraryNumber)



unique(sample_data(subsetG)$Cox)