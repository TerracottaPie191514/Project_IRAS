##  For Martijn Melissen - April 2023
### ImportData subset n=120 field study, n=6 per farm
#### n=120 METAGENOMIC and 16S

#### background:  https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html   #Introbackground
####              https://microbiome.github.io/OMA/    # chapter 10, clustering
#### Preprocces NG-tax:
####              https://www.frontiersin.org/articles/10.3389/fgene.2019.01366/full
####              https://f1000research.com/articles/5-1791

#### Load packages
library(phyloseq) # Data analysis and visualisation, also the basis of data object.
library(DT) # Interactive tables in html and markdown.
library(data.table) # Giving overview of data.
library(tidyverse) # Data handling and much more.
library(readxl) # Reading in excel files.
library(ape) # Phylogenetic package, used for creating random trees and as dependency for other packages.
library(magrittr) # Data handling, specifically assignment pipes.
library(microViz) # Both analysis and visualisation.
library(plyr) # to apply functions, transform data.
library(microbiome) # For data analysis and visualisation, reading phyloseq object.


### create phyloseq object
pseq <- read_phyloseq(otu.file= "ASV.biom1",
                      taxonomy.file = NULL,
                      metadata.file = "MetaData.csv",
                      type="biom", sep =";" )


treefile <- read_tree("all_asvTREE.tree")
ps <-merge_phyloseq(pseq, treefile)
ps # 180 samples

sort(sample_sums(ps))


### overview data
datatable(tax_table(ps))

### remove some contamination to filter out plant and eukaryote data that had chloroplast and mitochondrium bacterial dna
subset <- subset_taxa(ps, Domain !="NA")
subset <- subset_taxa(subset,Family !="f__Mitochondria=*")
subset <- subset_taxa(subset,Family !="f__Mitochondria")
subset <- subset_taxa(subset, Order !="o__Chloroplast")
subset <- subset_taxa(subset, Domain!="k__Archaea")

### remove taxa with zeros
subset <- prune_taxa(taxa_sums(subset) > 0, subset)

### subset phyloseq object n=120 metagenomics data
subset16S  <- subset_samples(subset,  Metagenomics == "yes" )    #n=120
subset16S <- prune_taxa(taxa_sums(subset16S) > 0, subset16S)
subset16S # 120 samples

# cleaning out all kinds of overlapping names from taxonomy table, removing ~*, =* and =<empty>, this helps to avoid problems with gene abundance later down the line
subset16S@tax_table = gsub("=\\*|~\\*|\\*|<empty>","",subset16S@tax_table)

# overview data
datatable(tax_table(subset16S))
rank_names(subset16S) # Shows classes and ARGs
sort(get_taxa_unique(subset16S, "Genus")) # Shows unique genera
sort(sample_sums(subset16S)) # Amount of unique taxa"per sample, the min is 46731 and max 393697, which is within a factor 10 difference
summary(sample_sums(subset16S)) # summary of the sampling depths
sample_variables(subset16S) # metadata variables

# Rewriting sampleIDs as sample_unique rownames to align with the other datasets

sample_names(subset16S) = sample_data(subset16S)$Sample_Unique
sample_names(subset16S)


# Stable "Farm2R1S1"  has the three lowest sampling depths of the dataset, the other nine samples are fairly average
subset16S %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% sample_sums() %>% sort()


# Amount of different taxa present.
sort(table(tax_table(subset16S)[, "Phylum"]))
sort(table(tax_table(subset16S)[, "Order"]))
sort(table(tax_table(subset16S)[, "Family"]))


# factorizing variables as not to create problems with visualization later down the line
sample_data(subset16S)$Cluster = as.factor(sample_data(subset16S)$Cluster)
sample_data(subset16S)$FlockSize = as.factor(sample_data(subset16S)$FlockSize)
sample_data(subset16S)$AgeParentStock = as.factor(sample_data(subset16S)$AgeParentStock)
sample_data(subset16S)$Age = as.factor(sample_data(subset16S)$Age)
sample_data(subset16S)$LibraryNumber = as.factor(sample_data(subset16S)$LibraryNumber)

# add stable column with shorter names
sample_data(subset16S)$FarmRoundStable = as.factor(sample_data(subset16S)$FarmRoundStable)
subset16S@sam_data$Stables = revalue(sample_data(subset16S)$FarmRoundStable, c("Farm1R1S1"="Stable1", "Farm1R1S2"="Stable2", "Farm2R1S1"="Stable3", "Farm2R1S2"="Stable4",
                                                                              "Farm2R2S1"="Stable5", "Farm2R2S2"="Stable6", "Farm3R1S1"="Stable7", "Farm3R1S2"="Stable8",
                                                                              "Farm4R1S1"="Stable9", "Farm4R1S2"="Stable10"))
# Shortening agent names
subset16S@sam_data$Cox[subset16S@sam_data$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
subset16S@sam_data$Cox[subset16S@sam_data$Cox == "narasin(monteban)"] = "Monteban"
subset16S@sam_data$Cox[subset16S@sam_data$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

# declutter R environment by removing objects that no longer serve a purpose
rm(pseq, ps, subset, treefile)
