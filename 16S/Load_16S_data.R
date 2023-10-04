##  For Martijn Melissen - April 2023
### ImportData subset n=120 field study, n=6 per farm
#### n=120 METAGENOMIC and 16S

#### background:  https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html   #Introbackground
####              https://microbiome.github.io/OMA/    # chapter 10, clustering
#### Preprocces NG-tax:
####              https://www.frontiersin.org/articles/10.3389/fgene.2019.01366/full
####              https://f1000research.com/articles/5-1791

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

# temp solution:

sample_names(subset16S) = sample_data(subset16S)$Sample_Unique
sample_names(subset16S)



# Stable "Farm2R1S1"  has the three lowest sampling depths of the dataset, the other nine samples are fairly average
subset16S %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% sample_sums() %>% sort()

# absolute abundances
plot_bar(subset16S, fill="Phylum", title = "Absolute abundances per sample")

# Amount of different taxa present.
sort(table(tax_table(subset16S)[, "Phylum"]))
sort(table(tax_table(subset16S)[, "Order"]))
sort(table(tax_table(subset16S)[, "Family"]))


# for plotting abundances of specific stables
subset16S %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% plot_bar(fill="Phylum")


# visualisation on Phyla, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- phyloseq::tax_glom(subset16S, "Phylum")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Phylum"]

psmelt(ps_prim) %>%
  ggplot(data = ., aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

# Check the amount of uniqe genera in samples which have and have not been treated with antibiotics
subset16S %>% ps_filter(AB == "no") %>% get_taxa_unique("Genus") # 93 different genera for non AB treated
subset16S %>% ps_filter(AB == "yes") %>% get_taxa_unique("Genus") # 74 different genera for AB treated
subset16S %>% get_taxa_unique("Genus") # 93 different genes in total, which are all present in non-treated


# Plots of relative abundances, fixing some genes that are clustered in the data twice, showing top 12 taxa and others are clustered

# relative abundance of Phyla shown by age and farm

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
  ggtitle("Relative abundance of genera by farm and age")



# Relative abundance for both stable and antibiotics used

subset16S %>% 
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

subset16S %>%
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


subset16S %>% 
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


# figures with other category

subset16S.rel <- microbiome::transform(subset16S, "compositional")

plot_composition(subset16S) + theme(legend.position = "bottom") +
  scale_fill_brewer("Phylum", palette = "Paired") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))


plot_composition(subset16S.rel, x.label = "Id") + theme(legend.position = "bottom") +
  scale_fill_brewer("Phylum", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))

plot_composition(subset16S.rel)

# factorizing variables as not to create problems with visualization later down the line
sample_data(subset16S)$Cluster = as.factor(sample_data(subset16S)$Cluster)
sample_data(subset16S)$FlockSize = as.factor(sample_data(subset16S)$FlockSize)
sample_data(subset16S)$AgeParentStock = as.factor(sample_data(subset16S)$AgeParentStock)
sample_data(subset16S)$Age = as.factor(sample_data(subset16S)$Age)
sample_data(subset16S)$LibraryNumber = as.factor(sample_data(subset16S)$LibraryNumber)

# declutter R environment by removing objects that no longer serve a purpose
rm(pseq, ps, ps_prim, subset, treefile)




# stupid plot

ps1.com <- subsetMG

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

subset16S.rel <- microbiome::transform(ps1.com.fam , "compositional")


plot_composition(subset16S.rel, x.label = "Id") + theme(legend.position = "bottom") +
  scale_fill_brewer("Phylum", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))

