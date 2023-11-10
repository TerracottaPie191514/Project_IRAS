### loading two subsets of metagenomic data into phyloseq format
externalps = import_biom("extern.biom")
internalps = import_biom("intern.biom")

externalps
internalps # internal has 102 less taxa, but this is unfiltered

# set Rank names
colnames(tax_table(externalps)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(internalps)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# filtering
externalps <- subset_taxa(externalps, Domain!="k__Archaea")
externalps <- subset_taxa(externalps, Domain!="k__Viruses")
externalps <- subset_taxa(externalps, Domain!="k__Eukaryota")

internalps <- subset_taxa(internalps, Domain!="k__Archaea")
internalps <- subset_taxa(internalps, Domain!="k__Viruses")
internalps <- subset_taxa(internalps, Domain!="k__Eukaryota")

externalps
internalps # now, there are 17 more taxa in internal
datatable(tax_table(externalps))

sample_names(internalps) = c("Firm_2_42", "Firm_2_47", "Firm_2_48", "Firm_2_49", "Firm_2_50", "Firm_2_51", "Firm_2_52", 
                             "Firm_2_56", "Firm_2_57", "Firm_2_58")

# Visualize in the same plot
dataset1 =  ps_filter(externalps)
dataset2 =  ps_filter(internalps)

dataset1 %<>% ps_mutate(dataset = "external")
dataset2 %<>% ps_mutate(dataset = "internal")

sample_names(dataset1) <- paste(sample_names(dataset1), "external", sep="_")
sample_names(dataset2) <- paste(sample_names(dataset2), "internal", sep="_")

combined <- phyloseq::merge_phyloseq(
  dataset1 %>% tax_fix(unknowns = c("g__")) %>%   tax_agg("Genus") %>% ps_get(),
  dataset2 %>% tax_fix(unknowns = c("g__")) %>% tax_agg("Genus") %>% ps_get()
)

combined %>% tax_fix(unknowns = c("g__")) %>%
  comp_barplot("Genus", facet_by = "dataset", n_taxa = 12, palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
               other_name = "Other Genera", merge_other = F, sample_order = "asis") +
  coord_flip() + ggtitle("External vs internal dataset")

# phylum
combined <- phyloseq::merge_phyloseq(
  dataset1  %>% tax_agg("Phylum") %>% ps_get(),
  dataset2  %>% tax_agg("Phylum") %>% ps_get()
)

combined %>% tax_fix(unknowns = c("g__")) %>%
  comp_barplot("Phylum", facet_by = "dataset", n_taxa = 4, palette = colorRampPalette(brewer.pal(5,"Accent"))(5),
               other_name = "Other Phyla", merge_other = F, sample_order = "asis") +
  coord_flip() + ggtitle("External vs internal dataset")

# procrustes


#Procrustes analyses

PCoA_BC_int = ordinate(internalps, "PCoA") 
PCoA_BC_ext = ordinate(externalps, "PCoA") 
procrustes = protest(PCoA_BC_int$vectors, PCoA_BC_ext$vectors)

plot_data <- data.frame(
  int_PC1 = procrustes$X[, 1],
  int_PC2 = procrustes$X[, 2],
  ext_PC1 = procrustes$Yrot[, 1],
  ext_PC2 = procrustes$Yrot[, 2]
)

sample_pairs = data.frame("int_sample"=rownames(PCoA_BC_int$vectors), "ext_sample"=rownames(PCoA_BC_ext$vectors)) 

ggplot(plot_data) +
  geom_point(aes(x=ext_PC1, y=ext_PC2), color = "green") +
  geom_point(aes(x=int_PC1, y=int_PC2), color = "blue") +
  geom_segment(aes(x=int_PC1,y=int_PC2,xend=ext_PC1,yend=ext_PC2),arrow=arrow(type = "closed", length=unit(0.2,"cm"))) +
  labs(title = "Procrustes Plot internal vs external server") # there is complete overlap, only one colour shows

# declutter R environment by removing objects that no longer serve a purpose
rm(externalps, internalps, sample_pairs, plot_data, procrustes, PCoA_BC_int, PCoA_BC_ext, dataset1, dataset2, combined) 
