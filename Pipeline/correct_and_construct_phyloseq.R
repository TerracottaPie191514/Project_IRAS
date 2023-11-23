
# count table has samples as row and ARGs as columns
counts <- read.csv( "results/Resfinder_mapped_AllSamples_Wide.tab",
                    sep = "\t", header = TRUE, check.names = FALSE, row.names = 1 )

# Taxonomy of each "OTU"
# Has ARG as row and tax+len as columns
ARGtable <- read.csv( "./Resfinder/ResfinderClustersFinal_Resfinder20200127_inclGeneLen+pheno_20201206d_TAXONOMY+ARGlen.tab",
                      header = TRUE, sep = "\t", row.names = 1 )

# Keep just "taxonomy" for the phyloseq
tax <- ARGtable %>% select(-ARGlength)


# Metadata (only samplenames for now
meta <- data.frame( "SampleID" <- rownames(counts) )
row.names(meta) <- meta$X.SampleID.....rownames.counts. # this way we retain the Sample_ID
#Only for ARGs in the dataset
args <- colnames(counts)

# Find gene length in ARG info table and correct each column of an ARG with its corresponding gene length to get fragments per kilobase (FPK)
fpk <- counts
for( arg in args ) {
  fpk[,arg] <- counts[,arg] / ARGtable[arg,"ARGlength"] * 1000
}

# Read in MetaPhlan total bacterial reads mapped to sample
mapped_reads <- read.csv( "results/total_mapped_reads_metaphlan.tab",
                          sep = "\t", header = FALSE, check.names = FALSE, row.names = 1 )
# Do the same for Kraken2 reads
mapped_reads_k2 <- read.csv( "results/total_mapped_reads_kraken2.tab",
                             sep = "\t", header = FALSE, check.names = FALSE, row.names = 1 )

# For looping over samples
samples <- rownames(mapped_reads)

# Repeat for kraken2

samples_k2 <- rownames(mapped_reads_k2)

# Dividing fpk by the scaling factor (sum of all fpk values in samples divided by a million) to get transcripts per kilobase million (TPM)
tpm <- counts
for( sample in samples ) {
  tpm[sample,] <- fpk[sample,] / (sum(fpk[sample,]) / 1000000)
}

# Use Metaphlan total mapped reads divided by a million as scaling factor, then divide counts by scaling factor to get fragments per millions (FPM)
fpm = counts
for( sample in samples ) {
  fpm[sample,] <- counts[sample,] / mapped_reads[sample,] * 1000000
}

# Correct for gene length as well to get fragments per kilobase million (FPKM)
fpkm = fpm
for( arg in args ) {
  fpkm[,arg] <- fpm[,arg] / ARGtable[arg,"ARGlength"] * 1000
}

# Now, we'll repeat creating fpkm for kraken2 data

fpm_k2 = counts
for( sample in samples_k2 ) {
  fpm_k2[sample,] <- counts[sample,] / mapped_reads_k2[sample,] * 1000000
}

# Correct for gene length as well to get fragments per kilobase million (FPKM)
fpkm_k2 = fpm_k2
for( arg in args ) {
  fpkm_k2[,arg] <- fpm_k2[,arg] / ARGtable[arg,"ARGlength"] * 1000
}



#Creating phyloseq objects and exporting the data

library(microbiome)
library(magrittr)
tax_mat <- as.matrix(tax)

fpkm %<>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) 

fpkm %<>% remove_rownames %>% column_to_rownames(var="name")

meta$X.SampleID.....rownames.counts.
OTU = otu_table(fpkm, taxa_are_rows=T)
TAX = tax_table(tax_mat)
taxa_names(TAX)
taxa_names(OTU)
pseq = phyloseq(OTU,TAX,sam_names=meta)

# Repeat for kraken2

fpkm_k2 %<>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) 

fpkm_k2 %<>% remove_rownames %>% column_to_rownames(var="name")
OTU_k2 = otu_table(fpkm_k2, taxa_are_rows=T)
pseq_k2 = phyloseq(OTU_k2,TAX,sam_names=meta)

tpm %<>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) 

tpm %<>% remove_rownames %>% column_to_rownames(var="name")

OTU2 = otu_table(tpm, taxa_are_rows=T)
pseq_tpm = phyloseq(OTU2,TAX,sam_names=meta)


OTU1 = as(otu_table(pseq), "matrix")
if(taxa_are_rows(pseq)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

saveRDS(pseq, "Phyloseq")
saveRDS(pseq_tpm, "Phyloseq_tpm")
saveRDS(pseq_k2, "./data/Phyloseq_k2")