library(tidyverse)

reads_k2 = read.table("total_mapped_reads_kraken2.tab", sep = "\t", header = FALSE)
reads_mp = read.table("total_mapped_reads_metaphlan.tab", sep = "\t", header = FALSE)
colnames(reads_k2) <- c("Sample", "Reads_k2")
colnames(reads_mp) <- c("Sample", "Reads_mp")


mapped_reads = left_join(reads_k2, reads_mp, by="Sample")

ggplot(mapped_reads, aes(x = Sample, y = Reads_k2, group = 1)) + geom_line()
ggplot(mapped_reads, aes(x = Sample, y = Reads_mp, group = 1)) + geom_line()

t.test(reads_k2$Reads_k2, y = reads_mp$Reads_mp, paired = TRUE)

# declutter R environment by removing objects that no longer serve a purpose
rm(reads_k2, reads_mp, mapped_reads) 