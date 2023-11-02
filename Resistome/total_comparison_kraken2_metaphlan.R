library(tidyverse)

reads_k2 = read.table("total_mapped_reads_kraken2.tab", sep = "\t", header = FALSE)
reads_mp = read.table("total_mapped_reads_metaphlan.tab", sep = "\t", header = FALSE)
colnames(reads_k2) <- c("Sample", "Reads_k2")
colnames(reads_mp) <- c("Sample", "Reads_mp")


mapped_reads = dplyr::left_join(reads_k2, reads_mp, by="Sample")

ggplot(mapped_reads, aes(x = Sample, y = Reads_k2, group = 1)) + geom_line()
ggplot(mapped_reads, aes(x = Sample, y = Reads_mp, group = 1)) + geom_line()

t.test(reads_k2$Reads_k2, y = reads_mp$Reads_mp, paired = TRUE)


mapped_reads$ratio = mapped_reads$Reads_k2 / mapped_reads$Reads_mp
print(mapped_reads$ratio)

ggplot(mapped_reads, aes(x = Sample, y = ratio)) + geom_point() + ylab("Ratio Kraken2/MetaPhlan4")

ggplot(mapped_reads, aes(x = ratio)) + geom_histogram(color="black", fill="white")

ggplot(mapped_reads, aes(x = ratio)) + geom_density()

# bland - altman plot
ggplot(mapped_reads, aes(x = Sample, y = ratio, label = Sample)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(mapped_reads$ratio), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(mapped_reads$ratio) - (1.96 * sd(mapped_reads$ratio)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(mapped_reads$ratio) + (1.96 * sd(mapped_reads$ratio)), colour = "red", size = 0.5) +
  ylab("Ratio Kraken2 / MetaPhlan4") 

# declutter R environment by removing objects that no longer serve a purpose
rm(reads_k2, reads_mp, mapped_reads) 
