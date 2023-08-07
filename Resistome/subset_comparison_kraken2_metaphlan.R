library(tidyverse)

subset_kraken2 = import_biom("FIRM-DNA-subset.json")
total_mapped_reads_kraken2 = sample_sums(subset_kraken2)
total_mapped_reads_kraken2 = data.frame(total_mapped_reads_kraken2)
rownames(total_mapped_reads_kraken2) = c("10_53", "9_39", "9_38", "14_30", "10_52")
total_mapped_reads_kraken2 <- tibble::rownames_to_column(total_mapped_reads_kraken2, "Sample")


total_mapped_reads_metaphlan <- data.frame (Sample  = c("10_53", "9_39", "9_38", "14_30", "10_52"),
                  total_mapped_reads_metaphlan = c("25900713", "8699809", "11619167", "12798500", "17231290")
)
total_mapped_reads_metaphlan$total_mapped_reads_metaphlan = as.numeric(total_mapped_reads_metaphlan$total_mapped_reads_metaphlan)


mapped_reads = left_join(total_mapped_reads_kraken2, total_mapped_reads_metaphlan, by="Sample")


ggplot(mapped_reads, aes(x = Sample, y = total_mapped_reads_kraken2, group = 1)) + geom_line()
ggplot(mapped_reads, aes(x = Sample, y = total_mapped_reads_metaphlan, group = 1)) + geom_line()

t.test(total_mapped_reads_metaphlan$total_mapped_reads_metaphlan, y = total_mapped_reads_kraken2$total_mapped_reads_kraken2, paired = TRUE)


total_mapped_reads_kraken2$total_mapped_reads_kraken2
# 6.F5.S2.CA.21.9.17 = 10_53  metaphlan: 25900713
# 6.F4.S2.CA.8.8.2017 = 9_39  metaphlan: 8699809
# 5.F4.S2.CA.8.8.2017 = 9_38  metaphlan: 11619167
# 2.F4.S2.CA.29.8.2017 = 14_30  metaphlan: 12798500
# 5.F5.S2.CA.21.9.17 = 10_52  metaphlan: 17231290