library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling


summary(sample_sums(subsetG))
otu_tab <- t(abundances(subsetG))
# rarefaction curve
p <- vegan::rarecurve(otu_tab,
                      step = 50, label = FALSE,
                      sample = min(rowSums(otu_tab),
                                   col = "blue", cex = 0.6))
# samples plateau so sufficient sequencing depth

#remove samples with lower sequencing depth?

set.seed(1337)

ps0.rar <- rarefy_even_depth(subsetG, sample.size = 46731)

# no data removed? function not advisable generally > ?rarefy_even_depth()


ps0.rar

plot_taxa_prevalence(ps0.rar, "Phylum")

hmp.div <- microbiome::alpha(ps0.rar, index = "all")

datatable(hmp.div)


hmp.meta <- meta(ps0.rar)
hmp.meta$sam_name <- rownames(hmp.meta)
hmp.div$sam_name <- rownames(hmp.div)
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")
colnames(div.df)

p <- ggboxplot(div.df,
               x = "Farm2",
               y = "diversity_shannon",
               fill = "Farm2",
               palette = "jco")


p <- p + rotate_x_text()
p


div.df2 <- div.df[, c("Farm2", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "diversity_fisher", "diversity_coverage")]

colnames(div.df2) <- c("Farm", "Inverse Simpson", "Gini-Simpson", "Shannon", "Fisher", "Coverage")

div_df_melt <- reshape2::melt(div.df2)


p <- ggboxplot(div_df_melt, x = "Farm", y = "value",
               fill = "Farm",
               palette = "jco",
               legend= "right",
               facet.by = "variable",
               scales = "free")

p <- p + rotate_x_text()
p <- p + rremove("x.text")

p


#ggsave("../Metataxonomic/Figures/Diversities.pdf", height = 4, width = 10)

#lev <- levels(div_df_melt$Variable)
lev = c("Farm1","Farm2","Farm3","Farm4")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
pval <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("****", "***", "**", "*", "n.s")
)

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")
  )
)


p2


library(picante)
ps0.rar.asvtab <- as.data.frame(ps0.rar@otu_table)
ps0.rar.tree <- ps0.rar@phy_tree
ps0.rar@phy_tree

df.pd <- pd(t(ps0.rar.asvtab), ps0.rar.tree,include.root=T)

datatable(df.pd)
hmp.meta$Phylogenetic_Diversity <- df.pd$PD

pd.plot <- ggboxplot(hmp.meta,
                     x = "Farm2",
                     y = "Phylogenetic_Diversity",
                     fill = "Farm2",
                     palette = "jco",
                     ylab = "Phylogenetic Diversity",
                     xlab = "Farm",
                     legend = "right"
)
pd.plot <- pd.plot + rotate_x_text()

pd.plot + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")
  )
)

#rarefy to equal library size or not?

lib.div <- microbiome::alpha(subsetG, index = "all")
lib.div2 <- richness(subsetG)
lib.div$ReadsPerSample <- sample_sums(subsetG)
lib.div$Observed <- lib.div2$observed
colnames(lib.div)
p1 <- ggscatter(lib.div, "diversity_shannon", "ReadsPerSample") +
  stat_cor(method = "pearson")
p2 <- ggscatter(lib.div, "diversity_inverse_simpson", "ReadsPerSample",
                add = "loess"
) +
  stat_cor(method = "pearson")
p3 <- ggscatter(lib.div, "Observed", "ReadsPerSample",
                add = "loess") +
  stat_cor(
    method = "pearson",
    label.x = 100,
    label.y = 50000
  )

ggarrange(p1, p2, p3, ncol = 2, nrow = 2)
