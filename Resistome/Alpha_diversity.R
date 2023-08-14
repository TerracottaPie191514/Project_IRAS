library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(picante)



# used the following guide: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html

otu_tab <- t(abundances(Rps))
otu_tab2 <- t(abundances(Rps_tpm))

# rarefaction curve
p <- vegan::rarecurve(otu_tab,
                      step = 50, label = FALSE,
                      sample = min(rowSums(otu_tab),
                                   col = "blue", cex = 0.6))
# virtually no samples are reaching a plateau so sequencing depth is not appropriate

p <- vegan::rarecurve(otu_tab2,
                      step = 50, label = FALSE,
                      sample = min(rowSums(otu_tab),
                                   col = "blue", cex = 0.6))
# rarefaction curves of TPM data are all converging towards the plateau, no rarefaction required

# rarefy to equal library size or not?

lib.div <- microbiome::alpha(Rps, index = "all")
lib.div2 <- richness(Rps)
lib.div$ReadsPerSample <- sample_sums(Rps)
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

# we can clearly see an increase in reads/sample when increasing abundance, so we require a rarefaction for FPKM data

set.seed(1337)

ps0.rar <- rarefy_even_depth(Rps, sample.size = 118) # we do not want to lose samples so lowest sample size is maintained

# function not advisable generally > ?rarefy_even_depth()

# create taxa prevalence plots, need to change "taxa" levels to actual taxa
colnames(ps0.rar@tax_table) = c("Phylum", "Order", "Class","Family") # Phylum = AMR_class_primary, Order = AMR_class_secondary, Class = ARGCluster90, Family = ID_Clust_Refsequence
plot_taxa_prevalence(ps0.rar, "Phylum")
pscopy = Rps
colnames(pscopy@tax_table) = c("Phylum", "Order", "Class","Family")
plot_taxa_prevalence(pscopy, "Phylum")
ps_tpmcopy = Rps_tpm
colnames(ps_tpmcopy@tax_table) = c("Phylum", "Order", "Class","Family")
plot_taxa_prevalence(ps_tpmcopy, "Phylum")


hmp.div <- microbiome::alpha(Rps, index = "all") # use ps0.rar if rarefied

datatable(hmp.div)


hmp.meta <- meta(Rps) # use ps0.rar if rarefied
hmp.meta$sam_name <- rownames(hmp.meta)
hmp.div$sam_name <- rownames(hmp.div)
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")
colnames(div.df)




#based on microbial agent
div.df$Cox[div.df$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
div.df$Cox[div.df$Cox == "narasin(monteban)"] = "Monteban"
div.df$Cox[div.df$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

div.df$Cox

ggboxplot(div.df,
               x = "Cox",
               y = "diversity_shannon",
               fill = "Cox",
               palette = "jco") + 
  rotate_x_text()



div.df2 <- div.df[, c("Cox", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "chao1", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Agent", "Inverse Simpson", "Gini-Simpson", "Shannon", "Chao1", "Coverage", "Pielou")


div_df_melt <- reshape2::melt(div.df2)


p <- ggboxplot(div_df_melt, x = "Agent", y = "value",
               fill = "Agent",
               palette = "jco",
               legend= "right",
               facet.by = "variable",
               scales = "free")

p <- p + rotate_x_text()
p <- p + rremove("x.text")

p


# with significance

#lev <- levels(div_df_melt$Variable)
lev = c("Maxiban","Sacox","Monteban","None")
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




hmp.meta$Cox[hmp.meta$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
hmp.meta$Cox[hmp.meta$Cox == "narasin(monteban)"] = "Monteban"
hmp.meta$Cox[hmp.meta$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

hmp.meta$Cox


pd.plot <- ggboxplot(hmp.meta,
                     x = "Cox",
                     y = "Phylogenetic_Diversity",
                     fill = "Cox",
                     palette = "jco",
                     ylab = "Phylogenetic Diversity",
                     xlab = "Antimicrobial",
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

# age / days

p <- ggboxplot(div.df,
               x = "Age",
               y = "diversity_shannon",
               fill = "Age",
               palette = "jco")


p <- p + rotate_x_text()
p



div.df2 <- div.df[, c("Age", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "chao1", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Age", "Inverse Simpson", "Gini-Simpson", "Shannon", "Chao1", "Coverage", "Pielou")


div_df_melt <- reshape2::melt(div.df2)


p <- ggboxplot(div_df_melt, x = "Age", y = "value",
               fill = "Age",
               palette = "jco",
               legend= "right",
               facet.by = "variable",
               scales = "free")

p <- p + rotate_x_text()
p <- p + rremove("x.text")

p


# with significance

#lev <- levels(div_df_melt$Variable)
lev = c("14","35")
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



pd.plot <- ggboxplot(hmp.meta,
                     x = "Age",
                     y = "Phylogenetic_Diversity",
                     fill = "Age",
                     palette = "jco",
                     ylab = "Phylogenetic Diversity",
                     xlab = "Age",
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


# farms / company

p <- ggboxplot(div.df,
               x = "Farm2",
               y = "diversity_shannon",
               fill = "Farm2",
               palette = "jco")


p <- p + rotate_x_text()
p


div.df2 <- div.df[, c("Farm2", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "chao1", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Farm", "Inverse Simpson", "Gini-Simpson", "Shannon", "Chao1", "Coverage", "Pielou")


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


# with significance

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

# based on AB

div.df2 <- div.df[, c("AB", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "chao1", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Antibiotics", "Inverse Simpson", "Gini-Simpson", "Shannon", "Chao1", "Coverage", "Pielou")


div_df_melt <- reshape2::melt(div.df2)


p <- ggboxplot(div_df_melt, x = "Antibiotics", y = "value",
               fill = "Antibiotics",
               palette = "jco",
               legend= "right",
               facet.by = "variable",
               scales = "free")

p <- p + rotate_x_text()
p <- p + rremove("x.text")

p


# with significance

#lev <- levels(div_df_melt$Variable)
lev = c("no","yes")
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



pd.plot <- ggboxplot(hmp.meta,
                     x = "AB",
                     y = "Phylogenetic_Diversity",
                     fill = "AB",
                     palette = "jco",
                     ylab = "Phylogenetic Diversity",
                     xlab = "Antibiotics Applied",
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


plot_richness(ps0.rar, x="Age", measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), color = "Age", nrow = 2)+
  geom_boxplot(alpha=0.6) + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

plot_richness(ps0.rar, x="Farm2", nrow = 2, color = "Farm2", title = "Alpha diversity metrics based on farm (rarefied)")+
  geom_boxplot(alpha=0.6) + theme_classic() +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.title.x = element_blank()) 
