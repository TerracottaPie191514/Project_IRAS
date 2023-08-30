library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(picante)
library(lme4)



# used the following guides: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html, https://rpubs.com/maddieSC/R_SOP_BRC_Oct_2019, https://rpubs.com/lconteville/713954

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
plot_taxa_prevalence(pscopy, "Phylum") # Sadly we can see entire phyla disappear, as well as many individual values, rarefaction is deemed not appropriate either way
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


# different way of plotting

plot_richness(ps0.rar, x="Age", measures=c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), color = "Age", nrow = 2)+
  geom_boxplot(alpha=0.6) + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

plot_richness(ps0.rar, x="Farm2", nrow = 2, color = "Farm2", title = "Alpha diversity metrics based on farm (rarefied)")+
  geom_boxplot(alpha=0.6) + theme_classic() +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.title.x = element_blank()) 


## Looking at significance

# Checking for normality

hist(lib.div$chao1, main="Chao richness", xlab="")
hist(lib.div$diversity_shannon, main="Shannon diversity", xlab="")
hist(lib.div$diversity_fisher, main="Fisher diversity", xlab="")
hist(lib.div$diversity_gini_simpson, main="Gini-Simpson diversity", xlab="")
hist(lib.div$evenness_simpson, main="Simpson evenness", xlab="")
#hist(1/lib.div$evenness_simpson, main="Inverse Simpson evenness", xlab="")
hist(lib.div$diversity_inverse_simpson, main="Inverse Simpson evenness", xlab="")
hist(lib.div$evenness_pielou, main="Pielou evenness", xlab="")


# If data is normally distributed we can use ANOVA / t-tests, if not we will use Kruskal-Wallis tests
# In this case, the data seems roughly normally distributed, we can use Shapiro-Wilk tests to test for normality for individual measures
shapiro.test(lib.div$chao1) # test deems it not normally distributed
shapiro.test(lib.div$diversity_shannon) # test deems this measure not normally distributed
shapiro.test(lib.div$diversity_fisher) # test deems this measure not normally distributed
shapiro.test(lib.div$diversity_gini_simpson) # test deems this measure not normally distributed
shapiro.test(lib.div$evenness_simpson) # test deems this measure normally distributed
shapiro.test(lib.div$diversity_inverse_simpson) # test deems this measure not normally distributed
shapiro.test(lib.div$evenness_pielou) # test deems this measure normally distributed

# Fairly small sample sizes however, and the shaprio-wilk test is not perfect, we will assume normality for all measures except for Gini-simpson diversity based on the graphs
# The variables that we are interested in are the Age, which Farm the samples are from, and whether antibiotics were applied, all of which are categorical variables.

# We will run ANOVAs for the normally distributed measures

# Age

# Normally distributed with only 2 levels, so we can use t-tests : 
# vis : #boxplot(lib.div$diversity_shannon ~  sample_data(Rps)$Age, ylab="Shannon's diversity")

t.test(lib.div$chao1 ~ sample_data(Rps)$Age)

t.test(lib.div$diversity_shannon ~ sample_data(Rps)$Age)

t.test(lib.div$diversity_fisher ~ sample_data(Rps)$Age)

t.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Age) # remove this later

t.test(lib.div$evenness_simpson ~ sample_data(Rps)$Age)

t.test(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$Age) # not significant

t.test(lib.div$evenness_pielou ~ sample_data(Rps)$Age) # not significant


# Non-normally distributed

wilcox.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Age) # gini simpson diversity does not seem to significantly differ across the different age groups

boxplot(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Age, ylab="Gini-Simpson's diversity") # the boxplots are quite similar so this is not


# Antibiotics

t.test(lib.div$chao1 ~ sample_data(Rps)$AB)

t.test(lib.div$diversity_shannon ~ sample_data(Rps)$AB)

t.test(lib.div$diversity_fisher ~ sample_data(Rps)$AB)

t.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$AB)

t.test(lib.div$evenness_simpson ~ sample_data(Rps)$AB)

t.test(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$AB)

t.test(lib.div$evenness_pielou ~ sample_data(Rps)$AB)

# Non-normally distributed

wilcox.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$AB) # none of the AB seem to significantly differ

# Farm has more than 2 levels, so we will use ANOVAs

aov.chao1.farm = aov(lib.div$chao1 ~ sample_data(Rps)$Farm2)
summary(aov.chao1.farm)
TukeyHSD(aov.chao1.farm)  # it seems that Farm 1 differs significantly from Farm 2 and 3 but not 4. Farm 2 differs from 1 and 4 but not 3, Farm 3 and 4 differ as well. ( 3 and 2 are similar, and 4 and 1 are similar)

aov.shannon.farm = aov(lib.div$diversity_shannon ~ sample_data(Rps)$Farm2)
summary(aov.shannon.farm)
TukeyHSD(aov.shannon.farm)

aov.fisher.farm = aov(lib.div$diversity_fisher ~ sample_data(Rps)$Farm2)
summary(aov.fisher.farm)
TukeyHSD(aov.fisher.farm)

aov.gini_simpson.farm = aov(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Farm2)
summary(aov.gini_simpson.farm)
TukeyHSD(aov.gini_simpson.farm)

aov.simpson.farm = aov(lib.div$evenness_simpson ~ sample_data(Rps)$Farm2)
summary(aov.simpson.farm)
TukeyHSD(aov.simpson.farm)

aov.inv_simpson.farm = aov(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$Farm2)
summary(aov.inv_simpson.farm)
TukeyHSD(aov.inv_simpson.farm)

aov.pielou.farm = aov(lib.div$evenness_pielou ~ sample_data(Rps)$Farm2)
summary(aov.pielou.farm)
TukeyHSD(aov.pielou.farm)

# Non-normally distributed

kruskal.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Farm2)

pairwise.wilcox.test(lib.div$diversity_gini_simpson, sample_data(Rps)$Farm2, p.adjust.method="fdr")

# In addition, it could be interesting to look at the concentration of DNA as a continuous variable


# Normally distributed

glm.chao1.age = glm(lib.div$chao1 ~ sample_data(Rps)$Conc...ng..µl.)
summary(glm.chao1.age)
plot(lib.div$chao1 ~ sample_data(Rps)$Conc...ng..µl.)
abline(glm.chao1.age)

glm.shannon.age = glm(lib.div$diversity_shannon ~ sample_data(Rps)$Conc...ng..µl.)
summary(glm.shannon.age)
plot(lib.div$diversity_shannon ~ sample_data(Rps)$Conc...ng..µl.)
abline(glm.shannon.age)

glm.fisher.age = glm(lib.div$diversity_fisher ~ sample_data(Rps)$Conc...ng..µl.)
summary(glm.fisher.age)
plot(lib.div$diversity_fisher ~ sample_data(Rps)$Conc...ng..µl.)
abline(glm.fisher.age)

glm.gini_simpson.age = glm(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Conc...ng..µl.)
summary(glm.gini_simpson.age)
plot(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Conc...ng..µl.)
abline(glm.gini_simpson.age)

glm.simpson.age = glm(lib.div$evenness_simpson ~ sample_data(Rps)$Conc...ng..µl.)
summary(glm.simpson.age)
plot(lib.div$evenness_simpson ~ sample_data(Rps)$Conc...ng..µl.)
abline(glm.simpson.age)

glm.inv_simpson.age = glm(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$Conc...ng..µl.)
summary(glm.inv_simpson.age)
plot(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$Conc...ng..µl.)
abline(glm.inv_simpson.age)

glm.pielou.age = glm(lib.div$evenness_pielou ~ sample_data(Rps)$Conc...ng..µl.)
summary(glm.pielou.age)
plot(lib.div$evenness_pielou ~ sample_data(Rps)$Conc...ng..µl.)
abline(glm.pielou.age)

# Non-normally distributed

gaussian.gini_simpson.conc = glm(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Conc...ng..µl., family="gaussian")
par(mfrow = c(1,2))
plot(gaussian.gini_simpson.conc, which=c(1,2))
qp.gini_simpson.conc = glm(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Conc...ng..µl., family="quasipoisson")
par(mfrow = c(1,2))
plot(qp.gini_simpson.conc, which=c(1,2)) # there are no huge differences between normal model (gaussian) and quasipoisson, we will use the latter regardless
summary(qp.gini_simpson.conc)
par(mfrow = c(1, 1))
#Plot
plot(log(lib.div$diversity_gini_simpson) ~ sample_data(Rps)$Conc...ng..µl., ylab="ln(Chao's richness)")
abline(qp.gini_simpson.conc)

# repeated measures, look at second guide to figure this out

rm.shannon.all = lmer(lib.div$diversity_shannon ~ sample_data(Rps)$AB + (1|sample_data(Rps)$Farm2))
summary(rm.shannon.all)

