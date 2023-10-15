library(data.table) # Alternative to data.frame
library(picante) # Used for calculating Phylogenetic diversities
library(lme4) # Repeated measures, add to report if used
library(QsRutils) # For the goods() function, to estimate coverage


# used the following guide: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html

otu_tab <- t(abundances(subset16S))
# rarefaction curve
vegan::rarecurve(otu_tab,
                      step = 50, label = FALSE,
                      sample = min(rowSums(otu_tab),
                                   col = "blue", cex = 0.6))
# samples plateau so sufficient sequencing depth

rarecurve(otu_tab, step=50)
abline(v=sample_sums(subset16S), lty='dotted', lwd=0.5)

summary(goods(otu_tab)) # there are no singletons in this data, already filtered out, means that richness estimates are probably unreliable or wrong?



#rarefy to equal library size or not?

lib.div <- microbiome::alpha(subset16S, index = "all")
lib.div2 <- richness(subset16S)
lib.div$ReadsPerSample <- sample_sums(subset16S)
lib.div$Observed <- lib.div2$observed
colnames(lib.div)
p1 = ggscatter(lib.div, "diversity_shannon", "ReadsPerSample", xlab = "Shannon diversity", add = "loess") +
  stat_cor(method = "pearson")
p2 = ggscatter(lib.div, "diversity_inverse_simpson", "ReadsPerSample",  xlab = "Inverse Simpson diversity", add = "loess") +
  stat_cor(method = "pearson")
p3 = ggscatter(lib.div, "observed", "ReadsPerSample",  xlab = "Observed", add = "loess") +
  stat_cor(method = "pearson")

df.pd <- pd(t(as.data.frame(subset16S@otu_table)), subset16S@phy_tree,include.root=T) # transposing for use in picante
lib.div$Phylogenetic_Diversity <- df.pd$PD

p4 = ggscatter(lib.div, "Phylogenetic_Diversity", "ReadsPerSample",  xlab = "Phylogenetic diversity", add = "loess") +
  stat_cor(method = "pearson")

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)


# remove samples with lower sequencing depth? not necessary for 16S dataset

#set.seed(1337)

#ps0.rar <- rarefy_even_depth(subset16S, sample.size = 46731)

# function not advisable generally > ?rarefy_even_depth()


#ps0.rar

#plot_taxa_prevalence(ps0.rar, "Phylum")

plot_taxa_prevalence(subset16S, "Phylum")

hmp.div <- microbiome::alpha(subset16S, index = "all") # use ps0.rar if rarefied

datatable(hmp.div)


hmp.meta <- meta(subset16S) # use ps0.rar if rarefied
hmp.meta$sam_name <- rownames(hmp.meta)
hmp.div$sam_name <- rownames(hmp.div)
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")
colnames(div.df)


#based on microbial agent
# Shortening names
div.df$Cox[div.df$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
div.df$Cox[div.df$Cox == "narasin(monteban)"] = "Monteban"
div.df$Cox[div.df$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

div.df2 <- div.df[, c("Cox", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "observed", "diversity_covhttp://127.0.0.1:35143/graphics/plot_zoom_png?width=1504&height=963erage", "evenness_pielou")]
colnames(div.df2) <- c("Agent", "Inverse Simpson", "Gini-Simpson", "Shannon", "Observed", "Coverage", "Pielou")


div_df_melt <- reshape2::melt(div.df2)

lev = c("Maxiban","Sacox","Monteban","None")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

ggboxplot(div_df_melt, x = "Agent", y = "value",
          fill = "Agent",
          palette = "lancet",
          legend= "right",
          facet.by = "variable",
          scales = "free",
          xlab = "Antimicrobial agent",
          title = "Alpha diversity metrics by microbial agent",
          outlier.shape = NA) + 
  rremove("x.text") + stat_compare_means(
    comparisons = L.pairs,
    label = "p.signif"
  ) + geom_jitter(size = 0.7, alpha = 0.9)

df.pd <- pd(t(as.data.frame(subset16S@otu_table)), subset16S@phy_tree,include.root=T) # transposing for use in picante
hmp.meta$Phylogenetic_Diversity <- df.pd$PD

# Shortening names
hmp.meta$Cox[hmp.meta$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
hmp.meta$Cox[hmp.meta$Cox == "narasin(monteban)"] = "Monteban"
hmp.meta$Cox[hmp.meta$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

ggboxplot(hmp.meta,
          x = "Cox",
          y = "Phylogenetic_Diversity",
          fill = "Cox",
          palette = "lancet",
          ylab = "Phylogenetic Diversity",
          xlab = "Antimicrobial agent",
          legend = "right",
          title = "Phylogenetic diversity by microbial agent",
          outlier.shape = NA) + rotate_x_text() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
  stat_compare_means(
    comparisons = L.pairs,
    label = "p.signif"
    ) + geom_jitter(size = 0.7, alpha = 0.9)
  


# age / days

div.df2 <- div.df[, c("Age", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "observed", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Age", "Inverse Simpson", "Gini-Simpson", "Shannon", "Observed", "Coverage", "Pielou")

div.df2$Age = as.factor(div.df2$Age)
div_df_melt <- reshape2::melt(div.df2)

ggboxplot(div_df_melt, x = "Age", y = "value",
          fill = "Age",
          palette = "lancet",
          legend= "right",
          facet.by = "variable",
          scales = "free",
          title = "Alpha diversity metrics by age",
          outlier.shape = NA) + 
  rremove("x.text") + stat_compare_means() + geom_jitter(size = 0.7, alpha = 0.9)


ggboxplot(hmp.meta,
          x = "Age",
          y = "Phylogenetic_Diversity",
          fill = "Age",
          palette = "lancet",
          ylab = "Phylogenetic Diversity",
          xlab = "Age",
          legend = "right",
          title = "Phylogenetic diversity by age",
          outlier.shape = NA) + rotate_x_text() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
  stat_compare_means(paired = TRUE) + geom_jitter(size = 0.7, alpha = 0.9)


# farms / company

div.df2 <- div.df[, c("Farm2", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "observed", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Farm", "Inverse Simpson", "Gini-Simpson", "Shannon", "Observed", "Coverage", "Pielou")

div_df_melt <- reshape2::melt(div.df2)

lev = c("Farm1","Farm2","Farm3","Farm4")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])


ggboxplot(div_df_melt, x = "Farm", y = "value",
          fill = "Farm",
          palette = "lancet",
          legend= "right",
          facet.by = "variable",
          scales = "free",
          order = lev,
          title = "Alpha diversity metrics by farm",
          outlier.shape = NA) + rotate_x_text() + rremove("x.text") + stat_compare_means(method = "wilcox.test",
            comparisons = L.pairs,
            label = "p.signif"
          ) + geom_jitter(size = 0.7, alpha = 0.9)


ggboxplot(hmp.meta,
          x = "Farm2",
          y = "Phylogenetic_Diversity",
          fill = "Farm2",
          palette = "lancet",
          ylab = "Phylogenetic Diversity",
          xlab = "Farm",
          legend = "right",
          title = "Phylogenetic diversity by farm",
          outlier.shape = NA) + rotate_x_text() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),  axis.title.x = element_blank()) +
  stat_compare_means(
    comparisons = L.pairs,
    label = "p.signif"
  ) + geom_jitter(size = 0.7, alpha = 0.9)


# based on AB

div.df2 <- div.df[, c("AB", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "observed", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("AB", "Inverse Simpson", "Gini-Simpson", "Shannon", "Observed", "Coverage", "Pielou")


div_df_melt <- reshape2::melt(div.df2)

ggboxplot(div_df_melt, x = "AB", y = "value",
          fill = "AB",
          palette = "lancet",
          legend= "right",
          facet.by = "variable",
          scales = "free",
          xlab = "Antibiotics used",
          title = "Alpha diversity metrics by antibiotic usage",
          outlier.shape = NA) + 
  rremove("x.text") + stat_compare_means() + geom_jitter(size = 0.7, alpha = 0.9)


ggboxplot(hmp.meta,
          x = "AB",
          y = "Phylogenetic_Diversity",
          fill = "AB",
          palette = "lancet",
          ylab = "Phylogenetic Diversity",
          xlab = "Antibiotics used",
          legend = "right",
          title = "Phylogenetic diversity by antibiotic usage",
          outlier.shape = NA) + rotate_x_text() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
  stat_compare_means() + geom_jitter(size = 0.7, alpha = 0.9)

## Looking at significance

# Checking for normality

hist(lib.div$observed, main="Observed richness", xlab="")
hist(lib.div$diversity_shannon, main="Shannon diversity", xlab="")
hist(lib.div$diversity_fisher, main="Fisher diversity", xlab="")
hist(lib.div$diversity_gini_simpson, main="Gini-Simpson diversity", xlab="")
hist(lib.div$diversity_inverse_simpson, main="Inverse Simpson evenness", xlab="")
hist(lib.div$evenness_pielou, main="Pielou evenness", xlab="")
hist(lib.div$diversity_coverage, main="Coverage diversity", xlab="")


# If data is normally distributed we can use ANOVA / t-tests, if not we will use Kruskal-Wallis tests
# In this case, the data seems roughly normally distributed for some metrics, we can use Shapiro-Wilk tests to test for normality for individual measures
shapiro.test(lib.div$observed) # test deems it  normally distributed p>0,05
shapiro.test(lib.div$diversity_shannon) # test deems this measure not normally distributed p<0,05
shapiro.test(lib.div$diversity_fisher) # test deems this measure normally distributed p>0,05
shapiro.test(lib.div$diversity_gini_simpson) # test deems this measure not normally distributed p<0,05
shapiro.test(lib.div$diversity_inverse_simpson) # test deems this measure not normally distributed p<0,05
shapiro.test(lib.div$evenness_pielou) # test deems this measure normally distributed p>0,05
shapiro.test(lib.div$diversity_coverage) # test deems this measure not normally distributed p<0,05

shapiro.test(lib.div$Phylogenetic_Diversity) # test deems this measure normally distributed p>0,05


# Based on shaprio-wilk tests we will assume normality for some measures 
# The variables that we are interested in are the Age, which Farm the samples are from, and whether antibiotics were applied, all of which are categorical variables.

# We will run ANOVAs for the normally distributed variables

# Age

# Normally distributed with only 2 levels, so we can use t-tests : 

t.test(lib.div$observed ~ sample_data(subset16S)$Age) # significant

t.test(lib.div$diversity_fisher ~ sample_data(subset16S)$Age)  # significant

t.test(lib.div$Phylogenetic_Diversity ~ sample_data(subset16S)$Age)  # significant


# Non-normally distributed

wilcox.test(lib.div$diversity_shannon ~ sample_data(subset16S)$Age) # shannon diversity seems to significantly differ across the different age groups

wilcox.test(lib.div$diversity_gini_simpson ~ sample_data(subset16S)$Age)  # significant

wilcox.test(lib.div$diversity_inverse_simpson ~ sample_data(subset16S)$Age)  # significant

wilcox.test(lib.div$evenness_pielou ~ sample_data(subset16S)$Age)  # significant

wilcox.test(lib.div$diversity_coverage ~ sample_data(subset16S)$Age)  # significant



# For age, the groups seems significantly different in all metrics except simpson evenness.

# Antibiotics

t.test(lib.div$observed ~ sample_data(subset16S)$AB) # significant

t.test(lib.div$diversity_fisher ~ sample_data(subset16S)$AB) # significant

t.test(lib.div$Phylogenetic_Diversity ~ sample_data(subset16S)$AB)  # significant


# Non-normally distributed

wilcox.test(lib.div$diversity_shannon ~ sample_data(subset16S)$AB) # shannon diversity does not seem to significantly differ across the different AB groups

wilcox.test(lib.div$diversity_gini_simpson ~ sample_data(subset16S)$AB) # not significant

wilcox.test(lib.div$diversity_inverse_simpson ~ sample_data(subset16S)$AB) # not significant

wilcox.test(lib.div$evenness_pielou ~ sample_data(subset16S)$AB) # not significant

wilcox.test(lib.div$diversity_coverage ~ sample_data(subset16S)$AB) # not significant



boxplot(lib.div$diversity_fisher ~ sample_data(subset16S)$AB, ylab="fisher") # the boxplots are quite similar so this is not unexpected

diversity_coverage

# used these functions to get means and sd per variable and alpha diversity metric
lib.div.ab = lib.div
lib.div.ab$AB = sample_data(subset16S)$AB

aggregate(lib.div.ab$observed, list(lib.div.ab$AB), FUN=mean) 
aggregate(lib.div.ab$observed, list(lib.div.ab$AB), FUN=sd) 



# AB does not seem to significantly differ in their alpha diversities except for observed and fisher diversity

# Farm has more than 2 levels, so we will use ANOVAs for normally distributed metrics

aov.observed.farm = aov(lib.div$observed ~ sample_data(subset16S)$Farm2)
summary(aov.observed.farm)
TukeyHSD(aov.observed.farm) # only not significant between 1 and 4 and 3 and 2

aov.fisher.farm = aov(lib.div$diversity_fisher ~ sample_data(subset16S)$Farm2)
summary(aov.fisher.farm)
TukeyHSD(aov.fisher.farm) # only not significant between 1 and 4 and 3 and 2


# Non-normally distributed

kruskal.test(lib.div$diversity_shannon ~ sample_data(subset16S)$Farm2) # shannon diversity seems to significantly differ across the different Farm2 groups
pairwise.wilcox.test(lib.div$diversity_shannon, sample_data(subset16S)$Farm2, p.adjust.method="fdr") # difference between 1 and 3 and 4 and all other farms

kruskal.test(lib.div$diversity_gini_simpson ~ sample_data(subset16S)$Farm2) # significant
pairwise.wilcox.test(lib.div$diversity_gini_simpson, sample_data(subset16S)$Farm2, p.adjust.method="fdr") # difference between farm 4 and 2 and 4 and 3


kruskal.test(lib.div$diversity_inverse_simpson ~ sample_data(subset16S)$Farm2) # not significant
pairwise.wilcox.test(lib.div$diversity_inverse_simpson, sample_data(subset16S)$Farm2, p.adjust.method="fdr") # difference between farm 4 and 2 and 4 and 3

kruskal.test(lib.div$evenness_pielou ~ sample_data(subset16S)$Farm2) # significant
pairwise.wilcox.test(lib.div$evenness_pielou, sample_data(subset16S)$Farm2, p.adjust.method="fdr") # difference between farm 4 and the other farms


# agent also has more than 2 levels, so we will use ANOVAs for normally distributed metrics

aov.observed.agent = aov(lib.div$observed ~ sample_data(subset16S)$Cox)
summary(aov.observed.agent)
TukeyHSD(aov.observed.agent) # only not significant between sacox and monteban

aov.fisher.agent = aov(lib.div$diversity_fisher ~ sample_data(subset16S)$Cox)
summary(aov.fisher.agent)
TukeyHSD(aov.fisher.agent) # only not significant maxiban & monteban and sacox & monteban


# Non-normally distributed

kruskal.test(lib.div$diversity_shannon ~ sample_data(subset16S)$Cox) # shannon diversity seems to significantly differ across the different Agent groups
pairwise.wilcox.test(lib.div$diversity_shannon, sample_data(subset16S)$Cox, p.adjust.method="fdr") # difference between none and sacox, diff between sacox and maxiban, no others

kruskal.test(lib.div$diversity_gini_simpson ~ sample_data(subset16S)$Cox) # significant
pairwise.wilcox.test(lib.div$diversity_gini_simpson, sample_data(subset16S)$Cox, p.adjust.method="fdr") # only diff between sacox and maxiban and sacox and none


kruskal.test(lib.div$diversity_inverse_simpson ~ sample_data(subset16S)$Cox) # not significant
pairwise.wilcox.test(lib.div$diversity_inverse_simpson, sample_data(subset16S)$Cox, p.adjust.method="fdr") # same as above

kruskal.test(lib.div$evenness_pielou ~ sample_data(subset16S)$Cox) # significant
pairwise.wilcox.test(lib.div$evenness_pielou, sample_data(subset16S)$Cox, p.adjust.method="fdr") # same as above




# Mixed models, variables might not be independent

aov.shannon.age_farm = aov(lib.div$diversity_shannon ~ sample_data(subset16S)$Age*sample_data(subset16S)$Farm2)
summary(aov.shannon.age_farm)

aov.shannon.age_farm = aov(lib.div$diversity_shannon ~ sample_data(subset16S)$Age+sample_data(subset16S)$Farm2)
summary(aov.shannon.age_farm)

aov.shannon.age_AB = aov(lib.div$diversity_shannon ~ sample_data(subset16S)$Age*sample_data(subset16S)$AB)
summary(aov.shannon.age_AB)

aov.shannon.age_AB = aov(lib.div$diversity_shannon ~ sample_data(subset16S)$Age+sample_data(subset16S)$AB)
summary(aov.shannon.age_AB)


aov.shannon.age_all = aov(lib.div$diversity_shannon ~ sample_data(subset16S)$Age*sample_data(subset16S)$AB*sample_data(subset16S)$Farm2*sample_data(subset16S)$Cox*sample_data(subset16S)$FarmRoundStable)
summary(aov.shannon.age_all)

# repeated measures, look at second guide to figure this out

rm.shannon.all = lmer(lib.div$diversity_shannon ~ sample_data(subset16S)$AB + (1|sample_data(subset16S)$Farm2))
summary(rm.shannon.all)

# declutter R environment by removing objects that no longer serve a purpose
rm(p1, p2, p3, ps0.rar, otu_tab, otu_tab2, L.pairs, pscopy, ps_tpmcopy, pval, lev, aov.observed.farm, aov.fisher.farm, aov.gini_simpson.farm, aov.inv_simpson.farm, aov.pielou.farm, aov.shannon.farm, div_df_melt, div.df, div.df2, glm.observed.age, glm.fisher.age, glm.gini_simpson.age, glm.inv_simpson.age, glm.shannon.age, glm.simpson.age, glm.pielou.age, gaussian.gini_simpson.conc, qp.gini_simpson.conc, df.pd, lib.div, lib.div2, hmp.div, hmp.meta)

rm(aov.shannon.age_farm, aov.shannon.age_AB, aov.shannon.age_all, aov.simpson.farm)
