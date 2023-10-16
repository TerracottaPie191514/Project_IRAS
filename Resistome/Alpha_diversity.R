library(data.table) # Alternative to data.frame
library(picante) # Used for calculating Phylogenetic diversities
library(lme4) # Repeated measures, add to report if used
library(QsRutils) # For the goods() function, to estimate coverage

# used the following guides: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html, https://rpubs.com/maddieSC/R_SOP_BRC_Oct_2019, https://rpubs.com/lconteville/713954

otu_tab <- t(abundances(Rps)) # can use veganotu() function for this
otu_tab2 <- t(abundances(Rps_tpm))
otu_tab3 <- t(abundances(Rps_mp))


# rarefaction curve
vegan::rarecurve(otu_tab,
                      step = 50, label = FALSE,
                      sample = min(rowSums(otu_tab),
                                   col = "blue", cex = 0.6))

# we can add lines to show sampling depths
rarecurve(otu_tab, step=50, ylab = "ARGs")
abline(v=sample_sums(Rps), lty='dotted', lwd=0.5)
  
rarecurve(otu_tab3, step=500000, ylab = "ARGs", label = FALSE)
abline(v=sample_sums(Rps), lty='dotted', lwd=0.5)


# virtually no samples are reaching a plateau so sequencing depth is not appropriate, undersampling for most of the dataset

vegan::rarecurve(otu_tab2,
                      step = 50, label = FALSE,
                      sample = min(rowSums(otu_tab),
                                   col = "blue", cex = 0.6))

# rarefaction curves of TPM data are all converging towards the plateau, no rarefaction required


# we use Good's coverage test to see the amount of singletons in the samples

summary(goods(otu_tab3)) # on average, 0.65% of the reads in the samples are singletons 

summary(goods(otu_tab2)) # there are no singletons in tpm

Rps %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% veganotu() %>% goods() %>% summary()
# for the stable Farm2R1S1, on average, 4.5% of the reads in the samples are singletons


# rarefy to equal library size or not?

lib.div <- microbiome::alpha(Rps, index = "all")
lib.div2 <- richness(Rps)
lib.div$ReadsPerSample <- sample_sums(Rps)
lib.div$chao1 <- lib.div2$chao1
colnames(lib.div)
p1 = ggscatter(lib.div, "diversity_shannon", "ReadsPerSample", xlab = "Shannon diversity", add = "loess") +
  stat_cor(method = "pearson")
p2 = ggscatter(lib.div, "diversity_inverse_simpson", "ReadsPerSample",  xlab = "Inverse Simpson diversity", add = "loess") +
  stat_cor(method = "pearson")
p3 = ggscatter(lib.div, "observed", "ReadsPerSample",  xlab = "Observed", add = "loess") +
  stat_cor(method = "pearson")

df.pd <- pd(t(as.data.frame(Rps@otu_table)), Rps@phy_tree,include.root=T) # transposing for use in picante
lib.div$Phylogenetic_Diversity <- df.pd$PD

p4 = ggscatter(lib.div, "Phylogenetic_Diversity", "ReadsPerSample",  xlab = "Phylogenetic diversity", add = "loess") +
  stat_cor(method = "pearson")

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

p4

# we can clearly see an increase in reads/sample when increasing abundance, so we require a rarefaction for FPKM data

set.seed(1337)

ps0.rar <- rarefy_even_depth(Rps, sample.size = 118) # we do not want to lose samples so lowest sample size is maintained, 492! OTUs are removed

ps0.rar <- srs_p(Rps) # we do not want to lose samples so lowest sample size is maintained, 492! OTUs are removed

# In the rarefaction curves, we can clearly see three outliers, with very large sample sizes and ARGs

# remove problematic samples
Rps %>% subset_samples(Sample_Unique != "10_1" & Sample_Unique != "10_2" & Sample_Unique != "10_3") 

sample_data(Rps_mp)$Sample_Unique = sample_names(Rps_mp)
sample_variables(Rps_mp)
Rps_mp %<>% subset_samples(Sample_Unique != "10_1" & Sample_Unique != "10_2" & Sample_Unique != "10_3") 

# function not advisable generally > ?rarefy_even_depth()

# In order to create taxa prevalence plots with these functions, we need to change our "taxa" levels to the names of actual taxa
# Phylum = AMR_class_primary, Order = AMR_class_secondary, Class = ARGCluster90, Family = ID_Clust_Refsequence
colnames(ps0.rar@tax_table) = c("Phylum", "Order", "Class","Family") 
plot_taxa_prevalence(ps0.rar, "Phylum")
pscopy = Rps
colnames(pscopy@tax_table) = c("Phylum", "Order", "Class","Family")
plot_taxa_prevalence(pscopy, "Phylum") # Sadly we can see entire phyla disappear, as well as many individual values, this rarefaction is deemed not appropriate either way
ps_tpmcopy = Rps_tpm
colnames(ps_tpmcopy@tax_table) = c("Phylum", "Order", "Class","Family")
plot_taxa_prevalence(ps_tpmcopy, "Phylum") # TPM data


# We will be using a different type of rarefaction

seed=1337
set.seed(seed)

ssum <- sample_sums(Rps)
ssum
ssum2 <- round(ssum/5000,0)
ssum2

nsamples <- length(sample_names(Rps)) # amount of samples

#first do 1st sample then build around it the rest
i=1
m.tmp <- subset_samples( Rps, sample_names(Rps) == sample_names(Rps)[i] )
m.tmp <- rarefy_even_depth( m.tmp, sample.size=ssum2[i], trimOTUs=FALSE, rngseed=seed )

Rps_rar <- m.tmp

for( i in 2:nsamples) {
  m.tmp <- subset_samples( Rps, sample_names(Rps) == sample_names(Rps)[i] )
  m.tmp <- rarefy_even_depth( m.tmp, sample.size=ssum2[i], trimOTUs=FALSE, rngseed=seed )
  Rps_rar <- merge_phyloseq(Rps_rar, m.tmp)
}


# specific variables ( niet wat er bedoeld werd, voor boxplots doen dit)

# non-AB
pscopy %>% ps_filter(AB == "no") %>% plot_taxa_prevalence("Phylum") + ggtitle("no")
# AB
pscopy %>% ps_filter(AB == "yes") %>% plot_taxa_prevalence("Phylum") + ggtitle("yes")

# age 14
pscopy %>% ps_filter(Age == "14") %>% plot_taxa_prevalence("Phylum") + ggtitle("14")
# age 35
pscopy %>% ps_filter(Age == "35") %>% plot_taxa_prevalence("Phylum") + ggtitle("35")

# farm1
pscopy %>% ps_filter(Farm2 == "Farm1") %>% plot_taxa_prevalence("Phylum") + ggtitle("farm1")
# farm2
pscopy %>% ps_filter(Farm2 == "Farm2") %>% plot_taxa_prevalence("Phylum") + ggtitle("farm2")
# farm3
pscopy %>% ps_filter(Farm2 == "Farm3") %>% plot_taxa_prevalence("Phylum") + ggtitle("farm3")
# farm4
pscopy %>% ps_filter(Farm2 == "Farm4") %>% plot_taxa_prevalence("Phylum") + ggtitle("farm4")




# calculating alpha diversity measures

hmp.div <- microbiome::alpha(Rps, index = "all") # use ps0.rar if rarefied

datatable(hmp.div)

hmp.meta <- meta(Rps)
hmp.meta$sam_name <- rownames(hmp.meta)
hmp.div$sam_name <- rownames(hmp.div)
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")
colnames(div.df)


#based on microbial agent
# Shortening names
div.df$Cox[div.df$Cox == "narasinandnicarbazin(maxiban)"] = "Maxiban"
div.df$Cox[div.df$Cox == "narasin(monteban)"] = "Monteban"
div.df$Cox[div.df$Cox == "salinomycin(Sacox120microGranulate)"] = "Sacox"

div.df2 <- div.df[, c("Cox", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "chao1", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Agent", "Inverse Simpson", "Gini-Simpson", "Shannon", "chao1", "Coverage", "Pielou")


div_df_melt <- reshape2::melt(div.df2)

lev = c("Maxiban","Sacox","Monteban","None")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

ggboxplot(div_df_melt, x = "Agent", y = "value",
               fill = "Agent",
               palette = "lancet",
               legend= "right",
               facet.by = "variable",
               scales = "free",
               title = "FPKM Alpha diversity metrics by microbial agent",
          outlier.shape = NA) + 
  rremove("x.text") + stat_compare_means(
    comparisons = L.pairs,
    method = "wilcox.test",
    label = "p.signif"
    ) + geom_jitter(size = 0.7, alpha = 0.9)

df.pd <- pd(t(as.data.frame(Rps@otu_table)), Rps@phy_tree,include.root=T) # transposing for use in picante
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
                     title = "FPKM phylogenetic diversity by microbial agent",
          outlier.shape = NA) + rotate_x_text() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")
  )
) + geom_jitter(size = 0.7, alpha = 0.9)


# age / days

div.df2 <- div.df[, c("Age", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "chao1", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Age", "Inverse Simpson", "Gini-Simpson", "Shannon", "chao1", "Coverage", "Pielou")

div.df2$Age = as.factor(div.df2$Age)
div_df_melt <- reshape2::melt(div.df2)

ggboxplot(div_df_melt, x = "Age", y = "value",
          fill = "Age",
          palette = "lancet",
          legend= "right",
          facet.by = "variable",
          scales = "free",
          title = "FPKM Alpha diversity metrics by age",
          outlier.shape = NA) + 
  rremove("x.text") + stat_compare_means(
    method = "wilcox.test",) + geom_jitter(size = 0.7, alpha = 0.9)


ggboxplot(hmp.meta,
          x = "Age",
          y = "Phylogenetic_Diversity",
          fill = "Age",
          palette = "lancet",
          ylab = "Phylogenetic Diversity",
          xlab = "Age",
          legend = "right",
          title = "FPKM phylogenetic diversity by age",
          outlier.shape = NA) + rotate_x_text() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
  stat_compare_means(method = "wilcox.test", paired = TRUE) + geom_jitter(size = 0.7, alpha = 0.9)
    

# farms / company

div.df2 <- div.df[, c("Farm2", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "chao1", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("Farm", "Inverse Simpson", "Gini-Simpson", "Shannon", "chao1", "Coverage", "Pielou")

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
          title = "FPKM Alpha diversity metrics by farm",
          outlier.shape = NA) + rotate_x_text() + rremove("x.text") +
  stat_compare_means(method = "wilcox.test",
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
          order = lev,
          title = "FPKM phylogenetic diversity by farm",
          outlier.shape = NA) + rotate_x_text() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
  stat_compare_means(
    comparisons = L.pairs,
    label = "p.signif"
    ) + geom_jitter(size = 0.7, alpha = 0.9)
  

# based on AB

div.df2 <- div.df[, c("AB", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "chao1", "diversity_coverage", "evenness_pielou")]
colnames(div.df2) <- c("AB", "Inverse Simpson", "Gini-Simpson", "Shannon", "chao1", "Coverage", "Pielou")


div_df_melt <- reshape2::melt(div.df2)

ggboxplot(div_df_melt, x = "AB", y = "value",
          fill = "AB",
          palette = "lancet",
          legend= "right",
          facet.by = "variable",
          scales = "free",
          title = "FPKM Alpha diversity metrics by antibiotic usage",
          outlier.shape = NA) + 
  rremove("x.text") + stat_compare_means(
    method = "wilcox.test") + geom_jitter(size = 0.7, alpha = 0.9)


ggboxplot(hmp.meta,
          x = "AB",
          y = "Phylogenetic_Diversity",
          fill = "AB",
          palette = "lancet",
          ylab = "Phylogenetic Diversity",
          xlab = "Antibiotics used",
          legend = "right",
          title = "FPKM phylogenetic diversity by antibiotic usage",
          outlier.shape = NA) + rotate_x_text() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
  stat_compare_means() + geom_jitter(size = 0.7, alpha = 0.9)



# alternative way of plotting

plot_richness(Rps, x="Age", measures=c("chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), color = "Age", nrow = 2)+
  geom_boxplot(alpha=0.6) + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

plot_richness(ps0.rar, x="Farm2", nrow = 2, color = "Farm2", title = "Alpha diversity metrics based on farm (rarefied)")+
  geom_boxplot(alpha=0.6) + theme_classic() +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.title.x = element_blank()) 


## Looking at significance

# Checking for normality

hist(lib.div$chao1, main="chao1 richness", xlab="")
hist(lib.div$diversity_shannon, main="Shannon diversity", xlab="")
hist(lib.div$diversity_fisher, main="Fisher diversity", xlab="")
hist(lib.div$diversity_gini_simpson, main="Gini-Simpson diversity", xlab="")
hist(lib.div$diversity_inverse_simpson, main="Inverse Simpson evenness", xlab="")
hist(lib.div$evenness_pielou, main="Pielou evenness", xlab="")
hist(lib.div$diversity_coverage, main="Coverage diversity", xlab="")




# If data is normally distributed we can use ANOVA / t-tests, if not we will use Kruskal-Wallis tests
# In this case, the data seems roughly normally distributed, we can use Shapiro-Wilk tests to test for normality for individual measures
shapiro.test(lib.div$chao1) # test deems it  normally distributed p>0,05
shapiro.test(lib.div$diversity_shannon) # test deems this measure not normally distributed p<0,05
shapiro.test(lib.div$diversity_fisher) # test deems this measure not normally distributed p<0,05
shapiro.test(lib.div$diversity_gini_simpson) # test deems this measure not normally distributed p<0,05
shapiro.test(lib.div$diversity_inverse_simpson) # test deems this measure not normally distributed p<0,05
shapiro.test(lib.div$evenness_pielou) # test deems this measure normally distributed p>0,05
shapiro.test(lib.div$diversity_coverage) # test deems this measure not normally distributed p>0,05


# Fairly small sample sizes however, and the shaprio-wilk test is not perfect, we will assume normality for all measures except for Shannon and Gini-simpson diversity based on the graphs
# The variables that we are interested in are the Age, which Farm the samples are from, and whether antibiotics were applied, all of which are categorical variables.

# We will run ANOVAs for the normally distributed variables

# Age

# Normally distributed with only 2 levels, so we can use t-tests : 

t.test(lib.div$chao1 ~ sample_data(Rps)$Age)

t.test(lib.div$diversity_shannon ~ sample_data(Rps)$Age) # shannon diversity does not seem to significantly differ across the different age groups

t.test(lib.div$diversity_fisher ~ sample_data(Rps)$Age)

t.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Age) 

t.test(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$Age) # not significant

t.test(lib.div$evenness_pielou ~ sample_data(Rps)$Age) # not significant

t.test(lib.div$evenness_pielou ~ sample_data(Rps)$Age) # not significant


# Non-normally distributed

wilcox.test(lib.div$chao1 ~ sample_data(Rps)$Age)

wilcox.test(lib.div$diversity_shannon ~ sample_data(Rps)$Age) # shannon diversity does not seem to significantly differ across the different age groups

wilcox.test(lib.div$diversity_fisher ~ sample_data(Rps)$Age)

wilcox.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Age) # remove this later

wilcox.test(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$Age) # not significant

wilcox.test(lib.div$evenness_pielou ~ sample_data(Rps)$Age) 


boxplot(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Age, ylab="Gini-Simpson's diversity") # the boxplots are quite similar so this is not unexpected

# For age, the groups seems significantly different in chao1, shannon, fisher diversity and simpson evenness, but not in gini-simpson, inverse simpson diversity, and pielou evenness.

# Antibiotics

t.test(lib.div$chao1 ~ sample_data(Rps)$AB) # not significant

t.test(lib.div$evenness_pielou ~ sample_data(Rps)$AB) # not significant

# used these functions to get means and sd per variable and alpha diversity metric
lib.div.ab = lib.div
lib.div.ab$AB = sample_data(Rps)$AB

aggregate(lib.div.ab$evenness_pielou, list(lib.div.ab$AB), FUN=mean) 
aggregate(lib.div.ab$evenness_pielou, list(lib.div.ab$AB), FUN=sd) 

# Non-normally distributed


wilcox.test(lib.div$diversity_shannon ~ sample_data(Rps)$AB) # shannon diversity does not seem to significantly differ across the different AB groups

wilcox.test(lib.div$diversity_fisher ~ sample_data(Rps)$AB)

wilcox.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$AB) # remove this later

wilcox.test(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$AB) # not significant

wilcox.test(lib.div$diversity_coverage ~ sample_data(Rps)$AB) 

boxplot(lib.div$evenness_pielou ~ sample_data(Rps)$AB, ylab="pielou") # the boxplots are quite similar so this is not unexpected

# AB does not seem to significantly differ in their alpha diversities 

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

aov.inv_simpson.farm = aov(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$Farm2)
summary(aov.inv_simpson.farm)
TukeyHSD(aov.inv_simpson.farm)

aov.pielou.farm = aov(lib.div$evenness_pielou ~ sample_data(Rps)$Farm2)
summary(aov.pielou.farm)
TukeyHSD(aov.pielou.farm)

# Non-normally distributed

kruskal.test(lib.div$chao1 ~ sample_data(Rps)$Farm2)
pairwise.wilcox.test(lib.div$chao1, sample_data(Rps)$Farm2, p.adjust.method="fdr")

kruskal.test(lib.div$diversity_shannon ~ sample_data(Rps)$Farm2) # shannon diversity does not seem to significantly differ across the different Farm2 groups
pairwise.wilcox.test(lib.div$diversity_shannon, sample_data(Rps)$Farm2, p.adjust.method="fdr")

kruskal.test(lib.div$diversity_fisher ~ sample_data(Rps)$Farm2)
pairwise.wilcox.test(lib.div$diversity_fisher, sample_data(Rps)$Farm2, p.adjust.method="fdr")

kruskal.test(lib.div$diversity_gini_simpson ~ sample_data(Rps)$Farm2) # remove this later
pairwise.wilcox.test(lib.div$diversity_gini_simpson, sample_data(Rps)$Farm2, p.adjust.method="fdr")

kruskal.test(lib.div$diversity_inverse_simpson ~ sample_data(Rps)$Farm2) # not significant
pairwise.wilcox.test(lib.div$diversity_inverse_simpson, sample_data(Rps)$Farm2, p.adjust.method="fdr")

kruskal.test(lib.div$evenness_pielou ~ sample_data(Rps)$Farm2) 
pairwise.wilcox.test(lib.div$evenness_pielou, sample_data(Rps)$Farm2, p.adjust.method="fdr")



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


# Mixed models, variables might not be independent

aov.shannon.age_farm = aov(lib.div$diversity_shannon ~ sample_data(Rps)$Age*sample_data(Rps)$Farm2)
summary(aov.shannon.age_farm)

aov.shannon.age_farm = aov(lib.div$diversity_shannon ~ sample_data(Rps)$Age+sample_data(Rps)$Farm2)
summary(aov.shannon.age_farm)

aov.shannon.age_AB = aov(lib.div$diversity_shannon ~ sample_data(Rps)$Age*sample_data(Rps)$AB)
summary(aov.shannon.age_AB)

aov.shannon.age_AB = aov(lib.div$diversity_shannon ~ sample_data(Rps)$Age+sample_data(Rps)$AB)
summary(aov.shannon.age_AB)


aov.shannon.age_all = aov(lib.div$diversity_shannon ~ sample_data(Rps)$Age*sample_data(Rps)$AB*sample_data(Rps)$Farm2*sample_data(Rps)$Cox*sample_data(Rps)$FarmRoundStable)
summary(aov.shannon.age_all)

# repeated measures, look at second guide to figure this out

rm.shannon.all = lmer(lib.div$diversity_shannon ~ sample_data(Rps)$AB + (1|sample_data(Rps)$Farm2))
summary(rm.shannon.all)

# declutter R environment by removing objects that no longer serve a purpose
rm(p1, p2, p3, ps0.rar, otu_tab, otu_tab2, L.pairs, pscopy, ps_tpmcopy, pval, lev, aov.chao1.farm, aov.fisher.farm, aov.gini_simpson.farm, aov.inv_simpson.farm, aov.pielou.farm, aov.shannon.farm, div_df_melt, div.df, div.df2, glm.chao1.age, glm.fisher.age, glm.gini_simpson.age, glm.inv_simpson.age, glm.shannon.age, glm.simpson.age, glm.pielou.age, gaussian.gini_simpson.conc, qp.gini_simpson.conc, df.pd, lib.div, lib.div2, hmp.div, hmp.meta)

rm(aov.shannon.age_farm, aov.shannon.age_AB, aov.shannon.age_all, aov.simpson.farm)


