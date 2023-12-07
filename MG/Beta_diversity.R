#### Load packages
library(scater) # For functions like plotReducedDim(), calculating dissimiilarity matrices etc. 
library(vegan) # used to run simper
library(nlme) # # for usage of llply(), to apply functions over lists

# Used the following guide : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

# single ordination visualisation

plot_ordination(subsetMG,  ordinate(subsetMG, method = "DPCoA", distance="bray"), "samples", color="Age", shape="AB")

# Different ordination methods based on BC dissimilarity (runs very slowly for MG)

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA", "DPCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(subsetMG, method=i, distance=dist)
  plot_ordination(subsetMG, ordi, "samples", color="Age", shape = "AB")
}, subsetMG, dist)

names(plist) <- ord_meths

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
ggplot(pdataframe, aes(Axis_1, Axis_2, color=Age, shape=AB)) + 
  geom_point(size=4) + 
  facet_wrap(~method, scales="free") +
  scale_fill_brewer(type="qual", palette="Set1") +
  scale_colour_brewer(type="qual", palette="Set1") +
  ggtitle("Different ordination methods for metagenomic data (Bray-Curtis)")

# PCoAs for different methods (fancy maken : + theme_classic() + scale_color_brewer("Farm2", palette = "Set2"))

# functionize plotting pcoa
plot_pcoa_ordination <- function(data, pcoa, var, title) {
  p <- plot_ordination(data, pcoa, color = var, shape = "AB") +
    geom_point(size = 3) +
    labs(title = title, color = var, shape = "Antibiotics used")
  
  return(p)
}

pcoa_bc = ordinate(subsetMG, "PCoA", "bray") 
pcoa_unifrac = ordinate(subsetMG, "PCoA", "unifrac") 
pcoa_wunifrac = ordinate(subsetMG, "PCoA", "wunifrac") 
pcoa_jsd = ordinate(subsetMG, "PCoA", "jsd") 
pcoa_jaccard = ordinate(subsetMG, "PCoA", "jaccard", binary=TRUE) 


plot_pcoa_ordination(subsetMG, pcoa_bc, "Age", "PCoA Bray Curtis")
plot_pcoa_ordination(subsetMG, pcoa_bc, "Farm2", "PCoA Bray Curtis")

plot_pcoa_ordination(subsetMG, pcoa_unifrac, "Age", "PCoA Unifrac")
plot_pcoa_ordination(subsetMG, pcoa_unifrac, "Farm2", "PCoA Unifrac")

plot_pcoa_ordination(subsetMG, pcoa_wunifrac, "Age", "PCoA Weighted Unifrac")
plot_pcoa_ordination(subsetMG, pcoa_wunifrac, "Farm2", "PCoA Weighted Unifrac")

plot_pcoa_ordination(subsetMG, pcoa_jsd, "Age", "PCoA Jensen-Shannon Divergence")
plot_pcoa_ordination(subsetMG, pcoa_jsd, "Farm2", "PCoA Jensen-Shannon Divergence")

plot_pcoa_ordination(subsetMG, pcoa_jaccard, "Age", "PCoA Jaccard")
plot_pcoa_ordination(subsetMG, pcoa_jaccard, "Farm2", "PCoA Jaccard")

plot_ordination(subsetMG, pcoa_jaccard, color = "Age", shape = "AB", label = "FarmRoundStable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard Age",color = "Age", shape = "Antibiotics used")



# BC plots for looking at percentage and total amount of bacterial reads mapped, removing 2_57 which is a giant outlier
subsetMG %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(subsetMG, "PCoA", "bray") , color = "ReadPerc", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC percentage",color = "ReadPerc") +
  scale_colour_viridis_c()

subsetMG %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(subsetMG, "PCoA", "bray") , color = "ReadTot", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC total",color = "ReadTot") +
  scale_colour_viridis_c()


pcoa_bc2 = subsetMG %>% subset_samples(Sample_Unique != "2_57") %>% ordinate("PCoA", "bray") 

pcoa_bc

subsetMG %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(pcoa_bc2, color = "ReadTot", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC concentration",color = "ReadTot") +
  scale_colour_viridis_c()

pcoa_bc3 = subsetMG_mp %>% subset_samples(Sample_Unique != "2_57") %>% ordinate("PCoA", "bray") 


subsetMG %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(pcoa_bc3, color = "ReadPerc", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC concentration",color = "ReadPerc") +
  scale_colour_viridis_c()


plot_scree(pcoa_jsd) #scree plots can be made for any of the PCoAs

# for changing specific labels etc

plot_ordination(subsetMG, pcoa_jaccard, color = "Age", shape = "AB", label = "FarmRoundStable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard Age",color = "Age", shape = "Antibiotics used")

# different way of plotting with scater and tses, this specifically is NMDS BC

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse <- transformCounts(tse, method = "relabundance")
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "BC", nmdsFUN = "monoMDS",
                exprs_values = "relabundance",
                keep_dist = TRUE)

plotReducedDim(tse, "BC", colour_by = "Age")

# PERMANOVAs

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse <- transformCounts(tse, method = "relabundance")

adonis2(t(assay(tse, "relabundance")) ~ Age, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ AB, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Cox, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Researcher, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ FeedProducent, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ LitterType, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ FeedType, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Gender, data = colData(tse), permutations = 9999) # NIET significant
adonis2(t(assay(tse, "relabundance")) ~ Stables, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ FlockSize, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Farm2, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ AgeParentStock, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ ReadPerc, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ ReadTot, data = colData(tse), permutations = 9999)


# Variances: Age: 0.045, AB: 0.0339, agent: 0.126, researcher: 0.0742, FP = 0.091, LT = 0.153, FT = 0.045, gender = 0.004,
# stable : 0.335, FS: 0.27, Farm: 0.232, APS: 0.261, readPerc: 0.133, readtot: 0.047
# Order: Stable>FS>APS>Farm>LT>ReadPerc>Cox>FP>LT>FP>Researcher>ReadTot>FT>Age>AB>Gender

# basically, composition seems to be different over every single variable, except for gender

# on genus level
tse_genus <- agglomerateByRank(tse, "Genus")
tse_genus <- transformCounts(tse_genus, method = "relabundance")

adonis2(t(assay(tse_genus, "relabundance")) ~ AB, data = colData(tse_genus), permutations = 9999) 
adonis2(t(assay(tse_genus, "relabundance")) ~ Cox, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ Researcher, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ FeedProducent, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ LitterType, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ FeedType, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ Gender, data = colData(tse_genus), permutations = 9999) # NIET significant
adonis2(t(assay(tse_genus, "relabundance")) ~ FarmRoundStable, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ FlockSize, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ Farm2, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ AgeParentStock, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ ReadTot, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ ReadPerc, data = colData(tse_genus), permutations = 9999)

# same results on genus level (and on phylum level, though p values become higher)
# for different ordination methods
ps1.rel <- microbiome::transform(subsetMG, "compositional")
metadf <- data.frame(sample_data(ps1.rel))

# alternative calculations
#otu <- abundances(ps1.rel)
#meta <- meta(ps1.rel)
#adonis2(t(otu) ~ Age, data = meta, permutations=9999, method = "bray")

#permanova = adonis(t(otu) ~ Age, data = meta, permutations=9999, method = "bray")
#permanova$aov.tab

unifrac.dist <- UniFrac(ps1.rel)

adonis2(unifrac.dist ~ Age, data = metadf)
adonis2(unifrac.dist ~ AB, data = metadf)
adonis2(unifrac.dist ~ Farm2, data = metadf)
adonis2(unifrac.dist ~ Cox, data = metadf)
adonis2(unifrac.dist ~ Researcher, data = metadf)
adonis2(unifrac.dist ~ LitterType, data = metadf)
adonis2(unifrac.dist ~ Gender, data = metadf) # not sign
adonis2(unifrac.dist ~ Stables, data = metadf)


# same patterns arise

wunifrac.dist <- UniFrac(ps1.rel, 
                         weighted = TRUE)

adonis2(wunifrac.dist ~ Age, data = metadf)
adonis2(wunifrac.dist ~ AB, data = metadf) 
adonis2(wunifrac.dist ~ Farm2, data = metadf)
adonis2(wunifrac.dist ~ Cox, data = metadf)
adonis2(wunifrac.dist ~ Researcher, data = metadf)
adonis2(wunifrac.dist ~ LitterType, data = metadf)
adonis2(wunifrac.dist ~ Gender, data = metadf) # not sign
adonis2(wunifrac.dist ~ Stables, data = metadf)


#  same patterns

jsd.dist <- phyloseq::distance(ps1.rel, "jsd")

adonis2(jsd.dist ~ Age, data = metadf)
adonis2(jsd.dist ~ AB, data = metadf) 
adonis2(jsd.dist ~ Farm2, data = metadf)
adonis2(jsd.dist ~ Cox, data = metadf)
adonis2(jsd.dist ~ Researcher, data = metadf)
adonis2(jsd.dist ~ LitterType, data = metadf)
adonis2(jsd.dist ~ Gender, data = metadf) # not sign
adonis2(jsd.dist ~ Stables, data = metadf)

# same is true for JSD

bray.dist <- phyloseq::distance(ps1.rel, "bray")

adonis2(bray.dist ~ Age, data = metadf)
adonis2(bray.dist ~ AB, data = metadf)
adonis2(bray.dist ~ Farm2, data = metadf)
adonis2(bray.dist ~ Cox, data = metadf)
adonis2(bray.dist ~ Researcher, data = metadf)
adonis2(bray.dist ~ LitterType, data = metadf)
adonis2(bray.dist ~ Gender, data = metadf) # not sign
adonis2(bray.dist ~ Stables, data = metadf)

# and BC

jaccard.dist <- phyloseq::distance(ps1.rel, "jaccard")

adonis2(jaccard.dist ~ Age, data = metadf)
adonis2(jaccard.dist ~ AB, data = metadf)
adonis2(jaccard.dist ~ Farm2, data = metadf)
adonis2(jaccard.dist ~ Cox, data = metadf)
adonis2(jaccard.dist ~ Researcher, data = metadf)
adonis2(jaccard.dist ~ LitterType, data = metadf)
adonis2(jaccard.dist ~ Gender, data = metadf) # not sign
adonis2(jaccard.dist ~ Stables, data = metadf)

# as well as jaccard

# PERMANOVA plots - Age

permanova_age <- adonis(t(assay(tse, "relabundance")) ~ Age, data = colData(tse), permutations = 9999)

coef <- coefficients(permanova_age)["Age1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "OTUs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial OTUs on age",
          rotate = TRUE,
          ggtheme = theme_minimal())

# AB

permanova_AB <- adonis(t(assay(tse, "relabundance")) ~ AB, data = colData(tse), permutations = 9999)

coef <- coefficients(permanova_AB)["AB1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "OTUs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial OTUs on AB",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Stable

permanova_stable <- adonis(t(assay(tse, "relabundance")) ~ Stables, data = colData(tse), permutations = 9999)

coef <- coefficients(permanova_stable)["Stables1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "OTUs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial OTUs on Stable",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Farm

permanova_farm <- adonis(t(assay(tse, "relabundance")) ~ Farm2, data = colData(tse), permutations = 9999)

coef <- coefficients(permanova_farm)["Farm21",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "OTUs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial OTUs on Farm",
          rotate = TRUE,
          ggtheme = theme_minimal())

# AB

permanova_agent <- adonis(t(assay(tse, "relabundance")) ~ Cox, data = colData(tse), permutations = 9999)

coef <- coefficients(permanova_agent)["Cox1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "OTUs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial OTUs on agent",
          rotate = TRUE,
          ggtheme = theme_minimal())

# same plots but for genera

#  Age

permanova_age <- adonis(t(assay(tse_genus, "relabundance")) ~ Age, data = colData(tse_genus), permutations = 9999)

coef <- coefficients(permanova_age)["Age1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "Genus",
          legend.title = "Genus contribution",
          title = "Impact of bacterial Genera on age",
          rotate = TRUE,
          ggtheme = theme_minimal())

# AB

permanova_AB <- adonis(t(assay(tse_genus, "relabundance")) ~ AB, data = colData(tse_genus), permutations = 9999)

coef <- coefficients(permanova_AB)["AB1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "Genus",
          legend.title = "Genus contribution",
          title = "Impact of bacterial Genera on AB",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Stable

permanova_stable <- adonis(t(assay(tse_genus, "relabundance")) ~ Stables, data = colData(tse_genus), permutations = 9999)

coef <- coefficients(permanova_stable)["Stables1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "Genus",
          legend.title = "Genus contribution",
          title = "Impact of bacterial Genera on Stable",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Farm

permanova_farm <- adonis(t(assay(tse_genus, "relabundance")) ~ Farm2, data = colData(tse_genus), permutations = 9999)

coef <- coefficients(permanova_farm)["Farm21",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "Genus",
          legend.title = "Genus contribution",
          title = "Impact of bacterial Genera on Farm",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Agent

permanova_agent <- adonis(t(assay(tse_genus, "relabundance")) ~ Cox, data = colData(tse_genus), permutations = 9999)

coef <- coefficients(permanova_agent)["Cox1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

df = data.frame(x = top.coef,
                y = factor(names(top.coef),
                           unique(names(top.coef))))

df$contr <- factor(ifelse(df$x < 0, "negative", "positive"), 
                   levels = c("negative", "positive"))

ggbarplot(df, x = "y", y = "x",
          fill = "contr",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          xlab = "Genus",
          legend.title = "Genus contribution",
          title = "Impact of bacterial Genera on agent",
          rotate = TRUE,
          ggtheme = theme_minimal())


# checking homogeneity condition - bray
# ANOVAs are performed on betadispers of our rel abund data to test whether groups are more variable than others

# Bray
ps.rel = microbiome::transform(subsetMG, "compositional")
meta = meta(ps.rel)
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Age))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$AB))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Farm2))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Stables))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Cox))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Researcher)) # not homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$LitterType)) # not homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Gender))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FlockSize))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$AgeParentStock)) 
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FeedProducent)) # not homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FeedType))

# ANOVAs do not work for ReadPerc, ReadToT, since they have different observations per group - continuous data
# We see that almost every variable is homogenously distributed, except for Researcher, LT and FP

# Jaccard
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Age))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$AB))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Farm2))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Stables))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Researcher)) # not homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$LitterType)) # not homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Gender))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FlockSize))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$AgeParentStock))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FeedProducent)) # not homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FeedType))

# same trends
# group variances are homogenous in most cases, so there are no differences in variances between groups -> anova / permanova condition met

# Tukey tests can be performed to see if and which groups differ in relation to variance

TukeyHSD(betadisper(vegdist(t(abundances(ps.rel))), meta$Farm2))


# different way of calculating homogeneity, permutation tests, null = no difference in dispersion between groups 
permutest(betadisper(vegdist(t(abundances(ps.rel))), meta$Age), pairwise = TRUE)

permutest(betadisper(unifrac.dist, metadf$Age), pairwise = TRUE) # looks like unifrac distances are homogenously dispersed for age
permutest(betadisper(unifrac.dist, metadf$AB), pairwise = TRUE) # not for AB though

permutest(betadisper(bray.dist, metadf$Age), pairwise = TRUE) # there are differences in P value with other method, but could be number of permutations

# SIMPER analyses

# We will automate simper with pre-existing scripts, sadly we cannot include all comparisons at once for it will cause the scripts to break

source("../Results/Scripts/Steinberger_scripts/simper_pretty.r")
source("../Results/Scripts/Steinberger_scripts/R_krusk.r")

#Age 

simper.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), interesting = c("Age"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "MG_age")

MG_age =  data.frame(read.csv("MG_age_clean_simper.csv"))

kruskal.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), csv = MG_age, interesting = c('Age'), output_name =  'MG_age')

KW_MG_age = data.frame(read.csv("MG_Age_krusk_simper.csv"))
KW_MG_age = KW_MG_age[KW_MG_age$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_MG_age = KW_MG_age[with(KW_MG_age, order(SIMPER, decreasing = TRUE)),]
KW_MG_age$OTU = as.factor(KW_MG_age$OTU)

KW_MG_age %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste("OTU =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined)

#AB
simper.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), interesting = c("AB"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "MG_AB")

MG_AB =  data.frame(read.csv("MG_AB_clean_simper.csv"))

kruskal.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), csv = MG_AB, interesting = c('AB'), output_name =  'MG_AB')

KW_MG_AB = data.frame(read.csv("MG_AB_krusk_simper.csv"))
KW_MG_AB = KW_MG_AB[KW_MG_AB$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_MG_AB = KW_MG_AB[with(KW_MG_AB, order(SIMPER, decreasing = TRUE)),]
KW_MG_AB$OTU = as.factor(KW_MG_AB$OTU)

KW_MG_AB %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste("OTU =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined)

#Farms
simper.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), interesting = c("Farm2"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "MG_Farm")

MG_Farm =  data.frame(read.csv("MG_Farm_clean_simper.csv"))

kruskal.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), csv = MG_Farm, interesting = c('Farm2'), output_name =  'MG_Farm')

KW_MG_Farm = data.frame(read.csv("MG_Farm_krusk_simper.csv"))
KW_MG_Farm = KW_MG_Farm[KW_MG_Farm$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_MG_Farm = KW_MG_Farm[with(KW_MG_Farm, order(SIMPER, decreasing = TRUE)),]
KW_MG_Farm$OTU = as.factor(KW_MG_Farm$OTU)

KW_MG_Farm %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste("OTU =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined)

# plots to look at specific ASVs (age)
abund = otu_table(subsetMG)/rowSums(otu_table(subsetMG))*100
boxplot(unlist(data.frame(abund["1624"])) ~ sample_data(subsetMG)$Age, ylab="% Relative abundance", main="OTU1")

# specific test
kruskal.test(unlist(data.frame(otu_table(subsetMG)["817"]), use.names = FALSE) ~ sample_data(subsetMG)$Age)

# declutter R environment by removing objects that no longer serve a purpose
rm(KW_MG_AB, KW_MG_age, KW_MG_Farm, MG_AB, MG_age, MG_Farm, permanova_AB, permanova_agent, permanova_farm, permanova_age, 
   permanova_stable, dist, ord_meths, pcoa_bc, pcoa_jaccard, pcoa_jsd, pcoa_unifrac, pcoa_wunifrac, simper.results, pdataframe, 
   metadf, abund, ps.rel, tse, tse_genus, bray.dist, jaccard.dist, jsd.dist, top.coef, unifrac.dist, wunifrac.dist,
   pcoa_bc2, plist, meta, otu, ps1.rel, coef)
