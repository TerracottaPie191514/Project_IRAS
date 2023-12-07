library(scater) # plotReducedDim
library(vegan) # used to run simper
library(nlme) # for usage of llply(), to apply functions over lists

# Used the following guides: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html, https://rfunctions.blogspot.com/2019/03/betadisper-and-adonis-homogeneity-of.html


# single ordination method
#plot_ordination(subset16S,  ordinate(subset16S, method = "DPCoA", distance="bray"), "samples", color="Age", shape="AB")

estimate_richness(subset16S) # no singletons

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA", "DPCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(subset16S, method=i, distance=dist)
  plot_ordination(subset16S, ordi, "samples", color="Age", shape = "AB")
}, subset16S, dist)

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
  ggtitle("Different ordination methods for 16S data (Bray-Curtis)")


# PCoAs for different methods (fancy maken : + theme_classic() + scale_color_brewer("Farm2", palette = "Set2"))

# functionize plotting pcoa
plot_pcoa_ordination <- function(data, pcoa, var, title) {
  p <- plot_ordination(data, pcoa, color = var, shape = "AB") +
    geom_point(size = 3) +
    labs(title = title, color = var, shape = "Antibiotics used")
  
  return(p)
}

pcoa_bc = ordinate(subset16S, "PCoA", "bray") 
pcoa_unifrac = ordinate(subset16S, "PCoA", "unifrac") 
pcoa_wunifrac = ordinate(subset16S, "PCoA", "wunifrac") 
pcoa_jsd = ordinate(subset16S, "PCoA", "jsd") 
pcoa_jaccard = ordinate(subset16S, "PCoA", "jaccard", binary=TRUE) 


plot_pcoa_ordination(subset16S, pcoa_bc, "Age", "PCoA Bray Curtis")
plot_pcoa_ordination(subset16S, pcoa_bc, "Farm2", "PCoA Bray Curtis")

# proper order of legend:
plot_ordination(subset16S, pcoa_bc, color = "Farm2", shape = "AB") +
  geom_point(size = 3) +
  labs(title = "PCoA Bray Curtis", color = "Farm", shape = "Antibiotics used")

plot_pcoa_ordination(subset16S, pcoa_unifrac, "Age", "PCoA Unifrac")
plot_pcoa_ordination(subset16S, pcoa_unifrac, "Farm2", "PCoA Unifrac")

plot_pcoa_ordination(subset16S, pcoa_wunifrac, "Age", "PCoA Weighted Unifrac")
plot_pcoa_ordination(subset16S, pcoa_wunifrac, "Farm2", "PCoA Weighted Unifrac")

plot_pcoa_ordination(subset16S, pcoa_jsd, "Age", "PCoA Jensen-Shannon Divergence")
plot_pcoa_ordination(subset16S, pcoa_jsd, "Farm2", "PCoA Jensen-Shannon Divergence")

plot_pcoa_ordination(subset16S, pcoa_jaccard, "Age", "PCoA Jaccard")
plot_pcoa_ordination(subset16S, pcoa_jaccard, "Farm2", "PCoA Jaccard")

#plot_ordination(subset16S, pcoa_jaccard, color = "Age", shape = "AB", label = "Stables") +  
#  geom_point(size = 3)  + labs(title = "PCoA Jaccard Age",color = "Age", shape = "Antibiotics used")


plot_scree(pcoa_bc) #scree plots can be made for any of the PCoAs, those that explain less than 10% of variance on first axis
plot_scree(pcoa_jaccard)

# NMDS

tse2 = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse2 %<>%  transformCounts( method = "relabundance")
tse2 %<>% runNMDS(FUN = vegan::vegdist, name = "BC", nmdsFUN = "monoMDS",
                  exprs_values = "relabundance",
                  keep_dist = TRUE)

tse2 %>% plotReducedDim("BC", colour_by = "Age") 

# PERMANOVAs

tse = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse <- transformCounts(tse, method = "relabundance")

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
adonis2(t(assay(tse, "relabundance")) ~ Age, data = colData(tse), permutations = 9999)

# variances: AB: 0.026, Cox: 0.102, Researcher: 0.06, FP : 0.067, LitterType: 0.061, FT :0.055, Gender: 0.007, 
# Stable: 0.167, FS: 0.1245, Farm 0.103, APS : 0.118, Age: 0.054
# Order: Stable>FS>APS>Farm>Cox>FP>LT>Researcher>FT>Age>AB>Gender

# Mixed models
adonis2(t(assay(tse, "relabundance")) ~ Stables * AB, data = colData(tse), permutations = 9999) 


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
adonis2(t(assay(tse_genus, "relabundance")) ~ Stables, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ FlockSize, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ Farm2, data = colData(tse_genus), permutations = 9999)
adonis2(t(assay(tse_genus, "relabundance")) ~ AgeParentStock, data = colData(tse_genus), permutations = 9999)

# same trends on genus level (and on phylum level, though p values become higher)


# for different ordination methods
ps1.rel <- microbiome::transform(subset16S, "compositional")
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
          xlab = "ASVs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial ASVs on age",
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
          xlab = "ASVs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial ASVs on AB",
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
          xlab = "ASVs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial ASVs on Stable",
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
          xlab = "ASVs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial ASVs on Farm",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Agent

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
          xlab = "ASVs",
          legend.title = "ASV contribution",
          title = "Impact of bacterial ASVs on agent",
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



permanova_age <- adonis(unifrac.dist ~ Age, data = metadf)
permanova_AB <- adonis2(unifrac.dist ~ AB, data = metadf)
permanova_farm <- adonis2(unifrac.dist ~ Farm2, data = metadf)
permanova_cox <- adonis2(unifrac.dist ~ Cox, data = metadf)
permanova_researcher <- adonis2(unifrac.dist ~ Researcher, data = metadf)
permanova_LitterType <- adonis2(unifrac.dist ~ LitterType, data = metadf)
permanova_cox <- adonis2(unifrac.dist ~ Cox, data = metadf)

# checking homogeneity condition - bray
# ANOVAs are performed on betadispers of our rel abund data to test whether groups are more variable than others

# Bray
ps.rel = microbiome::transform(subset16S, "compositional")
meta = meta(ps.rel)
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Age))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$AB))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Farm2))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Stables))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Cox))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Researcher))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$LitterType)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Gender)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FlockSize))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$AgeParentStock)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FeedProducent))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FeedType))

# Jaccard
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Age))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$AB))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Farm2))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Stables))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Researcher))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$LitterType)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Gender)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FlockSize))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$AgeParentStock)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FeedProducent))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FeedType))

# group variances are not homogenous in most cases, so there are differences in variances between groups -> bad for anova / permanova

# Tukey tests can be performed to see if and which groups differ in relation to variance

TukeyHSD(betadisper(vegdist(t(abundances(ps.rel))), meta$Farm2))


# different way of calculating homogeneity, permutation tests, null = no difference in dispersion between groups 
permutest(betadisper(vegdist(t(abundances(ps.rel))), meta$Age), pairwise = TRUE)

permutest(betadisper(unifrac.dist, metadf$Age), pairwise = TRUE) # looks like unifrac distances are homogenously dispersed for age
permutest(betadisper(unifrac.dist, metadf$AB), pairwise = TRUE) # not for AB though

permutest(betadisper(bray.dist, metadf$Age), pairwise = TRUE) # there are differences in P value with other method, but could be number of permutations

# SIMPER - we'll use MT as abbreviation for metataxonomics instead of 16s since R does not like its objects starting with numbers

source("../Results/Scripts/Steinberger_scripts/simper_pretty.r")
source("../Results/Scripts/Steinberger_scripts/R_krusk.r")

#Age 

simper.pretty(otu_table(subset16S), metrics = sample_data(subset16S), interesting = c("Age"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "MT_age")

MT_age =  data.frame(read.csv("MT_age_clean_simper.csv"))

kruskal.pretty(otu_table(subset16S), metrics = sample_data(subset16S), csv = MT_age, interesting = c('Age'), output_name =  'MT_age')

KW_MT_age = data.frame(read.csv("MT_Age_krusk_simper.csv"))
KW_MT_age = KW_MT_age[KW_MT_age$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_MT_age = KW_MT_age[with(KW_MT_age, order(SIMPER, decreasing = TRUE)),]
KW_MT_age$OTU = as.factor(KW_MT_age$OTU)

KW_MT_age %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste("ASV =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined) 

#AB
simper.pretty(otu_table(subset16S), metrics = sample_data(subset16S), interesting = c("AB"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "MT_AB")

MT_AB =  data.frame(read.csv("MT_AB_clean_simper.csv"))

kruskal.pretty(otu_table(subset16S), metrics = sample_data(subset16S), csv = MT_AB, interesting = c('AB'), output_name =  'MT_AB')

KW_MT_AB = data.frame(read.csv("MT_AB_krusk_simper.csv"))
KW_MT_AB = KW_MT_AB[KW_MT_AB$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_MT_AB = KW_MT_AB[with(KW_MT_AB, order(SIMPER, decreasing = TRUE)),]
KW_MT_AB$OTU = as.factor(KW_MT_AB$OTU)

KW_MT_AB %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>% 
  rowwise() %>% mutate(Combined = paste("ASV =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined)

#Farms - too many comparisons so maybe too extensive for report

simper.pretty(otu_table(subset16S), metrics = sample_data(subset16S), interesting = c("Farm2"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "MT_Farm")

MT_Farm =  data.frame(read.csv("MT_Farm_clean_simper.csv"))

kruskal.pretty(otu_table(subset16S), metrics = sample_data(subset16S), csv = MT_Farm, interesting = c('Farm2'), output_name =  'MT_Farm')

KW_MT_Farm = data.frame(read.csv("MT_Farm_krusk_simper.csv"))
KW_MT_Farm = KW_MT_Farm[KW_MT_Farm$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_MT_Farm = KW_MT_Farm[with(KW_MT_Farm, order(SIMPER, decreasing = TRUE)),]
KW_MT_Farm$OTU = as.factor(KW_MT_Farm$OTU)

KW_MT_Farm %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("Comparison", "SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste(Comparison, "ASV =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined) 

# plots to look at specific ASVs (age)
abund = otu_table(subset16S)/rowSums(otu_table(subset16S))*100
boxplot(unlist(data.frame(abund["224597762"])) ~ sample_data(subset16S)$Age, ylab="% Relative abundance", main="OTU1")

# specific test
kruskal.test(unlist(data.frame(otu_table(subset16S)["224597762"]), use.names = FALSE) ~ sample_data(subset16S)$Age)

# declutter R environment by removing objects that no longer serve a purpose
rm(KW_MT_AB, KW_MT_age, KW_MT_Farm, MT_AB, MT_age, MT_Farm, permanova_AB, permanova_agent, permanova_farm, permanova_age,
   permanova_stable, dist, ord_meths, pcoa_bc, pcoa_jaccard, pcoa_jsd, pcoa_unifrac, pcoa_wunifrac, simper.results,
   pdataframe, metadf, ps1.rel, df, tse, tse2, tse_genus, plist, bray.dist, unifrac.dist, wunifrac.dist, jsd.dist,
   jaccard.dist, coef, top.coef, meta, ps.rel, abund)