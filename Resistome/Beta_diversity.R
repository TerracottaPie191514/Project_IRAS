#library(microbiomeDataSets)
library(scater) # plotReducedDim
library(mia) # microbiome analysis package, making tse
library(vegan) # used to run simper
library(nlme) # for usage of llply(), to apply functions over lists

# Used the following guide : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

# Visualizing different kinds of ordination methods with BC distance matrix
#currently struggling to implement DPCoA for some reason, plot_ordination gives an error with DPCoA ordination
plot_ordination(Rps_scaled_copy,  ordinate(Rps_scaled_copy, method = "DPCoA", distance="bray"), "samples", color="Age", shape="AB")

estimate_richness(Rps)

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="Age", shape = "AB")
}, Rps, dist)

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
  ggtitle("Different ordination methods for resistomic data (Bray-Curtis)")

# Repeat for mp

plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="Age", shape = "AB")
}, Rps_mp, dist)

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
  ggtitle("Different ordination methods for resistomic data (Bray-Curtis)")

# Looking at taxa spread as well (irrelevant)
plot_ordination(Rps,  ordinate(Rps, method = "PCoA", distance="bray"), type = "split", color="Age", shape="AB")

# PCoAs for different methods, with Age and Farm as colors, and AB as shape

plot_pcoa_ordination <- function(data, pcoa, var, title) {
  p <- plot_ordination(data, pcoa, color = var, shape = "AB") +
    geom_point(size = 3) +
    labs(title = title, color = var, shape = "Antibiotics used")
  
  return(p)
}

pcoa_bc = ordinate(Rps, "PCoA", "bray")
pcoa_unifrac = ordinate(Rps, "PCoA", "unifrac") 
pcoa_wunifrac = ordinate(Rps, "PCoA", "wunifrac") 
pcoa_jsd = ordinate(Rps, "PCoA", "jsd") 
pcoa_jaccard = ordinate(Rps, "PCoA", "jaccard", binary=TRUE) 

plot_pcoa_ordination(Rps, pcoa_bc, "Age", "PCoA Bray Curtis")
plot_pcoa_ordination(Rps, pcoa_bc, "Farm2", "PCoA Bray Curtis")

plot_pcoa_ordination(Rps, pcoa_unifrac, "Age", "PCoA Unifrac")
plot_pcoa_ordination(Rps, pcoa_unifrac, "Farm2", "PCoA Unifrac")

plot_pcoa_ordination(Rps, pcoa_wunifrac, "Age", "PCoA Weighted Unifrac")
plot_pcoa_ordination(Rps, pcoa_wunifrac, "Farm2", "PCoA Weighted Unifrac")

plot_pcoa_ordination(Rps, pcoa_jsd, "Age", "PCoA Jensen-Shannon Divergence")
plot_pcoa_ordination(Rps, pcoa_jsd, "Farm2", "PCoA Jensen-Shannon Divergence")

plot_pcoa_ordination(Rps, pcoa_jaccard, "Age", "PCoA Jaccard")
plot_pcoa_ordination(Rps, pcoa_jaccard, "Farm2", "PCoA Jaccard")

# repeat for metaphlan
pcoa_bc = ordinate(Rps_mp, "PCoA", "bray")
pcoa_unifrac = ordinate(Rps_mp, "PCoA", "unifrac") 
pcoa_wunifrac = ordinate(Rps_mp, "PCoA", "wunifrac") 
pcoa_jsd = ordinate(Rps_mp, "PCoA", "jsd") 
pcoa_jaccard = ordinate(Rps_mp, "PCoA", "jaccard", binary=TRUE) 

plot_pcoa_ordination(Rps_mp, pcoa_bc, "Age", "PCoA Bray Curtis")
plot_pcoa_ordination(Rps_mp, pcoa_bc, "Farm2", "PCoA Bray Curtis")

plot_pcoa_ordination(Rps_mp, pcoa_unifrac, "Age", "PCoA Unifrac")
plot_pcoa_ordination(Rps_mp, pcoa_unifrac, "Farm2", "PCoA Unifrac")

plot_pcoa_ordination(Rps_mp, pcoa_wunifrac, "Age", "PCoA Weighted Unifrac")
plot_pcoa_ordination(Rps_mp, pcoa_wunifrac, "Farm2", "PCoA Weighted Unifrac")

plot_pcoa_ordination(Rps_mp, pcoa_jsd, "Age", "PCoA Jensen-Shannon Divergence")
plot_pcoa_ordination(Rps_mp, pcoa_jsd, "Farm2", "PCoA Jensen-Shannon Divergence")

plot_pcoa_ordination(Rps_mp, pcoa_jaccard, "Age", "PCoA Jaccard")
plot_pcoa_ordination(Rps_mp, pcoa_jaccard, "Farm2", "PCoA Jaccard")

#plot_ordination(Rps, pcoa_bc, type = "taxa", color = "AMR_class_primary") + 
#  geom_point(size = 3)  + labs(title = "PCoA primary AMR classes", color = "AMR_class_primary")
plot_scree(pcoa_jsd) #scree plots can be made for any of the PCoAs, and are made for those with less than 10% on first axis


plot_ordination(Rps_mp, pcoa_wunifrac, color = "Age", shape = "AB", label = "Sample_Unique") +
  geom_point(size = 3) # 2_57 outlier again
plot_ordination(Rps_mp, pcoa_jsd, color = "Age", shape = "AB", label = "Sample_Unique") +
  geom_point(size = 3) # 2_57 outlier again
 

#plot_ordination(Rps, pcoa_jsd, type = "split", color = "AMR_class_primary", shape = "AB") + 
#  geom_point(size = 3) + labs(title = "PCoA Bray Curtis Farms",color = "Farms", shape = "Antibiotics used")


# plot to look at concentration with a red/green gradient for metaphlan data
# because jaccard is hugely influenced by presence/absence, it will be impacted strongly by singletons
# whenever the data is not rounded, jaccard will not see any of our genes as being present as singletons
# case in point, the following is a plot of metaphlan data which is not rounded
pcoa_jaccard_mp = ordinate(Rps_mp, "PCoA", "jaccard", binary=TRUE) 

plot_ordination(Rps_mp, pcoa_jaccard_mp, color = "Conc...ng..µl.", shape = "AB", label = "Stable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard concentration",color = "Concentration", shape = "Antibiotics used") +
  scale_colour_gradient(low = "red", high = "green")
# there are also no differences to be found when comparing FPKM/ TPM or metaphlan and kraken2, this is a plot with kraken2 TPM data:
pcoa_jaccard_tpm = ordinate(Rps_tpm, "PCoA", "jaccard", binary=TRUE) 
plot_ordination(Rps_tpm, pcoa_jaccard_tpm, color = "Conc...ng..µl.", shape = "AB", label = "Stable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard concentration",color = "Concentration", shape = "Antibiotics used") +
  scale_colour_gradient(low = "red", high = "green")
# now, when rounding our data a completely different picture emerges
Rps_mp_rounded = Rps_mp
otu_table(Rps_mp_rounded) = otu_table(round(as((otu_table(Rps_mp)), "matrix")), taxa_are_rows(Rps_mp))
pcoa_jaccard_mp = ordinate(Rps_mp_rounded, "PCoA", "jaccard", binary=TRUE) 
plot_ordination(Rps_mp_rounded, pcoa_jaccard_mp, color = "Conc...ng..µl.", shape = "AB", label = "Stable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard concentration",color = "Concentration", shape = "Antibiotics used") +
  scale_colour_gradient(low = "red", high = "green")


pcoa_jaccard_k2 = ordinate(Rps_rounded, "PCoA", "jaccard", binary=TRUE) 
plot_ordination(Rps_rounded, pcoa_jaccard_k2, color = "Conc...ng..µl.", shape = "AB", label = "Stable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard concentration",color = "Concentration", shape = "Antibiotics used") +
  scale_colour_gradient(low = "red", high = "green")

Rps_k2_scaled = Rps
otu_table(Rps_k2_scaled) = otu_table(round(as((otu_table(Rps_k2_scaled)), "matrix")), taxa_are_rows(Rps_k2_scaled))
pcoa_jaccard_k2 = ordinate(Rps_k2_scaled, "PCoA", "jaccard", binary=TRUE) 
plot_ordination(Rps_k2_scaled, pcoa_jaccard_k2, color = "Conc...ng..µl.", shape = "AB", label = "Stable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard concentration",color = "Concentration", shape = "Antibiotics used") +
  scale_colour_gradient(low = "red", high = "green")

pcoa_jaccard_mp_scaled = ordinate(Rps_mp_scaled_copy, "PCoA", "jaccard") 

# BC plots for looking at percentage and total amount of bacterial reads mapped
Rps %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(Rps, "PCoA", "bray") , color = "ReadPerc", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC percentage",color = "ReadPerc") +
  scale_colour_viridis_c()

Rps %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(Rps, "PCoA", "bray") , color = "ReadTot", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC total",color = "ReadTot") +
  scale_colour_viridis_c()

# metaphlan
Rps_mp %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(Rps_mp, "PCoA", "bray") , color = "ReadPerc", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC percentage",color = "ReadPerc") +
  scale_colour_viridis_c()

Rps_mp %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(Rps_mp, "PCoA", "bray") , color = "ReadTot", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC total",color = "ReadTot") +
  scale_colour_viridis_c()


# Alternative way of plotting
plot_ordination(Rps, 
                pcoa_wunifrac, color="Farm2") + ggtitle("Weighted UniFrac") + geom_point(size = 2) + 
  theme_classic() + scale_color_brewer("Farm2", palette = "Set2")

# NMDS

tse2 = makeTreeSummarizedExperimentFromPhyloseq(Rps)
tse2 %<>%  transformCounts( method = "relabundance")
tse2 %<>% runNMDS(FUN = vegan::vegdist, name = "BC", nmdsFUN = "monoMDS",
                    exprs_values = "relabundance",
                    keep_dist = TRUE)

tse2 %>% plotReducedDim("BC", colour_by = "Age") 


# PERMANOVAs
tse = makeTreeSummarizedExperimentFromPhyloseq(Rps_copy)
tse <- transformCounts(tse, method = "relabundance")

adonis2(t(assay(tse, "relabundance")) ~ AB, data = colData(tse), permutations = 9999) # not significant
adonis2(t(assay(tse, "relabundance")) ~ Cox, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Researcher, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ FeedProducent, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ LitterType, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ FeedType, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Gender, data = colData(tse), permutations = 9999) # not significant
adonis2(t(assay(tse, "relabundance")) ~ Stables, data = colData(tse), permutations = 9999) 
adonis2(t(assay(tse, "relabundance")) ~ FlockSize, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Farm2, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ AgeParentStock, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Age, data = colData(tse), permutations = 9999)

# variances: AB: 0.0143, Cox: 0.1158, Researcher: 0.090, FP : 0.1224 , LitterType: 0.179, FT :0.030, Gender: 0.0067, 
# Stable: 0.293, FS: 0.224 , Farm 0.207, APS : 0.246, Age: 0.0302
# Order: Stable>APS>FS>Farm>LT>FP>Cox>Researcher>Age>FT>AB>Gender

# Mixed models
adonis2(t(assay(tse, "relabundance")) ~ Stables * AB, data = colData(tse), permutations = 9999) 


# basically, composition seems to be different over every single variable, except for gender

# on clust90 level
tse_clust90 <- agglomerateByRank(tse, "Order")
tse_clust90 <- transformCounts(tse_clust90, method = "relabundance")

adonis2(t(assay(tse_clust90, "relabundance")) ~ AB, data = colData(tse_clust90), permutations = 9999) # not sign
adonis2(t(assay(tse_clust90, "relabundance")) ~ Cox, data = colData(tse_clust90), permutations = 9999)
adonis2(t(assay(tse_clust90, "relabundance")) ~ Researcher, data = colData(tse_clust90), permutations = 9999)
adonis2(t(assay(tse_clust90, "relabundance")) ~ FeedProducent, data = colData(tse_clust90), permutations = 9999)
adonis2(t(assay(tse_clust90, "relabundance")) ~ LitterType, data = colData(tse_clust90), permutations = 9999)
adonis2(t(assay(tse_clust90, "relabundance")) ~ FeedType, data = colData(tse_clust90), permutations = 9999)
adonis2(t(assay(tse_clust90, "relabundance")) ~ Gender, data = colData(tse_clust90), permutations = 9999) # NIET significant
adonis2(t(assay(tse_clust90, "relabundance")) ~ Stables, data = colData(tse_clust90), permutations = 9999)
adonis2(t(assay(tse_clust90, "relabundance")) ~ FlockSize, data = colData(tse_clust90), permutations = 9999)
adonis2(t(assay(tse_clust90, "relabundance")) ~ Farm2, data = colData(tse_clust90), permutations = 9999)
adonis2(t(assay(tse_clust90, "relabundance")) ~ AgeParentStock, data = colData(tse_clust90), permutations = 9999)

# same trends on genus level


# for different ordination methods
ps1.rel <- microbiome::transform(Rps, "compositional")
metadf <- data.frame(sample_data(ps1.rel))

# alternative calculations
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


# same patterns arise, but AB is sign

wunifrac.dist <- UniFrac(ps1.rel, 
                         weighted = TRUE)

adonis2(wunifrac.dist ~ Age, data = metadf)
adonis2(wunifrac.dist ~ AB, data = metadf) # not sign
adonis2(wunifrac.dist ~ Farm2, data = metadf)
adonis2(wunifrac.dist ~ Cox, data = metadf)
adonis2(wunifrac.dist ~ Researcher, data = metadf)
adonis2(wunifrac.dist ~ LitterType, data = metadf)
adonis2(wunifrac.dist ~ Gender, data = metadf) # not sign
adonis2(wunifrac.dist ~ Stables, data = metadf)


#  sameish patterns

meta <- meta(ps1.rel)
adonis2(t(abundances(ps1.rel)) ~ Age, data = meta, permutations=9999, method = "bray")
adonis2(t(abundances(ps1.rel)) ~ AB, data = meta, permutations=9999, method = "bray") # not sign
adonis2(t(abundances(ps1.rel)) ~ Farm2, data = meta, permutations=9999, method = "bray")
adonis2(t(abundances(ps1.rel)) ~ Cox, data = meta, permutations=9999, method = "bray")
adonis2(t(abundances(ps1.rel)) ~ Researcher, data = meta, permutations=9999, method = "bray")
adonis2(t(abundances(ps1.rel)) ~ LitterType, data = meta, permutations=9999, method = "bray")
adonis2(t(abundances(ps1.rel)) ~ Gender, data = meta, permutations=9999, method = "bray") # not sign
adonis2(t(abundances(ps1.rel)) ~ Stables, data = meta, permutations=9999, method = "bray")


# and BC

adonis2(t(abundances(ps1.rel)) ~ Age, data = meta, permutations=9999, method = "jaccard")
adonis2(t(abundances(ps1.rel)) ~ AB, data = meta, permutations=9999, method = "jaccard") # not sign
adonis2(t(abundances(ps1.rel)) ~ Farm2, data = meta, permutations=9999, method = "jaccard")
adonis2(t(abundances(ps1.rel)) ~ Cox, data = meta, permutations=9999, method = "jaccard")
adonis2(t(abundances(ps1.rel)) ~ Researcher, data = meta, permutations=9999, method = "jaccard")
adonis2(t(abundances(ps1.rel)) ~ LitterType, data = meta, permutations=9999, method = "jaccard")
adonis2(t(abundances(ps1.rel)) ~ Gender, data = meta, permutations=9999, method = "jaccard") # not sign
adonis2(t(abundances(ps1.rel)) ~ Stables, data = meta, permutations=9999, method = "jaccard")

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
          xlab = "ARGs",
          legend.title = "ARG contribution",
          title = "Impact of bacterial ARGs on age",
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
          xlab = "ARGs",
          legend.title = "ARG contribution",
          title = "Impact of bacterial ARGs on AB",
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
          xlab = "ARGs",
          legend.title = "ARG contribution",
          title = "Impact of bacterial ARGs on Stable",
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
          xlab = "ARGs",
          legend.title = "ARG contribution",
          title = "Impact of bacterial ARGs on Farm",
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
          xlab = "ARGs",
          legend.title = "ARG contribution",
          title = "Impact of bacterial ARGs on agent",
          rotate = TRUE,
          ggtheme = theme_minimal())

# same plots but for argclust90

#  Age

permanova_age <- adonis(t(assay(tse_clust90, "relabundance")) ~ Age, data = colData(tse_clust90), permutations = 9999)

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
          xlab = "ARGClust90",
          legend.title = "ARGClust90 contribution",
          title = "Impact of ARGClust90 on age",
          rotate = TRUE,
          ggtheme = theme_minimal())

# AB

permanova_AB <- adonis(t(assay(tse_clust90, "relabundance")) ~ AB, data = colData(tse_clust90), permutations = 9999)

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
          xlab = "ARGClust90",
          legend.title = "ARGClust90 contribution",
          title = "Impact of ARGClust90 on AB",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Stable

permanova_stable <- adonis(t(assay(tse_clust90, "relabundance")) ~ Stables, data = colData(tse_clust90), permutations = 9999)

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
          xlab = "ARGClust90",
          legend.title = "ARGClust90 contribution",
          title = "Impact of ARGClust90 on Stable",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Farm

permanova_farm <- adonis(t(assay(tse_clust90, "relabundance")) ~ Farm2, data = colData(tse_clust90), permutations = 9999)

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
          xlab = "ARGClust90",
          legend.title = "ARGClust90 contribution",
          title = "Impact of ARGClust90 on Farm",
          rotate = TRUE,
          ggtheme = theme_minimal())

# Agent

permanova_agent <- adonis(t(assay(tse_clust90, "relabundance")) ~ Cox, data = colData(tse_clust90), permutations = 9999)

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
          xlab = "ARGClust90",
          legend.title = "ARGClust90 contribution",
          title = "Impact of ARGClust90 on agent",
          rotate = TRUE,
          ggtheme = theme_minimal())


# checking homogeneity condition - bray
# ANOVAs are performed on betadispers of our rel abund data to test whether groups are more variable than others

# Bray
ps.rel = microbiome::transform(Rps, "compositional")
meta = meta(ps.rel)
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Age)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$AB))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Farm2)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Stables))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Cox)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Researcher))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$LitterType)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$Gender)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FlockSize)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$AgeParentStock))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FeedProducent))
anova(betadisper(vegdist(t(abundances(ps.rel))), meta$FeedType)) # homogeneous

# Jaccard
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Age)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$AB))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Farm2)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Stables))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Researcher))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$LitterType)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$Gender)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FlockSize)) # homogeneous
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$AgeParentStock))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FeedProducent))
anova(betadisper(vegdist(t(abundances(ps.rel)), method = "jaccard"), meta$FeedType)) # homogeneous

# group variances are not homogenous in most cases, so there are differences in variances between groups -> bad for anova / permanova

# Tukey tests can be performed to see if and which groups differ in relation to variance

TukeyHSD(betadisper(vegdist(t(abundances(ps.rel))), meta$Farm2))


# different way of calculating homogeneity, permutation tests, null = no difference in dispersion between groups 
permutest(betadisper(vegdist(t(abundances(ps.rel))), meta$Age), pairwise = TRUE)

permutest(betadisper(unifrac.dist, metadf$Age), pairwise = TRUE) # looks like unifrac distances are homogenously dispersed for age
permutest(betadisper(unifrac.dist, metadf$AB), pairwise = TRUE) # not for AB though

# Simper analyses to see which species are most impactful to BC dissimilarity between groups
# Test for significance (OTU abundance will not be normally distributed so we will use kruskal wallis tests)

source("../Results/Scripts/Steinberger_scripts/simper_pretty.r")
source("../Results/Scripts/Steinberger_scripts/R_krusk.r")

#Age 
simper.pretty(otu_table(Rps), metrics = sample_data(Rps), interesting = c("Age"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "Rps_age")

Rps_age =  data.frame(read.csv("Rps_age_clean_simper.csv"))

kruskal.pretty(otu_table(Rps), metrics = sample_data(Rps), csv = Rps_age, interesting = c('Age'), output_name =  'Rps_age')

KW_Rps_age = data.frame(read.csv("Rps_Age_krusk_simper.csv"))
KW_Rps_age = KW_Rps_age[KW_Rps_age$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_Rps_age = KW_Rps_age[with(KW_Rps_age, order(SIMPER, decreasing = TRUE)),]
KW_Rps_age$OTU = as.factor(KW_Rps_age$OTU)

KW_Rps_age %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste("ARG =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined) 

#AB
simper.pretty(otu_table(Rps), metrics = sample_data(Rps), interesting = c("AB"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "Rps_AB")

Rps_AB =  data.frame(read.csv("Rps_AB_clean_simper.csv"))

kruskal.pretty(otu_table(Rps), metrics = sample_data(Rps), csv = Rps_AB, interesting = c('AB'), output_name =  'Rps_AB')

KW_Rps_AB = data.frame(read.csv("Rps_AB_krusk_simper.csv"))
KW_Rps_AB = KW_Rps_AB[KW_Rps_AB$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_Rps_AB = KW_Rps_AB[with(KW_Rps_AB, order(SIMPER, decreasing = TRUE)),]
KW_Rps_AB$OTU = as.factor(KW_Rps_AB$OTU)

KW_Rps_AB %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>% 
  rowwise() %>% mutate(Combined = paste("ARG =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined)

#Farms - too many comparisons so maybe too extensive for report

simper.pretty(otu_table(Rps), metrics = sample_data(Rps), interesting = c("Farm2"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "Rps_Farm")

Rps_Farm =  data.frame(read.csv("Rps_Farm_clean_simper.csv"))

kruskal.pretty(otu_table(Rps), metrics = sample_data(Rps), csv = Rps_Farm, interesting = c('Farm2'), output_name =  'Rps_Farm')

KW_Rps_Farm = data.frame(read.csv("Rps_Farm_krusk_simper.csv"))
KW_Rps_Farm = KW_Rps_Farm[KW_Rps_Farm$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_Rps_Farm = KW_Rps_Farm[with(KW_Rps_Farm, order(SIMPER, decreasing = TRUE)),]
KW_Rps_Farm$OTU = as.factor(KW_Rps_Farm$OTU)

KW_Rps_Farm %>% mutate_if(is.numeric, format.pval, 2) %>% dplyr::select("Comparison", "SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste(Comparison, "ARG =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined) 

# plots to look at specific ARGs (age)
abund = otu_table(Rps)/rowSums(otu_table(Rps))*100
boxplot(unlist(data.frame(abund["tet(44)_2_FN594949"])) ~ sample_data(Rps)$Age, ylab="% Relative abundance", main="ARG")

# specific test
kruskal.test(unlist(data.frame(otu_table(Rps)["tet(O/32/O)_5_FP929050"]), use.names = FALSE) ~ sample_data(Rps)$Age)

# declutter R environment by removing objects that no longer serve a purpose
rm(KW.results, dist, ord_meths, pcoa_bc, pcoa_jaccard, pcoa_jsd, pcoa_unifrac, pcoa_wunifrac, simper.results, pdataframe, metadf)
