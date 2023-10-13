#### Load packages
library(scater) # For functions like plotReducedDim(), calculating dissimiilarity matrices etc. 
library(mia) # microbiome analysis package, making tse
library(vegan) # used to run simper
library(plyr) # for llply, to apply functions
library(nlme) # 

# Used the following guide : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

# single ordination visualisation

plot_ordination(subsetMG,  ordinate(subsetMG, method = "DPCoA", distance="bray"), "samples", color="Age", shape="AB")

# Different ordination methods based on BC dissimilarity

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
  ggtitle("Different ordination methods for 16S data (Bray-Curtis)")

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

adonis2(t(assay(tse, "relabundance")) ~ AB, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Cox, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Researcher, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ FeedProducent, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ LitterType, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ FeedType, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Gender, data = colData(tse), permutations = 9999) # NIET significant
adonis2(t(assay(tse, "relabundance")) ~ FarmRoundStable, data = colData(tse), permutations = 9999) 
adonis2(t(assay(tse, "relabundance")) ~ FlockSize, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Farm2, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ AgeParentStock, data = colData(tse), permutations = 9999)

# basically, composition seems to be different over every single variable, except for gender

# on genus level
tse_genus <- agglomerateByRank(tse, "Genus")
tse_genus <- transformCounts(tse_genus, method = "relabundance")

adonis2(t(assay(tse, "relabundance")) ~ AB, data = colData(tse_genus), permutations = 9999)

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

# same results on genus level (and on phylum level, though p values become higher)

# for unifrac and wunifrac
ps1.rel <- microbiome::transform(subsetMG, "compositional")
otu <- abundances(ps1.rel)
meta <- meta(ps1.rel)

adonis2(t(otu) ~ Age, data = meta, permutations=9999, method = "bray")

adonis2(bray.dist ~ Age, data = metadf)

permanova = adonis(t(otu) ~ Age, data = meta, permutations=9999, method = "bray")
permanova$aov.tab

unifrac.dist <- UniFrac(ps1.rel)

adonis2(unifrac.dist ~ Age, data = metadf)
adonis2(unifrac.dist ~ AB, data = metadf)
adonis2(unifrac.dist ~ Farm2, data = metadf)
adonis2(unifrac.dist ~ Cox, data = metadf)
adonis2(unifrac.dist ~ Researcher, data = metadf)
adonis2(unifrac.dist ~ LitterType, data = metadf)
adonis2(unifrac.dist ~ Gender, data = metadf)
adonis2(unifrac.dist ~ FarmRoundStable, data = metadf)


# same patterns arise

wunifrac.dist <- UniFrac(ps1.rel, 
                         weighted = TRUE)

adonis2(wunifrac.dist ~ Age, data = metadf)
adonis2(wunifrac.dist ~ AB, data = metadf) # NOT significant
adonis2(wunifrac.dist ~ Farm2, data = metadf)
adonis2(wunifrac.dist ~ Cox, data = metadf)
adonis2(wunifrac.dist ~ Researcher, data = metadf)
adonis2(wunifrac.dist ~ LitterType, data = metadf)
adonis2(wunifrac.dist ~ Gender, data = metadf)
adonis2(wunifrac.dist ~ FarmRoundStable, data = metadf)


#  wunifrac also sees no significant difference between AB and non AB!

jsd.dist <- distance(ps1.rel, "jsd")

adonis2(jsd.dist ~ Age, data = metadf)
adonis2(jsd.dist ~ AB, data = metadf) # NOT significant
adonis2(jsd.dist ~ Farm2, data = metadf)
adonis2(jsd.dist ~ Cox, data = metadf)
adonis2(jsd.dist ~ Researcher, data = metadf)
adonis2(jsd.dist ~ LitterType, data = metadf)
adonis2(jsd.dist ~ Gender, data = metadf)
adonis2(jsd.dist ~ FarmRoundStable, data = metadf)

# same is true for JSD

bray.dist <- distance(ps1.rel, "bray")
bray.dist <- distance(subsetMG, "bray")


adonis2(bray.dist ~ Age, data = metadf)
adonis2(bray.dist ~ AB, data = metadf, permutations = 9999) # NOT significant
adonis2(bray.dist ~ Farm2, data = metadf)
adonis2(bray.dist ~ Cox, data = metadf)
adonis2(bray.dist ~ Researcher, data = metadf)
adonis2(bray.dist ~ LitterType, data = metadf)
adonis2(bray.dist ~ Gender, data = metadf)
adonis2(bray.dist ~ FarmRoundStable, data = metadf)

# and BC

jaccard.dist <- distance(ps1.rel, "jaccard")

adonis2(jaccard.dist ~ Age, data = metadf)
adonis2(jaccard.dist ~ AB, data = metadf) # NOT significant
adonis2(jaccard.dist ~ Farm2, data = metadf)
adonis2(jaccard.dist ~ Cox, data = metadf)
adonis2(jaccard.dist ~ Researcher, data = metadf)
adonis2(jaccard.dist ~ LitterType, data = metadf)
adonis2(jaccard.dist ~ Gender, data = metadf)
adonis2(jaccard.dist ~ FarmRoundStable, data = metadf)

# as well as jaccard..


adonis2(dist(otu_table(subsetMG), method='euclidean') ~ Age, 
       data=metadf)



metadf <- data.frame(sample_data(ps1.rel))



ps.disper <- betadisper(unifrac.dist, metadf$Age)
permutest(ps.disper, pairwise = TRUE)





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
head(KW_MG_age)

KW_MG_age %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste("ASV =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined) 

#AB
simper.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), interesting = c("AB"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "MG_AB")

MG_AB =  data.frame(read.csv("MG_AB_clean_simper.csv"))

kruskal.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), csv = MG_AB, interesting = c('AB'), output_name =  'MG_AB')

KW_MG_AB = data.frame(read.csv("MG_AB_krusk_simper.csv"))
KW_MG_AB = KW_MG_AB[KW_MG_AB$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_MG_AB = KW_MG_AB[with(KW_MG_AB, order(SIMPER, decreasing = TRUE)),]
head(KW_MG_AB)

KW_MG_AB %>% dplyr::select("SIMPER", "OTU", "fdr_krusk_p.val") %>% #round_df(3) %>%
  rowwise() %>% mutate(Combined = paste("ASV =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined) 



format.pval(KW_MG_AB$fdr_krusk_p.val, digits = 4) # rounding p-values to proper digits should be done

#Farms

simper.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), interesting = c("Farm2"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "MG_Farm")

MG_Farm =  data.frame(read.csv("MG_Farm_clean_simper.csv"))

kruskal.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), csv = MG_Farm, interesting = c('Farm2'), output_name =  'MG_Farm')

KW_MG_Farm = data.frame(read.csv("MG_Farm_krusk_simper.csv"))
KW_MG_Farm = KW_MG_Farm[KW_MG_Farm$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr
KW_MG_Farm = KW_MG_Farm[with(KW_MG_Farm, order(SIMPER, decreasing = TRUE)),]
head(KW_MG_Farm)

KW_MG_Farm %>% dplyr::select("Comparison", "SIMPER", "OTU", "fdr_krusk_p.val") %>%
  rowwise() %>% mutate(Combined = paste(Comparison, "ASV =", OTU, ", SIMPER =", SIMPER, ", p-value =", fdr_krusk_p.val)) %>% 
  dplyr::select(Combined) 

abund = otu_table(subsetMG)/rowSums(otu_table(subsetMG))*100


boxplot(unlist(data.frame(abund["1624"])) ~ sample_data(subsetMG)$Age, ylab="% Relative abundance", main="OTU1")


for (otu in KW.results$OTU) {
  print(otu)
  
}

kruskal.test(unlist(data.frame(otu_table(Rps)["tet(O/32/O)_5_FP929050"]), use.names = FALSE) ~ sample_data(Rps)$Age)

kruskal.test(unlist(data.frame(otu_table(Rps)["tet(O/W/32/O)_1_EF065523"]), use.names = FALSE) ~ sample_data(Rps)$Age)


# declutter R environment by removing objects that no longer serve a purpose
rm(KW.results, dist, ord_meths, pcoa_bc, pcoa_jaccard, pcoa_jsd, pcoa_unifrac, pcoa_wunifrac, simper.results, pdataframe, metadf)

