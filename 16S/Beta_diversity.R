#library(microbiomeDataSets)
library(scater) # plotReducedDim
library(mia) # microbiome analysis package, making tse
library(vegan) # used to run simper
library(nlme) # for usage of llply(), to apply functions over lists

# Used the following guide : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

plot_ordination(subset16S,  ordinate(subset16S, method = "DPCoA", distance="bray"), "samples", color="Age", shape="AB")

estimate_richness(subset16S)

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



filt.rar=data.frame(otu_table(subset16S))

dist_bc <- as.matrix(vegdist(filt.rar, method = "bray")) 

dist_bc[1:5, 1:5]
#sample_data(subset16S)$Age = as.factor(sample_data(subset16S)$Age)


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

plot_pcoa_ordination(subset16S, pcoa_unifrac, "Age", "PCoA Unifrac")
plot_pcoa_ordination(subset16S, pcoa_unifrac, "Farm2", "PCoA Unifrac")

plot_pcoa_ordination(subset16S, pcoa_wunifrac, "Age", "PCoA Weighted Unifrac")
plot_pcoa_ordination(subset16S, pcoa_wunifrac, "Farm2", "PCoA Weighted Unifrac")

plot_pcoa_ordination(subset16S, pcoa_jsd, "Age", "PCoA Jensen-Shannon Divergence")
plot_pcoa_ordination(subset16S, pcoa_jsd, "Farm2", "PCoA Jensen-Shannon Divergence")

plot_pcoa_ordination(subset16S, pcoa_jaccard, "Age", "PCoA Jaccard")
plot_pcoa_ordination(subset16S, pcoa_jaccard, "Farm2", "PCoA Jaccard")

plot_ordination(subset16S, pcoa_jaccard, color = "Age", shape = "AB", label = "FarmRoundStable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard Age",color = "Age", shape = "Antibiotics used")


plot_scree(pcoa_bc) #scree plots can be made for any of the PCoAs



# plots for looking at percentage and total amount of bacterial reads mapped


# Jaccard
Rps %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(Rps, "PCoA", "jaccard") , color = "ReadPerc", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard percentage",color = "ReadPerc") +
  scale_colour_viridis_c()

Rps %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(Rps, "PCoA", "jaccard") , color = "ReadTot", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard total",color = "ReadTot") +
  scale_colour_viridis_c()
# BC 
Rps %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(Rps, "PCoA", "bray") , color = "ReadPerc", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC percentage",color = "ReadPerc") +
  scale_colour_viridis_c()

Rps %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(ordinate(Rps, "PCoA", "bray") , color = "ReadTot", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC total",color = "ReadTot") +
  scale_colour_viridis_c()


pcoa_bc2 = Rps %>% subset_samples(Sample_Unique != "2_57") %>% ordinate("PCoA", "bray") 

pcoa_bc

Rps %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(pcoa_bc2, color = "ReadTot", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard concentration",color = "ReadTot") +
  scale_colour_viridis_c()

pcoa_bc3 = Rps_mp %>% subset_samples(Sample_Unique != "2_57") %>% ordinate("PCoA", "bray") 


Rps %>% subset_samples(Sample_Unique != "2_57") %>% plot_ordination(pcoa_bc3, color = "ReadPerc", label = "Sample_Unique") + 
  geom_point(size = 3)  + labs(title = "PCoA BC concentration",color = "ReadPerc") +
  scale_colour_viridis_c()


unwt.unifrac <- plot_ordination(subset16S, 
                                ordu.unwt.uni, color="Farm2") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
unwt.unifrac
ps1.rel <- microbiome::transform(subset16S2, "compositional")
ordu.wt.uni <- ordinate(ps1.rel , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(ps1.rel, 
                              ordu.wt.uni, color="AB") 
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() + scale_color_brewer("AB", palette = "Set2")
print(wt.unifrac)



metadf <- data.frame(sample_data(ps1.rel))

unifrac.dist <- UniFrac(ps1.rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

tse2 = makeTreeSummarizedExperimentFromPhyloseq(subset16S)
tse2 <- relAbundanceCounts(tse2)

tse2 <- transformCounts(tse2, method = "relabundance")
tse2 <- runNMDS(tse2, FUN = vegan::vegdist, name = "BC", nmdsFUN = "monoMDS",
                    exprs_values = "relabundance",
                    keep_dist = TRUE)

plotReducedDim(tse2, "BC", colour_by = "Age")


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
adonis2(t(assay(tse, "relabundance")) ~ FarmRoundStable, data = colData(tse), permutations = 9999) 
adonis2(t(assay(tse, "relabundance")) ~ FlockSize, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Farm2, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ AgeParentStock, data = colData(tse), permutations = 9999)
adonis2(t(assay(tse, "relabundance")) ~ Age, data = colData(tse), permutations = 9999)

# variances: AB: 0.026, Cox: 0.102, Researcher: 0.06, FP : 0.067, LitterType: 0.061, FT :0.055, Gender: 0.007, 
# Stable: 0.167, FS: 0.1245, Farm 0.103, APS : 0.118, Age: 0.054
# Order: Stable>FS>APS>Farm>Cox>FP>LT>Researcher>FT>Age>AB>Gender

adonis2(t(assay(tse, "relabundance")) ~ FarmRoundStable * Age, data = colData(tse), permutations = 9999) 


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
ps1.rel <- microbiome::transform(subset16S, "compositional")
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
bray.dist <- distance(subset16S, "bray")


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


adonis2(dist(otu_table(subset16S), method='euclidean') ~ Age, 
        data=metadf)




permanova_age <- adonis2(unifrac.dist ~ Age, data = metadf)
permanova_AB <- adonis2(unifrac.dist ~ AB, data = metadf)
permanova_farm <- adonis2(unifrac.dist ~ Farm2, data = metadf)
permanova_cox <- adonis2(unifrac.dist ~ Cox, data = metadf)
permanova_researcher <- adonis2(unifrac.dist ~ Researcher, data = metadf)
permanova_LitterType <- adonis2(unifrac.dist ~ LitterType, data = metadf)
permanova_cox <- adonis2(unifrac.dist ~ Cox, data = metadf)


ps.disper <- betadisper(unifrac.dist, metadf$Age)
permutest(ps.disper, pairwise = TRUE)


source("../Results/Scripts/Steinberger_scripts/simper_pretty.r")
source("../Results/Scripts/Steinberger_scripts/R_krusk.r")

simper.pretty(otu_table(subset16S), metrics = sample_data(Rps), interesting = c("Age", "AB", "Farm2"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "16S")

simper.results = data.frame(read.csv("16s_clean_simper.csv"))


simper.results = data.frame(read.csv("Rps_clean_simper_16s.csv"))

kruskal.pretty(otu_table(subset16S), metrics = sample_data(subset16S), csv = simper.results, interesting = c('Age'), output_name =  'Age')

kruskal.pretty(otu_table(subset16S), metrics = sample_data(subset16S), csv = simper.results, interesting = c('AB'), output_name =  'AB')


class(sample_data(subset16S))

KW.results = data.frame(read.csv("Age_krusk_simper.csv"))

KW.results = KW.results[KW.results$fdr_krusk_p.val < 0.05,] # filter out non-significant results, based on fdr

KW.results = KW.results[with(KW.results, order(OTU)),]
head(KW.results)

abund = otu_table(Rps)/rowSums(otu_table(Rps))*100


boxplot(unlist(data.frame(abund["tet(O/32/O)_5_FP929050"])) ~ sample_data(Rps)$Age, ylab="% Relative abundance", main="OTU1")


for (otu in KW.results$OTU) {
  print(otu)
  
}

kruskal.test(unlist(data.frame(otu_table(Rps)["tet(O/32/O)_5_FP929050"]), use.names = FALSE) ~ sample_data(Rps)$Age)

kruskal.test(unlist(data.frame(otu_table(Rps)["tet(O/W/32/O)_1_EF065523"]), use.names = FALSE) ~ sample_data(Rps)$Age)


# declutter R environment by removing objects that no longer serve a purpose
rm(KW.results, dist, ord_meths, pcoa_bc, pcoa_jaccard, pcoa_jsd, pcoa_unifrac, pcoa_wunifrac, simper.results, pdataframe, metadf)

