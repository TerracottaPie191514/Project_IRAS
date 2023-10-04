#library(microbiomeDataSets)
library(scater) # plotReducedDim
library(mia) # microbiome analysis package, making tse
library(vegan) # used to run simper
library(plyr) # for llply, to apply functions
library(nlme) # 

# Used the following guide : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

plot_ordination(subsetMG,  ordinate(subsetMG, method = "DPCoA", distance="bray"), "samples", color="Age", shape="AB")

estimate_richness(subsetMG)

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



filt.rar=data.frame(otu_table(subsetMG))

dist_bc <- as.matrix(vegdist(filt.rar, method = "bray")) 

dist_bc[1:5, 1:5]
#sample_data(subsetMG)$Age = as.factor(sample_data(subsetMG)$Age)


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


plot_scree(pcoa_bc) #scree plots can be made for any of the PCoAs


unwt.unifrac <- plot_ordination(subsetMG, 
                                ordu.unwt.uni, color="Farm2") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
unwt.unifrac
ps1.rel <- microbiome::transform(subsetMG2, "compositional")
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

tse2 = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse2 <- relAbundanceCounts(tse2)

tse2 <- transformCounts(tse2, method = "relabundance")
tse2 <- runNMDS(tse2, FUN = vegan::vegdist, name = "BC", nmdsFUN = "monoMDS",
                exprs_values = "relabundance",
                keep_dist = TRUE)

plotReducedDim(tse2, "BC", colour_by = "Age")


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

simper.pretty(otu_table(subsetMG), metrics = sample_data(Rps), interesting = c("Age", "AB", "Farm2"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "16S")

simper.results = data.frame(read.csv("16s_clean_simper.csv"))


simper.results = data.frame(read.csv("Rps_clean_simper_16s.csv"))

kruskal.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), csv = simper.results, interesting = c('Age'), output_name =  'Age')

kruskal.pretty(otu_table(subsetMG), metrics = sample_data(subsetMG), csv = simper.results, interesting = c('AB'), output_name =  'AB')


class(sample_data(subsetMG))

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

