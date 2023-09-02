library(microbiomeDataSets)
library(scater)
library(mia)
library(vegan)

# Used the following guide : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

# Visualizing different kinds of ordination methods

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
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
  scale_colour_brewer(type="qual", palette="Set1")


# PCoAs for different methods, with Age and Farm as colors, and AB as shape

# (fancy maken : + theme_classic() + scale_color_brewer("Farm2", palette = "Set2"))

pcoa_bc = ordinate(Rps, "PCoA", "bray") 

plot_ordination(Rps, pcoa_bc, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA Bray Curtis Age", color = "Age", shape = "Antibiotics used")

plot_ordination(Rps, pcoa_bc, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA Bray Curtis Farms",color = "Farms", shape = "Antibiotics used")

#plot_ordination(Rps, pcoa_bc, type = "taxa", color = "AMR_class_primary") + 
#  geom_point(size = 3)  + labs(title = "PCoA primary AMR classes", color = "AMR_class_primary")


plot_ordination(Rps, pcoa_bc, type = "split", color = "AMR_class_primary", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA Bray Curtis Farms",color = "Farms", shape = "Antibiotics used")

pcoa_unifrac = ordinate(Rps, "PCoA", "unifrac") 


plot_ordination(Rps, pcoa_unifrac, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA UniFrac Age",color = "Age", shape = "Antibiotics used")

plot_ordination(Rps, pcoa_unifrac, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA UniFrac Farms",color = "Farms", shape = "Antibiotics used")


pcoa_wunifrac = ordinate(Rps, "PCoA", "wunifrac") 


plot_ordination(Rps, pcoa_wunifrac, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA weighter UniFrac Age",color = "Age", shape = "Antibiotics used")

plot_ordination(Rps, pcoa_wunifrac, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA weighted Unifrac Farms",color = "Farms", shape = "Antibiotics used")

pcoa_jsd = ordinate(Rps, "PCoA", "jsd") 


plot_ordination(Rps, pcoa_jsd, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA Jensen-Shannon Divergence Age",color = "Age", shape = "Antibiotics used")

plot_ordination(Rps, pcoa_jsd, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA Jensen-Shannon Divergence Farms",color = "Farms", shape = "Antibiotics used")

pcoa_jaccard = ordinate(Rps, "PCoA", "jaccard", binary=TRUE) 

plot_ordination(Rps, pcoa_jaccard, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard Age",color = "Age", shape = "Antibiotics used")

plot_ordination(Rps, pcoa_jaccard, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA Jaccard Farms",color = "Farms", shape = "Antibiotics used")

# plot to look at concentration with a red/green gradient

plot_ordination(Rps, pcoa_jaccard, color = "Conc...ng..µl.", shape = "AB", label = "FarmRoundStable") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard concentration",color = "Conc...ng..µl.", shape = "Antibiotics used") +
  scale_colour_gradient(low = "red", high = "green")



#Rps.filtered <- core(Rps, detection = 10, prevalence = 0.05)

plot_scree(pcoa_jsd) #scree plots can be made for any of the PCoAs



# plots for AB where age = 35 (deprecated)
Rps@sam_data$Age==35
Rps2=subset_samples(Rps, Age == "35")


unwt.unifrac <- plot_ordination(Rps, 
                                ordu.unwt.uni, color="Farm2") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
unwt.unifrac
ps1.rel <- microbiome::transform(Rps2, "compositional")
ordu.wt.uni <- ordinate(ps1.rel , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(ps1.rel, 
                              ordu.wt.uni, color="AB") 
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() + scale_color_brewer("AB", palette = "Set2")
print(wt.unifrac)


ps1.rel <- microbiome::transform(Rps, "compositional")

metadf <- data.frame(sample_data(ps1.rel))


# waar komt dit vandaan? en het is weighted!?
unifrac.dist <- UniFrac(ps1.rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

tse2 = makeTreeSummarizedExperimentFromPhyloseq(Rps)
tse2 <- relAbundanceCounts(tse2)

tse2 <- transformCounts(tse2, method = "relabundance")
tse2 <- runNMDS(tse2, FUN = vegan::vegdist, name = "BC", nmdsFUN = "monoMDS",
                    exprs_values = "relabundance",
                    keep_dist = TRUE)

plotReducedDim(tse2, "BC", colour_by = "Age")

# Significance

psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Rps_veg = psotu2veg(Rps)
dist_bray = vegdist(Rps_veg, "bray")
betadisper(unifrac.dist,"Age")

adonis2(dist_bray ~ Age, data = metadf)


permanova_age <- adonis2(unifrac.dist ~ Age, data = metadf)
permanova_AB <- adonis2(unifrac.dist ~ AB, data = metadf)
permanova_farm <- adonis2(unifrac.dist ~ Farm2, data = metadf)
permanova_cox <- adonis2(unifrac.dist ~ Cox, data = metadf)
permanova_researcher <- adonis2(unifrac.dist ~ Researcher, data = metadf)
permanova_LitterType <- adonis2(unifrac.dist ~ LitterType, data = metadf)


ps.disper <- betadisper(unifrac.dist, metadf$Age)
permutest(ps.disper, pairwise = TRUE)


# Simper analyses to see which species are most impactful to BC dissimilarity between groups

# Transpose OTU table, and use farm as group structure
#simp_farm = simper(t(otu_table(Rps)), sample_data(Rps)$Farm2, permutations = 999)

# Test for significance (OTU abundance will not be normally distributed so we will use kruskal wallis tests)

source("../Results/Scripts/Steinberger_scripts/simper_pretty.r")
source("../Results/Scripts/Steinberger_scripts/R_krusk.r")

simper.pretty(otu_table(Rps), metrics = sample_data(Rps), interesting = c("Age", "AB", "Farm2"), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, output_name= "Rps")

simper.results = data.frame(read.csv("Age_clean_simper.csv"))

kruskal.pretty(otu_table(Rps), metrics = sample_data(Rps), csv = simper.results, interesting = c('Age'), output_name =  'Age')


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

kruskal.test(unlist(data.frame(otu_table(test)["tet(O/W/32/O)_1_EF065523"]), use.names = FALSE) ~ sample_data(Rps)$Age)
