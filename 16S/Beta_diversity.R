library(microbiomeDataSets)
library(scater)
library(mia)

# Used the following guide : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

filt.rar=data.frame(otu_table(subset_mg))

dist_bc <- as.matrix(vegdist(filt.rar, method = "bray")) 

dist_bc[1:5, 1:5]
#sample_data(subset_mg)$Age = as.factor(sample_data(subset_mg)$Age)


# PCoAs for different methods (fancy maken : + theme_classic() + scale_color_brewer("Farm2", palette = "Set2"))

pcoa_bc = ordinate(subset_mg, "PCoA", "bray") 

plot_ordination(subset_mg, pcoa_bc, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA Bray Curtis", color = "Age", shape = "Antibiotics used")

plot_ordination(subset_mg, pcoa_bc, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA Bray Curtis",color = "Farms", shape = "Antibiotics used")



pcoa_unifrac = ordinate(subset_mg, "PCoA", "unifrac") 


plot_ordination(subset_mg, pcoa_unifrac, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA UniFrac",color = "Age", shape = "Antibiotics used")

plot_ordination(subset_mg, pcoa_unifrac, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA UniFrac",color = "Farms", shape = "Antibiotics used")


pcoa_wunifrac = ordinate(subset_mg, "PCoA", "wunifrac") 


plot_ordination(subset_mg, pcoa_wunifrac, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA weighter UniFrac",color = "Age", shape = "Antibiotics used")

plot_ordination(subset_mg, pcoa_wunifrac, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA weighted Unifrac",color = "Farms", shape = "Antibiotics used")

pcoa_jsd = ordinate(subset_mg, "PCoA", "jsd") 


plot_ordination(subset_mg, pcoa_jsd, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA Jensen-Shannon Divergence",color = "Age", shape = "Antibiotics used")

plot_ordination(subset_mg, pcoa_jsd, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA Jensen-Shannon Divergence",color = "Farms", shape = "Antibiotics used")

pcoa_jaccard = ordinate(subset_mg, "PCoA", "jaccard", binary=TRUE) 

plot_ordination(subset_mg, pcoa_jaccard, color = "Age", shape = "AB") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard",color = "Age", shape = "Antibiotics used")

plot_ordination(subset_mg, pcoa_jaccard, color = "Farm2", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA Jaccard",color = "Farms", shape = "Antibiotics used")



#subset_mg.filtered <- core(subset_mg, detection = 10, prevalence = 0.05)

plot_scree(pcoa_bc) #scree plots can be made for any of the PCoAs


# plots for AB where age = 35 (depcrated)
subset_mg@sam_data$Age==35
subset_mg2=subset_samples(subset_mg, Age == "35")


unwt.unifrac <- plot_ordination(subset_mg, 
                                ordu.unwt.uni, color="Farm2") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
unwt.unifrac
ps1.rel <- microbiome::transform(subset_mg2, "compositional")
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

tse2 = makeTreeSummarizedExperimentFromPhyloseq(subset_mg)
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
