
filt.rar=data.frame(otu_table(subsetG))

dist_bc <- as.matrix(vegdist(filt.rar, method = "bray")) 

dist_bc[1:5, 1:5]

pcoa_bc = ordinate(subsetG, "PCoA", "bray") 


sample_data(subsetG)$Age = as.factor(sample_data(subsetG)$Age)
plot_ordination(subsetG, pcoa_bc, color = "Cluster") + 
  geom_point(size = 3) 



subsetG.filtered <- core(subsetG, detection = 10, prevalence = 0.05)
summarize_phyloseq(subsetG.filtered )
ordu.unwt.uni <- ordinate(subsetG, "PCoA", "unifrac", weighted=F)
barplot(ordu.unwt.uni$values$Eigenvalues[1:10])
plot_scree(ordu.unwt.uni)

unwt.unifrac <- plot_ordination(subsetG, 
                                ordu.unwt.uni, color="Farm2") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
unwt.unifrac

subsetG@sam_data$Age==35
subsetG2=subset_samples(subsetG, Age == "35")


rm(subsetG2)

sample_variables(subsetG)

sample_data(subsetG)$Age = as.factor(sample_data(subsetG)$Age)

sample_data(ps1.rel)$Cluster = as.factor(sample_data(ps1.rel)$Cluster)

unwt.unifrac <- plot_ordination(subsetG, 
                                ordu.unwt.uni, color="Farm2") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
unwt.unifrac
ps1.rel <- microbiome::transform(subsetG2, "compositional")
ordu.wt.uni <- ordinate(ps1.rel , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(ps1.rel, 
                              ordu.wt.uni, color="AB") 
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() + scale_color_brewer("AB", palette = "Set2")
print(wt.unifrac)


p <- plot_landscape(ps1.rel, 
                    "NMDS", 
                    "bray", 
                    col = "AB") +
  labs(title = paste("NMDS / Bray-Curtis"))   

p <- p + scale_color_brewer(palette = "Dark2")+ scale_fill_gradient(low = "#e0ecf4", high = "#6e016b") 
p



metadf <- data.frame(sample_data(ps1.rel))

unifrac.dist <- UniFrac(ps1.rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

tse2 = makeTreeSummarizedExperimentFromPhyloseq(subsetG)
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



permanova_age <- adonis(t(unifrac.dist) ~ Age, data = metadf, permutations = 9999)

permanova_age <- adonis(t(assay(tse2, "relabundance")) ~ Age, data = metadf, permutations = 9999)

permanova_age <- adonis(unifrac.dist ~ Age, data = metadf, permutations = 9999)

permanova_age$coef.sites

permanova_AB
permanova_farm
permanova_cox
permanova_researcher
permanova_LitterType

library(microbiomeDataSets)

coef <- coefficients(permanova_age)["Age1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

ggplot(data.frame(x = top.coef,
                  y = factor(names(top.coef),
                             unique(names(top.coef)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="",y="",title="Top Taxa") +
  theme_bw()

permanova_LitterType
permanova_age$aov.tab["Age","Pr(>F)"]

coef

permanova_age$coefficients


ps.disper <- betadisper(unifrac.dist, metadf$Age)
permutest(ps.disper, pairwise = TRUE)

assay(tse2)
