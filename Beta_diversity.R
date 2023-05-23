
data_otu <- read.table("data_loue_16S_nonnorm.txt", header = TRUE)



filt.rar=data.frame(otu_table(ps0.rar))

dist_bc <- as.matrix(vegdist(filt.rar, method = "bray")) 

dist_bc[1:5, 1:5]

pcoa_bc = ordinate(subsetG, "PCoA", "bray") 


sample_data(subsetG)$Age = as.factor(sample_data(subsetG)$Age)
plot_ordination(subsetG, pcoa_bc, color = "Cluster") + 
  geom_point(size = 3) 



ps0.rar.filtered <- core(ps0.rar, detection = 10, prevalence = 0.05)
summarize_phyloseq(ps0.rar.filtered )
ordu.unwt.uni <- ordinate(ps0.rar, "PCoA", "unifrac", weighted=F)
barplot(ordu.unwt.uni$values$Eigenvalues[1:10])
unwt.unifrac <- plot_ordination(ps0.rar, 
                                ordu.unwt.uni, color="Farm2") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
unwt.unifrac

sample_data(ps0.rar)$Age = as.factor(sample_data(ps0.rar)$Age)

sample_data(ps1.rel)$Cluster = as.factor(sample_data(ps1.rel)$Cluster)

unwt.unifrac <- plot_ordination(ps0.rar, 
                                ordu.unwt.uni, color="Farm2") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac") + geom_point(size = 2)
unwt.unifrac <- unwt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
unwt.unifrac
ps1.rel <- microbiome::transform(subsetG, "compositional")
ordu.wt.uni <- ordinate(ps1.rel , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(ps1.rel, 
                              ordu.wt.uni, color="Farm2") 
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() + scale_color_brewer("Farm2", palette = "Set2")
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

permanova <- adonis(unifrac.dist ~ Age, data = metadf)

summary(permanova)

ps.disper <- betadisper(unifrac.dist, metadf$Age)
permutest(ps.disper, pairwise = TRUE)
