#library(microbiomeDataSets)
library(scater) # plotReducedDim
library(mia) # microbiome analysis package, making tse
library(vegan) # used to run simper
library(plyr) # for llply, to apply functions
library(nlme) # for usage of llply(), to apply functions over lists

# Used the following guide : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html

# Visualizing different kinds of ordination methods with BC distance matrix
#currently struggling to implement DPCoA for some reason, plot_ordination gives an error with DPCoA ordination
plot_ordination(Rps_mp,  ordinate(Rps, method = "DPCoA", distance="bray"), "samples", color="Age", shape="AB")

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
  scale_colour_brewer(type="qual", palette="Set1")


#`Looking at taxa spread as well
plot_ordination(Rps,  ordinate(Rps, method = "DCA", distance="bray"), type = "split", color="Age", shape="AB")
plot_ordination(Rps,  ordinate(Rps, method = "CCA", distance="bray"), type = "split", color="Age", shape="AB")
plot_ordination(Rps,  ordinate(Rps, method = "RDA", distance="bray"), type = "split", color="Age", shape="AB")
plot_ordination(Rps,  ordinate(Rps, method = "NMDS", distance="bray"), type = "split", color="Age", shape="AB")
plot_ordination(Rps,  ordinate(Rps, method = "MDS", distance="bray"), type = "split", color="Age", shape="AB")
plot_ordination(Rps,  ordinate(Rps, method = "PCoA", distance="bray"), type = "split", color="Age", shape="AB")
`
# PCoAs for different methods, with Age and Farm as colors, and AB as shape

# (fancy maken : + theme_classic() + scale_color_brewer("Farm2", palette = "Set2"))

plot_pcoa_ordination <- function(data, pcoa, var, title) {
  p <- plot_ordination(data, pcoa, color = var, shape = "Antibiotics used") +
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


#plot_ordination(Rps, pcoa_bc, type = "taxa", color = "AMR_class_primary") + 
#  geom_point(size = 3)  + labs(title = "PCoA primary AMR classes", color = "AMR_class_primary")


plot_ordination(Rps, pcoa_bc, type = "split", color = "AMR_class_primary", shape = "AB") + 
  geom_point(size = 3) + labs(title = "PCoA Bray Curtis Farms",color = "Farms", shape = "Antibiotics used")

pcoa_unifrac = ordinate(Rps, "PCoA", "unifrac") 


# plot to look at concentration with a red/green gradient

plot_ordination(Rps, pcoa_jaccard, color = "Conc...ng..µl.", shape = "AB", label = "firm_id") + 
  geom_point(size = 3)  + labs(title = "PCoA Jaccard concentration",color = "Conc...ng..µl.", shape = "Antibiotics used") +
  scale_colour_gradient(low = "red", high = "green")

View(Rps@sam_data)

sample_data(Rps)['firm_id'] <- row.names(sample_data(Rps)) 


plot_scree(pcoa_jsd) #scree plots can be made for any of the PCoAs



# this is the same as above but looks different
plot_ordination(Rps, 
                pcoa_wunifrac, color="Farm2") + ggtitle("Unweighted UniFrac") + geom_point(size = 2) + 
  theme_classic() + scale_color_brewer("Farm2", palette = "Set2")


ps1.rel <- microbiome::transform(Rps, "compositional")

metadf <- data.frame(sample_data(ps1.rel))

# NMDS

tse2 = makeTreeSummarizedExperimentFromPhyloseq(Rps)
tse2 %<>% relAbundanceCounts()

tse2 %<>%  transformCounts( method = "relabundance")
tse2 %<>% runNMDS(FUN = vegan::vegdist, name = "BC", nmdsFUN = "monoMDS",
                    exprs_values = "relabundance",
                    keep_dist = TRUE)

tse2 %>% plotReducedDim("BC", colour_by = "Age") 

# Significance

psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


unifrac.dist <- UniFrac(ps1.rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

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


simper.results = data.frame(read.csv("Rps_clean_simper.csv"))

kruskal.pretty(otu_table(Rps), metrics = sample_data(Rps), csv = simper.results, interesting = c('Age'), output_name =  'Age')


kruskal.pretty(otu_table(Rps), metrics = sample_data(Rps), csv = simper.results, interesting = c('AB'), output_name =  'AB')


class(sample_data(Rps))

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

