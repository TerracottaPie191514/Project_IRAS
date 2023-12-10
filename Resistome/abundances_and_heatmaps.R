#### Load packages
library(ggpubr) # Publication quality figures, based on ggplot2.
library(RColorBrewer) # Color options.
library(microbiome) # For data analysis and visualisation.
library(microbiomeutilities) # Some utility tools for microbiome package.


# used the following guides:, https://microbiome.github.io/OMA/viz-chapter.html, https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/composition-plots.html#heatmaps, https://david-barnett.github.io/microViz/articles/web-only/compositions.html#sorting-the-barplot

# absolute abundances
plot_bar(Rps, fill="AMR_class_primary", title = "Absolute abundances per sample (FPKM)")
plot_bar(Rps_mp, fill="AMR_class_primary", title = "Absolute abundances per sample (FPKM) - Metaphlan")
plot_bar(Rps_tpm, fill="AMR_class_primary", title = "Absolute abundances per sample (TPM)")


# on average, for the three outliers, abundance is 
mean(sample_sums(Rps)[c("10_1","10_2","10_3")]) # 31928.03
# on average, the samples without the three outliers, abundance is
mean(sample_sums(Rps)[!sample_names(Rps) %in% c("10_1","10_2","10_3")]) # 6183.18
# so there are 5,1636x as much abundance in these samples

# MP:
# on average, for the three outliers, abundance is 
mean(sample_sums(Rps_mp)[c("10_1","10_2","10_3")]) # 78089.08
# on average, the samples without the three outliers, abundance is
mean(sample_sums(Rps_mp)[!sample_names(Rps_mp) %in% c("10_1","10_2","10_3")]) # 14828.26
# 5,266x as much abundance

# Amount of different AMR classes present.
sort(table(tax_table(Rps)[, "AMR_class_primary"]))
sort(table(tax_table(Rps)[, "ARGCluster90"]))

# for plotting abundances of specific stables
Rps %>% ps_filter(FarmRoundStable == c("Farm2R1S2")) %>% plot_bar(fill="AMR_class_primary")

# visualisation on AB at AMR primary class level, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- phyloseq::tax_glom(Rps_copy, "Phylum")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Phylum"]

ps_prim <- subset_taxa(Rps,Phylum %in% c("p__Actinomycetota", "p__Bacillota", "p__Bacteroidota", "p__Campylobacterota", "p__Pseudomonadota", "p__Verrucomicrobiota"))  %>% aggregate_top_taxa2("Phylum", top = 6) %>% phyloseq::tax_glom("Phylum") 


psmelt(ps_prim %>%  subset_taxa(!Phylum %in% c("Not determined", "Fosfomycin", "Quinolone", "Streptogramin"))) %>% # AB
  ggplot(data = ., aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  facet_wrap(~ Phylum, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "Primary AMR class") 

psmelt(ps_prim) %>% # Age
  ggplot(data = ., aes(x = Age, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  facet_wrap(~ Phylum, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "Primary AMR class") 

psmelt(ps_prim) %>% # Farm
  ggplot(data = ., aes(x = Farm2, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  facet_wrap(~ Phylum, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "Primary AMR class") 

psmelt(ps_prim) %>%  # Stable
  ggplot(data = ., aes(x = Stable, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  facet_wrap(~ Phylum, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "Primary AMR class") + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

psmelt(ps_prim) %>% # Agent
  ggplot(data = ., aes(x = Cox, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  facet_wrap(~ Phylum, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "Primary AMR class") 

# visualisation on AB at ARGclust90 level, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- Rps_copy %>% aggregate_top_taxa2("Order", top = 11) %>% phyloseq::tax_glom("Order")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Order"]

psmelt(ps_prim) %>% # AB
  ggplot(data = ., aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Order), height = 0, width = .2) +
  facet_wrap(~ Order, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "ARGCluster90") 

psmelt(ps_prim) %>% # Age
  ggplot(data = ., aes(x = Age, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Order), height = 0, width = .2) +
  facet_wrap(~ Order, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "ARGCluster90") 

psmelt(ps_prim) %>% # Farm
  ggplot(data = ., aes(x = Farm2, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Order), height = 0, width = .2) +
  facet_wrap(~ Order, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "ARGCluster90") 

psmelt(ps_prim) %>% # Stable
  ggplot(data = ., aes(x = Stables, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Order), height = 0, width = .2) +
  facet_wrap(~ Order, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "ARGCluster90") 

psmelt(ps_prim) %>% # Agent
  ggplot(data = ., aes(x = Cox, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Order), height = 0, width = .2) +
  facet_wrap(~ Order, scales = "free") +
  labs(x = "", y = "Abundance\n", color = "ARGCluster90") 


# relative abundances of primary AMR, split by AB (very ugly)
ps_rel_abund = transform_sample_counts(Rps, function(x){x / sum(x)})
plot_bar(ps_rel_abund, fill = "AMR_class_primary") +
  geom_bar(aes(color = AMR_class_primary, fill = AMR_class_primary), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ AB, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# visualisation on primary AMR classes, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- phyloseq::tax_glom(Rps, "AMR_class_primary")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "AMR_class_primary"]

psmelt(ps_prim) %>%
  ggplot(data = ., aes(x = FarmRoundStable, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

# Check the amount of uniqe ARGs in samples which have and have not been treated with antibiotics
Rps %>% ps_filter(AB == "no") %>% get_taxa_unique("ARGCluster90") # 163 different genes for non AB treated
Rps %>% ps_filter(AB == "yes") %>% get_taxa_unique("ARGCluster90") # 81 different genes for AB treated
Rps %>% get_taxa_unique("ARGCluster90") # 188 different genes in total, meaning 25 do not overlap


Rps_tpm %>% ps_filter(AB == "no") %>% get_taxa_unique("ARGCluster90") # 168 different genes for non AB treated
Rps_tpm %>% ps_filter(AB == "yes") %>% get_taxa_unique("ARGCluster90") # 97 different genes for AB treated
Rps_tpm %>% get_taxa_unique("ARGCluster90") # 186 different genes in total, meaning 18 do not overlap

# Plots of relative abundances, fixing some genes that are clustered in the data twice, showing top 12 taxa and others are clustered

# TPM and FPKM in the same plot, comparing non AB
dataset1 =  ps_filter(Rps)
dataset2 =  ps_filter(Rps_tpm)

dataset1 %<>% ps_mutate(dataset = "FPKM")
dataset2 %<>% ps_mutate(dataset = "TPM")

sample_names(dataset1) <- paste(sample_names(dataset1), "FPKM", sep="_")
sample_names(dataset2) <- paste(sample_names(dataset2), "TPM", sep="_")

combined <- phyloseq::merge_phyloseq(
  dataset1 %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>% tax_agg("ARGCluster90") %>% ps_get(),
  dataset2 %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>% tax_agg("ARGCluster90") %>% ps_get()
)

combined %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>%
  comp_barplot("ARGCluster90", facet_by = "dataset", n_taxa = 12, palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
               other_name = "Other cluster", merge_other = F, sample_order = "asis") +
  coord_flip() + ggtitle("FPKM vs TPM relative abundance of ARG clusters at 90% sequence identity")

# uitprobeersel

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>% ps_arrange(Age) %>% 
  comp_barplot(tax_level = "ARGCluster90", n_taxa = 12, sample_order = "asis", palette = colorRampPalette(brewer.pal(8,"Accent"))(13)) +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age", ., "days")),
    scales = "fixed") +
  coord_flip()


# MAKE PLOT OF K2 VS MP

dataset1 =  ps_filter(Rps)
dataset2 =  ps_filter(Rps_mp)

dataset1 %<>% ps_mutate(dataset = "k2")
dataset2 %<>% ps_mutate(dataset = "MP")

sample_names(dataset1) <- paste(sample_names(dataset1), "k2", sep="_")
sample_names(dataset2) <- paste(sample_names(dataset2), "MP", sep="_")

combined <- phyloseq::merge_phyloseq(
  dataset1 %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>% tax_agg("ARGCluster90") %>% ps_get(),
  dataset2 %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>% tax_agg("ARGCluster90") %>% ps_get()
)

# primary class
combined %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>%
  comp_barplot("AMR_class_primary", facet_by = "dataset", n_taxa = 14, palette = colorRampPalette(brewer.pal(8,"Accent"))(14),
               sample_order = "asis") +
  coord_flip() + ggtitle("Kraken2 vs Metaphlan relative abundance of primary AMR classes")

#Argclust90
combined %>% tax_fix(unknowns = c("blaOXA-493_clust", "cfr(B)", "dfrA16_clust", "dfrA7_dfrA17", "lnu(B)")) %>%
  comp_barplot("ARGCluster90", facet_by = "dataset", n_taxa = 12, palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
               other_name = "Other cluster", merge_other = F, sample_order = "asis") +
  coord_flip() + ggtitle("Kraken2 vs Metaphlan relative abundance of ARG clusters at 90% sequence identity")

# relative abundance of AMR primary class shown by age and farm

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Farm2") +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age = ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Farm", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of primary antimicrobial class by farm and age (FPKM)")

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Farm2") +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age = ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Farm", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of primary antimicrobial class by farm and age (TPM)")

# uitprobeersel (lelijk)
Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Farm2") +
  facet_grid(
    cols = vars(Age), rows = vars(AB),
    labeller = labeller(.cols = as_labeller(~ paste("Age = ", .))),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Farm", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of ARGs by farm and age (FPKM)")

# Relative abundance for both stable and antibiotics used

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "FarmRoundStable") +
  facet_wrap(
    facets = vars(AB), labeller = as_labeller(~ paste("Antiobotics used: ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Stable", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of primary antimicrobial class by stable and antibiotics used (FPKM)")

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "AMR_class_primary", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "FarmRoundStable") +
  facet_wrap(
    facets = vars(AB), labeller = as_labeller(~ paste("Antiobotics used: ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Stable", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of primary antimicrobial class by stable and antibiotics used (TPM)")


# Same plots but with ARGclust90

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "ARGCluster90", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Farm2") +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age = ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Farm", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of ARG by farm and age (FPKM)")

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "ARGCluster90", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(9),
    x = "Farm2") +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age = ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Farm", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of ARG by farm and age (TPM)")


Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "ARGCluster90", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
    x = "FarmRoundStable",
    n_taxa = 12, other_name = "Other ARG", merge_other = F) +
  facet_wrap(
    facets = vars(AB), labeller = as_labeller(~ paste("Antiobotics used: ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Stable", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of ARG by stable and antibiotics used (FPKM)")


Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>%
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "ARGCluster90", bar_width = 0.7, sample_order = "asis", 
    palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
    x = "FarmRoundStable",
    n_taxa = 12, other_name = "Other ARG", merge_other = F) +
  facet_wrap(
    facets = vars(AB), labeller = as_labeller(~ paste("Antiobotics used: ", .)),
    scales = "fixed"
  ) +
  coord_flip() +
  labs(x = "Stable", y = "Relative abundance") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) + 
  theme_bw() + 
  theme(panel.spacing.x = unit(6, "mm")) +
  ggtitle("Relative abundance of ARG by stable and antibiotics used (TPM)")


# 2_57 is a big outlier, for some reason barely has tet(W) and other common ARG clusters, and therefore a lot more other abundant ARG clusterss

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)")) %>% ps_filter(SampleID == "7.F1S2.21.09") %>% ps_calc_dominant(rank = "ARGCluster90") %>% 
  comp_barplot(tax_level = "ARGCluster90", n_taxa = 12, palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
               other_name = "Other ARG", merge_other = F) +
  coord_flip()

# rel abundance on primary AMR class level

tse = makeTreeSummarizedExperimentFromPhyloseq(Rps_copy)

tse <- transformCounts(tse, method = "relabundance")

tse_AMRClass <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)
tse_AMRClass <- transformCounts(tse,
                              assay.type = "counts",
                              method = "relabundance")
top_taxa <- getTopTaxa(tse_AMRClass,
                       top = 10,
                       assay.type = "relabundance")
AMRClass_renames <- lapply(rowData(tse_AMRClass)$Phylum,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_AMRClass)$Phylum <- as.character(AMRClass_renames)

# rel abundance figures, can order by specific taxa
miaViz::plotAbundance(tse_AMRClass,
                      assay.type = "relabundance",
                      rank = "Phylum",
                      order_rank_by = "abund")


tse_AMRClass$Farm2 = as.factor(tse_AMRClass$Farm2)
tse_AMRClass$AB = as.factor(tse_AMRClass$AB)

# Add AB plot on top

plots <- miaViz::plotAbundance(tse_AMRClass,
                               assay.type = "relabundance",
                               rank = "Phylum",
                               order_rank_by = "abund",
                               #                       order_sample_by = "o__Clostridiales",
                               order_sample_by = "AB",
                               features = "AB")

plots[[1]] <- plots[[1]] +
  theme(legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))
plots[[2]] <- plots[[2]] +
  theme(legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.direction = "vertical")

legend <- wrap_plots(as_ggplot(get_legend(plots[[1]])), as_ggplot(get_legend(plots[[2]])), ncol = 1) 
plots[[1]] <- plots[[1]] + theme(legend.position = "none")
plots[[2]] <- plots[[2]] + theme(legend.position = "none", axis.title.x=element_blank()) 

plot <- wrap_plots(plots[[2]], plots[[1]], ncol = 1, heights = c(2, 10))
wrap_plots(plot, legend, nrow = 1, widths = c(2, 1))


# heatmaps on phylum level

tse_AMRClass <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

tse_AMRClass <- transformCounts(tse_AMRClass, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)
tse_AMRClass <- transformCounts(tse_AMRClass, assay.type = "clr",
                              MARGIN = "features", 
                              method = "z", name = "clr_z")


top_taxa <- getTopTaxa(tse_AMRClass, top = 11)
tse_AMRClass <- tse_AMRClass[top_taxa, ]

# Phylum AB heatmap
tse_AMRClass@metadata$anno_colors$AB = c(yes = "darkred",no ="darkblue")


sechm(tse_AMRClass,
      features = rownames(tse_AMRClass),
      assayName = "clr",
      do.scale = TRUE,
      top_annotation = c("AB"), 
      gaps_at = "AB",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)



# heatmap with AB and stable

sechm(tse_AMRClass,
      features = rownames(tse_AMRClass),
      assayName = "clr",
      do.scale = TRUE,
      top_annotation = c("AB"), 
      gaps_at = "Stables",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)

tse_AMRClass@metadata$anno_colors$Stable = (brewer.pal(n=10, name = "Set3"))



sechm(tse_AMRClass,
      features = rownames(tse_AMRClass),
      assayName = "clr",
      do.scale = TRUE,
      top_annotation = c("AB"), 
      gaps_at = "AB",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)

# Phylum heatmap
mat <- assay(tse_AMRClass, "clr_z")

pheatmap(mat)

# Phylum heatmap hierarchal clustering with AB

# Clustering both samples and features hierarchically 

taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

taxa_tree # based on this three, we'll create two clusters

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

taxa_clusters <- cutree(tree = taxa_hclust, k = 2) # 2 clusters based on tree figure

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters

rowData(tse_AMRClass)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_AMRClass))), ]

# Prints taxa and their clusters
rowData(tse_AMRClass)$clusters


sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

# to view the tree, run
sample_tree

# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 2)) # 2 clusters based on methods in Clustering.R script

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 

# Order data based on 
tse_AMRClass <- tse_AMRClass[ , rownames(sample_data)]

# Add sample type data
sample_data$sample_types <- colData(tse_AMRClass)$AB

sample_data


breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
#colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1) replaced with viridis pallette

pheatmap(mat, annotation_row = taxa_clusters,
         annotation_col = sample_data,
         breaks = breaks,
         color = colorRampPalette(viridis(256))(length(breaks)-1))


# heatmaps on ARG level

tse = makeTreeSummarizedExperimentFromPhyloseq(Rps)
tse <- transformCounts(tse, method = "relabundance")
tse <- transformCounts(tse, MARGIN = "samples", method = "clr", assay.type = "counts", pseudocount=1)
tse <- transformCounts(tse, assay.type = "clr",
                       MARGIN = "features", 
                       method = "z", name = "clr_z")
top_taxa <- getTopTaxa(tse, top = 20)
tse <- tse[top_taxa, ]

# ARG heatmap AB
tse@metadata$anno_colors$AB = c(yes = "darkred",no ="darkblue")

sechm(tse, 
      features = rownames(tse), 
      assayName = "clr", 
      do.scale = TRUE, 
      top_annotation = c("AB"), 
      gaps_at = "AB",
      hmcols = viridis(256),
      cluster_cols = TRUE, cluster_rows = TRUE)

# ARG heatmap
mat <- assay(tse, "clr_z")

pheatmap(mat)

# ARG heatmap hierarchal clustering with AB

# Clustering both samples and features hierarchically 

taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

taxa_tree # based on this three, we'll create two clusters

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

taxa_clusters <- cutree(tree = taxa_hclust, k = 2) # 2 clusters based on methods in Clustering.R script

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters

rowData(tse_AMRClass)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_AMRClass))), ]

# Prints taxa and their clusters
rowData(tse_AMRClass)$clusters


sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

# to view the tree, run
sample_tree

# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 2))

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 

# Order data based on 
tse_AMRClass <- tse_AMRClass[ , rownames(sample_data)]

# Add sample type data
sample_data$sample_types <- colData(tse_AMRClass)$AB

sample_data


breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )

pheatmap(mat, annotation_row = taxa_clusters,
         annotation_col = sample_data,
         breaks = breaks,
         color = colorRampPalette(viridis(256))(length(breaks)-1))


# declutter environment
rm(ps_prim, tse, tse_phylum, taxmat, taxic, taxa_tree, taxa_hclust, taxa_clusters, subset16S.rel, sample_tree,
   sample_hclust, sample_data, ps1.com, ps1.com.fam, plots, plot, phylum_renamed, mat, legend, guide_italics,
   breaks, new.tax, sample_clusters, samples_ordered, taxa_ordered, top_taxa)
