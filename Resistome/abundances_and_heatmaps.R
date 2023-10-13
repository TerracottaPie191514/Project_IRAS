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

# Amount of different AMR classes present.
sort(table(tax_table(Rps)[, "AMR_class_primary"]))
sort(table(tax_table(Rps)[, "ARGCluster90"]))

# for plotting abundances of specific stables
Rps %>% ps_filter(FarmRoundStable == c("Farm2R1S2")) %>% plot_bar(fill="AMR_class_primary")


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
               other_name = "Other ARG", merge_other = F, sample_order = "asis") +
  coord_flip() + ggtitle("FPKM vs TPM relative abundance of ARGs")

# uitprobeersel

Rps_tpm %>% tax_fix(unknowns = c("lnu(B)","cfr(B)","blaOXA-493_clust", "dfrA16_clust", "dfrA7_dfrA17")) %>% ps_arrange(Age) %>% 
  comp_barplot(tax_level = "ARGCluster90", n_taxa = 12, sample_order = "asis", palette = colorRampPalette(brewer.pal(8,"Accent"))(13)) +
  facet_wrap(
    facets = vars(Age), labeller = as_labeller(~ paste("Age", ., "days")),
    scales = "fixed") +
  coord_flip()

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


# Same plots but with ARG

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


# 2_57 is a big outlier, for some reason barely has tet(W) and other common ARGs, and therefore a lot more other abundant ARGs

Rps %>% tax_fix(unknowns = c("lnu(B)","cfr(B)")) %>% ps_filter(SampleID == "7.F1S2.21.09") %>% ps_calc_dominant(rank = "ARGCluster90") %>% 
  comp_barplot(tax_level = "ARGCluster90", n_taxa = 12, palette = colorRampPalette(brewer.pal(8,"Accent"))(13),
               other_name = "Other ARG", merge_other = F) +
  coord_flip()


# declutter environment
rm(ps_prim, ps_rel_abund, dataset1, dataset2, combined)
