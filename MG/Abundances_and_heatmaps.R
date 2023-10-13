#### Load packages
library(ggpubr) # Publication quality figures, based on ggplot2.
library(RColorBrewer) # Color options.
library(microbiome) # For data analysis and visualisation.
library(microbiomeutilities) # Some utility tools for microbiome package.


# used the following guides:, https://microbiome.github.io/OMA/viz-chapter.html, https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/composition-plots.html#heatmaps, https://david-barnett.github.io/microViz/articles/web-only/compositions.html#sorting-the-barplot


# absolute abundances, since there are a lot of phyla (43), we will only include the top 5 phyla 

subsetMG %>% aggregate_top_taxa2("Phylum", top = 5) %>% plot_bar(fill="Phylum", title = "Absolute abundances per sample")

# for plotting abundances of specific stables
subsetMG %>% aggregate_top_taxa2("Phylum", top = 5) %>% ps_filter(FarmRoundStable == c("Farm2R1S1")) %>% plot_bar(fill="Phylum")


# visualisation on Phyla, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- subsetMG %>% aggregate_top_taxa2("Phylum", top = 11) %>% phyloseq::tax_glom("Phylum")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Phylum"]

psmelt(ps_prim) %>% 
  ggplot(aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 

# visualisation on Genera, more data for samples which have not been treated with AB, but also many more samples in this group (18 vs 102)

ps_prim <- subsetMG %>% aggregate_top_taxa2("Genus", top = 22) %>% phyloseq::tax_glom("Genus")
taxa_names(ps_prim) <- phyloseq::tax_table(ps_prim)[, "Genus"]

psmelt(ps_prim) %>% 
  ggplot(aes(x = AB, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") 



# Plots of relative abundances, fixing some genes that are clustered in the data twice, showing top 12 taxa and others are clustered

# relative abundance of Phyla shown by age and farm

subsetMG %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "Phylum", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of genera by farm and age")



# Relative abundance for both stable and antibiotics used

subsetMG %>% 
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "Phylum", bar_width = 0.7, sample_order = "asis", 
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
  ggtitle("Relative abundance of genera by stable and antibiotics used")


# Same plots but with Order

subsetMG %>% tax_fix(unknowns = c("o__")) %>%
  ps_arrange(Farm2) %>%
  ps_mutate(
    Farm2 = factor(Farm2, rev(unique(Farm2)))
  ) %>%
  comp_barplot(
    tax_level = "Order", bar_width = 0.7, sample_order = "asis", 
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


subsetMG %>% tax_fix(unknowns = c("o__")) %>% 
  ps_arrange(FarmRoundStable) %>%
  ps_mutate(
    FarmRoundStable = factor(FarmRoundStable, rev(unique(FarmRoundStable)))
  ) %>%
  comp_barplot(
    tax_level = "Order", bar_width = 0.7, sample_order = "asis", 
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


# relabundance with other category
subsetMG.rel <- subsetMG  %>% aggregate_top_taxa2("Phylum", top = 5) %>% microbiome::transform("compositional")

plot_composition(subsetMG) + theme(legend.position = "bottom") +
  scale_fill_brewer("Phylum", palette = "Paired") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))


plot_composition(subsetMG.rel, x.label = "Id") + theme(legend.position = "bottom") +
  scale_fill_brewer("Phylum", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90))  + 
  ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun)


# declutter R environment by removing objects that no longer serve a purpose
rm(ps_prim)