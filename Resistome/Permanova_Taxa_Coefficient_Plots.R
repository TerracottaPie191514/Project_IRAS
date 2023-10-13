
# Used the following guides: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/multivariate-comparisons-of-microbial-community-composition.html, https://microbiome.github.io/course_2022_radboud/beta-diversity-demo.html

tse = makeTreeSummarizedExperimentFromPhyloseq(subsetMG)
tse <- transformCounts(tse, method = "relabundance")

tse_genus <- agglomerateByRank(tse, "Genus")
tse_genus <- transformCounts(tse_genus, method = "relabundance")


permanova_age <- adonis2(t(assay(tse_genus, "relabundance")) ~ Age, data = colData(tse_genus), permutations = 9999)
permanova_age
permanova_AB <- adonis(t(assay(tse_genus, "relabundance")) ~ AB, data = colData(tse_genus), permutations = 9999)
permanova_cox <- adonis(t(assay(tse_genus, "relabundance")) ~ Cox, data = colData(tse_genus), permutations = 9999)
permanova_researcher <- adonis(t(assay(tse_genus, "relabundance")) ~ Researcher, data = colData(tse_genus), permutations = 9999)
permanova_feedpr <- adonis(t(assay(tse_genus, "relabundance")) ~ FeedProducent, data = colData(tse_genus), permutations = 9999)
permanova_LT <- adonis(t(assay(tse_genus, "relabundance")) ~ LitterType, data = colData(tse_genus), permutations = 9999)
permanova_FT <- adonis(t(assay(tse_genus, "relabundance")) ~ FeedType, data = colData(tse_genus), permutations = 9999)
permanova_Gender <- adonis(t(assay(tse_genus, "relabundance")) ~ Gender, data = colData(tse_genus), permutations = 9999)
permanova_Stable <- adonis(t(assay(tse_genus, "relabundance")) ~ Stable, data = colData(tse_genus), permutations = 9999)
permanova_FS <- adonis(t(assay(tse_genus, "relabundance")) ~ FlockSize, data = colData(tse_genus), permutations = 9999)
permanova_farm <- adonis(t(assay(tse_genus, "relabundance")) ~ Farm2, data = colData(tse_genus), permutations = 9999)
permanova_AP <- adonis(t(assay(tse_genus, "relabundance")) ~ AgeParentStock, data = colData(tse_genus), permutations = 9999)




permanova_feedpr$aov.tab
coef <- coefficients(permanova_feedpr)["FeedProducent1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

ggplot(data.frame(x = top.coef,
                  y = factor(names(top.coef),
                             unique(names(top.coef)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="",y="",title="Producent, P=1e-04, R^2 = 0.06593 ") +
  theme_bw()
