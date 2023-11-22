#Procrustes analyses
copy16S = subset16S
# get the samples in the same order
sample_names(copy16S) = sample_names(subsetMG)

PCoA_BC_16s = ordinate(copy16S, "PCoA")
PCoA_BC_MG = ordinate(subsetMG, "PCoA")
procrustes = protest(PCoA_BC_16s$vectors, PCoA_BC_MG$vectors)

plot_data <- data.frame(
  MT_PC1 = procrustes$X[, 1],
  MT_PC2 = procrustes$X[, 2],
  MG_PC1 = procrustes$Yrot[, 1],
  MG_PC2 = procrustes$Yrot[, 2])

# with arrows pointing from MG to 16S
ggplot(plot_data) +
  geom_point(aes(x=MT_PC1, y=MT_PC2), color = "blue") +
  geom_segment(aes(x=MT_PC1,y=MT_PC2,xend=MG_PC1,yend=MG_PC2),arrow=arrow(type = "closed", length=unit(0.2,"cm"))) +
  labs(title = "Procrustes Plot Metataxonomic vs metagenomic", x = "PC1", y = "PC2") + 
  scale_color_manual(values = c("16S" = "black", "MG" = "blue"))

# plot with both points
#ggplot(plot_data) +
#  geom_point(aes(x=MT_PC1, y=MT_PC2), color = "green") +
#  geom_point(aes(x=MG_PC1, y=MG_PC2), color = "blue") +
#  geom_segment(aes(x = MT_PC1, xend = MG_PC1, y = MT_PC2, yend = MG_PC2), linetype = "solid") +
#  labs(title = "Procrustes Plot Metataxonomic vs metagenomic")

# Resistome vs MG
copyRps = Rps
# get the samples in the same order
sample_names(copyRps) = sample_names(subsetMG)

PCoA_BC_Rps = ordinate(copyRps, "PCoA") 
PCoA_BC_MG = ordinate(subsetMG, "PCoA") 
procrustes = protest(PCoA_BC_Rps$vectors, PCoA_BC_MG$vectors)

plot_data <- data.frame(
  R_PC1 = procrustes$X[, 1],
  R_PC2 = procrustes$X[, 2],
  MG_PC1 = procrustes$Yrot[, 1],
  MG_PC2 = procrustes$Yrot[, 2])

# resistome (k2) points to MG
ggplot(plot_data) +
  geom_point(aes(x=R_PC1, y=R_PC2), color = "blue") +
  geom_segment(aes(x=R_PC1,y=R_PC2,xend=MG_PC1,yend=MG_PC2),arrow=arrow(type = "closed", length=unit(0.2,"cm"))) +
  scale_color_manual(values = c("16S" = "black", "MG" = "blue")) +
  guides(color = guide_legend(title = "Data Type")) +
  labs(title = "Procrustes Plot Resistome vs Metagenome", x = "PC1", y = "PC2") 

# adds labels to see if samples line up
ggplot(plot_data) +
  geom_point(aes(x=R_PC1, y=R_PC2)) +
  geom_point(aes(x=MG_PC1, y=MG_PC2), color = "blue")+
  geom_segment(aes(x=R_PC1,y=R_PC2,xend=MG_PC1,yend=MG_PC2),arrow=arrow(length=unit(0.2,"cm"))) + 
  geom_text(aes(x = MG_PC1, y = MG_PC2, label = rownames(plot_data))) +
  geom_text(aes(x = R_PC1, y = R_PC2, label = rownames(plot_data))) +
  labs(title = "Procrustes Plot")

# MP vs k2
PCoA_BC_MP = ordinate(Rps_mp, "PCoA") 
PCoA_BC_k2 = ordinate(Rps, "PCoA") 
procrustes = protest(PCoA_BC_MP$vectors, PCoA_BC_k2$vectors)

plot_data <- data.frame(
  MP_PC1 = procrustes$X[, 1],
  MP_PC2 = procrustes$X[, 2],
  k2_PC1 = procrustes$Yrot[, 1],
  k2_PC2 = procrustes$Yrot[, 2])

ggplot(plot_data) +
  geom_point(aes(x=MP_PC1, y=MP_PC2), color = "blue") +
  geom_segment(aes(x=MP_PC1,y=MP_PC2,xend=k2_PC1,yend=k2_PC2),arrow=arrow(type = "closed", length=unit(0.2,"cm"))) +
  labs(title = "Procrustes Plot MetaPhlAn vs Kraken 2", x = "PC1", y = "PC2") + 
  scale_color_manual(values = c("16S" = "black", "MG" = "blue")) +
  guides(color = guide_legend(title = "Data Type"))

  # declutter R environment by removing objects that no longer serve a purpose
rm(PCoA_BC_16s, PCoA_BC_MG, PCoA_BC_Rps, PCoA_BC_MP, PCoA_BC_k2, plot_data, procrustes, copy16S, copyRps) 
