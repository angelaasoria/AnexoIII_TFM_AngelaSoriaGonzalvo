library(Matrix)
library(igraph)
library(scran)
library(pheatmap)
library(scater)
library(patchwork)
library(bluster)
library(SingleCellExperiment)
library(MetBrewer)

#----------------------- Load data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
sce.dim.reduc <- readRDS('sce.dim.reduc.rds')

# Generate palettes
palette <- met.brewer('Renoir', n = 14)

#----------------------- Two steps clustering ----------------------------------

center <- round(dim(sce.dim.reduc)[2] / 4)

set.seed(100)
# Para k = 15
cluster.15 <- clusterCells(sce.dim.reduc, use.dimred="PCA", full=TRUE,
                           BLUSPARAM=TwoStepParam(first = KmeansParam(centers = center),
                                                  second = NNGraphParam(cluster.fun = 'louvain',
                                                                        k = 15)))
sce.dim.reduc$cluster15 <- cluster.15$clusters
cluster.15 <- sce.dim.reduc$cluster15
table(sce.dim.reduc$cluster15)

# Para k = 25
cluster.25 <- clusterCells(sce.dim.reduc, use.dimred="PCA", full=TRUE,
                           BLUSPARAM=TwoStepParam(first = KmeansParam(centers = center),
                                                  second = NNGraphParam(cluster.fun = 'louvain',
                                                                        k = 25)))
sce.dim.reduc$cluster25 <- cluster.25$clusters
cluster.25 <- sce.dim.reduc$cluster25
table(sce.dim.reduc$cluster25)

# Para k = 40
cluster.40 <- clusterCells(sce.dim.reduc, use.dimred="PCA", full=TRUE,
                           BLUSPARAM=TwoStepParam(first = KmeansParam(centers = center),
                                                  second = NNGraphParam(cluster.fun = 'louvain',
                                                                        k = 40)))
sce.dim.reduc$cluster40 <- cluster.40$clusters
cluster.40 <- sce.dim.reduc$cluster40
table(sce.dim.reduc$cluster40)


#----------------------- Option 2: Louvain method ------------------------------

set.seed(100)
# Para k = 15
cluster.15.lou <- clusterCells(sce.dim.reduc, use.dimred="PCA", full=TRUE,
                           BLUSPARAM = NNGraphParam(cluster.fun = 'louvain',k = 15))
sce.dim.reduc$cluster15.lou <- cluster.15.lou$clusters
cluster.15.lou <- sce.dim.reduc$cluster15.lou
table(sce.dim.reduc$cluster15.lou)

# Para k = 40
cluster.40.lou <- clusterCells(sce.dim.reduc, use.dimred="PCA", full=TRUE,
                           BLUSPARAM= NNGraphParam(cluster.fun = 'louvain',k = 40))
sce.dim.reduc$cluster40.lou <- cluster.40.lou$clusters
cluster.40.lou <- sce.dim.reduc$cluster40.lou
table(sce.dim.reduc$cluster40.lou)


#----------------------- Cluster plots -----------------------------------------

# PCA
PCA15 <- plotReducedDim(sce.dim.reduc, dimred="PCA",colour_by="cluster15", text_by="cluster15")+
          xlab("PCA1")+
          ylab("PCA2")+
          ggtitle('k = 15')+
          scale_color_manual(values = palette, name = "Grupos")

# tSNE
TSNE15 <-plotReducedDim(sce.dim.reduc, dimred="TSNE.seed.100.per.90",colour_by="cluster15", text_by="cluster15")+
          xlab("tSNE1")+
          ylab("tSNE2")+
          ggtitle('k = 15')+
          scale_color_manual(values = palette, name = "Grupos")

# UMAP
UMAP15 <-plotReducedDim(sce.dim.reduc, dimred="UMAP.seed.100.n300",colour_by="cluster15", text_by="cluster15")+
          xlab("UMAP1")+
          ylab("UMAP2")+
          ggtitle('k = 15')+
          scale_color_manual(values = palette, name = "Grupos")


# PCA
PCA25 <- plotReducedDim(sce.dim.reduc, dimred="PCA",colour_by="cluster25", text_by="cluster25")+
  xlab("PCA1")+
  ylab("PCA2")+
  ggtitle('k = 25')+
  scale_color_manual(values = palette, name = "Grupos")

# tSNE
TSNE25 <-plotReducedDim(sce.dim.reduc, dimred="TSNE.seed.100.per.90",colour_by="cluster25", text_by="cluster25")+
  xlab("tSNE1")+
  ylab("tSNE2")+
  ggtitle('k = 25')+
  scale_color_manual(values = palette, name = "Grupos")

# UMAP
UMAP25 <-plotReducedDim(sce.dim.reduc, dimred="UMAP.seed.100.n300",colour_by="cluster25", text_by="cluster25")+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle('k = 25')+
  scale_color_manual(values = palette, name = "Grupos")


# PCA
PCA40 <- plotReducedDim(sce.dim.reduc, dimred="PCA",colour_by="cluster40", text_by="cluster40")+
          xlab("PCA1")+
          ylab("PCA2")+
          ggtitle('k = 40')+
          scale_color_manual(values = palette, name = "Grupos")

# tSNE
TSNE40 <- plotReducedDim(sce.dim.reduc, dimred="TSNE.seed.100.per.90",colour_by="cluster40", text_by="cluster40")+
            xlab("tSNE1")+
            ylab("tSNE2")+
            ggtitle('k = 40')+
            scale_color_manual(values = palette, name = "Grupos")

# UMAP
UMAP40 <- plotReducedDim(sce.dim.reduc, dimred="UMAP.seed.100.n300",colour_by="cluster40", text_by="cluster40")+
            xlab("UMAP1")+
            ylab("UMAP2")+
            ggtitle('k = 40')+
            scale_color_manual(values = palette, name = "Grupos")

comb.pcas.two <- PCA15 + PCA40
comb.tsnes.two <- TSNE15 + TSNE40
comb.umaps.two <- UMAP15 + UMAP40

jpeg(file = './figures/5.clustering/PCA_k15k40.jpeg', width=10, height=6, units="in", res=300)
print(comb.pcas.two)
dev.off()
svg(file = './figures/5.clustering/PCA_k15k40.svg', width=10, height=6)
print(comb.pcas.two)
dev.off()

jpeg(file = './figures/5.clustering/TSNE_k15k40.jpeg', width=10, height=6, units="in", res=300)
print(comb.tsnes.two)
dev.off()
svg(file = './figures/5.clustering/TSNE_k15k40.svg', width=10, height=6)
print(comb.tsnes.two)
dev.off()

jpeg(file = './figures/5.clustering/UMAP_k15k40.jpeg', width=10, height=6, units="in", res=300)
print(comb.umaps.two)
dev.off()
svg(file = './figures/5.clustering/UMAP_k15k40.svg', width=10, height=6)
print(comb.umaps.two)
dev.off()

#----------------------- Clustering evaluation ---------------------------------

##### k = 15 #####
# 1. Cluster purity
purity.15 <- neighborPurity(reducedDim(sce.dim.reduc, 'PCA'), 
                            clusters = sce.dim.reduc$cluster15)
purity.df.15 <- as.data.frame(purity.15)
purity.df.15$maximum <- factor(purity.df.15$maximum)
purity.df.15$cluster <- factor(sce.dim.reduc$cluster15)
plot.purity.15 <- ggplot(purity.df.15, aes(x = cluster, 
                                           y = purity, colour = maximum))+
  ggbeeswarm::geom_quasirandom(method = 'smiley')

# 2. Cluster width
width.15 <- approxSilhouette(reducedDim(sce.dim.reduc, 'PCA'), 
                          clusters = sce.dim.reduc$cluster15)
width.df.15 <- as.data.frame(width.15)
width.df.15$closest <- factor(ifelse(width.df.15$width > 0, 
                                  sce.dim.reduc$cluster15,
                                  width.df.15$other))
width.df.15$cluster <- factor(sce.dim.reduc$cluster15)
plot.width.15 <- ggplot(width.df.15, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")+
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = -0.5)

# 3. RMSD per cluster
rmsd.15 <- clusterRMSD(reducedDim(sce.dim.reduc, 'PCA'), sce.dim.reduc$cluster15)
rmsd.df.15 <- data.frame(rmsd = rmsd.15, cluster = names(rmsd.15))
plot.rmsd.15 <- ggplot(rmsd.df.15, aes(x=cluster, y=rmsd, fill = cluster))+ 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 5) +
  geom_hline(yintercept = 10) 

# 4. RMSD vs size factors
by.clust.15 <- split(sizeFactors(sce.dim.reduc), sce.dim.reduc$cluster15)
sf.by.clust.15 <- vapply(by.clust.15, mean, 0)
point.df.15 <- data.frame(rmsd = rmsd.15, sf = sf.by.clust.15, label = names(rmsd.15))
pointplot.15 <- ggplot(point.df.15, aes(rmsd,sf)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = label)) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 2, linetype = 'dashed', color = "#a8554e") + 
  geom_vline(xintercept = 5) +
  geom_vline(xintercept = 10, linetype = 'dashed', color = "#a8554e")

##### k = 40 #####
# 1. Cluster purity
purity.40 <- neighborPurity(reducedDim(sce.dim.reduc, 'PCA'), 
                            clusters = sce.dim.reduc$cluster40)
purity.df.40 <- as.data.frame(purity.40)
purity.df.40$maximum <- factor(purity.df.40$maximum)
purity.df.40$cluster <- factor(sce.dim.reduc$cluster40)
plot.purity.40 <- ggplot(purity.df.40, aes(x = cluster, y = purity, colour = maximum))+
  ggbeeswarm::geom_quasirandom(method = 'smiley')

# 2. Cluster width
width.40 <- approxSilhouette(reducedDim(sce.dim.reduc, 'PCA'), 
                             clusters = sce.dim.reduc$cluster40)
width.df.40 <- as.data.frame(width.40)
width.df.40$closest <- factor(ifelse(width.df.40$width > 0, 
                                     sce.dim.reduc$cluster40,
                                     width.df.40$other))
width.df.40$cluster <- factor(sce.dim.reduc$cluster40)
plot.width.40 <- ggplot(width.df.40, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")+
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = -0.5)

# 3. RMSD per cluster
rmsd.40 <- clusterRMSD(reducedDim(sce.dim.reduc, 'PCA'), sce.dim.reduc$cluster40)
rmsd.df.40 <- data.frame(rmsd = rmsd.40, cluster = names(rmsd.40))
plot.rmsd.40 <- ggplot(rmsd.df.40, aes(x=cluster, y=rmsd, fill = cluster))+ 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 5) +
  geom_hline(yintercept = 10) 

# 4. RMSD vs size factors
by.clust.40 <- split(sizeFactors(sce.dim.reduc), sce.dim.reduc$cluster40)
sf.by.clust.40 <- vapply(by.clust.40, mean, 0)
point.df.40 <- data.frame(rmsd = rmsd.40, sf = sf.by.clust.40, label = names(rmsd.40))
pointplot.40 <- ggplot(point.df.40, aes(rmsd,sf)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = label)) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 2, linetype = 'dashed', color = "#a8554e") + 
  geom_vline(xintercept = 5) +
  geom_vline(xintercept = 10, linetype = 'dashed', color = "#a8554e")

plot.purity <- plot.purity.15 + plot.purity.40
plot.width <- plot.width.15 + plot.width.40
plot.rmsd <- plot.rmsd.15 + plot.rmsd.40
plot.pointplot <- pointplot.15 + pointplot.40

jpeg(file = './figures/5.clustering/purity_k15k40.jpeg', width=10, height=6, units="in", res=300)
print(plot.purity)
dev.off()
svg(file = './figures/5.clustering/purity_k15k40.svg', width=10, height=6)
print(plot.purity)
dev.off()

jpeg(file = './figures/5.clustering/width_k15k40.jpeg', width=10, height=6, units="in", res=300)
print(plot.width)
dev.off()
svg(file = './figures/5.clustering/width_k15k40.svg', width=10, height=6)
print(plot.width)
dev.off()

jpeg(file = './figures/5.clustering/rmsd_k15k40.jpeg', width=10, height=6, units="in", res=300)
print(plot.rmsd)
dev.off()
svg(file = './figures/5.clustering/rmsd_k15k40.svg', width=10, height=6)
print(plot.rmsd)
dev.off()

jpeg(file = './figures/5.clustering/pointplot_k15k40.jpeg', width=10, height=6, units="in", res=300)
print(plot.pointplot)
dev.off()
svg(file = './figures/5.clustering/pointplot_k15k40.svg', width=10, height=6)
print(plot.pointplot)
dev.off()

multipanel.eval <- plot.purity / plot.pointplot
multipanel.clust <- comb.pcas.two / comb.tsnes.two / comb.umaps.two
multipanel.eval
#----------------------- Clustering evaluation - LOUVAIN -----------------------

##### k = 15 #####
# 1. Cluster purity
purity.15.lou <- neighborPurity(reducedDim(sce.dim.reduc, 'PCA'), 
                            clusters = sce.dim.reduc$cluster15.lou)
purity.df.15.lou <- as.data.frame(purity.15.lou)
purity.df.15.lou$maximum <- factor(purity.df.15.lou$maximum)
purity.df.15.lou$cluster <- factor(sce.dim.reduc$cluster15.lou)
plot.purity.15.lou <- ggplot(purity.df.15.lou, aes(x = cluster, 
                                           y = purity, colour = maximum))+
  ggbeeswarm::geom_quasirandom(method = 'smiley')

# 2. Cluster width
width.15.lou <- approxSilhouette(reducedDim(sce.dim.reduc, 'PCA'), 
                             clusters = sce.dim.reduc$cluster15.lou)
width.df.15.lou <- as.data.frame(width.15.lou)
width.df.15.lou$closest <- factor(ifelse(width.df.15.lou$width > 0, 
                                     sce.dim.reduc$cluster15.lou,
                                     width.df.15.lou$other))
width.df.15.lou$cluster <- factor(sce.dim.reduc$cluster15.lou)
plot.width.15.lou <- ggplot(width.df.15.lou, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")+
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = -0.5)

# 3. RMSD per cluster
rmsd.15.lou <- clusterRMSD(reducedDim(sce.dim.reduc, 'PCA'), sce.dim.reduc$cluster15.lou)
rmsd.df.15.lou <- data.frame(rmsd = rmsd.15.lou, cluster = names(rmsd.15.lou))
plot.rmsd.15.lou <- ggplot(rmsd.df.15.lou, aes(x=cluster, y=rmsd, fill = cluster))+ 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 5) +
  geom_hline(yintercept = 10) 

# 4. RMSD vs size factors
by.clust.15.lou <- split(sizeFactors(sce.dim.reduc), sce.dim.reduc$cluster15.lou)
sf.by.clust.15.lou <- vapply(by.clust.15.lou, mean, 0)
point.df.15.lou <- data.frame(rmsd = rmsd.15.lou, sf = sf.by.clust.15.lou, label = names(rmsd.15.lou))
pointplot.15.lou <- ggplot(point.df.15.lou, aes(rmsd,sf)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = label)) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 2, linetype = 'dashed', color = "#a8554e") + 
  geom_vline(xintercept = 5) +
  geom_vline(xintercept = 10, linetype = 'dashed', color = "#a8554e")


##### k = 40 #####
# 1. Cluster purity
purity.40.lou <- neighborPurity(reducedDim(sce.dim.reduc, 'PCA'), 
                                clusters = sce.dim.reduc$cluster40.lou)
purity.df.40.lou <- as.data.frame(purity.40.lou)
purity.df.40.lou$maximum <- factor(purity.df.40.lou$maximum)
purity.df.40.lou$cluster <- factor(sce.dim.reduc$cluster40.lou)
plot.purity.40.lou <- ggplot(purity.df.40.lou, aes(x = cluster, 
                                                   y = purity, colour = maximum))+
  ggbeeswarm::geom_quasirandom(method = 'smiley')

# 2. Cluster width
width.40.lou <- approxSilhouette(reducedDim(sce.dim.reduc, 'PCA'), 
                                 clusters = sce.dim.reduc$cluster40.lou)
width.df.40.lou <- as.data.frame(width.40.lou)
width.df.40.lou$closest <- factor(ifelse(width.df.40.lou$width > 0, 
                                         sce.dim.reduc$cluster40.lou,
                                         width.df.40.lou$other))
width.df.40.lou$cluster <- factor(sce.dim.reduc$cluster40.lou)
plot.width.40.lou <- ggplot(width.df.40.lou, aes(x=cluster, y=width, colour=closest)) +
  ggbeeswarm::geom_quasirandom(method="smiley")+
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = -0.5)

# 3. RMSD per cluster
rmsd.40.lou <- clusterRMSD(reducedDim(sce.dim.reduc, 'PCA'), sce.dim.reduc$cluster40.lou)
rmsd.df.40.lou <- data.frame(rmsd = rmsd.40.lou, cluster = names(rmsd.40.lou))
plot.rmsd.40.lou <- ggplot(rmsd.df.40.lou, aes(x=cluster, y=rmsd, fill = cluster))+ 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 5) +
  geom_hline(yintercept = 10) 

# 4. RMSD vs size factors
by.clust.40.lou <- split(sizeFactors(sce.dim.reduc), sce.dim.reduc$cluster40.lou)
sf.by.clust.40.lou <- vapply(by.clust.40.lou, mean, 0)
point.df.40.lou <- data.frame(rmsd = rmsd.40.lou, sf = sf.by.clust.40.lou, label = names(rmsd.40.lou))
pointplot.40.lou <- ggplot(point.df.40.lou, aes(rmsd,sf)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = label)) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 2, linetype = 'dashed', color = "#a8554e") + 
  geom_vline(xintercept = 5) +
  geom_vline(xintercept = 10, linetype = 'dashed', color = "#a8554e")

plot.purity.lou <- plot.purity.15.lou + plot.purity.40.lou
plot.width.lou <- plot.width.15.lou + plot.width.40.lou
plot.rmsd.lou <- plot.rmsd.15.lou + plot.rmsd.40.lou
plot.pointplot.lou <- pointplot.15.lou + pointplot.40.lou

jpeg(file = './figures/5.clustering/purity_k15k40_lou.jpeg', width=10, height=6, units="in", res=300)
print(plot.purity.lou)
dev.off()
svg(file = './figures/5.clustering/purity_k15k40_lou.svg', width=10, height=6)
print(plot.purity.lou)
dev.off()

jpeg(file = './figures/5.clustering/width_k15k40_lou.jpeg', width=10, height=6, units="in", res=300)
print(plot.width.lou)
dev.off()
svg(file = './figures/5.clustering/width_k15k40_lou.svg', width=10, height=6)
print(plot.width.lou)
dev.off()

jpeg(file = './figures/5.clustering/rmsd_k15k40_lou.jpeg', width=10, height=6, units="in", res=300)
print(plot.rmsd.lou)
dev.off()
svg(file = './figures/5.clustering/rmsd_k15k40_lou.svg', width=10, height=6)
print(plot.rmsd.lou)
dev.off()

jpeg(file = './figures/5.clustering/pointplot_k15k40_lou.jpeg', width=10, height=6, units="in", res=300)
print(plot.pointplot.lou)
dev.off()
svg(file = './figures/5.clustering/pointplot_k15k40_lou.svg', width=10, height=6)
print(plot.pointplot.lou)
dev.off()

#----------------------- Save data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
saveRDS(sce.dim.reduc, 'sce.cluster.rds')
