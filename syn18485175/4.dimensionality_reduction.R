library(scran)
library(PCAtools)
library(scater)
library(SingleCellExperiment)
library(ggthemes)
library(MetBrewer)
library(patchwork)

#----------------------- Load data --------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/results')
sce.hvg <- readRDS('./sce.hvg.rds')

# Generate palettes
palette <- met.brewer('Renoir')
palette2 <- c(palette[4], palette[9])
palette4 <- c(palette[4], palette[8], palette[9], palette[11])

#----------------------- Principal component analysis --------------------------

sce.hvg <- runPCA(sce.hvg, exprs_values = 'logcounts')
percent.var <- attr(SingleCellExperiment::reducedDim(sce.hvg), 'percentVar')
elbow.point <- PCAtools::findElbowPoint(percent.var)
# Number of PCAs chosen by elbow method
elbow.point

jpeg(file = './figures/4.dimensionality_reduction/variability.percentage.jpeg', width=10, height=6, units="in", res=300)
plot(percent.var, xlab = 'Componente principal', ylab = 'Porcentaje de variabilidad explicada (%)')
abline(v = elbow.point, col = "#a8554e")
dev.off()
svg(file = './figures/4.dimensionality_reduction/variability.percentage.svg', width=10, height=6)
plot(percent.var, xlab = 'Componente principal', ylab = 'Porcentaje de variabilidad explicada (%)')
abline(v = elbow.point, col = "#a8554e")
dev.off()

#----------------------- Plot PCA ----------------------------------------------

reducedDim(sce.hvg, 'PCA') <- reducedDim(sce.hvg, 'PCA')[ ,1:elbow.point]
set.seed(100)

# By sex
plot.pca.sex <- plotReducedDim(sce.hvg, dimred="PCA", colour_by="sex", point_alpha = 0.3)+
  xlab("PCA1")+
  ylab("PCA2")+
  scale_color_manual(values = palette2, name = "Sexo", labels = c('Hombre', 'Mujer'))

# By condition
plot.pca.condition <- plotReducedDim(sce.hvg, dimred="PCA", colour_by="condition", point_alpha = 0.3)+
  xlab("PCA1")+
  ylab("PCA2")+
  scale_color_manual(values = palette2, name = "Condición", labels = c('AD', 'Control'))

# By condition and sex
plot.pca.groupsex <- plotReducedDim(sce.hvg, dimred="PCA", colour_by="group.sex", point_alpha = 0.3)+
  xlab("PCA1")+
  ylab("PCA2")+
  scale_color_manual(values = palette4, name = "Condición y sexo", 
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

comb.plots.pca <- plot.pca.groupsex + plot.pca.condition / plot.pca.sex

jpeg(file = './figures/4.dimensionality_reduction/pca.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.pca)
dev.off()
svg(file = './figures/4.dimensionality_reduction/pca.svg', width=10, height=6)
print(comb.plots.pca)
dev.off()

## Elbow components
# By sex
plot.elbow.pca.sex <- plotReducedDim(sce.hvg, dimred="PCA", colour_by="sex", 
                                     ncomponents = elbow.point)+
  scale_color_manual(values = palette2, name = "Sexo", labels = c('Hombre', 'Mujer'))

# By condition
plot.elbow.pca.condition <- plotReducedDim(sce.hvg, dimred="PCA", colour_by="condition", 
                                           ncomponents = elbow.point)+
  scale_color_manual(values = palette2, name = "Condición", labels = c('AD', 'Control'))

# By condition and sex
plot.elbow.pca.groupsex <-plotReducedDim(sce.hvg, dimred="PCA", colour_by="group.sex", 
                                         ncomponents = elbow.point)+
  scale_color_manual(values = palette4, name = "Condición y sexo", 
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

jpeg(file = './figures/4.dimensionality_reduction/pca.elbow.sex.jpeg', width=10, height=6, units="in", res=300)
print(plot.elbow.pca.sex)
dev.off()
svg(file = './figures/4.dimensionality_reduction/pca.elbow.sex.svg', width=10, height=6)
print(plot.elbow.pca.sex)
dev.off()

jpeg(file = './figures/4.dimensionality_reduction/pca.elbow.condition.jpeg', width=10, height=6, units="in", res=300)
print(plot.elbow.pca.condition)
dev.off()
svg(file = './figures/4.dimensionality_reduction/pca.elbow.condition.svg', width=10, height=6)
print(plot.elbow.pca.condition)
dev.off()

jpeg(file = './figures/4.dimensionality_reduction/pca.elbow.groupsex.jpeg', width=10, height=6, units="in", res=300)
print(plot.elbow.pca.groupsex)
dev.off()
svg(file = './figures/4.dimensionality_reduction/pca.elbow.groupsex.svg', width=10, height=6)
print(plot.elbow.pca.groupsex)
dev.off()

#----------------------- Plot UMAP ---------------------------------------------

# Neighbors 300
set.seed(100)
sce.hvg <- runUMAP(sce.hvg, dimred = 'PCA', n_neighbors = 300)
reducedDim(sce.hvg, 'UMAP.seed.100.n300') <- reducedDim(sce.hvg, 'UMAP')

# By sex
plot.umap.sex.n300 <- plotReducedDim(sce.hvg, dimred = 'UMAP.seed.100.n300', 
                                     colour_by = "sex")+
  xlab('UMAP1')+
  ylab('UMAP2')+
  scale_color_manual(values = palette2, 
                     name = "Sexo")

# By condition
plot.umap.condition.n300 <- plotReducedDim(sce.hvg, dimred = 'UMAP.seed.100.n300', 
                                           colour_by = "condition")+
  xlab('UMAP1')+
  ylab('UMAP2')+
  scale_color_manual(values = palette2, 
                     name = "Condición", labels = c('AD', 'Control'))

# By condition and sex
plot.umap.groupsex.n300 <- plotReducedDim(sce.hvg, dimred = 'UMAP.seed.100.n300', 
                                          colour_by = "group.sex")+
  xlab('UMAP1')+
  ylab('UMAP2')+
  scale_color_manual(values = palette4, name = "Condición", 
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre','Control_Mujer'))

comb.plots.umap.n300 <- plot.umap.groupsex.n300 + plot.umap.condition.n300 / plot.umap.sex.n300 

jpeg(file = './figures/4.dimensionality_reduction/umap.n300.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.umap.n300)
dev.off()
svg(file = './figures/4.dimensionality_reduction/umap.n300.svg', width=10, height=6)
print(comb.plots.umap.n300)
dev.off()

#----------------------- tSNE --------------------------------------------------

set.seed(100)
# Perplexity 5
sce.hvg <- scater::runTSNE(sce.hvg, dimred="PCA", perplexity=5)
out5 <- plotReducedDim(sce.hvg, dimred="TSNE", colour_by="group.sex") + 
  ggtitle("perplexity = 5")+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette4, name = "Condición y sexo", 
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre','Control_Mujer'))

# Perplexity 20
sce.hvg <- scater::runTSNE(sce.hvg, dimred="PCA", perplexity=20)
out20 <- plotReducedDim(sce.hvg, dimred="TSNE", colour_by="group.sex") + 
  ggtitle("perplexity = 20")+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette4, name = "Condición y sexo", 
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre','Control_Mujer'))

# Perplexity 80
sce.hvg <- scater::runTSNE(sce.hvg, dimred="PCA", perplexity=80)
out80 <- plotReducedDim(sce.hvg, dimred="TSNE", colour_by="group.sex") + 
  ggtitle("perplexity = 80")+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette4, name = "Condición y sexo", 
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre','Control_Mujer'))

# Perplexity 120
sce.hvg <- scater::runTSNE(sce.hvg, dimred="PCA", perplexity=120)
out120 <- plotReducedDim(sce.hvg, dimred="TSNE", colour_by="group.sex") + 
  ggtitle("perplexity = 120")+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette4, name = "Condición y sexo", 
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre','Control_Mujer'))

comb.plots.tsne <- (out5 + out20) / (out80 + out120)

jpeg(file = './figures/4.dimensionality_reduction/tsne.perplexities.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.tsne)
dev.off()
svg(file = './figures/4.dimensionality_reduction/tsne.perplexities.svg', width=10, height=6)
print(comb.plots.tsne)
dev.off()

#----------------------- Plot tSNE ---------------------------------------------

# Chosen perplexity:90
set.seed(100)
sce.hvg <- scater::runTSNE(sce.hvg, dimred="PCA", perplexity = 90)
reducedDim(sce.hvg, "TSNE.seed.100.per.90") <- reducedDim(sce.hvg, "TSNE")

# By condition
plot.tsne.condition <- plotReducedDim(sce.hvg, dimred="TSNE.seed.100.per.90", colour_by="condition", point_alpha=0.3)+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette2, name = "Condición")

# By sex
plot.tsne.sex <- plotReducedDim(sce.hvg, dimred="TSNE.seed.100.per.90", colour_by="sex", point_alpha=0.3)+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette2, name = "Sexo")

# By group and sex
plot.tsne.groupsex <- plotReducedDim(sce.hvg, dimred="TSNE.seed.100.per.90", colour_by="group.sex", point_alpha=0.3)+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette4, name = "Condición y sexo")

comb.plots.tsne <- plot.tsne.groupsex + plot.tsne.condition / plot.tsne.sex 

jpeg(file = './figures/4.dimensionality_reduction/tsne.per90.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.tsne)
dev.off()
svg(file = './figures/4.dimensionality_reduction/tsne.per90.svg', width=10, height=6)
print(comb.plots.tsne)
dev.off()

#----------------------- Save ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/results')
saveRDS(sce.hvg, './sce.dim.reduc.rds')


