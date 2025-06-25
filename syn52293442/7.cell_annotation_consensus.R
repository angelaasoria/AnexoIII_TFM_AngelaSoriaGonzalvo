library(scater)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(openxlsx)

### Script for visualizing the results of the final cell annotation and for
#   the dataset processing for creating sce.annot.

pacman::p_load(SingleCellExperiment, scater, scran, dplyr, BiocParallel, 
               scDblFinder, AnnotationHub, tidyverse, patchwork, ggvenn, broom,
               kableExtra,HDF5Array, viridis, DropletUtils,
               devtools, celda, ggvenn, bluster, cluster, ggpubr, SCINA, SingleR,CHETAH)

#----------------------- Load data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn52293442_pfc/results/')
sce.annot <- readRDS('./sce.annot.rds')
sce.norm <- readRDS('./sce.norm.rds')

# Generate palettes
palette <- met.brewer('Renoir')
palette2 <- c(palette[4], palette[9])
palette4 <- c(palette[4], palette[8], palette[9], palette[11])
palette7 <- c(palette[4], palette[5], palette[6], palette[8], palette[9],
              palette[11], palette[12])

#----------------------- Creation of consensus results -------------------------

# Transform intermediate nodes in unknown 
sce.annot$CHETAH.labels <- str_replace(sce.annot$CHETAH.labels, "Node1", "unknown")
sce.annot$CHETAH.labels <- str_replace(sce.annot$CHETAH.labels, "Node2", "unknown")
sce.annot$CHETAH.labels <- str_replace(sce.annot$CHETAH.labels, "Node3", "unknown")
sce.annot$CHETAH.labels <- str_replace(sce.annot$CHETAH.labels, "Node4", "unknown")
sce.annot$CHETAH.labels <- str_replace(sce.annot$CHETAH.labels, "Unassigned", "unknown")
sce.annot$SCINA.labels <- str_replace(sce.annot$SCINA.labels, "endotelio", "unknown")
sce.annot$SCINA.labels <- str_replace(sce.annot$SCINA.labels, "fibroblasto", "unknown")

# Combine results from 3 tools
all.results <- data.frame(sce.annot$SCINA.labels,
                          sce.annot$SingleR.labels,
                          sce.annot$CHETAH.labels)
# Make tsv file from consensus
all.results$cell.type <- apply(all.results, 1, function(x) names(which.max(table(as.character(x)))))
all.results$n.support <- apply(all.results[, -dim(all.results)[2]], 1, function(x) max(table(as.character(x))))
all.results$cell.type[which(all.results$n.support == 1)] = 'unknown'
apply(all.results, 2, table)
write.table(all.results, file = './matrix.annot.tsv', sep = '\t', row.names = FALSE)

#----------------------- Creation of annotated SCE -----------------------------

sce.annot$cell.type <- all.results$cell.type
sce.norm$cell.type <- all.results$cell.type
saveRDS(sce.norm, file = './sce.norm.rds')

#----------------------- Before filtering --------------------------------------

tab <- table(label = sce.annot$cell.type, cluster = sce.annot$cluster40)
y <- as.data.frame(tab)

# Cell types by cluster
celltype.BF <- ggplot(y, aes(x=cluster, y=Freq, fill=label)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = palette7) +
  ggtitle("Número de células por tipo celular por cluster") +
  xlab("Cluster") + 
  ylab("Número de células")
ggsave(filename = './figures/7.cell_annotation_consensus/celltype_BF.jpeg', 
       plot = celltype.BF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/celltype_BF.svg', 
       plot = celltype.BF, width = 10, height = 6, device = 'svg')

# Number of cells by cluster
ncells.BF <- ggplot(y, aes(x=label, y=Freq, fill=label)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = palette7) +
  ggtitle("Número de células por tipo celular") +
  xlab("Tipo celular") + 
  ylab("Número de células")
ggsave(filename = './figures/7.cell_annotation_consensus/ncells_BF.jpeg', 
       plot = ncells.BF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/ncells_BF.svg', 
       plot = ncells.BF, width = 10, height = 6, device = 'svg')

# Heatmap
jpeg(file = './figures/7.cell_annotation_consensus/heatmap_BF.jpeg', width=10, height=6, units="in", res=300)
pheatmap::pheatmap(log2(tab+10), cluster_rows = F, cluster_cols = F, 
                   color=colorRampPalette(c("white", "blue"))(101))
dev.off()
svg(file = './figures/7.cell_annotation_consensus/heatmap_BF.svg', width=10, height=6)
pheatmap::pheatmap(log2(tab+10), cluster_rows = F, cluster_cols = F, 
                   color=colorRampPalette(c("white", "blue"))(101))
dev.off()

#----------------------- Cell filtering ----------------------------------------

# Select unknown cells to remove them
unknown.remove <- grep('unknown', sce.annot$cell.type)
# Get which clusters these cells belong to
cluster.unknown <- sce.annot@colData[unknown.remove, 'cluster40']
# Number of cells to be removed from each cluster
removed <- table(cluster.unknown)
# Total of cells in each cluster
total <- table(sce.annot$cluster40)
# Cells kept in each cluster after removing unknown
kept <- total - removed

df.removed <- data.frame(removed = removed)
df.removed$type <- 'removed'
colnames(df.removed) <- c('Cluster', 'Freq', 'Type')

df.kept <- data.frame(kept = kept)
df.kept$type <- 'kept'
colnames(df.kept) <- c('Cluster', 'Freq', 'Type')

# Combine the rows of the two dataframes
df.combined <- as.data.frame(rbind(df.removed, df.kept))

eliminados.BF <- ggplot(df.combined, aes(x=Cluster, y=Freq, fill=Type)) +
                        geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
                        ggtitle("Células consenso") +
                        xlab("Cluster") + 
                        ylab("Número de células")+
                        scale_fill_manual(values = c('skyblue', 'tomato3'))
ggsave(filename = './figures/7.cell_annotation_consensus/eliminados_BF.jpeg', 
       plot = eliminados.BF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/eliminados_BF.svg', 
       plot = eliminados.BF, width = 10, height = 6, device = 'svg')

#tSNE final
tsne.BF <- plotReducedDim((sce.annot), dimred="TSNE.seed.100.per.90", colour_by="cell.type", text_by= 'cluster40')+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +  
  scale_fill_manual(values = palette7, aesthetics="colour")
ggsave(filename = './figures/7.cell_annotation_consensus/tsne_BF.jpeg', 
       plot = tsne.BF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/tsne_BF.svg', 
       plot = tsne.BF, width = 10, height = 6, device = 'svg')

#PCA final
pca.BF <- plotReducedDim((sce.annot), dimred="PCA", colour_by="cell.type", text_by= 'cluster40')+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +  
  scale_fill_manual(values = palette7, aesthetics="colour")
ggsave(filename = './figures/7.cell_annotation_consensus/pca_BF.jpeg', 
       plot = pca.BF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/pca_BF.svg', 
       plot = pca.BF, width = 10, height = 6, device = 'svg')

#UMAP final
umap.BF <- plotReducedDim((sce.annot), dimred="UMAP.seed.100.n300", colour_by="cell.type", text_by= 'cluster40')+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +  
  scale_fill_manual(values = palette7, aesthetics="colour")
ggsave(filename = './figures/7.cell_annotation_consensus/umap_BF.jpeg', 
       plot = umap.BF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/umap_BF.svg', 
       plot = umap.BF, width = 10, height = 6, device = 'svg')

# Check that the unknown are not condition-exclusive

umap.data <- as.data.frame(reducedDim(sce.annot, 'UMAP.seed.100.n300'))
umap.data$group.sex <- sce.annot$group.sex
umap.data$cell.type <- sce.annot$cell.type

disease.cell.unknown <- ggplot(umap.data, aes(x = V1, y = V2, color = cell.type)) +
  geom_point() +
  facet_wrap(. ~ group.sex) +
  theme(legend.key.size = unit(0.25, 'cm')) +
  scale_fill_manual(values = palette7, aesthetics = 'colour')
ggsave(filename = './figures/7.cell_annotation_consensus/disease_cell_unknown.jpeg', 
       plot = disease.cell.unknown, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/disease_cell_unknown.svg', 
       plot = disease.cell.unknown, width = 10, height = 6, device = 'svg')

# Filtering
# Cell distribution along the clusters
table(sce.annot$cluster40, sce.annot$cell.type)
sce.annot <- sce.annot[ ,!sce.annot$cell.type == 'unknown']

# See what cells are outliers and remove them manually
table(sce.annot$cluster40, sce.annot$cell.type)
sce.annot <- sce.annot[ , !(sce.annot$cell.type != 'Oligodendrocitos' & sce.annot$cluster40 == '2')]
sce.annot <- sce.annot[ , !(sce.annot$cell.type != 'OPC' & sce.annot$cluster40 == '3')]
sce.annot <- sce.annot[ , !(sce.annot$cell.type != 'Astrocitos' & sce.annot$cluster40 == '6')]
sce.annot <- sce.annot[ , !(sce.annot$cell.type != 'Microglía' & sce.annot$cluster40 == '7')]
sce.annot <- sce.annot[ , !(sce.annot$cell.type != 'Inhibidoras' & sce.annot$cluster40 == '9')]
sce.annot <- sce.annot[ , !(sce.annot$cell.type != 'Excitadoras' & sce.annot$cluster40 == '11')]


as.data.frame(table(sce.annot$cluster40, sce.annot$cell.type)) %>% 
  group_by(Var2)  %>% 
  summarise(cond_disp = sum(Freq))

#----------------------- After filtering --------------------------------------

tab <- table(label = sce.annot$cell.type, cluster = sce.annot$cluster40)
y <- as.data.frame(tab)

# Heatmap
jpeg(file = './figures/7.cell_annotation_consensus/heatmap_AF.jpeg', width=10, height=6, units="in", res=300)
pheatmap::pheatmap(log2(tab+10), cluster_rows = F, cluster_cols = F, 
                   color=colorRampPalette(c("white", "blue"))(101))
dev.off()
svg(file = './figures/7.cell_annotation_consensus/heatmap_AF.svg', width=10, height=6)
pheatmap::pheatmap(log2(tab+10), cluster_rows = F, cluster_cols = F, 
                   color=colorRampPalette(c("white", "blue"))(101))
dev.off()

# Select unknown cells to remove them
unknown.remove <- grep('unknown', sce.annot$cell.type)
# Get which clusters these cells belong to
cluster.unknown <- sce.annot@colData[unknown.remove, 'cluster40']
# Number of cells to be removed from each cluster
removed <- table(cluster.unknown)
# Total of cells in each cluster
total <- table(sce.annot$cluster40)
# Cells kept in each cluster after removing unknown
kept <- total - removed

df.removed <- data.frame(removed = removed)
df.removed$type <- 'removed'
colnames(df.removed) <- c('Cluster', 'Freq', 'Type')

df.kept <- data.frame(kept = kept)
df.kept$type <- 'kept'
colnames(df.kept) <- c('Cluster', 'Freq', 'Type')

# Combine the rows of the two dataframes
df.combined <- as.data.frame(rbind(df.removed, df.kept))

#tSNE final cell type
tsne.AF.celltype <- plotReducedDim((sce.annot), dimred="TSNE.seed.100.per.90", colour_by="cell.type", text_by= 'cluster40')+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +  
  scale_fill_manual(values = palette7, aesthetics="colour")
jpeg(file = './figures/7.cell_annotation_consensus/tsne_AF_celltype.jpeg', width=10, height=6, units="in", res=300)
print(tsne.AF.celltype)
dev.off()
svg(file = './figures/7.cell_annotation_consensus/tsne_AF_celltype.svg', width=10, height=6)
print(tsne.AF.celltype)
dev.off()

#tSNE final group and sex
tsne.AF.groupsex <- plotReducedDim((sce.annot), dimred="TSNE.seed.100.per.90", colour_by="group.sex", text_by= 'cluster40')+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +  
  scale_fill_manual(values = palette4, aesthetics="colour")
jpeg(file = './figures/7.cell_annotation_consensus/tsne_AF_groupsex.jpeg', width=10, height=6, units="in", res=300)
print(tsne.AF.groupsex)
dev.off()
svg(file = './figures/7.cell_annotation_consensus/tsne_AF_groupsex.svg', width=10, height=6)
print(tsne.AF.groupsex)
dev.off()


#PCA final
pca.AF <- plotReducedDim((sce.annot), dimred="PCA", colour_by="cell.type", text_by= 'cluster40')+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +  
  scale_fill_manual(values = palette7, aesthetics="colour")
jpeg(file = './figures/7.cell_annotation_consensus/pca_AF.jpeg', width=10, height=6, units="in", res=300)
print(pca.AF)
dev.off()
svg(file = './figures/7.cell_annotation_consensus/pca_AF.svg', width=10, height=6)
print(pca.AF)
dev.off()

#UMAP final
umap.AF <- plotReducedDim((sce.annot), dimred="UMAP.seed.100.n300", colour_by="cell.type", text_by= 'cluster40')+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +  
  scale_fill_manual(values = palette7, aesthetics="colour")
jpeg(file = './figures/7.cell_annotation_consensus/umap_AF.jpeg', width=10, height=6, units="in", res=300)
print(umap.AF)
dev.off()
svg(file = './figures/7.cell_annotation_consensus/umap_AF.svg', width=10, height=6)
print(umap.AF)
dev.off()


# More plots

# AD vs control per cluster
x = as.data.frame(table(Condition = sce.annot$condition, Cluster = sce.annot$cluster40))
ADvsCtrl.AF <- ggplot(x, aes(x=Cluster, y=Freq, fill=Condition)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = palette2) + 
  ggtitle("Enfermos vs controles por cluster") +
  xlab("Cluster") + 
  ylab("Número de células")
ggsave(filename = './figures/7.cell_annotation_consensus/ADvsCtrl_AF.jpeg', 
       plot = ADvsCtrl.AF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/ADvsCtrl_AF.svg', 
       plot = ADvsCtrl.AF, width = 10, height = 6, device = 'svg')

# Groups per cluster
x = as.data.frame(table(Condition = sce.annot$group.sex, Cluster = sce.annot$cluster40))
cluster.groups.AF <- ggplot(x, aes(x=Cluster, y=Freq, fill=Condition)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = palette4) + 
  ggtitle("Grupos por cluster") +
  xlab("Cluster") + 
  ylab("Número de células")
ggsave(filename = './figures/7.cell_annotation_consensus/cluster_groups_AF.jpeg', 
       plot = cluster.groups.AF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/cluster_groups_AF.svg', 
       plot = cluster.groups.AF, width = 10, height = 6, device = 'svg')

# Cell types per group
x = as.data.frame(table(Condition = sce.annot$group.sex, Cluster = sce.annot$cell.type))
cluster.groups.celltype.AF <- ggplot(x, aes(x=Cluster, y=Freq, fill=Condition)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = palette4) + 
  ggtitle("Tipos celulares por grupo") +
  xlab("Grupo") + 
  ylab("Número de células")
ggsave(filename = './figures/7.cell_annotation_consensus/cluster_groups_celltype_AF.jpeg', 
       plot = cluster.groups.celltype.AF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/cluster_groups_celltype_AF.svg', 
       plot = cluster.groups.celltype.AF, width = 10, height = 6, device = 'svg')


#
counts <- as.data.frame(colData(sce.annot)) %>% 
  group_by(cell.type, group.sex)  %>% 
  dplyr::summarize(count = n())

counts$patata = paste0(counts$cell.type,"_",counts$group.sex)

barplot.AF <- ggplot(counts, aes(x=patata, y=count, fill = cell.type, label = count)) + 
  geom_bar(stat = "identity")+
  xlab('Grupo') + 
  ylab('Conteos')+
  theme(axis.text.x = 
          element_text(angle = 90, vjust = 1, hjust = 1)) +
  labs(fill = 'Tipo celular')+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  coord_flip() 
ggsave(filename = './figures/7.cell_annotation_consensus/barplot_AF.jpeg', 
       plot = barplot.AF, width = 10, height = 6, units = 'in', dpi = 300,
       device = 'jpeg')
ggsave(filename = './figures/7.cell_annotation_consensus/barplot_AF.svg', 
       plot = barplot.AF, width = 10, height = 6, device = 'svg')


#----------------------- Save data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn52293442_pfc/results/')
saveRDS(sce.annot, './sce.annot.clean.rds')

x = as.data.frame(table(sce.annot$cell.type))
y = as.vector(x$Freq)
names(y) = x$Var1
z = (y*100)/sum(y)

proportions = data.frame(cell.type= names(z),
                         percentage = z,
                         total = y)

write.xlsx(proportions, file = './cell.abundance.xlsx', rowNames = TRUE)




