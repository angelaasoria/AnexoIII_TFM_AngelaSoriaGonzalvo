# Classification of cell types. Three tools are used:
# SingleR: Allows cell-level detection and is based on reference datasets.
#         Does not allow marking cells as “Unknown”.
#         https://bioconductor.org/packages/release/bioc/html/SingleR.html
# CHETAH: Allows detection at the cellular level and is based on reference datasets.
#         Does allow marking cells as “Unknown”.
#         https://www.bioconductor.org/packages/release/bioc/html/CHETAH.html
# SCINA: Enables detection at the cellular level and is based on marker genes.
#         Yes allows marking cells as “Unknown”.
#         https://github.com/jcao89757/SCINA

library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(BiocParallel)
library(scDblFinder)
library(AnnotationHub)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(broom)
library(kableExtra)
library(HDF5Array)
library(DropletUtils)
library(devtools)
library(celda)
library(bluster)
library(ggpubr)
library(SCINA)
library(SingleR)
library(CHETAH)
library(openxlsx)
library(Seurat)
library(HGNChelper)
library(cowplot)
library(MetBrewer)

#----------------------- Load data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
sce.norm <- readRDS('sce.norm.rds')
sce.cluster <- readRDS('sce.cluster.rds')

#random.rows <- sample(nrow(counts(sce.norm)), size = 100, replace = FALSE)
#sce.norm <- sce.norm[random.rows, ]
sce.matrix <- assay(sce.norm, 'logcounts')
rownames(sce.matrix) <- rowData(sce.norm)$gene.symbol

# Load reference datasets
sce.downsample <- readRDS('/clinicfs/projects/i63/tfm_hipocampo/notacion_celular/human_atlas_dataset/sce_allan_red_downsample.rds')
sce.downsample$tipos_celulares <- recode(
  sce.downsample$tipos_celulares,
  "Astro" = "Astrocitos",
  "Exc" = "Excitadoras",
  "Inh" = "Inhibidoras",
  "Micro" = "Microglía",
  "Oligo" = "Oligodendrocitos",
  "OPC" = "OPC"
)

# Generate palettes
palette <- met.brewer('Renoir', n = 11)
palette[10] <- palette[11]
palette7 <- c(palette[4], palette[5], palette[6], palette[8], palette[9],
              palette[11], palette[12])

#----------------------- Load marker genes -------------------------------------

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# CELLMARKER
cellmarker.db.exp <- read.csv(file = '/clinicfs/projects/i63/tfm_hipocampo/notacion_celular/genes_marcadores/cellmarker_experiment.csv')
cellmarker.db <- cellmarker.db.exp

cellmarker.db <- cellmarker.db[cellmarker.db$Cancer == 'Normal cell', ]
cellmarker.db <- cellmarker.db[cellmarker.db$Tissue.Type %in% 
                                 c('Cortex', 'Prefrontal cortex', 'Brain'), ]
# Create columns with each cell type and their markers
cellmarker.brain <- cellmarker.db %>%
  group_by(`Cell.name`) %>%
  summarize(geneset = list(`Cell.marker`))

##### OPTION A #####
gm.cellmarker <- setNames(cellmarker.brain$geneset, cellmarker.brain$Cell.name)

gm.cellmarker.mc <- list(
  "Astrocitos" = c(
    gm.cellmarker[["Astrocyte"]],
    gm.cellmarker[["Fibrous astrocyte"]],
    gm.cellmarker[["Mature astrocyte"]]
  ),
  "Oligodendrocitos" = c(
    gm.cellmarker[["Oligodendrocyte"]],
    gm.cellmarker[["Mature oligodendrocyte"]],
    gm.cellmarker[["Oligodendrocyte-like cell"]]
  ),
  "Microglía" = c(
    gm.cellmarker[["Microglial cell"]],
    gm.cellmarker[["Homeostatic microglial cell"]],
    gm.cellmarker[["Activated microglial cell"]],
    gm.cellmarker[["Disease-associated microglial cell"]],
    gm.cellmarker[["Demyelinating microglial cell"]],
    gm.cellmarker[["Microglia-like cell"]]
  ),
  "OPCs" = c(
    gm.cellmarker[["Oligodendrocyte precursor cell"]],
    gm.cellmarker[["Oligodendrocyte progenitor cell"]],
    gm.cellmarker[["Oligodendrocyte-like cell"]]
  ),
  "Neuronas excitadoras" = c(
    gm.cellmarker[["Excitatory neuron"]],
    gm.cellmarker[["Glutamatergic neuron"]],
    gm.cellmarker[["Deep layer neuron"]],
    gm.cellmarker[["Upper layer cortical neuron"]]
  ),
  "Neuronas inhibidoras" = c(
    gm.cellmarker[["Inhibitory neuron"]],
    gm.cellmarker[["GABAergic neuron"]],
    gm.cellmarker[["Interneuron"]],
    gm.cellmarker[["Cholinergic neuron"]],
    gm.cellmarker[["Dopaminergic neuron"]],
    gm.cellmarker[["Serotonergic neuron"]]
  )
)


##### OPTION B: SCTYPE #####
sctype.db <- gene_sets_prepare('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx',
                               'Brain')[[1]]
gm.sctype <- sctype.db

gm.sctype.mc <- list(
  "Astrocitos" = gm.sctype[["Astrocytes"]],
  "Oligodendrocitos" = c(gm.sctype[["Oligodendrocytes"]], gm.sctype[["Myelinating Schwann cells"]]),
  "Microglia" = gm.sctype[["Microglial cells"]],
  "OPCs" = gm.sctype[["Oligodendrocyte precursor cells"]],
  "Neuronas excitadoras" = gm.sctype[["Glutamatergic neurons"]],
  "Neuronas inhibidoras" = c(
    gm.sctype[["GABAergic neurons"]],
    gm.sctype[["Cholinergic neurons"]],
    gm.sctype[["Dopaminergic neurons"]],
    gm.sctype[["Serotonergic neurons"]]
  )
)

##### OPTION C: MIX #####
intersection.list <- lapply(names(gm.sctype.mc), function(cell_type){
  intersect(gm.sctype.mc[[cell_type]], gm.sctype.mc[[cell_type]])
})
names(intersection.list) <- names(gm.sctype.mc)
names(intersection.list) <- c('Astrocitos', 'Oligodendrocitos', 'Microglía', 'OPC',
                              'Excitadoras', 'Inhibidoras')
excitadoras <- unlist(c(intersection.list['Excitadoras']))
excitadoras <- intersect(excitadoras, rowData(sce.norm)$gene.symbol)

inhibidoras <- unlist(c(intersection.list['Inhibidoras']))
inhibidoras <- intersect(inhibidoras, rowData(sce.norm)$gene.symbol)
inhibidoras

astrocitos <- unlist(c(intersection.list['Astrocitos']))
astrocitos <- intersect(astrocitos, rowData(sce.norm)$gene.symbol)

oligodendrocitos <- unlist(c(intersection.list['Oligodendrocitos']))
oligodendrocitos <- intersect(oligodendrocitos, rowData(sce.norm)$gene.symbol)

opc <- unlist(c(intersection.list['OPC']))
opc <- intersect(opc, rowData(sce.norm)$gene.symbol)

microglia <- unlist(c(intersection.list['Microglía']))
microglia <- intersect(microglia, rowData(sce.norm)$gene.symbol)

gm.manual.2 <- list(excitadoras, inhibidoras, astrocitos, oligodendrocitos, opc, microglia)
names(gm.manual.2) <- c('Excitadoras', 'Inhibidoras', 'Astrocitos', 'Oligodendrocitos', 'OPC', 'Microglía')

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
saveRDS(gm.manual.2, './gm.manual.2.rds')

##### OPTION D: MANUAL #####
gm.manual <- readRDS('/clinicfs/projects/i63/tfm_hipocampo/notacion_celular/genes_marcadores/genes_marcadores.rds')

excitadoras <- unlist(c(gm.manual['Ex']))
excitadoras <- intersect(excitadoras, rowData(sce.norm)$gene.symbol)

inhibidoras <- unlist(c(gm.manual['Inh']))
inhibidoras <- intersect(inhibidoras, rowData(sce.norm)$gene.symbol)

astrocitos <- unlist(c(gm.manual['Astrocitos'], gm.sctype['Astrocytes']))
astrocitos <- intersect(astrocitos, rowData(sce.norm)$gene.symbol)

oligodendrocitos <- unlist(c(gm.manual['Oligodendrocitos'], gm.sctype['Oligodendrocytes']))
oligodendrocitos <- intersect(oligodendrocitos, rowData(sce.norm)$gene.symbol)

opc <- unlist(c(gm.manual['OPC']))
opc <- intersect(opc, rowData(sce.norm)$gene.symbol)

microglia <- unlist(c(gm.manual['Microglia'], gm.sctype['Microglial cells']))
microglia <- intersect(microglia, rowData(sce.norm)$gene.symbol)

endotelio <- unlist(c(gm.manual['Endotelio'], gm.sctype['Endothelial cells']))
endotelio <- intersect(endotelio, rowData(sce.norm)$gene.symbol)

fibroblasto <- unlist(c(gm.manual['Fibroblastos']))
fibroblasto <- intersect(fibroblasto, rowData(sce.norm)$gene.symbol)

gm.manual.mix.complete <- list(excitadoras, inhibidoras, astrocitos, oligodendrocitos, 
                               opc, microglia, endotelio, fibroblasto)
names(gm.manual.mix.complete) <- c('Excitadoras', 'Inhibidoras', 'Astrocitos', 
                                   'Oligodendrocitos', 'OPC', 'Microglia', 
                                   'Endotelio', 'Fibroblasto')

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
saveRDS(gm.manual.mix.complete, 'gm.manual.mix.rds')

gm.manual.3 <- lapply(names(gm.manual.2), function(cell_type) {
  c(gm.manual.2[[cell_type]], gm.manual.mix.complete[[cell_type]])
})
names(gm.manual.3) = c("Excitadoras", "Inhibidoras", "Astrocitos", "Oligodendrocitos","OPC","Microglía")
setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
saveRDS(gm.manual.3, './gm.manual.3.rds')

# Chosen clustering
sce.norm$cluster40 <- sce.cluster$cluster40


#----------------------- Cell type identification SCINA ------------------------
#--------------------------- Cell level, Unknown: YES --------------------------
# With SCINA we will analyse using the marker genes obtained from ‘The human 
# protein atlas’, which are ‘gm.manual’. In addition we are going to use a second
# set of marker genes called ‘gm.manual.combi’ which are the ones obtained by
# combining two extra databases: CellMarker and sc-type.

# In addition, with regard to the marker genes, as the ‘Unknow’ option is 
# available, we will use only the major cell types of interest, discarding cell 
# types such as ‘Neuronal death’ or ‘endothelium’.

result.scina <- SCINA::SCINA(exp = sce.matrix,
                             signatures = gm.manual.3,
                             allow_unknown = TRUE)

colData(sce.norm)$SCINA.labels <- result.scina$cell_labels
reducedDim(sce.norm, "TSNE.seed.100.per.90") <- reducedDim(sce.cluster, "TSNE.seed.100.per.90")
reducedDim(sce.norm, "UMAP.seed.100.n300") <- reducedDim(sce.cluster, "UMAP.seed.100.n300")
reducedDim(sce.norm, "PCA") <- reducedDim(sce.cluster, "PCA")

# Plot
PCA.scina <- plotReducedDim(sce.norm, dimred="PCA",colour_by="SCINA.labels", text_by="cluster40")+
  xlab("PCA1")+
  ylab("PCA2")+
  scale_color_manual(values = palette7, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/PCA_scina.jpeg', width=10, height=6, units="in", res=300)
print(PCA.scina)
dev.off()
svg(file = './figures/6.cell_annotation/PCA_scina.svg', width=10, height=6)
print(PCA.scina)
dev.off()

TSNE.scina <- plotReducedDim(sce.norm, dimred="TSNE.seed.100.per.90",colour_by="SCINA.labels", text_by="cluster40")+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette7, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/TSNE_scina.jpeg', width=10, height=6, units="in", res=300)
print(TSNE.scina)
dev.off()
svg(file = './figures/6.cell_annotation/TSNE_scina.svg', width=10, height=6)
print(TSNE.scina)
dev.off()

UMAP.scina <- plotReducedDim(sce.norm, dimred="UMAP.seed.100.n300",colour_by="SCINA.labels", text_by="cluster40")+
  xlab("UMAP1")+
  ylab("UMAP2")+
  scale_color_manual(values = palette7, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/UMAP_scina.jpeg', width=10, height=6, units="in", res=300)
print(UMAP.scina)
dev.off()
svg(file = './figures/6.cell_annotation/UMAP_scina.svg', width=10, height=6)
print(UMAP.scina)
dev.off()

#----------------------- Cell type identification SINGLER ----------------------
#--------------------------- Cell level, Unknown: NO --------------------------
# With SINGLER we perform cell type classification using a reference dataset 
# (sce.downsample) with predefined cell types ('tipos_celulares').

result.singler <- SingleR::SingleR(test = sce.matrix,
                                   ref = sce.downsample,
                                   labels = colData(sce.downsample)$tipos_celulares)
colData(sce.norm)$SingleR.labels <- result.singler$labels

# Plot
pca.singler <- plotReducedDim(sce.norm, dimred="PCA",colour_by="SingleR.labels", 
                              text_by="cluster40")+
  xlab("PCA1")+
  ylab("PCA2")+
  scale_color_manual(values = palette7, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/PCA_SingleR.jpeg', width=10, height=6, units="in", res=300)
print(pca.singler)
dev.off()
svg(file = './figures/6.cell_annotation/PCA_SingleR.svg', width=10, height=6)
print(pca.singler)
dev.off()

tsne.singler <- plotReducedDim(sce.norm, dimred="TSNE.seed.100.per.90",
                               colour_by="SingleR.labels", text_by="cluster40")+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette7, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/tSNE_SingleR.jpeg', width=10, height=6, units="in", res=300)
print(tsne.singler)
dev.off()
svg(file = './figures/6.cell_annotation/tSNE_SingleR.svg', width=10, height=6)
print(tsne.singler)
dev.off()

umap.singler <- plotReducedDim(sce.norm, dimred="UMAP.seed.100.n300",
                               colour_by="SingleR.labels", text_by="cluster40")+
  xlab("UMAP1")+
  ylab("UMAP2")+
  scale_color_manual(values = palette7, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/UMAP_SingleR.jpeg', width=10, height=6, units="in", res=300)
print(umap.singler)
dev.off()
svg(file = './figures/6.cell_annotation/UMAP_SingleR.svg', width=10, height=6)
print(umap.singler)
dev.off()

#----------------------- Cell type identification CHETAH ----------------------
#--------------------------- Cell level, Unknown: YES --------------------------
# We perform cell type classification using a reference dataset (sce.downsample) 
# with predefined cell types ('tipos_celulares'). CHETAH compares gene expression 
# profiles between our data and the reference to assign cell identities.

sce.chetah <- sce.norm
rownames(sce.chetah) <- rowData(sce.norm)$gene.symbol
valid.genes <- !(is.na(rownames(sce.chetah)) | rownames(sce.chetah) == "")
sce.chetah <- sce.chetah[valid.genes, ]

result.chetah <- CHETAH::CHETAHclassifier(input = sce.chetah,
                                          ref_cells = sce.downsample,
                                          ref_ct = 'tipos_celulares')

colData(sce.norm)$CHETAH.labels <- result.chetah$celltype_CHETAH

# Plot
pca.chetah <- plotReducedDim(sce.norm, dimred="PCA",colour_by="CHETAH.labels", 
                             text_by="cluster40")+
  xlab("PCA1")+
  ylab("PCA2")+
  scale_color_manual(values = palette, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/PCA_CHETAH.jpeg', width=10, height=6, units="in", res=300)
print(pca.chetah)
dev.off()
svg(file = './figures/6.cell_annotation/PCA_CHETAH.svg', width=10, height=6)
print(pca.chetah)
dev.off()


tsne.chetah <- plotReducedDim(sce.norm, dimred="TSNE.seed.100.per.90",
                              colour_by="CHETAH.labels", text_by="cluster40")+
  xlab("tSNE1")+
  ylab("tSNE2")+
  scale_color_manual(values = palette, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/TSNE_CHETAH.jpeg', width=10, height=6, units="in", res=300)
print(tsne.chetah)
dev.off()
svg(file = './figures/6.cell_annotation/TSNE_CHETAH.svg', width=10, height=6)
print(tsne.chetah)
dev.off()


umap.chetah <- plotReducedDim(sce.norm, dimred="UMAP.seed.100.n300",
                              colour_by="CHETAH.labels", text_by="cluster40")+
  xlab("UMAP1")+
  ylab("UMAP2")+
  scale_color_manual(values = palette, name = "Tipo celular")
jpeg(file = './figures/6.cell_annotation/UMAP_CHETAH.jpeg', width=10, height=6, units="in", res=300)
print(umap.chetah)
dev.off()
svg(file = './figures/6.cell_annotation/UMAP_CHETAH.svg', width=10, height=6)
print(umap.chetah)
dev.off()

#----------------------- Save data ---------------------------------------------

saveRDS(sce.norm, 'sce.annot.rds')
