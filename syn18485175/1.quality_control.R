library(SingleCellExperiment)
library(scater)
library(scran)
library(Matrix)
library(scDblFinder)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(MetBrewer)

#----------------------- Load data ---------------------------------------------

# Load sce
setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/results/')
sce.raw <- readRDS("sce.raw.rds")

coldata <- as.data.frame(colData(sce.raw))
rowdata  <- as.data.frame(rowData(sce.raw))

# Generate palettes
palette <- met.brewer('Renoir')
palette2 <- c(palette[4], palette[9])
palette4 <- c(palette[4], palette[8], palette[9], palette[11])
palette5 <- c(palette[4], palette[6], palette[8], palette[9], palette[11])

#----------------------- Doublet identification -------------------------------

dim(sce.raw)
sce.raw <- sce.raw[Matrix::rowSums(counts(sce.raw) > 0) > 0, ]
# Removed non detected genes
dim(sce.raw) 

# Identify doublets per patient
sce.raw <- scDblFinder::scDblFinder(sce.raw, samples = 'PatientID')

# Library size. Number of genes counts per cells
sce.raw$nCount <- Matrix::colSums(counts(sce.raw))

# Number of genes per cell
sce.raw$nFeatures <- Matrix::colSums(counts(sce.raw)>0)

# Parameters of mitochondrial, ribosomal and ERCC features
is.mito <- grepl("^MT-", rowData(sce.raw)$gene.symbol)
mito.genes <- rowData(sce.raw)$gene.symbol[is.mito]
mito.genes <- rownames(sce.raw)[rowData(sce.raw)$gene.symbol %in% mito.genes]

is.ribo <- grepl("^RP[SL]", rowData(sce.raw)$gene.symbol)
ribo.genes <- rowData(sce.raw)$gene.symbol[is.ribo]
ribo.genes <- rownames(sce.raw)[rowData(sce.raw)$gene.symbol %in% ribo.genes]

is.ercc <- grepl("^ERCC", rowData(sce.raw)$gene.symbol)
ercc.genes <- rowData(sce.raw)$gene.symbol[is.ercc]
ercc.genes <- rownames(sce.raw)[rowData(sce.raw)$gene.symbol %in% ercc.genes]

sce.raw <- scater::addPerCellQC(sce.raw,
                                subsets = list(ERCC = ercc.genes,
                                               RIBO = ribo.genes,
                                               MITO = mito.genes))
length(c(mito.genes, ribo.genes, ercc.genes))


#----------------------- QC plots before filtering------------------------------

# Density plot number for number of genes
jpeg(file = './figures/1.quality_control/nFeatures_density_BF.jpeg', width=10, height=6, units="in", res=300)
plot(density(sce.raw$nFeatures), 
     xlab = 'Número de genes',
     ylab = 'Densidad')
dev.off()
svg('./figures/1.quality_control/nFeatures_density_BF.svg', width = 10, height = 6)
plot(density(sce.raw$nFeatures), 
     xlab = 'Número de genes',
     ylab = 'Densidad')
dev.off()


# Number of genes for each variable
plot.nFeatures.condition.BF <- plotColData(
  sce.raw, x = 'condition', y = 'nFeatures', colour_by = 'condition') +
  xlab('Condición') +
  ylab('Número de genes') +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_colour_manual(values = palette2) +
  scale_x_discrete(labels = c('AD', 'Control'))

plot.nFeatures.sex.BF <- plotColData(
  sce.raw, x = 'sex', y = 'nFeatures', colour_by = 'sex') +
  xlab('Sexo') +
  ylab('Número de genes') +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_colour_manual(values = palette2)

plot.nFeatures.groupsex.BF <- plotColData(
  sce.raw, x = 'group.sex', y = 'nFeatures', colour_by = 'group.sex') +
  xlab('Condición y sexo') +
  ylab('Número de genes') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_colour_manual(values = palette4)

plot.nFeatures.patient.BF <- plotColData(
  sce.raw, x = 'PatientID', y = 'nFeatures', colour_by = 'group.sex') +
  xlab('Paciente') +
  ylab('Número de genes') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_colour_manual(values = palette4, name = 'Condición y sexo',
                      labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

comb.plots.nFeatures.BF <- (plot.nFeatures.condition.BF + plot.nFeatures.sex.BF + plot.nFeatures.groupsex.BF) / plot.nFeatures.patient.BF

ggsave('./figures/1.quality_control/nFeatures_BF.jpeg', plot = comb.plots.nFeatures.BF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/nFeatures_BF.svg', plot = comb.plots.nFeatures.BF, 
       width = 10, height = 6, units = 'in', device = 'svg')


# Density plot for number of counts
jpeg(file = './figures/1.quality_control/nCounts_density_BF.jpeg', width=10, height=6, units="in", res=300)
plot(density(sce.raw$nCount), 
     xlab = 'Tamaño de librería',
     ylab = 'Densidad')
dev.off()
svg('./figures/1.quality_control/nCounts_density_BF.svg', width = 10, height = 6)
plot(density(sce.raw$nCount), 
     xlab = 'Tamaño de librería', 
     ylab = 'Densidad')
dev.off()


# Number of counts for each variable
plot.nCount.condition.BF <- plotColData(
  sce.raw, x = 'condition', y = 'nCount', colour_by = 'condition') +
  xlab('Condición') +
  ylab('Tamaño de librería') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_colour_manual(values = palette2) +
  scale_x_discrete(labels = c('AD', 'Control'))

plot.nCount.sex.BF <- plotColData(
  sce.raw, x = 'sex', y = 'nCount', colour_by = 'sex') +
  xlab('Sexo') +
  ylab('Tamaño de librería') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_colour_manual(values = palette2)

plot.nCount.groupsex.BF <- plotColData(
  sce.raw, x = 'group.sex', y = 'nCount', colour_by = 'group.sex') +
  xlab('Condición y sexo') +
  ylab('Tamaño de librería') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_colour_manual(values = palette4)

plot.nCount.patient.BF <- plotColData(
  sce.raw, x = 'PatientID', y = 'nCount', colour_by = 'group.sex') +
  xlab('Paciente') +
  ylab('Tamaño de librería') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_colour_manual(values = palette4, name = 'Condición y sexo',
                      labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

comb.plots.nCount.BF <- (plot.nCount.condition.BF + plot.nCount.sex.BF + plot.nCount.groupsex.BF) / plot.nCount.patient.BF

ggsave('./figures/1.quality_control/nCounts_BF.jpeg', plot = comb.plots.nCount.BF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/nCounts_BF.svg', plot = comb.plots.nCount.BF, 
       width = 10, height = 6, units = 'in', device = 'svg')


# Percentage of mitochondrial genes for each variable
plot.mito.condition.BF <- plotColData(
  sce.raw, x = 'condition', y = 'subsets_MITO_percent', colour_by = 'condition') +
  xlab('Condición') +
  ylab('Fracción de genes mitocondriales') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2) +
  scale_x_discrete(labels = c('AD', 'Control'))

plot.mito.sex.BF <- plotColData(
  sce.raw, x = 'sex', y = 'subsets_MITO_percent', colour_by = 'sex') +
  xlab('Sexo') +
  ylab('Fracción de genes mitocondriales') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2)

plot.mito.groupsex.BF <- plotColData(
  sce.raw, x = 'group.sex', y = 'subsets_MITO_percent', colour_by = 'group.sex') +
  xlab('Condición y sexo') +
  ylab('Fracción de genes mitocondriales') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4)

plot.mito.patient.BF <- plotColData(
  sce.raw, x = 'PatientID', y = 'subsets_MITO_percent', colour_by = 'group.sex') +
  xlab('Paciente') +
  ylab('Fracción de genes mitocondriales') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4, name = 'Condición y sexo')

comb.plots.mito.BF <- (plot.mito.condition.BF + plot.mito.sex.BF + plot.mito.groupsex.BF) / plot.mito.patient.BF

ggsave('./figures/1.quality_control/mito_BF.jpeg', plot = comb.plots.mito.BF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/mito_BF.svg', plot = comb.plots.mito.BF, 
       width = 10, height = 6, units = 'in', device = 'svg')


# Number of counts vs number of genes for each cell by each variable
plot.nCountVSnFeatures.condition.BF <- plotColData(
  sce.raw, x = 'nCount', y = 'nFeatures', colour_by = 'condition') +
  xlab('Tamaño de librería') +
  ylab('Número de genes') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2,
                     labels = c('AD', 'Control'), 
                     name = "condition") +
  guides(color = guide_legend(override.aes = list(size = 3)))

plot.nCountVSnFeatures.sex.BF <- plotColData(
  sce.raw, x = 'nCount', y = 'nFeatures', colour_by = 'sex') +
  xlab('Tamaño de librería') +
  ylab('Número de genes') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2, 
                     name = "sex") +
  guides(color = guide_legend(override.aes = list(size = 3)))

plot.nCountVSnFeatures.groupsex.BF <- plotColData(
  sce.raw, x = 'nCount', y = 'nFeatures', colour_by = 'group.sex') +
  xlab('Tamaño de librería') +
  ylab('Número de genes') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4,
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'), 
                     name = "Condición y sexo") +
  guides(color = guide_legend(override.aes = list(size = 3)))

plot.nCountVSnFeatures.patient.BF <- plotColData(
  sce.raw, x = 'nCount', y = 'nFeatures', colour_by = 'group.sex') +
  xlab('Tamaño de librería') +
  ylab('Número de genes') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4,
                     name = "Condición y sexo") +
  guides(color = guide_legend(override.aes = list(size = 3)))

comb.plots.nCountVSnFeatures.BF <- (plot.nCountVSnFeatures.condition.BF + plot.nCountVSnFeatures.sex.BF + plot.nCountVSnFeatures.groupsex.BF) / plot.nCountVSnFeatures.patient.BF

ggsave('./figures/1.quality_control/nCountVSnFeatures_BF.jpeg', plot = comb.plots.nCountVSnFeatures.BF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/nCountVSnFeatures_BF.svg', plot = comb.plots.nCountVSnFeatures.BF, 
       width = 10, height = 6, units = 'in', device = 'svg')


# Number of cells
plot.ncells.condition <- ggplot(coldata, aes(x = condition, fill = condition)) +
  geom_bar() +
  scale_fill_manual(values = palette2) +
  labs(
    x = 'Condición',
    y = 'Número de células') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = 'none')

plot.ncells.sex <- ggplot(coldata, aes(x = sex, fill = condition)) +
  geom_bar() +
  scale_fill_manual(values = palette2) +
  labs(
    x = 'Sexo',
    y = 'Número de células') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = 'none')

plot.ncells.groupsex <- ggplot(coldata, aes(x = group.sex, fill = group.sex)) +
  geom_bar() +
  scale_fill_manual(values = palette4) +
  labs(
    x = 'Condición y sexo',
    y = 'Número de células') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = 'none')

plot.ncells.patient <- ggplot(coldata, aes(x = PatientID, fill = group.sex)) +
  geom_bar() +
  scale_fill_manual(values = palette4) +
  labs(
    x = 'Paciente',
    y = 'Número de células') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

comb.plots.ncells.BF <- (plot.ncells.condition + plot.ncells.sex + plot.ncells.groupsex) / plot.ncells.patient

ggsave('./figures/1.quality_control/ncells_BF.jpeg', plot = comb.plots.ncells.BF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/ncells_BF.svg', plot = comb.plots.ncells.BF, 
       width = 10, height = 6, units = 'in', device = 'svg')

#----------------------- Cells QC -----------------------------------------------

sce.qc <- sce.raw
coldata <- as.data.frame(colData(sce.qc))

# Correct by library size, in the paper <200
discard.count <- sce.qc$sum < 200
sce.qc$discard.count <- discard.count
print("Células eliminadas por tamaño librería:")
table(discard.count)

# Correct by mitochondrial reads.
discard.mito <- sce.qc$subsets_MITO_percent > 5
sce.qc$discard.mito <- discard.mito
print("Células eliminadas por proporción mitocondrial:")
table(discard.mito)

# Correct doublets
discard.doublet <- sce.qc@colData$scDblFinder.class == 'doublet'
sce.qc$discard.doublet <- discard.doublet
print("Células eliminadas por dobletes:")
table(discard.doublet)

# Total discarded cells
discard = discard.count | discard.mito | discard.doublet
sce.qc$discard <- discard
coldata <- as.data.frame(colData(sce.qc))
print("Células descartadas:")
table(discard)

#----------------------- Genes QC ---------------------------------------------

pseudogenes <- rowData(sce.raw)$gene.symbol[grep('P$', rowData(sce.raw)$gene.symbol)]
# Processed pseudogenes
pseudogenes1 <- rowData(sce.raw)$gene.symbol[grep("P[0-9]+$",rowData(sce.raw)$gene.symbol)]
#tRNA
trna <- rowData(sce.raw)$gene.symbol[grep('^TR.-', rowData(sce.raw)$gene.symbol)]
# Small nuclear RNA
snrna <- rowData(sce.raw)$gene.symbol[grep('^RNU', rowData(sce.raw)$gene.symbol)]
# Small nucleolar RNA
snorna <- rowData(sce.raw)$gene.symbol[grep("^SNORD",rowData(sce.raw)$gene.symbol)] # SNORD# for “small nucleolar RNA, C/D box” genes
snorna1 <- rowData(sce.raw)$gene.symbol[grep("^SNORA",rowData(sce.raw)$gene.symbol)] # SNORA# for “small nucleolar RNA, H/ACA box” genes
snorna2 <- rowData(sce.raw)$gene.symbol[grep("^SCARNA",rowData(sce.raw)$gene.symbol)] # SCARNA# for “small Cajal body‐specific RNA” genes
#Ribosomal RNA
rrna <- rowData(sce.raw)$gene.symbol[grep('^RNA[0-9]', rowData(sce.raw)$gene.symbol)]
# (pseudo?)genes not described https://www.biostars.org/p/51456/
rprna <- rowData(sce.raw)$gene.symbol[grep("^RP[0-9]",rowData(sce.raw)$gene.symbol)]
# Mitochondrial genes
mito.genes <- rowData(sce.raw)$gene.symbol[grepl("^MT-", rowData(sce.raw)$gene.symbol)]
# Ribosomal genes
ribo.genes <- rowData(sce.raw)$gene.symbol[grepl("^RP[SL]", rowData(sce.raw)$gene.symbol)]
# ERCC genes
ercc.genes <- rowData(sce.raw)$gene.symbol[grepl("^ERCC", rowData(sce.raw)$gene.symbol)]


# Discard genes that are expressed in less than 2 cells, following the paper
print("Genes a descartar por aparecer en menos de 2 celulas")
selected.features <-  nexprs(sce.raw,byrow = TRUE,detection_limit = 1) > 2
table(selected.features) # FALSE is discarded

genes.discard <- c(pseudogenes, pseudogenes1, trna, snrna, snorna, snorna1, 
                   snorna2, rrna, rprna, mito.genes, ribo.genes, ercc.genes)
print('Genes a descartar:')
length(genes.discard)

not.selected.pseudogenes <- rowData(sce.raw)$gene.symbol %in% genes.discard
print('Genes a descartar finales:')
# Discard the features from genes.discard and the false ones in selected.features
table(not.selected.pseudogenes | selected.features == FALSE) 

nkeep <- !not.selected.pseudogenes & selected.features

#----------------------- Final dataset -----------------------------------------

# We keep the rows/features from nkeep and the columns/cells that are NOT discard
sce.filt <- sce.qc[nkeep, !discard]
dim(sce.qc)
dim(sce.filt)
setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/results')
saveRDS(sce.filt, file = "./sce.filt.rds")

#----------------------- QC plots after filtering------------------------------

coldata.filt <- as.data.frame(colData(sce.filt))
coldata$discarded <- !(rownames(coldata) %in% colnames(sce.filt))
table(coldata$discarded)


# Number of cells
plot.ncells.condition.AF <- ggplot(coldata.filt, aes(x = condition, fill = condition)) +
  geom_bar() +
  scale_fill_manual(values = palette2) +
  labs(
    x = 'Condición',
    y = 'Número de células',
    fill = 'Condición') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

plot.ncells.sex.AF <- ggplot(coldata.filt, aes(x = sex, fill = sex)) +
  geom_bar() +
  scale_fill_manual(values = palette2) +
  labs(
    x = 'Sexo',
    y = 'Número de células',
    fill = 'Sexo') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

plot.ncells.groupsex.AF <- ggplot(coldata.filt, aes(x = group.sex, fill = group.sex)) +
  geom_bar() +
  scale_fill_manual(values = palette4) +
  labs(
    x = 'Condición y sexo',
    y = 'Número de células',
    fill = 'Condición y sexo') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    legend.position = 'none')

plot.ncells.patient.AF <- ggplot(coldata.filt, aes(x = PatientID, fill = group.sex)) +
  geom_bar() +
  scale_fill_manual(values = palette4) +
  labs(
    x = 'Paciente',
    y = 'Número de células',
    fill = 'Condición y sexo') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


comb.plots.ncells.AF <- (plot.ncells.condition.AF + plot.ncells.sex.AF + plot.ncells.groupsex.AF) / plot.ncells.patient.AF

ggsave('./figures/1.quality_control/ncells_AF.jpeg', plot = comb.plots.ncells.AF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/ncells_AF.svg', plot = comb.plots.ncells.AF, 
       width = 10, height = 6, units = 'in', device = 'svg')

# Density plot number for number of genes
jpeg('./figures/1.quality_control/nFeatures_density_AF.jpeg', 
     width = 10, height = 6, units = "in", res = 300)
plot(density(sce.filt$nFeatures), 
     xlab = 'Número de genes', 
     ylab = 'Densidad')
dev.off()

svg('./figures/1.quality_control/nFeatures_density_AF.svg', 
    width = 10, height = 6)
plot(density(sce.filt$nFeatures), 
     xlab = 'Número de genes', 
     ylab = 'Densidad')
dev.off()

# Number of genes for each variable
plot.nFeatures.condition.AF <- plotColData(
  sce.filt, x = 'condition', y = 'nFeatures', colour_by = 'condition') +
  xlab('Condición') +
  ylab('Número de genes') +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2) +
  scale_x_discrete(labels = c('AD', 'Control'))

plot.nFeatures.sex.AF <- plotColData(
  sce.filt, x = 'sex', y = 'nFeatures', colour_by = 'sex') +
  xlab('Sexo') +
  ylab('Número de genes') +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2)

plot.nFeatures.groupsex.AF <- plotColData(
  sce.filt, x = 'group.sex', y = 'nFeatures', colour_by = 'group.sex') +
  xlab('Condición y sexo') +
  ylab('Número de genes') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4)

plot.nFeatures.patient.AF <- plotColData(
  sce.filt, x = 'PatientID', y = 'nFeatures', colour_by = 'group.sex') +
  xlab('Paciente') +
  ylab('Número de genes') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4, name = 'Condición y sexo')

comb.plots.nFeatures.AF <- (plot.nFeatures.condition.AF + plot.nFeatures.sex.AF + plot.nFeatures.groupsex.AF) / plot.nFeatures.patient.AF

ggsave('./figures/1.quality_control/nFeatures_AF.jpeg', plot = comb.plots.nFeatures.AF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/nFeatures_AF.svg', plot = comb.plots.nFeatures.AF, 
       width = 10, height = 6, units = 'in', device = 'svg')

# Density plot for number of counts
jpeg('./figures/1.quality_control/nCounts_density_AF.jpeg', 
     width = 10, height = 6, units = "in", res = 300)
plot(density(sce.filt$nCount), 
     xlab = 'Tamaño de librería', 
     ylab = 'Densidad')
dev.off()

svg('./figures/1.quality_control/nCounts_density_AF.svg', 
    width = 10, height = 6)
plot(density(sce.filt$nCount), 
     xlab = 'Tamaño de librería', 
     ylab = 'Densidad')
dev.off()


# Number of counts for each variable
plot.nCount.condition.AF <- plotColData(
  sce.filt, x = 'condition', y = 'nCount', colour_by = 'condition') +
  xlab('Condición') +
  ylab('Tamaño de librería') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2) +
  scale_x_discrete(labels = c('AD', 'Control'))

plot.nCount.sex.AF <- plotColData(
  sce.filt, x = 'sex', y = 'nCount', colour_by = 'sex') +
  xlab('Sexo') +
  ylab('Tamaño de librería') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2)

plot.nCount.groupsex.AF <- plotColData(
  sce.filt, x = 'group.sex', y = 'nCount', colour_by = 'group.sex') +
  xlab('Condición y sexo') +
  ylab('Tamaño de librería') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4)

plot.nCount.patient.AF <- plotColData(
  sce.filt, x = 'PatientID', y = 'nCount', colour_by = 'group.sex') +
  xlab('Paciente') +
  ylab('Tamaño de librería') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4, name = 'Condición y sexo')

comb.plots.nCount.AF <- (plot.nCount.condition.AF + plot.nCount.sex.AF + plot.nCount.groupsex.AF) / plot.nCount.patient.AF

ggsave('./figures/1.quality_control/nCounts_AF.jpeg', plot = comb.plots.nCount.AF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/nCounts_AF.svg', plot = comb.plots.nCount.AF, 
       width = 10, height = 6, units = 'in', device = 'svg')

# Percentage of mitochondrial genes for each variable
plot.mito.condition.AF <- plotColData(
  sce.filt, x = 'condition', y = 'subsets_MITO_percent', colour_by = 'condition') +
  xlab('Condición') +
  ylab('Fracción de genes mitocondriales') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2) +
  scale_x_discrete(labels = c('AD', 'Control'))

plot.mito.sex.AF <- plotColData(
  sce.filt, x = 'sex', y = 'subsets_MITO_percent', colour_by = 'sex') +
  xlab('Sexo') +
  ylab('Fracción de genes mitocondriales') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2)

plot.mito.groupsex.AF <- plotColData(
  sce.filt, x = 'group.sex', y = 'subsets_MITO_percent', colour_by = 'group.sex') +
  xlab('Condición y sexo') +
  ylab('Fracción de genes mitocondriales') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4)

plot.mito.patient.AF <- plotColData(
  sce.filt, x = 'PatientID', y = 'subsets_MITO_percent', colour_by = 'group.sex') +
  xlab('Paciente') +
  ylab('Fracción de genes mitocondriales') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4, name = 'Condición y sexo')

comb.plots.mito.AF <- (plot.mito.condition.AF + plot.mito.sex.AF + plot.mito.groupsex.AF) / plot.mito.patient.AF

ggsave('./figures/1.quality_control/mito_AF.jpeg', plot = comb.plots.mito.AF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/mito_AF.svg', plot = comb.plots.mito.AF, 
       width = 10, height = 6, units = 'in', device = 'svg')


# Number of counts vs number of genes for each cell by each variable

plot.nCountVSnFeatures.condition.AF <- plotColData(
  sce.filt, x = 'nCount', y = 'nFeatures', colour_by = 'condition') +
  xlab('Tamaño de librería') +
  ylab('Número de genes') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2,
                     labels = c('AD', 'Control'), 
                     name = "Condición") +
  guides(color = guide_legend(override.aes = list(size = 3)))

plot.nCountVSnFeatures.sex.AF <- plotColData(
  sce.filt, x = 'nCount', y = 'nFeatures', colour_by = 'sex') +
  xlab('Tamaño de librería') +
  ylab('Número de genes') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette2,
                     name = "Sexo") +
  guides(color = guide_legend(override.aes = list(size = 3)))

plot.nCountVSnFeatures.groupsex.AF <- plotColData(
  sce.filt, x = 'nCount', y = 'nFeatures', colour_by = 'group.sex') +
  xlab('Tamaño de librería') +
  ylab('Número de genes') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4,
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'), 
                     name = "Condición y sexo") +
  guides(color = guide_legend(override.aes = list(size = 3)))

plot.nCountVSnFeatures.patient.AF <- plotColData(
  sce.filt, x = 'nCount', y = 'nFeatures', colour_by = 'group.sex') +
  xlab('Tamaño de librería') +
  ylab('Número de genes') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = palette4,
                     name = "Condición y sexo") +
  guides(color = guide_legend(override.aes = list(size = 3)))

comb.plots.nCountVSnFeatures.AF <- (plot.nCountVSnFeatures.condition.AF + plot.nCountVSnFeatures.sex.AF + plot.nCountVSnFeatures.groupsex.AF) / plot.nCountVSnFeatures.patient.AF

ggsave('./figures/1.quality_control/nCountVSnFeatures_AF.jpeg', plot = comb.plots.nCountVSnFeatures.AF, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/nCountVSnFeatures_AF.svg', plot = comb.plots.nCountVSnFeatures.AF, 
       width = 10, height = 6, units = 'in', device = 'svg')


# Discard plot, number of cells
discard.plot <- ggplot(coldata, aes(x = factor(PatientID), fill = discard)) +
  geom_bar(position = position_stack(reverse = TRUE))+
  scale_fill_manual(values = c("FALSE" = 'skyblue', "TRUE" = "tomato3"), name = 'Descartadas') +
  labs(
    x = 'Paciente',
    y = 'Número de células') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave('./figures/1.quality_control/discard_AF.jpeg', plot = discard.plot, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/discard_AF.svg', plot = discard.plot, 
       width = 10, height = 6, units = 'in', device = 'svg')

# Discard plot with colors depending on the discarding reason
coldata_plot <- coldata %>%
  mutate(
    status = case_when(
      discard == FALSE ~ "Conservado",
      # Discarded cells
      discard.mito == TRUE ~ "Descartado por % mitocondrial",
      discard.count == TRUE ~ "Descartado por tamaño de librería",
      discard.doublet == TRUE ~ "Descartado por doblete"
    ),
    status = factor(status, levels = c(
      "Conservado",
      "Descartado por % mitocondrial",
      "Descartado por tamaño de librería",
      "Descartado por doblete"
    ))
  )

discard.plot.colors <- ggplot(coldata_plot, aes(x = factor(PatientID), fill = status)) +
  geom_bar(position = position_stack(reverse =TRUE), color = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "Conservado" = palette5[1],  
      "Descartado por % mitocondrial" = palette5[2],
      "Descartado por tamaño de librería" = palette5[4],
      "Descartado por doblete" = palette5[5]
    ),
    name = "Estado"
  ) +
  labs(
    x = "Paciente",
    y = "Número de células"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
discard.plot.colors
ggsave('./figures/1.quality_control/discard_AF_colors.jpeg', plot = discard.plot.colors, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./figures/1.quality_control/discard_AF_colors.svg', plot = discard.plot.colors, 
       width = 10, height = 6, units = 'in', device = 'svg')
