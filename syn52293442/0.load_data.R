library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(S4Vectors)
library(Matrix)

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn52293442_pfc')

#----------------------- Crate SingleCellExperiment ----------------------------

cortex <- readRDS('data/Prefrontal_cortex.rds')
class(cortex)
sce <- as.SingleCellExperiment(cortex)
coldata <- as.data.frame(colData(sce))
dim(sce)

#----------------------- Generate metadata -------------------------------------

metadata <- read.table("metadata/Supplementary_Table_1_sample_metadata.txt", sep = "\t", header = T)
metadata <- metadata[metadata$region == "PFC",]

metadata2 <- read.table("metadata/MIT_ROSMAP_Multiomics_assay_multiome_metadata.csv", sep = ",", header = T)
metadata3 <- read.table("metadata/MIT_ROSMAP_Multiomics_assay_RNAseq_metadata.csv", sep = ",", header = T)
metadata4 <- read.table("metadata/MIT_ROSMAP_Multiomics_assay_snATACseq_metadata.csv", sep = ",", header = T)
metadata5 <- read.table("metadata/MIT_ROSMAP_Multiomics_assay_snRNAseq_metadata.csv", sep = ",", header = T)
metadata6 <- read.table("metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv", sep = ",", header = T)
metadata7 <- read.table("metadata/MIT_ROSMAP_Multiomics_individual_metadata.csv", sep = ",", header = T)
rosmap_metadata <- read.table("metadata/ROSMAP_clinical.csv", sep = ",", header = T)

metadata <- merge(metadata, unique(metadata7[ ,c(1,19)]),
                  by.x = 'subject', by.y = 'subject', all.x = TRUE)

metadata <- merge(metadata, rosmap_metadata, 
                  by.x = 'individualID', by.y = 'individualID')

write.table(metadata, file = "metadata/final_metadata.txt")

#----------------------- colData -----------------------------------------------

coldata <- as.data.frame(colData(sce))
coldata <- merge(coldata, metadata,
                 by.x = 'projid', by.y = 'projid')

# Create 'sex' column in coldata
coldata <- coldata %>%
  mutate(sex = recode(msex.x,
                      '0' = 'Mujer',
                      '1' = 'Hombre'))

# Create 'condition' column
coldata <- coldata %>%
  mutate(condition = recode(pathAD,
                            'non-AD' = 'Control',
                            'AD' = 'AD'))
# Creo columna de covariable
coldata <- coldata %>%
  mutate(group.sex = paste(condition, sex, sep = '_'))

colData(sce) <- DataFrame(coldata)

#----------------------- RowData -----------------------------------------------
sce <- sce.raw
# In this case, the original rowdata is in genesymbol, so we get ensembl id

rownames(sce)
mart <- biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
gene.info <-biomaRt::getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
                           filters = 'hgnc_symbol',
                           values = rownames(sce),
                           mart = mart)

row.annot <- DataFrame(gene.symbol = rownames(sce),
                       ensembl.id = gene.info$ensembl_gene_id[match(rownames(sce), gene.info$hgnc_symbol)],
                       row.names = rownames(sce))
rowData(sce) <- DataFrame(row.annot)

#----------------------- Save data ---------------------------------------------

saveRDS(sce, 'results/sce.raw.rds')

