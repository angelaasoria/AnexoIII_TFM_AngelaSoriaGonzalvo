library(readr)
library(data.table)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(biomaRt)
library(stringr)

#----------------------- Load data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/data')

# Función para leer y formatear datos
load_counts <- function(expr_file, barcode_file, features_file){
  # Leer matriz de conteo
  counts <- readMM(expr_file)
  # Leer nombres de genes (features)
  features <- read.table(features_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  # Leer metadatos de barcodes (columnas)
  barcodes <- read.table(barcode_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Asignar nombres a las filas y columnas
  rownames(counts) <- features$V1
  colnames(counts) <- barcodes$TAG
  
  return(counts)
}

counts_total <- load_counts('notfiltered_count_matrix.mtx',
                            'notfiltered_column_metadata.txt',
                            'notfiltered_gene_row_names.txt')


#----------------------- ColData -----------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/metadata')

#metadata <- read.table('./metadatos.finales.txt')

sce <- SingleCellExperiment(assays = list(counts = counts_total))

# Generar colData con 2 columnas: barcodes y patient id
colData(sce)$PatientID <- sub(".*\\.", "", colnames(sce))
colData(sce)$barcode <- sub("\\..*", "", colnames(sce))
coldata <- as.data.frame(colData(sce))

## Código borja generar_metadatos
setwd("/clinicfs/projects/i63/tfm_hipocampo/syn18485175/")

metadatos = read.table("metadata/snRNAseqPFC_BA10_biospecimen_metadata.csv", sep = ",", header = T)
metadatos2 = read.table("metadata/snRNAseqPFC_BA10_assay_scRNAseq_metadata.csv", sep = ",", header = T)
metadatos3 <- read.table("metadata/snRNAseqPFC_BA10_id_mapping.csv", sep = ",", header = T)
metadatos4 <- read.table("metadata/snRNAseqPFC_BA10_Sample_key.csv", sep = ",", header = T)
rosmap_metadata = read.table("metadata/ROSMAP_clinical.csv", sep = ",", header = T)
rosmap_metadata <- rosmap_metadata[rosmap_metadata$projid %in% metadatos4$projid,]

metadatos <-  merge(metadatos, rosmap_metadata, 
                    by.x = "projid", 
                    by.y = "projid")

metadatos_pacientes <- readxl::read_xlsx("metadata/41586_2019_1195_MOESM3_ESM.xlsx", sheet = 2)
metadatos3 <- unique(metadatos3[,c(2,3)])

metadatos_pacientes <-  merge(metadatos_pacientes, metadatos3, 
                              by.x = "Subject", 
                              by.y = "Subject")

metadatos_pacientes <- metadatos_pacientes[ ,c(1,2,6,7,8,10,11)]
metadatos <-  merge(metadatos, metadatos_pacientes, 
                    by.x = "projid", 
                    by.y = "projid")
metadatos$PatientID <- sub(".*ROS", "", metadatos$Subject)

write.table(metadatos, file = "metadata/metadatos.finales.txt")

# Create coldata including metadata
coldata <- merge(coldata, metadatos,
                 by.x = 'PatientID',
                 by.y = 'PatientID',
                 all.x = TRUE)

# Create 'sex' column in coldata
coldata <- coldata %>%
  mutate(sex = recode(msex,
                       '0' = 'Mujer',
                       '1' = 'Hombre'))

# Create 'condition' column
coldata <- coldata %>%
  mutate(condition = recode(`pathologic diagnosis of AD`,
                       'NO' = 'Control',
                       'YES' = 'AD'))

# Creo columna de covariable
coldata <- coldata %>%
  mutate(group.sex = paste(condition, sex, sep = '_'))

colData(sce) <- DataFrame(coldata)

#----------------------- RowData -----------------------------------------------

# Sacar gene symbol
rownames(sce)
mart <- biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
gene.info <-biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = rownames(sce),
                           mart = mart)

row.annot <- DataFrame(ensembl.id = rownames(sce),
                       gene.symbol = gene.info$hgnc_symbol[match(rownames(sce), gene.info$ensembl_gene_id)],
                       row.names = rownames(sce))
rowData(sce) <- DataFrame(row.annot)

#----------------------- Save data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/results')
saveRDS(sce, "sce.raw.rds")
