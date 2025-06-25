library(readr)
library(data.table)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(biomaRt)
library(stringr)

#----------------------- Load data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/data')

# Función para leer y formatear datos
load_counts <- function(expr_file, barcode_file, features_file, sample_name) {
  # Leer datos de expresión
  counts <- readMM(expr_file)
  
  # Leer barcodes
  barcodes <- data.table::fread(barcode_file, header = FALSE)$V1
  barcodes <- paste0(barcodes, "_", sample_name)
  
  features <- data.table::fread(features_file, header = FALSE)$V1
  
  rownames(counts) <- features
  colnames(counts) <- barcodes
  
  return(counts)
}

samples <- c(paste0('AD', c(1,2,3,5,7,8,9,10,11,12,13)), 
             paste0('C', c(1,2,3,4,5,6,7,8,9,11,12)), 
             paste0('P', c(2,3,5,6,7,9,10,11,12,13)))

list_counts <- list()
for (sample_name in samples) {
  list_counts[[sample_name]] <- load_counts(
    paste0(sample_name, "_matrix.mtx.gz"),
    paste0(sample_name, "_barcodes.tsv.gz"),
    paste0(sample_name, "_features.tsv.gz"),
    sample_name
  )
}

counts_total <- do.call(cbind, list_counts)

#----------------------- ColData -----------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/metadata')

sce <- SingleCellExperiment(assays = list(counts = counts_total))

# Generar colData con 2 columnas: barcodes y patient id
colData(sce)$PatientID <- sub(".*_", "", colnames(sce))
colData(sce)$barcode <- sub("_.*", "", colnames(sce))
coldata <- as.data.frame(colData(sce))


## Código borja generar_metadatos
metadatos <- read.table("snRNAseqAD_TREM2_assay_scRNAseq_metadata.csv", sep = ",", header = T)
metadatos2 <- read.table("snRNAseqAD_TREM2_biospecimen_metadata.csv", sep = ",", header = T)
metadatos3 <- readxl::read_xlsx("41591_2019_695_MOESM2_ESM.xlsx", sheet = 13)
colnames(metadatos3) <- metadatos3[4,]
metadatos3 <- metadatos3[-c(1:4),]

metadatos <- merge(metadatos, metadatos2, 
                   by = "specimenID")  
metadatos_prueba <-  merge(metadatos, metadatos3, 
                           by.x = "specimenID", 
                           by.y = "Sample ID in snRNA-seq")

rosmap_metadata <- read.table("ROSMAP_clinical.csv", sep = ",", header = T)
rosmap_metadata <- rosmap_metadata[rosmap_metadata$individualID %in% metadatos$individualID,]

metadatos <-  merge(metadatos_prueba, rosmap_metadata, 
                    by.x = "individualID", 
                    by.y = "individualID")

coldata <- merge(coldata, metadatos,
                 by.x = 'PatientID',
                 by.y = 'specimenID',
                 all.x = TRUE)

# Creo columna de sex
coldata <- coldata %>%
  mutate(sex = recode(Sex,
                       'female' = 'Mujer',
                       'male' = 'Hombre'))
# Creo columna de condicion
coldata <- coldata %>%
  mutate(
    condition = case_when(
      str_detect(coldata$`Sample Source`, "AD - Rush|TREM2 - R62H - Rush") ~ "AD",
      str_detect(coldata$`Sample Source`, "Control - Rush") ~ "Control"
    )
  )

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

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
saveRDS(sce, "sce.raw.rds")


