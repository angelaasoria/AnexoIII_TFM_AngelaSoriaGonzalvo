library(readr)
library(data.table)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(biomaRt)
library(stringr)
library(GEOquery)
library(openxlsx)

#----------------------- Load data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/GSE157827/data')

# Get files from GEO database and prepare them
geo <- GEOquery::getGEOSuppFiles(GEO = 'GSE157827')
untar('./GSE157827/GSE157827_RAW.tar')

# Rename removing the GSM*
files <- list.files(pattern = "GSM.*", recursive = TRUE)
new.names <- gsub("GSM[0-9]+_", "", files)
file.rename(from = files, to = new.names)

# Function for reading data
load_counts <- function(expr_file, barcode_file, features_file, sample_name) {
  # Read expression data
  counts <- readMM(expr_file)
  
  # Read barcodes
  barcodes <- data.table::fread(barcode_file, header = FALSE)$V1
  barcodes <- paste0(barcodes, "_", sample_name)
  
  features <- data.table::fread(features_file, header = FALSE)$V1
  
  rownames(counts) <- features
  colnames(counts) <- barcodes
  
  return(counts)
}

samples <- c(paste0('AD', c(1,2,4,5,6,8,9,10,13,19,20,21)), 
             paste0('NC', c(3,7,11,12,14,15,16,17,18)))

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

setwd('/clinicfs/projects/i63/tfm_hipocampo/GSE157827/metadata')

sce <- SingleCellExperiment(assays = list(counts = counts_total))

# Generar colData con 2 columnas: barcodes y patient id
colData(sce)$PatientID <- sub(".*_", "", colnames(sce))
colData(sce)$barcode <- sub("_.*", "", colnames(sce))
coldata <- as.data.frame(colData(sce))

metadata <- read.xlsx("pnas.2008762117.sd01.xlsx", 1)
colnames(metadata)[colnames(metadata) == 'ID'] <- 'PatientID'

coldata <- colData(sce)
coldata <- merge(coldata, metadata,
                 by.x = 'PatientID',
                 by.y = 'PatientID',
                 all.x = TRUE)
coldata <- as.data.frame(coldata)

# Create sex column# Creacoldatate sex column
coldata <- coldata %>%
  mutate(sex = recode(SEX,
                       'F' = 'Mujer',
                       'M' = 'Hombre'))
# Create condition column
coldata <- coldata %>%
  mutate(
    condition = case_when(
      str_detect(coldata$CONDITION, "AD") ~ "AD",
      str_detect(coldata$CONDITION, "NC") ~ "Control"
    )
  )

# Create covariable column
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

setwd('/clinicfs/projects/i63/tfm_hipocampo/GSE157827/results')
saveRDS(sce, "sce.raw.rds")
