
pacman::p_load(SingleCellExperiment, scater, scran, dplyr, BiocParallel, 
               scDblFinder, AnnotationHub, tidyverse, patchwork, ggvenn, broom,
               kableExtra,HDF5Array, MetBrewer, DropletUtils,
               devtools, celda, ggvenn, bluster, cluster, ggpubr, MAST, limma, edgeR)


#--------------------------- Load data -----------------------------------------

setwd("/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results")
sce.annot.clean <- readRDS("./sce.annot.clean.rds")

#------------------------- Create functions ------------------------------------

# MAST requires normalization by TPM
## Calculate log2(TPM + 1) counts
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

## Extract results from lrtest
extract_results = function(zz) {
  
  results = data.frame(rownames(zz))
  colnames(results) = "gene_id"
  
  # Selection of
  results[,"lambda"] = data.frame(zz[,3,1]) #Statistic
  results[,"p.value"] = data.frame(zz[,3,3]) #p.value hurdle
  results[,"p.adjusted"] = p.adjust(results[,"p.value"], "BH") #p.adjusted BH
  
  return(results)
}

## logFC standard calculation
LFC.standard = function(data, results, case, control) {
  
  # Calculate cpms and select data of each group. Sum pseudocount factor (1)
  cpm_case = cpm(counts(data))[results[,"gene_id"],data$group.sex==case]  + 1
  cpm_control = cpm(counts(data))[results[,"gene_id"],data$group.sex==control]  + 1
  
  # Calculate logFC
  LFC = as.data.frame(log2(rowMeans(cpm_case)/rowMeans(cpm_control)))
  colnames(LFC) = "logFC.standard"
  LFC$gene_id = rownames(cpm_case)
  
  return(LFC)
}


## logFC MAST calculation
LFC.MAST = function(zlm, control, case, results) {
  coefnames = colnames(zlm@LMlike@modelMatrix)
  
  # Constrast 0
  contrast0 = setNames(rep(0, length(coefnames)), coefnames)
  contrast0[control] = 1
  contrast1=data.frame("case"=setNames(rep(0, length(coefnames)), coefnames))
  contrast1[case,] = 1
  
  
  control=Hypothesis(control, colnames(zlm@LMlike@modelMatrix))
  case=Hypothesis(case, colnames(zlm@LMlike@modelMatrix))
  
  fc=getLogFC(zlm, contrast0, contrast1)
  sum =summary(zlm, logFC=fc)$datatable
  sum=sum[contrast == "case" & component=='logFC',-c(2,3)]
  sum$SE <- (sum[, "ci.hi"] - sum[, "ci.lo"])/ 3.92
  colnames(sum)[c(1,4)]=c("gene_id","logFC")
  sum[is.na(sum)]=0
  
  results = merge(results, sum, by ="gene_id", sort=F)
  
  return(results)
}


#------------------------- Select cell type ------------------------------------

cell.types <- unique(sce.annot.clean$cell.type)

for (i in 1:length(cell.types)){
  
  cell.type = cell.types[i]
  print(paste0("Analizando ",cell.type))
  
  sce.dge = sce.annot.clean[ ,sce.annot.clean$cell.type == cell.type]
  
  dir.create(paste0("./figures/8.MAST/", cell.type,"/"), showWarnings=F)
  path <- paste0("./figures/8.MAST/", cell.type,"/")
  
  
  #------------------------- Get length of each transcript -------------------

  ## "We find that MAST performs best with log-transformed, scale-normalized data that has been thresholded, such as log2(transcripts per million+1)"
  # MAST requiere las TPM para corregir, por lo que nos vamos a descargar la longitud de cada transcrito
  print("Genes iniciales")
  dim(sce.dge)
  
  rowData <- as.data.frame(rowData(sce.dge)) 
  rowData(sce.dge) <- DataFrame(rowData)
  
  rownames(sce.dge) <- rowData(sce.dge)[,'ensembl.id']
  x <- as.data.frame(rowData(sce.dge))
  
  ## Load the length
  transcript <- read.delim("/clinicfs/projects/i63/tfm_hipocampo/MAST/transcript_length_borja.txt", 
                           header = TRUE, sep = "\t")
  
  ## Remove duplicates
  transcript <- transcript[!duplicated(transcript$Gene.stable.ID),]
  ## Add names to rownames
  rownames(transcript) <- transcript$Gene.stable.ID

  ## Keep only the common genes between sce and our database
  # Keep common genes in SCE
  x <- rownames(rowData(sce.dge)) %in% rownames(transcript) 
  sce.dge  <- sce.dge[x, ]
  
  # Keep common genes in DATABASE
  x <- rownames(transcript) %in% rownames(rowData(sce.dge))  
  transcript <- transcript[x, ]
  
  ## Check thats correct
  ensembl.id <- rowData(sce.dge )[,'ensembl.id']
  transcript <- transcript[match(ensembl.id, transcript$Gene.stable.ID),]
  all.equal(rowData(sce.dge)[,'ensembl.id'], transcript$Gene.stable.ID)
  
  ## Add length info of each transcript to sce.dge
  rowData(sce.dge)[,"gene.start"] <- transcript$Gene.start..bp.
  rowData(sce.dge)[,"gene.end"] <- transcript$Gene.end..bp.
  rowData(sce.dge)[,"gene.length"] <- transcript$Gene.end..bp. - transcript$Gene.start..bp. + 1
  
  x <- as.data.frame(rowData(sce.dge))
  patata <- merge(x, transcript, by.x = "ensembl.id", by.y = "Gene.stable.ID")
  rowData(sce.dge) <- DataFrame(patata)
  
  print("Genes que se incluirán en el analisis DEG")
  dim(sce.dge)
  
  #------------------------- Normalize by TPM ----------------------------------
  
  ### Compute log2(TPM + 1) values
  # TPM is very similar to RPKM and FPKM. The only difference is the order of operations. Here’s how you calculate TPM:
  ## Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
  ## Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
  ## Divide the RPK values by the “per million” scaling factor. This gives you TPM.
  
  library(Matrix)
  library(DelayedArray)
  
  assay(sce.dge, "tpm") <- tpm(counts = assay(sce.dge, "counts"), len = rowData(sce.dge)[,"gene.length"])
  sce.dge@assays@data[["log2tpm"]] <- log1p(assay(sce.dge, "tpm")) / log(2)
  
  
  ## Scaling the genes
  # Scale the counts and create a new variable that is saved in sce.
  # For scaling or standarizing, first we first subtract the total mean from each 
  # value and then divide each value by the standard deviation.
  ## Count for each column how many counts are not 0
  cdr2_cells = colSums(assay(sce.dge, "log2tpm")>0)
  colData(sce.dge)$cngeneson = scale(cdr2_cells)  
  
  #------------------------- Modelling with MAST -------------------------------
  
  print(paste0("Modelando para ", cell.type))
  
  # 1.- Build the MAST object
  sce.assay_cells <- SceToSingleCellAssay(sce.dge , check_sanity = FALSE)
  rowData(sce.assay_cells)$primerid = rownames(sce.assay_cells)
  
  # 2.- Prepare the Hurdle model and adjust it to the data
  zlmCond = zlm(~0+group.sex+cngeneson+PatientID, sce.assay_cells, exprs_values = "log2tpm", method = "bayesglm", ebayes = TRUE, parallel = TRUE)
  # Save the model
  saveRDS(zlmCond, file = paste0("./figures/8.MAST/", cell.type,"/zlmCond.RDS"))
  # Save dataset
  saveRDS(sce.dge, file = paste0("./figures/8.MAST/", cell.type,"/sce.dge.RDS"))
}



