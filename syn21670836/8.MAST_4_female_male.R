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

### MAIN ANALYSIS ###

#------------------------- Select cell type ------------------------------------

cell.types <- unique(sce.annot.clean$cell.type)

for (i in 1:length(cell.types)){
  
  cell.type = cell.types[i]
  print(paste0("Analizando ",cell.type))
  
  sce.dge = sce.annot.clean[ ,sce.annot.clean$cell.type == cell.type]
  
  dir.create(paste0("./figures/8.MAST/", cell.type,"/"), showWarnings=F)
  path <- paste0("./figures/8.MAST/", cell.type,"/")
  
  # Load model
  zlmCond <- readRDS(paste0("./figures/8.MAST/", cell.type,"/zlmCond.RDS"))
  # Load dataset
  sce.dge <- readRDS(paste0("./figures/8.MAST/", cell.type,"/sce.dge.RDS"))
  
  #------------------------- Hypotesis contrast --------------------------------
  
  ### (ADMujer - ControlMujer) - (ADHombre - ControlHombre)  ###
  
  # Test hypothesis
  zlmCond_cells <- lrTest(zlmCond, Hypothesis("(group.sexAD_Mujer - group.sexControl_Mujer) - (group.sexAD_Hombre - group.sexControl_Hombre) "))
  saveRDS(zlmCond_cells, file = paste0("./figures/8.MAST/", cell.type,"/zlmCond_cells_FM.RDS"))
  
  # Load data
  zlmCond <- readRDS(paste0(path, "zlmCond.RDS"))
  sce.dge <- readRDS(paste0(path, "sce.dge.RDS"))
  
  # Load models for each contrast
  zlmCond_cells_FM <- readRDS(paste0(path, "zlmCond_cells_FM.RDS"))
  results_2 <- extract_results(zlmCond_cells_FM)

  zlmCond_cellsF <- readRDS(paste0(path, "zlmCond_cells_F.RDS"))
  resultsF <- extract_results(zlmCond_cellsF)
  
  zlmCond_cellsM <- readRDS(paste0(path, "zlmCond_cells_M.RDS"))
  resultsM <- extract_results(zlmCond_cellsM)
  
  LFC.sF <- LFC.standard(data = sce.dge, results = resultsF, case = "AD_Mujer", control = "Control_Mujer")
  LFC.sM <- LFC.standard(data = sce.dge, results = resultsM, case = "AD_Hombre", control = "Control_Hombre")
  LFC.mFM = LFC.MAST(zlm = zlmCond, control = c("group.sexControl_Mujer","group.sexAD_Hombre"), 
                     case =  c("group.sexControl_Hombre","group.sexAD_Mujer"), 
                     results = results_2)
  
  # Calculate logFC with CPM
  logFC_FM <- merge(LFC.sF, LFC.sM, by ="gene_id")
  logFC_FM$logFC.cpm <- logFC_FM[,"logFC.standard.x"] - logFC_FM[,"logFC.standard.y"]
  logFC_FM <- logFC_FM[,c("gene_id","logFC.cpm")]
  
  results <- LFC.mFM
  results$logFC.CPM <- logFC_FM$logFC.cpm

  results$gene.symbol <- rowData(sce.dge)[,"gene.symbol"]
  results <- results[order(results$p.adjusted), ]
  
  write.table(results, file=paste0(path, cell.type ,"res_FM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  write.table(results[results$p.adjusted < 0.05,], file=paste0(path, cell.type,"res_FM_sig.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  
  print("MAST para sex completado")
}

