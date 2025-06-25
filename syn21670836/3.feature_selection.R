library(scran)

#----------------------- Load data --------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
sce.norm <- readRDS('./sce.norm.rds')

library(scran)

#----------------------- Quantifying variance of log-counts --------------------

sce.norm$block = paste0(sce.norm$libraryBatch, "_", sce.norm$sequencingBatch)
# Modelling gene variation
dec <- modelGeneVar(sce.norm, assay.type = "logcounts", block=sce.norm$block)
blocked.stats <- dec$per.block

par(mar = c(5, 5, 2, 1))
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  curfit <- metadata(current)
  
  jpeg(file = paste0('./figures/3.feature_selection/modelgenevar_overfitting_',i,'.jpeg'), 
       width=10, height=6, units="in", res=300)
  plot(current$mean, current$total, main = i,
       xlab = "Media de expresión log-normalizada", 
       ylab = "Varianza de expresión log-normalizada")
  curve(curfit$trend(x), col = "#a8554e", add = TRUE, lwd = 2)
  dev.off()
  svg(file = paste0('./figures/3.feature_selection/modelgenevar_overfitting_',i,'.svg'), 
      width=10, height=6)
  plot(current$mean, current$total, main = i,
       xlab = "Media de expresión log-normalizada", 
       ylab = "Varianza de expresión log-normalizada")
  curve(curfit$trend(x), col = "#a8554e", add = TRUE, lwd = 2)
  dev.off()
}


# To avoid overfitting, we add the parameter density.weights = FALSE
dec <- modelGeneVar(sce.norm, assay.type = "logcounts", 
                    block=sce.norm$block, density.weights = FALSE)
blocked.stats <- dec$per.block

par(mar = c(5, 5, 2, 1))
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  curfit <- metadata(current)
  
  jpeg(file = paste0('./figures/3.feature_selection/modelgenevar_goodfit_',i,'.jpeg'), 
       width=10, height=6, units="in", res=300)
  plot(current$mean, current$total, main = i,
       xlab = "Media de expresión log-normalizada", 
       ylab = "Varianza de expresión log-normalizada") 
  curve(curfit$trend(x), col = "#a8554e", add = TRUE, lwd = 2)
  dev.off()
  svg(file = paste0('./figures/3.feature_selection/modelgenevar_goodfit_',i,'.svg'), 
      width=10, height=6)
  plot(current$mean, current$total, main = i,
       xlab = "Media de expresión log-normalizada", 
       ylab = "Varianza de expresión log-normalizada") 
  curve(curfit$trend(x), col = "#a8554e", add = TRUE, lwd = 2)
  dev.off()
}


# Take a look at most variable genes
dec[order(dec$bio, decreasing=TRUE),]

#----------------------- Get HVGs ------------------------------------------

tops <- getTopHVGs(dec, n = 3000)
sce.hvg<- sce.norm[tops, ]


# Save the original as an alternative experiment of the chosen
altExp(sce.hvg, 'original') <- sce.norm
altExpNames(sce.hvg)

# If we need to recover original data:
# sce.original <- altExp(sce.hvg, 'original', withColData=TRUE)

#----------------------- Save data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn21670836/results')
saveRDS(sce.hvg, './sce.hvg.rds')
