library(scran)

#----------------------- Load data --------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/results')
sce.norm <- readRDS('./sce.norm.rds')

#----------------------- Quantifying variance of log-counts --------------------

# Modelling gene variation
dec <- modelGeneVar(sce.norm, assay.type = "logcounts")
curfit <- metadata(dec)

par(mar = c(5, 5, 2, 1))
jpeg(file = './figures/3.feature_selection/modelgenevar_overfitting.jpeg', 
      width=10, height=6, units="in", res=300)
plot(dec$mean, dec$total,
      xlab = "Media de expresión log-normalizada", 
      ylab = "Varianza de expresión log-normalizada")
curve(curfit$trend(x), col = "#a8554e", add = TRUE, lwd = 2)
dev.off()
svg(file = './figures/3.feature_selection/modelgenevar_overfitting.svg', 
    width=10, height=6)
plot(dec$mean, dec$total,
      xlab = "Media de expresión log-normalizada", 
      ylab = "Varianza de expresión log-normalizada")
curve(curfit$trend(x), col = "#a8554e", add = TRUE, lwd = 2)
dev.off()

# To avoid overfitting, we add the parameter density.weights = FALSE
dec <- modelGeneVar(sce.norm, assay.type = "logcounts", density.weights = FALSE)
curfit <- metadata(dec)

par(mar = c(5, 5, 2, 1))
jpeg(file = './figures/3.feature_selection/modelgenevar_goodfit.jpeg', 
      width=10, height=6, units="in", res=300)
plot(dec$mean, dec$total, main = i,
      xlab = "Media de expresión log-normalizada", 
      ylab = "Varianza de expresión log-normalizada") 
curve(curfit$trend(x), col = "#a8554e", add = TRUE, lwd = 2)
dev.off()
svg(file = './figures/3.feature_selection/modelgenevar_goodfit.svg', 
    width=10, height=6)
plot(dec$mean, dec$total, main = i,
      xlab = "Media de expresión log-normalizada", 
      ylab = "Varianza de expresión log-normalizada") 
curve(curfit$trend(x), col = "#a8554e", add = TRUE, lwd = 2)
dev.off()


# Take a look at most variable genes
dec[order(dec$bio, decreasing=TRUE),]

#----------------------- Get HVGs ----------------------------------------------

tops <- getTopHVGs(dec, n = 3000)
sce.hvg <- sce.norm[tops, ]

# Save the original as an alternative experiment of the chosen
altExp(sce.hvg, 'original') <- sce.norm
altExpNames(sce.hvg)

# If we need to recover original data:
# sce.original <- altExp(sce.hvg, 'original', withColData=TRUE)

#----------------------- Save data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn18485175/results')
saveRDS(sce.hvg, './sce.hvg.rds')


