library(ggplot2)
library(scales)
library(ggthemes)
library(scuttle)
library(scran)
library(MetBrewer)
library(patchwork)

#----------------------- Load data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn52293442_pfc/results/')
sce.filt.n <- readRDS('./sce.filt.rds')

# Generate palettes
palette <- met.brewer('Renoir')
palette2 <- c(palette[4], palette[9])
palette4 <- c(palette[4], palette[8], palette[9], palette[11])

#----------------------- Library size factors ----------------------------------

# Normalize per library size factor each cell
lib.sce.filt <- librarySizeFactors(sce.filt.n) 
sce.filt.n@colData$library.size.factors <- lib.sce.filt
summary(lib.sce.filt)

#----------------------- Plot not normalized data ------------------------------

jpeg(file = './figures/2.normalization/factorfreq.BN.jpeg', width=10, height=6, units="in", res=300)
hist(lib.sce.filt, xlab="Factor de normalización", ylab="Frecuencia", col="grey80")
dev.off()
svg(file = './figures/2.normalization/factorfreq.BN.svg', width=10, height=6)
hist(lib.sce.filt, xlab="Factor de normalización", ylab="Frecuencia", col="grey80")
dev.off()


jpeg(file = './figures/2.normalization/factorfreq.log.BN.jpeg', width=10, height=6, units="in", res=300)
hist(log10(lib.sce.filt), xlab="log10(Factor de normalización)", ylab="Frecuencia", col="grey80")
dev.off()
svg(file = './figures/2.normalization/factorfreq.log.BN.svg', width=10, height=6)
hist(log10(lib.sce.filt), xlab="log10(Factor de normalización)", ylab="Frecuencia", col="grey80")
dev.off()


lib.size <- colSums(counts(sce.filt.n))
jpeg(file = './figures/2.normalization/sizefactor.BN.jpeg', width=10, height=6, units="in", res=300)
plot(lib.size, lib.sce.filt, xlab="Tamaño de la librería", ylab="Factor de normalización")
dev.off()
svg(file = './figures/2.normalization/sizefactor.BN.svg', width=10, height=6)
plot(lib.size, lib.sce.filt, xlab="Tamaño de la librería", ylab="Factor de normalización")
dev.off()

#----------------------- Normalization by deconvolution ------------------------

set.seed(100)
quick.clusters <- scran::quickCluster(sce.filt.n, method = "igraph")
sce.filt.n <- computeSumFactors(sce.filt.n, clusters = quick.clusters)
sce.norm <- logNormCounts(sce.filt.n, 
                          assay.type = "counts",
                          log=TRUE, 
                          pseudo.count=1)

print("Summary normalización por deconvolución")
summary(sizeFactors(sce.norm))
sce.norm@colData$deconvolution.size.factors <- sizeFactors(sce.norm)

# Plot after normalizing
jpeg(file = './figures/2.normalization/factorfreq.AN.jpeg', width=10, height=6, units="in", res=300)
hist(sizeFactors(sce.norm), xlab="Factor de normalización",ylab="Frecuencia", col="grey80")
dev.off()
svg(file = './figures/2.normalization/factorfreq.AN.svg', width=10, height=6)
hist(sizeFactors(sce.norm), xlab="Factor de normalización",ylab="Frecuencia", col="grey80")
dev.off()


jpeg(file = './figures/2.normalization/factorfreq.log.AN.jpeg', width=10, height=6, units="in", res=300)
hist(log10(sizeFactors(sce.norm)), xlab="log10(Factor de normalización)", ylab="Frecuencia", col="grey80")
dev.off()
svg(file = './figures/2.normalization/factorfreq.log.AN.svg', width=10, height=6)
hist(log10(sizeFactors(sce.norm)), xlab="log10(Factor de normalización)", ylab="Frecuencia", col="grey80")
dev.off()


jpeg(file = './figures/2.normalization/sizefactor.AN.jpeg', width=10, height=6, units="in", res=300)
plot(lib.size, sizeFactors(sce.norm), xlab="Tamaño de la librería", ylab="Factor de normalización")
dev.off()
svg(file = './figures/2.normalization/sizefactor.AN.svg', width=10, height=6)
plot(lib.size, sizeFactors(sce.norm), xlab="Tamaño de la librería", ylab="Factor de normalización")
dev.off()

#----------------------- Plot normalized data ----------------------------------

sce.norm@colData$nCountAF <- colSums(counts(sce.filt.n))
sce.norm@colData$nCountnorm = Matrix::colSums(logcounts(sce.norm))
coldata.norm <- as.data.frame(sce.norm@colData)

# Before and after normalization

# By condition
plot.size.condition.BN <- ggplot(coldata.norm, aes(x=condition, y=nCountAF, fill=condition)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Condición")+
  ylab("Tamaño de la librería sin normalizar")+
  scale_fill_manual(values = palette2, name= "Condición", labels = c('AD', 'Control'))

plot.size.condition.AN <- ggplot(coldata.norm, aes(x=condition, y=nCountnorm, fill=condition)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Condición")+
  ylab("Tamaño de la librería log2-normalizado")+
  scale_fill_manual(values = palette2, name= "Condición", labels = c('AD', 'Control'))

comb.plots.size.condition <- plot.size.condition.BN + plot.size.condition.AN

jpeg(file = './figures/2.normalization/size.condition.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.size.condition)
dev.off()
svg(file = './figures/2.normalization/size.condition.svg', width=10, height=6)
print(comb.plots.size.condition)
dev.off()

# By sex
plot.size.sex.BN <- ggplot(coldata.norm, aes(x=sex, y=nCountAF, fill=sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Sexo")+
  ylab("Tamaño de la librería sin normalizar")+
  scale_fill_manual(values = palette2, name= "Sexo", labels = c('AD', 'Control'))

plot.size.sex.AN <- ggplot(coldata.norm, aes(x=sex, y=nCountnorm, fill=sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Sexo")+
  ylab("Tamaño de la librería log2-normalizado")+
  scale_fill_manual(values = palette2, name= "Sexo", labels = c('AD', 'Control'))

comb.plots.size.sex <- plot.size.sex.BN + plot.size.sex.AN

jpeg(file = './figures/2.normalization/size.sex.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.size.sex)
dev.off()
svg(file = './figures/2.normalization/size.sex.svg', width=10, height=6)
print(comb.plots.size.sex)
dev.off()

# By condition and sex

plot.size.groupsex.BN <- ggplot(coldata.norm, aes(x=group.sex, y=nCountAF, fill=group.sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Condición y sexo")+
  ylab("Tamaño de la librería sin normalizar")+
  scale_fill_manual(values = palette4, name= "Condición y sexo", 
                    labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

plot.size.groupsex.AN <- ggplot(coldata.norm, aes(x=group.sex, y=nCountnorm, fill=group.sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Condición y sexo")+
  ylab("Tamaño de la librería log2-normalizado")+
  scale_fill_manual(values = palette4, name= "Condición y sexo", 
                    labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

comb.plots.size.groupsex <- plot.size.groupsex.BN + plot.size.groupsex.AN

jpeg(file = './figures/2.normalization/size.groupsex.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.size.groupsex)
dev.off()
svg(file = './figures/2.normalization/size.groupsex.svg', width=10, height=6)
print(comb.plots.size.groupsex)
dev.off()


#By patient
plot.size.patient.BN <- ggplot(coldata.norm, aes(x=projid, y=nCountAF, fill=group.sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Paciente")+
  ylab("Tamaño de la librería sin normalizar")+
  scale_fill_manual(values = palette4, name = 'Condición y sexo',
                    labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

plot.size.patient.AN <- ggplot(coldata.norm, aes(x=projid, y=nCountnorm, fill=group.sex)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Paciente")+
  ylab("Tamaño de la librería log2-normalizado")+
  scale_fill_manual(values = palette4, name = 'Condición y sexo',
                    labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

comb.plots.size.patient <- plot.size.patient.BN + plot.size.patient.AN

jpeg(file = './figures/2.normalization/size.patient.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.size.patient)
dev.off()
svg(file = './figures/2.normalization/size.patient.svg', width=10, height=6)
print(comb.plots.size.groupsex)
dev.off()

#----------------------- Methods comparison ------------------------------------

methods.condition <- ggplot(coldata.norm, aes(
  x=library.size.factors, y=deconvolution.size.factors, color=condition)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="black") +
  theme_classic() +
  xlab("Factor de normalización por tamaño de librería")+
  ylab("Factor de normalización por deconvolución")+
  scale_color_manual(values = palette, name = "Condición", labels = c('AD', 'Control'))

methods.groupsex <- ggplot(coldata.norm, aes(
  x = library.size.factors, y = deconvolution.size.factors, color = group.sex)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="black") +
  theme_classic() +
  xlab("Factor de normalización por tamaño de librería")+
  ylab("Factor de normalización por deconvolución")+
  scale_color_manual(values = palette, name = "Condición y sexo", 
                     labels = c('AD_Hombre', 'AD_Mujer', 'Control_Hombre', 'Control_Mujer'))

comb.plots.methods.comparison <- methods.condition + methods.groupsex

jpeg(file = './figures/2.normalization/methods.comparison.jpeg', width=10, height=6, units="in", res=300)
print(comb.plots.methods.comparison)
dev.off()
svg(file = './figures/2.normalization/methods.comparison.svg', width=10, height=6)
print(comb.plots.methods.comparison)
dev.off()

#----------------------- Save data ---------------------------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/syn52293442_pfc/results')
saveRDS(sce.norm, file = "./sce.norm.rds")