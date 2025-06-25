library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)


#---------------------- Load data and do enrichment ----------------------------

setwd('/clinicfs/projects/i63/tfm_hipocampo/metaanalysis_angela/dorsolateral prefrontal cortex')
dir <- getwd()
cell.types <- list.dirs(path = dir, full.names = FALSE, recursive = FALSE)
results_GO <- list()
results_KEGG <- list()

for (cell.type in cell.types) {
  contrasts <- list.dirs(path = file.path(dir, cell.type), full.names = FALSE, recursive = FALSE)
  for (contrast in contrasts) {
    dir.tsv <- file.path(dir, cell.type, contrast, 'sig.genes.tsv')
    
    if (file.exists(dir.tsv)) {
      tryCatch({
        # Read data
        df <- read.delim(dir.tsv, header = TRUE, sep = '\t')
        original_gene_list <- df$logFC
        names(original_gene_list) <- df$ENSEMBL
        gene_list <- na.omit(original_gene_list)
        gene_list <- sort(gene_list, decreasing = TRUE)
        
        #---------------------- GO ---------------------------------------------
        
        gseG <- gseGO(
          geneList = gene_list,
          ont = "BP",
          keyType = "ENSEMBL",
          OrgDb = 'org.Hs.eg.db',
          minGSSize = 3,
          maxGSSize = 800,
          pvalueCutoff = 0.05,
          pAdjustMethod = "none"
        )
        
        # Save GO
        if (nrow(gseG) > 0) {
          results_GO[[paste0(cell.type, "_", contrast)]] <- gseG
        }
        
        #---------------------- KEGG -------------------------------------------
        
        ids <- bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
        dedup_ids <- ids[!duplicated(ids$ENSEMBL), ]
        df2 <- df[df$ENSEMBL %in% dedup_ids$ENSEMBL, ]
        df2$ENTREZID <- dedup_ids$ENTREZID
        kegg_gene_list <- df2$logFC
        names(kegg_gene_list) <- df2$ENTREZID
        kegg_gene_list <- na.omit(sort(kegg_gene_list, decreasing = TRUE))
        
        gseK <- gseKEGG(
          geneList = kegg_gene_list,
          organism = "hsa",
          keyType = "ncbi-geneid",
          minGSSize = 3,
          maxGSSize = 800,
          pvalueCutoff = 0.05,
          pAdjustMethod = "none"
        )
        
        # Save KEGG
        if (nrow(gseK) > 0) {
          results_KEGG[[paste0(cell.type, "_", contrast)]] <- gseK
        }
        
      }, error = function(e) {
        message(paste("Error en", cell.type, contrast, ":", e$message))
      })
    }
  }
}

#---------------------- Manual result evaluation -------------------------------

# Select manually by GO ID

# Astrocitos
selected_terms<-results_GO$Astrocitos_FM@result[c('GO:0043087','GO:0002682','GO:0030099',
                                                  'GO:0007178','GO:0006954','GO:0006979',
                                                  'GO:0007005','GO:0007612','GO:0048863',
                                                  'GO:0060070','GO:0034620','GO:0002250',
                                                  'GO:0048736','GO:0014070'),]
filtered_enrich <- results_GO$Astrocitos_FM
filtered_enrich@result <- selected_terms
filtered_enrich@geneSets <- filtered_enrich@geneSets[selected_terms$ID]
results_GO$Astrocitos_FM <- filtered_enrich

# Excitadoras
results_GO$Excitadoras_FM[,c('Description', 'ID')]

selected_terms<-results_GO$Excitadoras_FM@result[c('GO:0007166','GO:0006812','GO:0009615',
                                                   'GO:0048584','GO:0043410','GO:0007265',
                                                   'GO:0007264','GO:0002285','GO:0045667',
                                                   'GO:0031589','GO:0061077','GO:0051084',
                                                   'GO:0042026','GO:0006458','GO:0016032'),]

filtered_enrich <- results_GO$Excitadoras_FM
filtered_enrich@result <- selected_terms
filtered_enrich@geneSets <- filtered_enrich@geneSets[selected_terms$ID]
results_GO$Excitadoras_FM <- filtered_enrich

# Inhibidoras
results_GO$Inhibidoras_FM[,c('Description', 'ID')]

selected_terms<-results_GO$Inhibidoras_FM@result[c('GO:0016055','GO:0007265','GO:0007166',
                                                   'GO:0007165','GO:0008219','GO:0012501',
                                                   'GO:0002218','GO:0002221','GO:0007249',
                                                   'GO:0010033','GO:0009605','GO:0032103',
                                                   'GO:0009617','GO:0044283','GO:0022411'),]
filtered_enrich <- results_GO$Inhibidoras_FM
filtered_enrich@result <- selected_terms
filtered_enrich@geneSets <- filtered_enrich@geneSets[selected_terms$ID]
results_GO$Inhibidoras_FM <- filtered_enrich

# Microglia
results_GO$Microglía_FM[,c('Description', 'ID')]
selected_terms<-results_GO$Microglía_FM@result[c('GO:0006950','GO:0001890','GO:0007417',
                                                  'GO:0033554', 'GO:0050793', 'GO:0045444', 
                                                  'GO:0051171', 'GO:0010468', 'GO:0031323', 
                                                  'GO:0045861', 'GO:0032502', 'GO:0048856', 
                                                  'GO:0051239', 'GO:0036211', 'GO:1901360'),]
filtered_enrich <- results_GO$Microglía_FM
filtered_enrich@result <- selected_terms
filtered_enrich@geneSets <- filtered_enrich@geneSets[selected_terms$ID]
results_GO$Microglía_FM <- filtered_enrich

# Oligodendrocitos
results_GO$Oligodendrocitos_FM[,c('Description', 'ID')] # We keep all of them

# OPC
results_GO$OPC_FM[,c('Description', 'ID')]
selected_terms<-results_GO$OPC_FM@result[c('GO:0033554','GO:0006950','GO:0002376',
                                           'GO:0072593','GO:0006796','GO:0016310',
                                           'GO:0006508','GO:0006629','GO:0036211',
                                           'GO:0035556','GO:0007267','GO:0010628',
                                           'GO:1901701','GO:0022607','GO:0065009'),]
                                              
filtered_enrich <- results_GO$Microglía_FM
filtered_enrich@result <- selected_terms
filtered_enrich@geneSets <- filtered_enrich@geneSets[selected_terms$ID]
results_GO$Microglía_FM <- filtered_enrich

#---------------------- Create plots -----------------------

for (result_name in names(results_GO)) {
  cell.type <- strsplit(result_name, "_")[[1]][1]
  contrast <- sub(paste0(cell.type, "_"), "", result_name)
  
  output.dir <- file.path(dir, cell.type, contrast, 'functional_characterization')
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  
  gseG <- results_GO[[result_name]]
  
  if (!is.null(gseG) && nrow(gseG) > 0){
    tryCatch({
    # Dotplot GO
    dotplotGO <- dotplot(gseG, showCategory = 9, split = ".sign", title = cell.type) + 
      facet_grid(.~.sign) +
      theme(axis.text.y = element_text(size = 8)) +
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')

    ggsave(file.path(output.dir, 'dotplot.GO.jpeg'), dotplotGO, width = 10, height = 6, dpi = 300)
    ggsave(file.path(output.dir, 'dotplot.GO.svg'), dotplotGO, width = 10, height = 6)
    
    # Enrichmap GO
    gseG <- pairwise_termsim(gseG)
    enrichmapGO <- emapplot(gseG, showCategory = 10, cex_label_category = 0.6, title = cell.type)+
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
  
    jpeg(file.path(output.dir, 'enrichmap.GO.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(enrichmapGO)
    dev.off()
    svg(file.path(output.dir, 'enrichmap.GO.svg'), 
        width = 10, height = 6)
    print(enrichmapGO)
    dev.off()
    
    # Ridge plot GO
    ridgeplotGO <- ridgeplot(gseG) + labs(x = "Perfil de enrequecimiento") +
      theme(
        axis.text.y = element_text(size = 6)) +
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
    
    jpeg(file.path(output.dir, 'ridgeplot.GO.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(ridgeplotGO)
    dev.off()
    svg(file.path(output.dir, 'ridgeplot.GO.svg'), 
        width = 10, height = 6)
    print(ridgeplotGO)
    dev.off()
    
    #GSEA Plot
    plotverde <- gseaplot(gseG, by = "all", title = gseG$Description[1], geneSetID = 1)
    jpeg(file.path(output.dir, 'plotverde.GO.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(plotverde)
    dev.off()
    svg(file.path(output.dir, 'plotverde.GO.svg'), 
        width = 10, height = 6)
    print(plotverde)
    dev.off()
    }, error = function(e) message("Dotplot falló para ", cell.type, ": ", e$message))
  }
  
  gseK <- results_KEGG[[result_name]]
  
  if (!is.null(gseK) && nrow(gseK) > 0){
    # Dotplot KEGG
    dotplotKEGG <- dotplot(gseK, showCategory = 10, split = ".sign", title = cell.type) + 
      facet_grid(.~.sign) +
      theme(axis.text.y = element_text(size = 8))+
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
    
    ggsave(file.path(output.dir, 'dotplot.KEGG.jpeg'), dotplotKEGG, width = 10, height = 6, dpi = 300)
    ggsave(file.path(output.dir, 'dotplot.KEGG.svg'), dotplotKEGG, width = 10, height = 6)
    
    # Enrichmap KEGG
    gseK <- pairwise_termsim(gseK)
    enrichmapKEGG <- emapplot(gseK, showCategory = 10, cex_label_category = 0.6, title = cell.type)+
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
    
    jpeg(file.path(output.dir, 'enrichmap.KEGG.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(enrichmapKEGG)
    dev.off()
    svg(file.path(output.dir, 'enrichmap.KEGG.svg'), 
        width = 10, height = 6)
    print(enrichmapKEGG)
    dev.off()
    
    # Ridge plot KEGG
    ridgeplotKEGG <- ridgeplot(gseK) + labs(x = "Perfil de enrequecimiento") +
      theme(
        axis.text.y = element_text(size = 6)) +
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
    
    jpeg(file.path(output.dir, 'ridgeplot.KEGG.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(ridgeplotKEGG)
    dev.off()
    svg(file.path(output.dir, 'ridgeplot.KEGG.svg'), 
        width = 10, height = 6)
    print(ridgeplotKEGG)
    dev.off()
    
    #GSEA Plot
    plotverde <- gseaplot(gseK, by = "all", title = gseK$Description[1], geneSetID = 1)
    jpeg(file.path(output.dir, 'plotverde.KEGG.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(plotverde)
    dev.off()
    svg(file.path(output.dir, 'plotverde.KEGG.svg'), 
        width = 10, height = 6)
    print(plotverde)
    dev.off()
  }
}

