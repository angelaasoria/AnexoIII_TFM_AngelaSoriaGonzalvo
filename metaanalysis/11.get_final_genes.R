library(tidyverse)
library(UpSetR)
library(MetBrewer)

setwd('/clinicfs/projects/i63/tfm_hipocampo/metaanalysis_angela/dorsolateral prefrontal cortex')

palette <- met.brewer('Renoir')
palette6 <- c(palette[4], palette[6], palette[8], palette[9], palette[11], palette[12])
#-------------------- Plot logFC density ---------------------------------------

# Create dataframe from all significant genes to see logFC distribution
dir <- getwd()
df.logfc <- data.frame()

cell.types <- list.dirs(
  path = dir,
  full.names = FALSE,
  recursive = FALSE
)

for (cell.type in cell.types) {
  dir.tsv <- file.path(dir, cell.type, 'FM', 'sig.genes.tsv')
  if (file.exists(dir.tsv)) {
    tryCatch({
      data <- read_tsv(dir.tsv, show_col_types = FALSE)
      data$cell.type <- cell.type
      data.processed <- data[ ,c('ENSEMBL', 'SYMBOL', 'logFC', 'cell.type')]
      df.logfc <- rbind(df.logfc, data.processed)
      
      message(paste('Procesado:', dir.tsv))
    }, error = function(e) {
      message(paste('Error procesando', dir.tsv, ':', e$message))
    })
  } else {
    message(paste('No se encontró archivo para', cell.type))
  }
}

# Density plot
density <- ggplot(data=df.logfc, aes(x=abs(logFC), group=cell.type, fill=cell.type)) +
                labs(fill = 'Tipo celular') +
                geom_vline(xintercept = 0.25) +
                geom_density(adjust=1.5, alpha=.4)
ggsave('./density_logfc.jpeg', plot = density, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./density_logfc.svg', plot = density, 
       width = 10, height = 6, units = 'in', device = 'svg')


grid.density <- ggplot(data=df.logfc, aes(x=abs(logFC), group=cell.type, fill=cell.type)) +
                    geom_density(adjust=1.5) +
                    facet_wrap(~cell.type) +
                    geom_vline(xintercept = 0.25) +
                    theme(
                      legend.position='none',
                      panel.spacing = unit(0.1, 'lines'),
                      axis.ticks.x=element_blank()
                    )
ggsave('./grid_density_logfc.jpeg', plot = grid.density, 
       width = 10, height = 6, units = 'in', device = 'jpeg', dpi = 300)
ggsave('./grid_density_logfc.svg', plot = grid.density, 
       width = 10, height = 6, units = 'in', device = 'svg')

#--------------------- Get genes with |logFC| > 0.25 ---------------------------

for (cell.type in cell.types) {
  contrasts <- list.dirs(
    path = file.path(dir, cell.type),
    full.names = FALSE,
    recursive = FALSE
  )
  for (contrast in contrasts) {
    dir.tsv <- file.path(dir, cell.type, contrast, 'sig.genes.tsv')
    output.dir <- file.path(dir, cell.type, contrast, 'filtered.sig.genes.tsv')
    
    if (file.exists(dir.tsv)) {
      tryCatch({
        data <- read_tsv(dir.tsv, show_col_types = FALSE)
        data.filtered <- data %>% 
          filter(abs(logFC) > 0.25 & !is.na(SYMBOL))
        
        if (nrow(data.filtered) > 0) {
          write_tsv(data.filtered, output.dir)
          message(paste('Archivo filtrado creado:', output.dir))
        } else {
          message(paste('No hay genes con |logFC| > 0.25 en:', dir.tsv))
        }
        
      }, error = function(e) {
        message(paste('Error procesando', dir.tsv, ':', e$message))
      })
    } else {
      message(paste('No se encontró archivo para', cell.type, 'y contraste', contrast))
    }
  }
}

#--------------------- Upset plot ----------------------------------------------

# Get common genes between each contrast
for (cell.type in cell.types) {
  contrast.genes <- list()
  contrasts <- list.dirs(
    path = file.path(dir, cell.type),
    full.names = FALSE,
    recursive = FALSE
  )
  for (contrast in contrasts) {
    dir.tsv <- file.path(dir, cell.type, contrast, 'filtered.sig.genes.tsv')
    if (file.exists(dir.tsv)) {
      tryCatch({
        data <- read_tsv(dir.tsv, show_col_types = FALSE)
        # Make a contrast for the upregulated in male (+) and another one for the upr. in female (-)
        contrast.genes[[paste0(contrast, " (up)")]] <- data %>% 
          filter(logFC > 0) %>% 
          pull(SYMBOL)
        
        contrast.genes[[paste0(contrast, " (down)")]] <- data %>% 
          filter(logFC < 0) %>% 
          pull(SYMBOL)
        
      }, error = function(e) {
        message(paste('Error procesando', dir.tsv, ':', e$message))
      })
    } else {
      message(paste('No se encontró archivo para', cell.type, 'y contraste', contrast))
    }
  }

  if (length(contrast.genes) > 0) {
    contrast.genes <- contrast.genes[order(names(contrast.genes))]
    contrast.order <- sort(names(contrast.genes))
    all.genes <- unique(unlist(contrast.genes))
    binary.matrix <- sapply(contrast.genes, function(x) as.integer(all.genes %in% x))
    rownames(binary.matrix) <- all.genes
  
    upset.plot <- UpSetR::upset(
      as.data.frame(binary.matrix),
      nsets = length(contrast.genes),
      nintersects = NA,
      order.by = "freq",
      mainbar.y.label = "Número de genes",
      sets.x.label = "Genes por contraste",
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      line.size = 1,
      point.size = 2.8, 
      sets.bar.color = palette6,
      keep.order = TRUE,
      sets = contrast.order
    )

    jpeg(paste0(dir, '/upsetplot_', cell.type, '.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(upset.plot)
    dev.off()
    svg(paste0(dir, '/upsetplot_', cell.type, '.svg'), 
        width = 10, height = 6)
    print(upset.plot)
    dev.off()
    
  }
}

