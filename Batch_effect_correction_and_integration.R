# Cargar los paquetes necesarios
library(tidyverse)
library(AnnotationDbi)
library(GEOquery)
library(edgeR)
library(sva)

# Leer datos de expresión del conjunto GSE145645
GSE145645 <- read.table("GSE145645_XLGBM.txt", header = T)

# Obtener información de GEO
GSE145645_geo_info <- getGEO("GSE145645", GSEMatrix = TRUE)
phenoData_GSE145645 <- pData(GSE145645_geo_info[[1]])

# Seleccionar y renombrar columnas del dataframe de phenoData
phenoData_GSE145645 <- phenoData_GSE145645[, c("geo_accession", "description.1", "disease state:ch1")]
colnames(phenoData_GSE145645) <- c("SampleID", "SampleName", "Group")

# Agregar columna de Batch
phenoData_GSE145645$Batch <- "GSE145645"

# Ajustar los nombres de filas y columnas en el conjunto de datos GSE145645
rownames(GSE145645) <- GSE145645$EnsemblID
GSE145645$EnsemblID <- NULL
GSE145645$Gene <- NULL
GSE145645$Category <- NULL

# Filtrar datos y renombrar columnas
phenoData_GSE145645 <- phenoData_GSE145645[phenoData_GSE145645$Group != "Normal", ]
colnames(GSE145645)[colnames(GSE145645) == "GBM11v2"] <- "GBM11"
GSE145645 <- GSE145645[, colnames(GSE145645) %in% phenoData_GSE145645$SampleName]

# Leer datos de expresión para TCGA-LGG
TCGA_LGG <- read.table("TCGA-LGG.htseq_counts.tsv", header = T)
rownames(TCGA_LGG) <- TCGA_LGG$Ensembl_ID
TCGA_LGG$Ensembl_ID <- NULL

# Crear phenoData para TCGA-LGG
phenoData_TCGA_LGG <- data.frame(colnames(TCGA_LGG), colnames(TCGA_LGG))
colnames(phenoData_TCGA_LGG) <- c("SampleID", "SampleName")
phenoData_TCGA_LGG$Group <- "LGG"
phenoData_TCGA_LGG$Batch <- "TCGA_LGG"

# Leer datos de expresión para TCGA-GBM
TCGA_GBM <- read.table("TCGA-GBM.htseq_counts.tsv", header = T)
rownames(TCGA_GBM) <- TCGA_GBM$Ensembl_ID
TCGA_GBM$Ensembl_ID <- NULL

# Crear phenoData para TCGA-GBM
phenoData_TCGA_GBM <- data.frame(colnames(TCGA_GBM), colnames(TCGA_GBM))
colnames(phenoData_TCGA_GBM) <- c("SampleID", "SampleName")
phenoData_TCGA_GBM$Group <- "GBM"
phenoData_TCGA_GBM$Batch <- "TCGA_GBM"

# Leer datos de expresión para GSE47774
GSE47774 <- read.table("GSE47774_SEQC_ILM_BGI.txt", header = T)

# Cargar paquetes necesarios para mapeo de IDs
library(org.Hs.eg.db)

# Mapear IDs de transcriptos a ENSEMBL
GSE47774$ENSEMBL <- mapIds(org.Hs.eg.db,
  keys = GSE47774$TranscriptID,
  column = "ENSEMBL",
  keytype = "REFSEQ",
  multiVals = "first"
)

# Identificar y manejar valores duplicados
duplicates <- GSE47774$ENSEMBL[duplicated(GSE47774$ENSEMBL)]
GSE47774_unique <- GSE47774[!duplicated(GSE47774$ENSEMBL), ]
GSE47774$ENSEMBL <- make.unique(GSE47774$ENSEMBL)

# Eliminar filas con valores NA
GSE47774 <- drop_na(GSE47774)

# Ajustar los nombres de filas en el conjunto de datos GSE47774
rownames(GSE47774) <- GSE47774$ENSEMBL
GSE47774$ENSEMBL <- NULL
GSE47774$TranscriptID <- NULL

# Obtener información de GEO para GSE47774
GSE47774_geo_info <- getGEO("GSE47774", GSEMatrix = TRUE)
phenoData_GSE47774 <- pData(GSE47774_geo_info[[1]])

# Seleccionar y renombrar columnas del dataframe de phenoData
phenoData_GSE47774 <- phenoData_GSE47774[, c("geo_accession", "title", "characteristics_ch1")]
colnames(phenoData_GSE47774) <- c("SampleID", "SampleName", "Group")

# Agregar columna de Batch
phenoData_GSE47774$Batch <- "GSE47774"

# Filtrar datos para el grupo específico
phenoData_GSE47774 <- phenoData_GSE47774[phenoData_GSE47774$Group == "seqc sample: B", ]
GSE47774 <- GSE47774[, colnames(GSE47774) %in% phenoData_GSE47774$SampleName]

# Combinar datos de phenoData
Metadata_experiments <- rbind(
  phenoData_GSE145645,
  phenoData_TCGA_GBM,
  phenoData_TCGA_LGG,
  phenoData_GSE47774
)

# Normalización de datos de expresión
# Asegurarse de que las matrices tengan el mismo conjunto de genes
common_genes <- Reduce(intersect, list(rownames(GSE145645), rownames(TCGA_GBM), rownames(TCGA_LGG), rownames(GSE47774)))

# Filtrar matrices para que solo contengan genes comunes
GSE145645_filtered <- GSE145645[common_genes, ]
TCGA_GBM_filtered <- TCGA_GBM[common_genes, ]
TCGA_LGG_filtered <- TCGA_LGG[common_genes, ]
GSE47774_filtered <- GSE47774[common_genes, ]

# Combinar matrices filtradas
combined_expr <- cbind(GSE145645_filtered, TCGA_GBM_filtered, TCGA_LGG_filtered, GSE47774_filtered)

# Crear un objeto DGEList para la normalización
dge <- DGEList(counts = combined_expr)

# Calcular factores de normalización TMM
dge <- calcNormFactors(dge, method = "TMM")

# Extraer datos normalizados
expr_tmm <- cpm(dge, log = TRUE)

# Crear un vector de batch (asumiendo que hay cuatro batches)
batch <- factor(c(
  rep("GSE145645", ncol(GSE145645_filtered)),
  rep("TCGA_GBM", ncol(TCGA_GBM_filtered)),
  rep("TCGA_LGG", ncol(TCGA_LGG_filtered)),
  rep("GSE47774", ncol(GSE47774_filtered))
))

# Crear un diseño de modelo (solo un intercepto en este caso)
design <- model.matrix(~1, data = data.frame(batch = batch))

# Aplicar SVA para ajustar el efecto de batch
modcombat <- ComBat(dat = expr_tmm, batch = batch, mod = design)

# Función para PCA y gráfico
plot_pca <- function(data, batch, title) {
  pca <- prcomp(t(data), scale. = TRUE)
  pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Batch = batch)
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3) +
    labs(title = title) +
    theme_bw()
}

# PCA antes de la corrección
pca_before <- plot_pca(expr_tmm, batch, "PCA Before Batch Effect Correction")

# PCA después de la corrección
pca_after <- plot_pca(modcombat, batch, "PCA After Batch Effect Correction")

# Mostrar gráficos
plot_grid(pca_after, pca_before, ncol = 2)

# Función para UMAP y gráfico
plot_umap <- function(data, batch, title) {
  umap_result <- umap(t(data))
  umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], Batch = batch)
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Batch)) +
    geom_point(size = 3) +
    labs(title = title) +
    theme_bw()
}

# UMAP antes de la corrección
umap_before <- plot_umap(expr_tmm, batch, "UMAP Before Batch Effect Correction")

# UMAP después de la corrección
umap_after <- plot_umap(modcombat, batch, "UMAP After Batch Effect Correction")

# Mostrar gráficos
print(umap_before)
print(umap_after)

plot_grid(umap_before, umap_after, ncol = 2)
