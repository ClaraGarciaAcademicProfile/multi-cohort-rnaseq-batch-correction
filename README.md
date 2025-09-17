# Análisis de expresión génica en gliomas

Análisis comparativo de datos de expresión génica entre diferentes tipos de gliomas usando datasets públicos de GEO y TCGA.

## Datasets utilizados

- **GSE145645**: Datos de glioblastoma (counts)
- **TCGA-GBM**: Datos de glioblastoma multiforme del proyecto TCGA (counts)
- **TCGA-LGG**: Datos de gliomas de bajo grado del proyecto TCGA  (counts)
- **GSE47774**: Dataset del consorcio SEQC (counts)

## Metodología

El análisis incluye:

1. Carga y procesamiento de datos de múltiples fuentes
2. Normalización TMM (Trimmed Mean of M-values) 
3. Corrección de efectos de batch usando ComBat
4. Visualización mediante PCA y UMAP

## Estructura del código

- `analysis.R`: Script principal con todo el pipeline de análisis
- Los datos deben colocarse en el directorio raíz del proyecto

## Requisitos
```r
library(tidyverse)
library(AnnotationDbi)
library(GEOquery)
library(edgeR)
library(sva)
library(org.Hs.eg.db)
library(umap)
library(cowplot)

## Instalación de dependencias
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("AnnotationDbi", "GEOquery", "edgeR", "sva", "org.Hs.eg.db"))

# CRAN packages
install.packages(c("tidyverse", "umap", "cowplot"))
```
## Uso

- Coloca los archivos de datos en el directorio del proyecto
- Ejecuta el script analysis.R
- Los gráficos PCA y UMAP se generarán automáticamente

## Archivos de datos esperados (de ejemplo)

- GSE145645_XLGBM.txt
- TCGA-LGG.htseq_counts.tsv
- TCGA-GBM.htseq_counts.tsv
- GSE47774_SEQC_ILM_BGI.txt

## Resultados
El análisis genera visualizaciones que muestran el efecto de la corrección de batch sobre la estructura de los datos, para evaluar si las diferencias observadas se deben a efectos técnicos o biológicos reales.