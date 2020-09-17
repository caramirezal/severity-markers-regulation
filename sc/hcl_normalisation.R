library(Seurat)
library(dplyr)
library(Matrix)
library(rocket)
library(viridis)
library(extrafont)
library(ggplot2)

setwd('~/sc/severity-markers-regulation/')

## Gene set visualization in Human Cell Landscape
## Reading Human Cell Landscape data
hcl <- readRDS('data/hcl_normalized_adult_seu.rds')

## Subsetting to adult tissue
hcl$'adult' <- grepl('adult', tolower(hcl$orig.ident)) &
        ! ( grepl('fetal', tolower(hcl$celltype)))
hcl <- subset(hcl, adult == TRUE)

###############################################################################
## Correlation of SPINT2 and TMPRSS2 in tissues
celltypes <- c('Epithelial cell',
               'Endothelial cell',
               'Enterocyte progenitor',
               'Epithelial cell (intermediated)',
               'AT2 cell',
               'Goblet cell',
               'Stratified epithelial cell',
               'Kidney intercalated cell',
               'Loop of Henle',
               'Enterocyte',
               'Enterocyte progenitor',
               'Gastric chief cell',
               'Gastric endocrine cell',
               'Hepatocyte/Endodermal cell',
               'Fibroblast',
               'Macrophage',
               'M2 Macrophage',
               'Stromal cell',
               'Smooth muscle cell',
               'Endothelial cell (APC)',
               'Endothelial cell (endothelial to mesenchymal transition)',
               'Sinusoidal endothelial cell',
               'Ventricle cardiomyocyte',
               'Proximal tubule progenitor',
               'Pancreas exocrine cel',
               'Dendritic cell',
               'hESC')

####################################################################
## Re-normalisation
hcl <- subset(hcl, celltype %in% celltypes)

## SCT transform
hcl <- SCTransform(hcl, variable.features.n = 12000)

saveRDS(hcl, 'data/hcl_sct_norm.rds')



