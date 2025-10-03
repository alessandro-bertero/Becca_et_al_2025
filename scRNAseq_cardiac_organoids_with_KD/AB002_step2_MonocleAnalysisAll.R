# Set up the environment with required packages
library(monocle3)
library(tidyr)
library(stringr)
library(grid)
library(viridis)
library(ggthemes)
library(dplyr)
library(tidyverse)
library(ggsignif)
library(umap)
library(heatmap3)
library(plyr)
library(edgeR)
library(compareGroups)
library(dbscan)
library(MAST)
library(RColorBrewer)
library(Seurat)
library(htmlwidgets)
library(Matrix)
library(ggplot2)
library(reshape2)

# Session options
options(stringsAsFactors = FALSE)
set.seed(12345)

# Set up the ggplot default params
theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size=15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle=45, hjust = 1)))


# Set up filename prefixes and output folders with the data
dir.create(file.path("/home/shared_folder/Output", format(Sys.Date(), "%y%m%d")))
dir.create(file.path("/home/shared_folder/Output", format(Sys.Date(), "%y%m%d"), "Monocle"))
monocle_prefix = paste0("/home/shared_folder/Output/", format(Sys.Date(), "%y%m%d"), "/Monocle/")

#palettes
okabe_pal=c("#E69F00","#56B4E9","#009E73","#f0E442","#0072B2","#D55E00","#CC79A7","#000000")
pal = c("#000000","#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff", "#920000","#924900","#db6d00","#24ff24","#ffff6d")

#load data
load("/home/shared_folder/RData/seurat_filtered2.RData")
load("/home/shared_folder/RData/seurat_meta_norm.RData")
genes_data = as.data.frame(Features(data_seurat))
rm(data_seurat)
colnames(genes_data) = "gene_short_name"
rownames(genes_data) = genes_data$gene_short_name

cds <- new_cell_data_set(counts_post,
                         cell_metadata = meta_norm,
                         gene_metadata = genes_data)
rm(counts_post)

#preprocess cds
cds <- preprocess_cds(cds, num_dim = 100)

#dimension reduction
cds <- reduce_dimension(cds, reduction_method = "UMAP", umap.min_dist = 0.2)
cds_batch <- reduce_dimension(cds_batch, reduction_method = "UMAP", umap.min_dist = 0.2)

#PCA plots
plot_pc_variance_explained(cds)
ggsave(filename = paste0(monocle_prefix, "PCA_VarExp.pdf"),
       width = 12, height = 8)

p1 = plot_cells(cds, reduction_method = "PCA", color_cells_by = "nCount_RNA")
p2 = plot_cells(cds, reduction_method = "PCA", color_cells_by = "nFeature_RNA")
p1 + p2
ggsave(filename = paste0(monocle_prefix, "PCA_CountsFeature.pdf"),
       width = 20, height = 8)

p1 = plot_cells(cds, reduction_method = "PCA", color_cells_by = "percent.mt")
p2 = plot_cells(cds, reduction_method = "PCA", color_cells_by = "percent.ribo")
p1 + p2
ggsave(filename = paste0(monocle_prefix, "PCA_MtRibo.pdf"),
       width = 20, height = 8)

PCA = data.frame(p1$data$sample_name, libID = p1$data$sample_name, PC1 = p1$data$data_dim_1, PC2 = p1$data$data_dim_2, row.names = 1)
colData(cds)$PC1 = PCA$PC1
colData(cds)$PC2 = PCA$PC2

p1 = plot_cells(cds, reduction_method = "PCA", color_cells_by = "replicate", label_cell_groups = FALSE)
p2 = plot_cells(cds, reduction_method = "PCA", color_cells_by = "day", label_cell_groups = FALSE)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "PCA_replicates.pdf"),
       width = 20, height = 8)

#clustering
cds <- cluster_cells(cds, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)
cds_batch <- cluster_cells(cds_batch, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)

p1 = plot_cells(cds, norm_method = "log", reduction_method = "UMAP", color_cells_by = "replicate", label_cell_groups = FALSE)
p2 = plot_cells(cds, norm_method = "log", reduction_method = "UMAP", color_cells_by = "day", label_cell_groups = FALSE)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_replicates.pdf"),
       width = 20, height = 12)

p1 = plot_cells(cds, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nCount_RNA")
p2 = plot_cells(cds, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nFeature_RNA")
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_QC_counts.pdf"),
       width = 20, height = 12)
p1 = plot_cells(cds, norm_method = "log", reduction_method = "UMAP", color_cells_by = "percent.mt")
p2 = plot_cells(cds, norm_method = "log", reduction_method = "UMAP", color_cells_by = "percent.ribo")
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_QC_mtRibo.pdf"),
       width = 20, height = 12)
p1 = plot_cells(cds, norm_method = "log", reduction_method = "UMAP", color_cells_by = "PC1")
p2 = plot_cells(cds, norm_method = "log", reduction_method = "UMAP", color_cells_by = "PC2")
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_QC_PCA.pdf"),
       width = 20, height = 12)
UMAP = data.frame(p1$data$sample_name, libID = p1$data$sample_name, UMAP1 = p1$data$data_dim_1, UMAP2 = p1$data$data_dim_2, monocle_cluster = p1$data$cell_color, row.names = 1)
colData(cds)$UMAP1 = UMAP$UMAP1
colData(cds)$UMAP2 = UMAP$UMAP2
colData(cds)$monocle_cluster = UMAP$monocle_cluster

#filtering
cds_subset <- cds[, 
                  colData(cds) %>%
                    subset(
                      PC1 > -50
                    ) %>%
                    row.names
]

geneInfo <- read.table("/home/shared_folder/geneInfo.tab", quote="\"", comment.char="", skip = 1)
colnames(geneInfo) = c("geneID", "gene_short_name", "gene_type")
geneInfo_pc = geneInfo %>% dplyr::filter(gene_type == "protein_coding")
cds_subset <- cds_subset[
  cds_subset@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% geneInfo_pc$gene_short_name,
]

p1 = plot_cells(cds_subset, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nCount_RNA")
p2 = plot_cells(cds_subset, norm_method = "log", reduction_method = "UMAP", color_cells_by = "day", label_cell_groups = FALSE)
p3 = plot_cells(cds_subset, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE)
p1 + p2 + p3
ggsave(filename = paste0(monocle_prefix, "UMAP_subset_old_clusters.pdf"),
       width = 24, height = 8)

#time trajectory
cds_subset <- reduce_dimension(cds_subset, reduction_method = "UMAP", umap.min_dist = 0.2)
cds_subset <- cluster_cells(cds_subset, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)

p1 = plot_cells(cds_subset, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nCount_RNA")
p2 = plot_cells(cds_subset, norm_method = "log", reduction_method = "UMAP", color_cells_by = "day", label_cell_groups = FALSE)
p3 = plot_cells(cds_subset, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE)
p1 + p2 + p3
ggsave(filename = paste0(monocle_prefix, "UMAP_subset_new_clusters.pdf"),
       width = 24, height = 8)

cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset, norm_method = "log", reduction_method = "UMAP", genes = c("TNNT2", "TTN", "TBXT", "EOMES"))
ggsave(filename = paste0(monocle_prefix, "UMAP_genes.pdf"),
       width = 20, height = 20)

save(cds_subset,
     file = paste0("/home/shared_folder/RData/cds_filtered_all.RData"))




