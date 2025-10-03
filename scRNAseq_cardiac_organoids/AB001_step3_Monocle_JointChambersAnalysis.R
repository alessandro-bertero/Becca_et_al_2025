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
library(ensembldb)
library(topGO)
library(limma)
library(org.Hs.eg.db)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(dqrng)
library(clusterProfiler)
library(enrichplot)

# Session options
options(stringsAsFactors = FALSE)

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
load("/home/shared_folder/RData/cds_chambers_annotated.RData")

p1 = plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
df = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
df = df %>% dplyr::mutate(cluster_name = ifelse(cluster == "1", "mature_cardio_AT",
                                                ifelse(cluster == "2", "cardiac_fibroblasts_AT",
                                                       ifelse(cluster == "3", "proliferating_progenitors_AT",
                                                              ifelse(cluster == "4", "endothelial_cells_AT",
                                                                     ifelse(cluster == "5", "advanced_cardio_AT", ""))))))
colData(cds_AT)$cluster_name = df$cluster_name

p1 = plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
fibroblasts = choose_cells(cds_RV)
endothelial = choose_cells(cds_RV)
df = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
df = df %>% mutate(cluster = ifelse(libID %in% endothelial@colData@rownames, "7",
                                    ifelse(libID %in% fibroblasts@colData@rownames, "6", cluster)))
df = df %>% dplyr::mutate(cluster_name = ifelse(cluster == "1", "mature_cardio_RV",
                                                ifelse(cluster == "2", "immature_cardio_RV",
                                                       ifelse(cluster == "3", "proliferating_progenitors_RV",
                                                              ifelse(cluster == "4", "undifferentiated_cells_RV",
                                                                     ifelse(cluster == "5", "advanced_cardio_RV",
                                                                            ifelse(cluster == "6", "cardiac_fibroblasts_RV",
                                                                                   ifelse(cluster == "7", "endothelial_cells_RV", ""))))))))
colData(cds_RV)$cluster_name = df$cluster_name

p1 = plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
df = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
endothelial = choose_cells(cds_LV)
df = df %>% mutate(cluster = ifelse(libID %in%endothelial@colData@rownames, "5", cluster))
df = df %>% dplyr::mutate(cluster_name = ifelse(cluster == "1", "mature_cardio_LV",
                                                ifelse(cluster == "2", "cardiac_fibroblasts_LV",
                                                       ifelse(cluster == "3", "proliferating_progenitors_LV",
                                                              ifelse(cluster == "4", "epithelial_cells_LV", 
                                                                     ifelse(cluster == "5", "endothelial_cells_LV", ""))))))
colData(cds_LV)$cluster_name = df$cluster_name

cds_all <- combine_cds(list(cds_AT, cds_RV, cds_LV))
cds_all <- preprocess_cds(cds_all)
cds_all <- reduce_dimension(cds_all, reduction_method = "UMAP", umap.min_dist = 0.05)
cds_all <- cluster_cells(cds_all, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)

chamber = c()
identity = c()
meta = cds_all@colData
for(i in meta$cluster){
  chamber = append(chamber, strsplit(i, "_")[[1]][3])
  identity = append(identity, paste0(strsplit(i, "_")[[1]][1], "_", strsplit(i, "_")[[1]][2]))
}
p1 = plot_cells(cds_all, norm_method = "log", reduction_method = "UMAP")
meta$UMAP1 = p1[["data"]][["data_dim_1"]]
meta$UMAP2 = p1[["data"]][["data_dim_2"]]
meta$chamber = chamber
meta$identity = identity
colData(cds_all)$identity = identity
colData(cds_all)$chamber = chamber
plot_cells(cds_all, norm_method = "log", reduction_method = "UMAP", color_cells_by = "chamber",  label_cell_groups = FALSE, cell_size = 1, alpha = 0.5) + scale_color_manual(values = c("#CC79A7", "#0072B2", "#D55E00"))
ggsave(paste0(monocle_prefix, "UMAP_chamber.pdf"))
plot_cells(cds_all, norm_method = "log", reduction_method = "UMAP", color_cells_by = "identity",  label_cell_groups = FALSE, cell_size = 1, alpha = 0.5) + scale_color_manual(values =  c("#920000","#1CBE4F","#AA0DFE","#ff6db6","#E69F00","#F6222E","#56B4E9","grey"))
ggsave(paste0(monocle_prefix, "UMAP_identity.pdf"))

markers = c("COL3A1", "ERG", "MKI67", "TNNT2", "TTN", "KRT19", "POU5F1", "NR2F2", "IRX4", "IRX1")
for(i in markers){
  plot_cells(cds_all, genes = i, norm_method = "log", label_cell_groups = FALSE, cell_size = 1)
  ggsave(paste0(monocle_prefix, "UMAP_", i, ".pdf"))
}
meta_all = cds_all@colData
save(cds_all, meta_all,
     file = paste0("/home/shared_folder/RData/cds_joined.RData"))

##cluster transfering
#GSE263193 data
system2("wget", "https://zenodo.org/records/10932845/files/scRoadmap_CMDiff_ShinyApp_iSEE.rds?download=1")
scRoadmap <- readRDS("scRoadmap_CMDiff_ShinyApp_iSEE.rds")
meta = as.data.frame(scRoadmap@colData)
seurat_ref <- as.Seurat(scRoadmap)
genes_AB001 = rowData(cds_ref)$gene_short_name

genes_ref = as.data.frame(seurat_ref@assays[["RNA"]]@counts@Dimnames[[1]])
colnames(genes_ref) = "gene_short_name"
rownames(genes_ref) = genes_ref$gene_short_name
cds_ref <- new_cell_data_set(
  seurat_ref@assays$RNA@counts,
  cell_metadata = meta,
  gene_metadata = genes_ref
)
cds_subset = cds_ref[rowData(cds_ref)$gene_short_name %in% genes_AB001]

counts <- assay(cds_subset, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(cds_subset) <- log2(t(t(counts)/size.factors) + 1)
seurat_ref_filtered = as.Seurat(cds_subset, data = NULL)
seurat_ref_filtered <- FindVariableFeatures(seurat_ref_filtered, selection.method = "vst", nfeatures = 8000)
seurat_ref_filtered = ScaleData(seurat_ref_filtered)
seurat_ref_filtered = RunPCA(seurat_ref_filtered, npcs = 30)
DimPlot(seurat_ref_filtered, group.by = "ident", label = F)
ggsave(paste0(monocle_prefix, "clusters_ref_filtered.pdf"), width = 10, height = 10)
genes_ref = seurat_ref_filtered@assays[["originalexp"]]@counts@Dimnames[[1]]

#AB001 data
cds_subset_AB001 = cds_all[rowData(cds_all)$gene_short_name %in% genes_ref]
counts <- assay(cds_subset_AB001, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(cds_subset_AB001) <- log2(t(t(counts)/size.factors) + 1)
seurat_subset_AB001 <- as.Seurat(cds_subset_AB001, data = NULL)

seurat_subset_AB001 <- FindVariableFeatures(seurat_subset_AB001, selection.method = "vst", nfeatures = 8000)
seurat_subset_AB001 = ScaleData(seurat_subset_AB001)
seurat_subset_AB001 = RunPCA(seurat_subset_AB001 , npcs = 30)

#anchors
anchors = FindTransferAnchors(reference = seurat_ref_filtered, query = seurat_subset_AB001)
predictions_AB001 <- TransferData(anchorset = anchors, 
                                  refdata = seurat_ref_filtered$ident,
                                  reference = seurat_ref_filtered, query = seurat_subset_AB001)
DimPlot(predictions_AB001, group.by = "predicted.id", label = F)
ggsave(paste0(monocle_prefix, "predicted_GSE263193.pdf"), width = 10, height = 10)

#quantification
meta_predicted = predictions_AB001@meta.data
meta_predicted$n = rep(c(1), times = length(rownames(meta_predicted)))
score = data.frame(predicted.id = unique(meta_predicted$predicted.id))
perc = data.frame(predicted.id = unique(meta_predicted$predicted.id))
for(i in unique(meta_predicted$identity)){
  meta_i = meta_predicted %>% dplyr::filter(identity == i)
  meta_i = meta_i %>% dplyr::group_by(predicted.id) %>% dplyr::summarise(sum = sum(n), score = sum(predicted.id.score))
  meta_i = meta_i %>% dplyr::mutate(score = score/sum*100) %>% dplyr::mutate(sum = sum/sum(sum)*100)
  meta_score = meta_i %>% dplyr::select(c("predicted.id", "score"))
  colnames(meta_score) = c("predicted.id", i)
  meta_perc = meta_i %>% dplyr::select(c("predicted.id", "sum"))
  colnames(meta_perc) = c("predicted.id", i)
  score = score %>% dplyr::left_join(meta_score, by = "predicted.id")
  perc = perc %>% dplyr::left_join(meta_perc, by = "predicted.id")
}

melt = melt(perc)
scores = melt(score)
colnames(melt) = c("predicted.id", "cluster", "perc")
colnames(scores) = c("predicted.id", "cluster", "score")
melt = melt %>% dplyr::left_join(scores, by = c("predicted.id", "cluster")) %>% dplyr::arrange(cluster)
melt[is.na(melt)] <- 0
melt$n = as.numeric(melt$cluster)

ggplot(melt, aes(x = cluster, y = perc, fill = predicted.id)) + geom_bar(position="fill", stat="identity") + scale_fill_colorblind()
ggsave(paste0(monocle_prefix, "distribution_perc_GSE263193.pdf"), width = 10, height = 10)

meta_prediction = predictions_AB001@meta.data
save(meta_prediction,
     file = paste0("/home/shared_folder/RData/meta_cluster_transfering.RData"))