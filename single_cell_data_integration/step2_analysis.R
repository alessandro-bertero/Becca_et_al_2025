library(dplyr)
library(Seurat)
library(monocle3)
library(ggplotify)
library(ggplot2)

load("/home/shared_folder/preprocessed.RData")
pal = c("#3283FE", "#FA0087", "#009E73", "#FBE426", "#56B4E9", "#FEAF16", "#DEA0FD", "#1CBE4F", "#F6222E", "#1CFFCE", "#325A9B", "#AA0DFE","#D55E00", "#2ED9FF", "#f0E442", "#1C8356", "#0072B2", "#CC79A7")

dir.create("/home/shared_folder/plots")
meta = rbind(meta1, meta2, meta3)
meta = meta %>% mutate(timepoint = ifelse(timepoint == "embryonic day 3", "day 3",
                                     ifelse(timepoint == "embryonic day 4", "day 4",
                                     ifelse(timepoint == "embryonic day 5", "day 5",
                                     ifelse(timepoint == "embryonic day 6", "day 6",
                                     ifelse(timepoint == "embryonic day 7", "day 7",
                                     ifelse(timepoint == "day16", "day 16",
                                     ifelse(timepoint == "5W", "week 5",
                                     ifelse(timepoint == "6W", "week 6",
                                     ifelse(timepoint == "7W", "week 7",
                                     ifelse(timepoint == "9W", "week 9",
                                     ifelse(timepoint == "10W", "week 10",
                                     ifelse(timepoint == "13W", "week 13",
                                     ifelse(timepoint == "15W", "week 15",
                                     ifelse(timepoint == "17W", "week 17",
                                     ifelse(timepoint == "20W", "week 20",
                                     ifelse(timepoint == "22W", "week 22",
                                     ifelse(timepoint == "23W", "week 23",
                                     ifelse(timepoint == "24W", "week 24", "")))))))))))))))))))

meta = meta %>% mutate(phenotype = ifelse(phenotype == "not applicable", "undifferentiated",
                                   ifelse(phenotype == "Advanced Mesoderm", "advanced mesoderm",
                                   ifelse(phenotype == "Primitive Streak", "primitive streak",
                                   ifelse(phenotype == "Emergent Mesoderm", "emergent mesoderm",
                                   ifelse(phenotype == "Nascent Mesoderm", "nascent mesoderm",
                                   ifelse(phenotype == "LA", "left atrium",
                                   ifelse(phenotype == "LV", "left ventricle",
                                   ifelse(phenotype == "RA", "right atrium",
                                   ifelse(phenotype == "RV", "right ventricle", "")))))))))))
rownames(meta) = meta$cell
genes = as.data.frame(rownames(counts))
colnames(genes) = "gene_short_name"
rownames(genes) = genes$gene_short_name

counts = cbind(counts1, counts2, counts3)

#RV
dir.create("/home/shared_folder/plots/RV")
meta_RV = meta_seurat %>% dplyr::filter(!phenotype %in% c("left atrium", "left ventricle", "right atrium"))
counts_RV = counts_seurat %>% dplyr::select(meta_RV$cell)
cds_RV <- new_cell_data_set(as.matrix(counts_RV),
                         cell_metadata = meta_RV,
                         gene_metadata = genes)
cds_RV <- preprocess_cds(cds_RV, method = "PCA", num_dim = 50, norm_method = "log")
cds_RV <- align_cds(cds_RV, alignment_group = "stage")
cds_RV <- reduce_dimension(cds_RV, umap.min_dist = 1, reduction_method = "UMAP")
plot_cells(cds_RV, color_cells_by = "stage", cell_size = 0.5)
ggsave("/home/shared_folder/plots/RV/UMAP_stages_RV.pdf", width = 8, height = 8)

cds_RV <- cluster_cells(cds_RV)
cds_RV <- learn_graph(cds_RV)
cds_RV <- order_cells(cds_RV)

plot_cells(cds_RV, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/RV/UMAP_pseudotime_RV.pdf", width = 8, height = 8)
plot_cells(cds_RV, genes = c("CTCF", "GATA4"), show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/RV/UMAP_genes_RV.pdf", width = 12, height = 8)
plot_cells(cds_RV, color_cells_by = "phenotype", show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/RV/UMAP_phenotype_RV.pdf", width = 8, height = 8)
plot_cells(cds_RV, color_cells_by = "timepoint", show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/RV/UMAP_timepoint_RV.pdf", width = 8, height = 8)

cds_subset_RV <- cds_RV[rowData(cds_RV)$gene_short_name %in% c("CTCF", "GATA4", "POU5F1", "TBXT", "TTN", "TNNI3", "EOMES"),]

plot_genes_in_pseudotime(
  cds_subset_RV,
  color_cells_by = "pseudotime",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/RV/genes_pseudotime_RV.pdf", width = 12, height = 12)

plot_genes_in_pseudotime(
  cds_subset_RV,
  color_cells_by = "stage",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/RV/genes_stage_RV.pdf", width = 12, height = 12)


plot_genes_in_pseudotime(
  cds_subset_RV,
  color_cells_by = "phenotype",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/RV/genes_phenotype_RV.pdf", width = 12, height = 12)


plot_genes_in_pseudotime(
  cds_subset_RV,
  color_cells_by = "timepoint",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/RV/genes_timepoint_RV.pdf", width = 12, height = 12)


#RA
dir.create("/home/shared_folder/plots/RA")
meta_RA = meta_seurat %>% dplyr::filter(!phenotype %in% c("left atrium", "left ventricle", "right ventricle"))
counts_RA = counts_seurat %>% dplyr::select(meta_RA$cell)
cds_RA <- new_cell_data_set(as.matrix(counts_RA),
                            cell_metadata = meta_RA,
                            gene_metadata = genes)
cds_RA <- preprocess_cds(cds_RA, method = "PCA", num_dim = 50, norm_method = "log")
cds_RA <- align_cds(cds_RA, alignment_group = "stage")
cds_RA <- reduce_dimension(cds_RA, umap.min_dist = 1, reduction_method = "UMAP")
plot_cells(cds_RA, color_cells_by = "stage", cell_size = 0.5)
ggsave("/home/shared_folder/plots/RA/UMAP_stages_RA.pdf", width = 8, height = 8)

cds_RA <- cluster_cells(cds_RA)
cds_RA <- learn_graph(cds_RA)
cds_RA <- order_cells(cds_RA)

plot_cells(cds_RA, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/RA/UMAP_pseudotime_RA.pdf", width = 8, height = 8)
plot_cells(cds_RA, genes = c("CTCF", "GATA4"), show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/RA/UMAP_genes_RA.pdf", width = 12, height = 8)
plot_cells(cds_RA, color_cells_by = "phenotype", show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/RA/UMAP_phenotype_RA.pdf", width = 8, height = 8)
plot_cells(cds_RA, color_cells_by = "timepoint", show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/RA/UMAP_timepoint_RA.pdf", width = 8, height = 8)

cds_subset_RA <- cds_RA[rowData(cds_RA)$gene_short_name %in% c("CTCF", "GATA4", "POU5F1", "TBXT", "TTN", "TNNI3", "EOMES"),]

plot_genes_in_pseudotime(
  cds_subset_RA,
  color_cells_by = "pseudotime",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/RA/genes_pseudotime_RA.pdf", width = 12, height = 12)

plot_genes_in_pseudotime(
  cds_subset_RA,
  color_cells_by = "stage",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/RA/genes_stage_RA.pdf", width = 12, height = 12)


plot_genes_in_pseudotime(
  cds_subset_RA,
  color_cells_by = "phenotype",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/RA/genes_phenotype_RA.pdf", width = 12, height = 12)


plot_genes_in_pseudotime(
  cds_subset_RA,
  color_cells_by = "timepoint",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/RA/genes_timepoint_RA.pdf", width = 12, height = 12)

#LV
dir.create("/home/shared_folder/plots/LV")
meta_LV = meta_seurat %>% dplyr::filter(!phenotype %in% c("right atrium", "left atrium", "right ventricle"))
counts_LV = counts_seurat %>% dplyr::select(meta_LV$cell)
cds_LV <- new_cell_data_set(as.matrix(counts_LV),
                            cell_metadata = meta_LV,
                            gene_metadata = genes)
cds_LV <- preprocess_cds(cds_LV, method = "PCA", num_dim = 50, norm_method = "log")
cds_LV <- align_cds(cds_LV, alignment_group = "stage")
cds_LV <- reduce_dimension(cds_LV, umap.min_dist = 1, reduction_method = "UMAP")
plot_cells(cds_LV, color_cells_by = "stage", cell_size = 1, label_cell_groups = FALSE)
ggsave("/home/shared_folder/plots/LV/UMAP_stages_LV.pdf", width = 8, height = 8)

cds_LV <- cluster_cells(cds_LV)
cds_LV <- learn_graph(cds_LV)
cds_LV <- order_cells(cds_LV)

plot_cells(cds_LV, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 1, label_cell_groups = FALSE)
ggsave("/home/shared_folder/plots/LV/UMAP_pseudotime_LV.pdf", width = 8, height = 8)
plot_cells(cds_LV, genes = c("CTCF", "GATA4"), show_trajectory_graph = FALSE, cell_size = 1, label_cell_groups = FALSE)
ggsave("/home/shared_folder/plots/LV/UMAP_genes_LV.pdf", width = 12, height = 8)
plot_cells(cds_LV, color_cells_by = "phenotype", show_trajectory_graph = FALSE, cell_size = 1, label_cell_groups = FALSE)
ggsave("/home/shared_folder/plots/LV/UMAP_phenotype_LV.pdf", width = 8, height = 8)
plot_cells(cds_LV, color_cells_by = "timepoint", show_trajectory_graph = FALSE, cell_size = 1, label_cell_groups = FALSE)
ggsave("/home/shared_folder/plots/LV/UMAP_timepoint_LV.pdf", width = 8, height = 8)

cds_subset_LV <- cds_LV[rowData(cds_LV)$gene_short_name %in% c("CTCF", "GATA4", "POU5F1", "TBXT", "TTN", "TNNI3", "EOMES"),]

plot_genes_in_pseudotime(
  cds_subset_LV,
  color_cells_by = "pseudotime",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/LV/genes_pseudotime_LV.pdf", width = 12, height = 12)

plot_genes_in_pseudotime(
  cds_subset_LV,
  color_cells_by = "stage",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/LV/genes_stage_LV.pdf", width = 12, height = 12)


plot_genes_in_pseudotime(
  cds_subset_LV,
  color_cells_by = "phenotype",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/LV/genes_phenotype_LV.pdf", width = 12, height = 12)


plot_genes_in_pseudotime(
  cds_subset_LV,
  color_cells_by = "timepoint",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/LV/genes_timepoint_LV.pdf", width = 12, height = 12)

#LA
dir.create("/home/shared_folder/plots/LA")
meta_LA = meta_seurat %>% dplyr::filter(!phenotype %in% c("right atrium", "left ventricle", "right ventricle"))
counts_LA = counts_seurat %>% dplyr::select(meta_LA$cell)
cds_LA <- new_cell_data_set(as.matrix(counts_LA),
                            cell_metadata = meta_LA,
                            gene_metadata = genes)
cds_LA <- preprocess_cds(cds_LA, method = "PCA", num_dim = 50, norm_method = "log")
cds_LA <- align_cds(cds_LA, alignment_group = "stage")
cds_LA <- reduce_dimension(cds_LA, umap.min_dist = 1, reduction_method = "UMAP")
plot_cells(cds_LA, color_cells_by = "stage", cell_size = 0.5)
ggsave("/home/shared_folder/plots/LA/UMAP_stages_LA.pdf", width = 8, height = 8)

cds_LA <- cluster_cells(cds_LA)
cds_LA <- learn_graph(cds_LA)
cds_LA <- order_cells(cds_LA)

plot_cells(cds_LA, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/LA/UMAP_pseudotime_LA.pdf", width = 8, height = 8)
plot_cells(cds_LA, genes = c("CTCF", "GATA4"), show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/LA/UMAP_genes_LA.pdf", width = 12, height = 8)
plot_cells(cds_LA, color_cells_by = "phenotype", show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/LA/UMAP_phenotype_LA.pdf", width = 8, height = 8)
plot_cells(cds_LA, color_cells_by = "timepoint", show_trajectory_graph = FALSE, cell_size = 0.5)
ggsave("/home/shared_folder/plots/LA/UMAP_timepoint_LA.pdf", width = 8, height = 8)

cds_subset_LA <- cds_LA[rowData(cds_LA)$gene_short_name %in% c("CTCF", "GATA4", "POU5F1", "TBXT", "TTN", "TNNI3", "EOMES"),]

plot_genes_in_pseudotime(
  cds_subset_LA,
  color_cells_by = "pseudotime",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/LA/genes_pseudotime_LA.pdf", width = 12, height = 12)

plot_genes_in_pseudotime(
  cds_subset_LA,
  color_cells_by = "stage",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/LA/genes_stage_LA.pdf", width = 12, height = 12)


plot_genes_in_pseudotime(
  cds_subset_LA,
  color_cells_by = "phenotype",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/LA/genes_phenotype_LA.pdf", width = 12, height = 12)


plot_genes_in_pseudotime(
  cds_subset_LA,
  color_cells_by = "timepoint",
  min_expr = 0.5
)
ggsave("/home/shared_folder/plots/LA/genes_timepoint_LA.pdf", width = 12, height = 12)
