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

##chamber separation
#AT
cds_AT <- cds[, 
              colData(cds) %>%
                subset(
                  cardiac_chamber == "AT"
                ) %>%
                row.names
]
cds_AT <- reduce_dimension(cds_AT, reduction_method = "UMAP", umap.min_dist = 0.05)
cds_AT <- cluster_cells(cds_AT, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)

p1 = plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nCount_RNA", cell_size = 0.7)
p2 = plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_clusters_AT.pdf"),
       width = 20, height = 8)

p1 = plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 0.7)
p2 = plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_clusters_Phase_AT.pdf"),
       width = 20, height = 8)

cds_AT <- learn_graph(cds_AT)
plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", genes = c("TNNT2", "TTN", "EOMES", "POU5F1"), cell_size = 0.7)
ggsave(filename = paste0(monocle_prefix, "UMAP_diff_markers_AT.pdf"),
       width = 20, height = 20)

cds_AT <- order_cells(cds_AT)
p1 = plot_cells(cds_AT,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5,
                cell_size = 0.7)
p1
ggsave(filename = paste0(monocle_prefix, "UMAP_pseudotime_AT.pdf"),
       width = 12, height = 12)

pseudotime = data.frame(p1$data$sample_name, libID = p1$data$sample_name, pseudotime = p1$data$cell_color, row.names = 1)
colData(cds_AT)$pseudotime = pseudotime$pseudotime

#RV
cds_RV <- cds[, 
              colData(cds) %>%
                subset(
                  cardiac_chamber == "RV"
                ) %>%
                row.names
]
cds_RV <- reduce_dimension(cds_RV, reduction_method = "UMAP", umap.min_dist = 0.05)
cds_RV <- cluster_cells(cds_RV, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)

p1 = plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nCount_RNA", cell_size = 0.7)
p2 = plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_clusters_RV.pdf"),
       width = 20, height = 8)

p1 = plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 0.7)
p2 = plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_clusters_Phase_RV.pdf"),
       width = 20, height = 8)

cds_RV <- learn_graph(cds_RV)
plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", genes = c("TNNT2", "TTN", "EOMES", "POU5F1"), cell_size = 0.7)
ggsave(filename = paste0(monocle_prefix, "UMAP_diff_markers_RV.pdf"),
       width = 20, height = 20)

cds_RV <- order_cells(cds_RV)
p1 = plot_cells(cds_RV,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5,
                cell_size = 0.7)
p1
ggsave(filename = paste0(monocle_prefix, "UMAP_pseudotime_RV.pdf"),
       width = 12, height = 12)

pseudotime = data.frame(p1$data$sample_name, libID = p1$data$sample_name, pseudotime = p1$data$cell_color, row.names = 1)
colData(cds_RV)$pseudotime = pseudotime$pseudotime

#LV
cds_LV <- cds[, 
              colData(cds) %>%
                subset(
                  cardiac_chamber == "LV"
                ) %>%
                row.names
]
cds_LV <- reduce_dimension(cds_LV, reduction_method = "UMAP", umap.min_dist = 0.05)
cds_LV <- cluster_cells(cds_LV, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)

p1 = plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nCount_RNA", cell_size = 0.7)
p2 = plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_clusters_LV.pdf"),
       width = 20, height = 8)

p1 = plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 0.7)
p2 = plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE, cell_size = 0.7)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_clusters_Phase_LV.pdf"),
       width = 20, height = 8)

cds_LV <- learn_graph(cds_LV, use_partition = FALSE)
plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", genes = c("TNNT2", "TTN", "EOMES", "POU5F1"), cell_size = 0.7)
ggsave(filename = paste0(monocle_prefix, "UMAP_diff_markers_LV.pdf"),
       width = 20, height = 20)

cds_LV <- order_cells(cds_LV)
p1 = plot_cells(cds_LV,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=1.5,
                cell_size = 0.7)
p1
ggsave(filename = paste0(monocle_prefix, "UMAP_pseudotime_LV.pdf"),
       width = 12, height = 12)

pseudotime = data.frame(p1$data$sample_name, libID = p1$data$sample_name, pseudotime = p1$data$cell_color, row.names = 1)
colData(cds_LV)$pseudotime = pseudotime$pseudotime

#save
save(cds_AT, cds_RV, cds_LV,
     file = paste0("/home/shared_folder/RData/cds_chambers.RData"))

##cluster classification
#AT
marker_cluster_AT <- top_markers(cds_AT, group_cells_by="cluster", 
                                 reference_cells=1000, cores=8)
top_marker_cluster_AT <- marker_cluster_AT %>%
  dplyr::filter(fraction_expressing >= 0.10) %>%
  dplyr::filter(marker_test_q_value < 0.05) %>%
  dplyr::group_by(cell_group) %>%
  dplyr::top_n(10, pseudo_R2) %>%
  dplyr::arrange(cell_group)

AT_cluster_list = unique(top_marker_cluster_AT %>% pull(gene_id))

plot_genes_by_group(cds_AT,
                    AT_cluster_list,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)
ggsave(filename = paste0(monocle_prefix, "UMAP_AT_cluster_markers.pdf"),
       width = 10, height = 20)

for(i in top_marker_cluster_AT$cell_group){
  GO = ensembldb::select(org.Hs.eg.db, keys = (top_marker_cluster_AT %>% dplyr::filter(cell_group == i))$gene_id, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  go = goana(GO$ENTREZID, species = "Hs", convert = T)
  res = topGO(go, ontology = "BP", number = 30, truncate.term = 50)
  ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
    labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of BP of cluster ", i, " genes"))
  ggsave(filename = paste0(monocle_prefix, "GO_AT_cluster ", i, "_markers.pdf"),
         width = 8, height = 12)
}

generate_garnett_marker_file(top_marker_cluster_AT, file="/home/shared_folder/marker_file_AT.txt")

#RV
marker_cluster_RV <- top_markers(cds_RV, group_cells_by="cluster", 
                                 reference_cells=1000, cores=8)
top_marker_cluster_RV <- marker_cluster_RV %>%
  dplyr::filter(fraction_expressing >= 0.10) %>%
  dplyr::filter(marker_test_q_value < 0.05) %>%
  dplyr::group_by(cell_group) %>%
  dplyr::top_n(10, pseudo_R2) %>%
  dplyr::arrange(cell_group)

RV_cluster_list = unique(top_marker_cluster_RV %>% pull(gene_id))

plot_genes_by_group(cds_RV,
                    RV_cluster_list,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)
ggsave(filename = paste0(monocle_prefix, "UMAP_RV_cluster_markers.pdf"),
       width = 10, height = 20)

for(i in top_marker_cluster_RV$cell_group){
  GO = ensembldb::select(org.Hs.eg.db, keys = (top_marker_cluster_RV %>% dplyr::filter(cell_group == i))$gene_id, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  go = goana(GO$ENTREZID, species = "Hs", convert = T)
  res = topGO(go, ontology = "BP", number = 30, truncate.term = 50)
  ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
    labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of BP of cluster ", i, " genes"))
  ggsave(filename = paste0(monocle_prefix, "GO_RV_cluster ", i, "_markers.pdf"),
         width = 8, height = 12)
}

generate_garnett_marker_file(top_marker_cluster_RV, file="/home/shared_folder/marker_file_RV.txt")

#LV
marker_cluster_LV <- top_markers(cds_LV, group_cells_by="cluster", 
                                 reference_cells=1000, cores=8)
top_marker_cluster_LV <- marker_cluster_LV %>%
  dplyr::filter(fraction_expressing >= 0.10) %>%
  dplyr::filter(marker_test_q_value < 0.05) %>%
  dplyr::group_by(cell_group) %>%
  dplyr::top_n(10, pseudo_R2) %>%
  dplyr::arrange(cell_group)

LV_cluster_list = unique(top_marker_cluster_LV %>% pull(gene_id))

plot_genes_by_group(cds_LV,
                    LV_cluster_list,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)
ggsave(filename = paste0(monocle_prefix, "UMAP_LV_cluster_markers.pdf"),
       width = 10, height = 20)

for(i in top_marker_cluster_LV$cell_group){
  GO = ensembldb::select(org.Hs.eg.db, keys = (top_marker_cluster_LV %>% dplyr::filter(cell_group == i))$gene_id, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  go = goana(GO$ENTREZID, species = "Hs", convert = T)
  res = topGO(go, ontology = "BP", number = 30, truncate.term = 50)
  ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
    labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of BP of cluster ", i, " genes"))
  ggsave(filename = paste0(monocle_prefix, "GO_LV_cluster ", i, "_markers.pdf"),
         width = 8, height = 12)
}

generate_garnett_marker_file(top_marker_cluster_LV, file="/home/shared_folder/marker_file_LV.txt")

#save
save(cds_AT, cds_RV, cds_LV, top_marker_cluster_AT, top_marker_cluster_RV, top_marker_cluster_LV,
     file = paste0(monocle_prefix, "_cds_chambers_annotated.RData"))


##differential gene expressiona anlysis
GSEA = c("GO:0061337", "GO:0086003", "GO:0055007", "GO:0010657", "GO:0051450", "GO:0030029", "GO:0000278")

#AT
p1 = plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", cell_size = 0.7)
p1
ggsave(filename = paste0(monocle_prefix, "clusters_AT.pdf"),
       width = 8, height = 8)

clusters = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
colData(cds_AT)$cluster_AT = clusters$cluster

for(i in unique(clusters$cluster)){
  clusters = clusters %>% dplyr::mutate(cfr = ifelse(cluster == i, i, 0))
  colData(cds_AT)$cfr = clusters$cfr
  gene_fits <- fit_models(cds_AT, expression_family = "poisson", model_formula_str = "~cfr")
  fit_coefs <- coefficient_table(gene_fits)
  fit_coefs = fit_coefs %>% dplyr::filter(term != "(Intercept)")
  diff = fit_coefs %>% dplyr::filter (q_value < 0.05) %>%
    dplyr::select(gene_short_name, term, q_value, estimate, test_val)
  
  diff_pos = diff %>% dplyr::filter(estimate > 0) %>% dplyr::arrange(abs(test_val))
  top_pos = diff_pos[c(1:12),]
  cds_subset <- cds_AT[
    cds_AT@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_pos$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_AT", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_AT_cluster_", i, "_pos.pdf"),
         width = 10, height = 24)
  plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", genes = c(diff_pos$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_AT_cluster_", i, "_pos.pdf"),
         width = 20, height = 8)
  
  diff_neg = diff %>% dplyr::filter(estimate < 0) %>% dplyr::arrange(abs(test_val))
  top_neg = diff_neg[c(1:12),]
  cds_subset <- cds_AT[
    cds_AT@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_neg$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_AT", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_AT_cluster_", i, "_neg.pdf"),
         width = 10, height = 24)
  plot_cells(cds_AT, norm_method = "log", reduction_method = "UMAP", genes = c(diff_neg$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_AT_cluster_", i, "_neg.pdf"),
         width = 20, height = 8)
  
  diff = as.data.frame(diff) %>% dplyr::mutate(color = ifelse(estimate > 0, "up", ifelse(estimate < 0, "down", "other")))
  ggplot(diff, aes(x = estimate, y = -log10(q_value))) + geom_point(aes(color = color, size = 0.5, alpha = 0.5)) +
    scale_colour_manual(values = c("#FA0087", "#3283FE")) +
    ggrepel::geom_text_repel(aes(label = gene_short_name), max.overlaps = 20, size = 2, hjust = 1.2)
  ggsave(filename = paste0(monocle_prefix, "volcano_plot_cluster_", i, "_AT.pdf"),
         width = 20, height = 8)
  
  genes = diff$estimate
  names(genes) = diff$gene_short_name
  genes = sort(genes, decreasing = TRUE)
  gse = gseGO(geneList = genes, 
              ont = "BP", 
              keyType = "SYMBOL",
              pvalueCutoff = 1, 
              OrgDb = org.Hs.eg.db, 
              pAdjustMethod = "BH",
              by = "fgsea")
  res_gsea = gse@result
  res_gsea$n = c(1:length(rownames(res_gsea)))
  write.csv(res_gsea, file = paste0(monocle_prefix, "gsea_AT_cluster_", i, ".csv"))
  dotplot(gse, showCategory = 30, split = ".sign") + facet_grid(.~.sign)
  ggsave(filename = paste0(monocle_prefix, "dotplot_AT_cluster_", i, ".pdf"), width = 14, height = 20)
  
  for(j in GSEA){if(j %in% res_gsea$ID){
    res_j = res_gsea %>% dplyr::filter(ID == j);
    n = res_j$n;
    gseaplot(gse, by ="runningScore", color.line = "#009E73", title = paste0(gse$Description[n], " p = ", gse$p.adjust[n]), geneSetID = n) + theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7));
    ggsave(filename = paste0(monocle_prefix, "dotplot_AT_cluster_", i, "_", gse$Description[n], ".pdf"), width = 14, height = 8);
  }}
  
  write.csv(as.data.frame(fit_coefs) %>% dplyr::select(-c(model, model_summary, status)), paste0(monocle_prefix, "dge_AT_cluster_", i, ".csv"))
}


#RV
p1 = plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", cell_size = 0.7)
p1
ggsave(filename = paste0(monocle_prefix, "clusters_RV.pdf"),
       width = 8, height = 8)

clusters = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
colData(cds_RV)$cluster_RV = clusters$cluster

for(i in unique(clusters$cluster)){
  clusters = clusters %>% dplyr::mutate(cfr = ifelse(cluster == i, i, 0))
  colData(cds_RV)$cfr = clusters$cfr
  gene_fits <- fit_models(cds_RV, expression_family = "poisson", model_formula_str = "~cfr")
  fit_coefs <- coefficient_table(gene_fits)
  fit_coefs = fit_coefs %>% dplyr::filter(term != "(Intercept)")
  diff = fit_coefs %>% dplyr::filter (q_value < 0.05) %>%
    dplyr::select(gene_short_name, term, q_value, estimate, test_val)
  
  diff_pos = diff %>% dplyr::filter(estimate > 0) %>% dplyr::arrange(abs(test_val))
  top_pos = diff_pos[c(1:12),]
  cds_subset <- cds_RV[
    cds_RV@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_pos$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_RV", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_RV_cluster_", i, "_pos.pdf"),
         width = 10, height = 24)
  plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", genes = c(diff_pos$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_RV_cluster_", i, "_pos.pdf"),
         width = 20, height = 8)
  
  diff_neg = diff %>% dplyr::filter(estimate < 0) %>% dplyr::arrange(abs(test_val))
  top_neg = diff_neg[c(1:12),]
  cds_subset <- cds_RV[
    cds_RV@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_neg$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_RV", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_RV_cluster_", i, "_neg.pdf"),
         width = 10, height = 24)
  plot_cells(cds_RV, norm_method = "log", reduction_method = "UMAP", genes = c(diff_neg$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_RV_cluster_", i, "_neg.pdf"),
         width = 20, height = 8)
  
  diff = as.data.frame(diff) %>% dplyr::mutate(color = ifelse(estimate > 0, "up", ifelse(estimate < 0, "down", "other")))
  ggplot(diff, aes(x = estimate, y = -log10(q_value))) + geom_point(aes(color = color, size = 0.5, alpha = 0.5)) +
    scale_colour_manual(values = c("#FA0087", "#3283FE")) +
    ggrepel::geom_text_repel(aes(label = gene_short_name), max.overlaps = 20, size = 2, hjust = 1.2)
  ggsave(filename = paste0(monocle_prefix, "volcano_plot_cluster_", i, "_RV.pdf"),
         width = 20, height = 8)
  
  genes = diff$estimate
  names(genes) = diff$gene_short_name
  genes = sort(genes, decreasing = TRUE)
  gse = gseGO(geneList = genes, 
              ont = "BP", 
              keyType = "SYMBOL",
              pvalueCutoff = 1, 
              OrgDb = org.Hs.eg.db, 
              pAdjustMethod = "BH",
              by = "fgsea")
  res_gsea = gse@result
  res_gsea$n = c(1:length(rownames(res_gsea)))
  write.csv(res_gsea, file = paste0(monocle_prefix, "gsea_RV_cluster_", i, ".csv"))
  dotplot(gse, showCategory = 30, split = ".sign") + facet_grid(.~.sign)
  ggsave(filename = paste0(monocle_prefix, "dotplot_RV_cluster_", i, ".pdf"), width = 14, height = 20)
  
  for(j in GSEA){if(j %in% res_gsea$ID){
    res_j = res_gsea %>% dplyr::filter(ID == j);
    n = res_j$n;
    gseaplot(gse, by ="runningScore", color.line = "#009E73", title = paste0(gse$Description[n], " p = ", gse$p.adjust[n]), geneSetID = n) + theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7));
    ggsave(filename = paste0(monocle_prefix, "dotplot_RV_cluster_", i, "_", gse$Description[n], ".pdf"), width = 14, height = 8);
  }}
  
  write.csv(as.data.frame(fit_coefs) %>% dplyr::select(-c(model, model_summary, status)), paste0(monocle_prefix, "dge_RV_cluster_", i, ".csv"))
}

#LV
p1 = plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", cell_size = 0.7)
p1
ggsave(filename = paste0(monocle_prefix, "clusters_LV.pdf"),
       width = 8, height = 8)

clusters = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
colData(cds_LV)$cluster_LV = clusters$cluster

for(i in unique(clusters$cluster)){
  clusters = clusters %>% dplyr::mutate(cfr = ifelse(cluster == i, i, 0))
  colData(cds_LV)$cfr = clusters$cfr
  gene_fits <- fit_models(cds_LV, expression_family = "poisson", model_formula_str = "~cfr")
  fit_coefs <- coefficient_table(gene_fits)
  fit_coefs = fit_coefs %>% dplyr::filter(term != "(Intercept)")
  diff = fit_coefs %>% dplyr::filter (q_value < 0.05) %>%
    dplyr::select(gene_short_name, term, q_value, estimate, test_val)
  
  diff_pos = diff %>% dplyr::filter(estimate > 0) %>% dplyr::arrange(abs(test_val))
  top_pos = diff_pos[c(1:12),]
  cds_subset <- cds_LV[
    cds_LV@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_pos$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_LV", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_LV_cluster_", i, "_pos.pdf"),
         width = 10, height = 24)
  plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", genes = c(diff_pos$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_LV_cluster_", i, "_pos.pdf"),
         width = 20, height = 8)
  
  diff_neg = diff %>% dplyr::filter(estimate < 0) %>% dplyr::arrange(abs(test_val))
  top_neg = diff_neg[c(1:12),]
  cds_subset <- cds_LV[
    cds_LV@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_neg$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_LV", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_LV_cluster_", i, "_neg.pdf"),
         width = 10, height = 24)
  plot_cells(cds_LV, norm_method = "log", reduction_method = "UMAP", genes = c(diff_neg$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_LV_cluster_", i, "_neg.pdf"),
         width = 20, height = 8)
  
  diff = as.data.frame(diff) %>% dplyr::mutate(color = ifelse(estimate > 0, "up", ifelse(estimate < 0, "down", "other")))
  ggplot(diff, aes(x = estimate, y = -log10(q_value))) + geom_point(aes(color = color, size = 0.5, alpha = 0.5)) +
    scale_colour_manual(values = c("#FA0087", "#3283FE")) +
    ggrepel::geom_text_repel(aes(label = gene_short_name), max.overlaps = 20, size = 2, hjust = 1.2)
  ggsave(filename = paste0(monocle_prefix, "volcano_plot_cluster_", i, "_LV.pdf"),
         width = 20, height = 8)
  
  genes = diff$estimate
  names(genes) = diff$gene_short_name
  genes = sort(genes, decreasing = TRUE)
  gse = gseGO(geneList = genes, 
              ont = "BP", 
              keyType = "SYMBOL",
              pvalueCutoff = 1, 
              OrgDb = org.Hs.eg.db, 
              pAdjustMethod = "BH",
              by = "fgsea")
  res_gsea = gse@result
  res_gsea$n = c(1:length(rownames(res_gsea)))
  write.csv(res_gsea, file = paste0(monocle_prefix, "gsea_LV_cluster_", i, ".csv"))
  dotplot(gse, showCategory = 30, split = ".sign") + facet_grid(.~.sign)
  ggsave(filename = paste0(monocle_prefix, "dotplot_LV_cluster_", i, ".pdf"), width = 14, height = 20)
  
  for(j in GSEA){if(j %in% res_gsea$ID){
    res_j = res_gsea %>% dplyr::filter(ID == j);
    n = res_j$n;
    gseaplot(gse, by ="runningScore", color.line = "#009E73", title = paste0(gse$Description[n], " p = ", gse$p.adjust[n]), geneSetID = n) + theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7));
    ggsave(filename = paste0(monocle_prefix, "dotplot_LV_cluster_", i, "_", gse$Description[n], ".pdf"), width = 14, height = 8);
  }}
  
  write.csv(as.data.frame(fit_coefs) %>% dplyr::select(-c(model, model_summary, status)), paste0(monocle_prefix, "dge_LV_cluster_", i, ".csv"))
}




