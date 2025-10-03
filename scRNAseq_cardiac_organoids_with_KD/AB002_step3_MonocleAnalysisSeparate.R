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
pal = c("#3283FE", "#FA0087", "#009E73", "#FBE426", "#56B4E9", "#FEAF16", "#DEA0FD", "#1CBE4F", "#F6222E", "#1CFFCE", "#325A9B", "#AA0DFE","#D55E00", "#2ED9FF", "#f0E442", "#1C8356", "#0072B2", "#CC79A7")
samples = c("#FFB6DB", "#FA0087", "#B6DBFF", "#3283FE", "#808080", "#000000")

#load data
load("/home/shared_folder/RData/cds_filtered_all.RData")

#day separation
cds_day4 <- cds_subset[, 
                       colData(cds_subset) %>%
                         subset(
                           day == "D4"
                         ) %>%
                         row.names
]
cds_day4 <- reduce_dimension(cds_day4, reduction_method = "UMAP", umap.min_dist = 0.2)
cds_day4 <- cluster_cells(cds_day4, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)

p1 = plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nCount_RNA")
p2 = plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_clusters_day4.pdf"),
       width = 20, height = 8)

cds_day4 <- learn_graph(cds_day4)
plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = c("TNNT2", "TTN", "TBXT", "EOMES"))
ggsave(filename = paste0(monocle_prefix, "UMAP_genes_day4.pdf"),
       width = 20, height = 20)

cds_day7 <- cds_subset[, 
                       colData(cds_subset) %>%
                         subset(
                           day == "D7"
                         ) %>%
                         row.names
]
cds_day7 <- reduce_dimension(cds_day7, reduction_method = "UMAP", umap.min_dist = 0.2)
cds_day7 <- cluster_cells(cds_day7, reduction_method = "UMAP", random_seed = 1234, partition_qval = 0.05)

p1 = plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "nCount_RNA")
p2 = plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", label_cell_groups = FALSE)
p1 + p2
ggsave(filename = paste0(monocle_prefix, "UMAP_clusters_day7.pdf"),
       width = 20, height = 8)

cds_day7 <- learn_graph(cds_day7)
plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = c("TNNT2", "TTN", "TBXT", "EOMES"))
ggsave(filename = paste0(monocle_prefix, "UMAP_genes_day7.pdf"),
       width = 20, height = 20)

##cluster classification
#day4
marker_cluster_day4 <- top_markers(cds_day4, group_cells_by="cluster", 
                                   reference_cells=1000, cores=8)
top_marker_cluster_day4 <- marker_cluster_day4 %>%
  dplyr::filter(fraction_expressing >= 0.10) %>%
  dplyr::filter(marker_test_q_value < 0.05) %>%
  dplyr::group_by(cell_group) %>%
  dplyr::top_n(10, pseudo_R2) %>%
  dplyr::arrange(cell_group)

day4_cluster_list = unique(top_marker_cluster_day4 %>% pull(gene_id))

plot_genes_by_group(cds_day4,
                    day4_cluster_list,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)
ggsave(filename = paste0(monocle_prefix, "UMAP_day4_cluster_markers.pdf"),
       width = 10, height = 20)

for(i in top_marker_cluster_day4$cell_group){
  GO = ensembldb::select(org.Hs.eg.db, keys = (top_marker_cluster_day4 %>% dplyr::filter(cell_group == i))$gene_id, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  go = goana(GO$ENTREZID, species = "Hs", convert = T)
  res = topGO(go, ontology = "BP", number = 30, truncate.term = 50)
  ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
    labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of BP of cluster ", i, " genes"))
  ggsave(filename = paste0(monocle_prefix, "GO_day4_cluster ", i, "_markers.pdf"),
         width = 8, height = 12)
}

generate_garnett_marker_file(top_marker_cluster_day4, file="/home/shared_folder/marker_file_day4.txt")

#day7
marker_cluster_day7 <- top_markers(cds_day7, group_cells_by="cluster", 
                                   reference_cells=1000, cores=8)
top_marker_cluster_day7 <- marker_cluster_day7 %>%
  dplyr::filter(fraction_expressing >= 0.10) %>%
  dplyr::filter(marker_test_q_value < 0.05) %>%
  dplyr::group_by(cell_group) %>%
  dplyr::top_n(10, pseudo_R2) %>%
  dplyr::arrange(cell_group)

day7_cluster_list = unique(top_marker_cluster_day7 %>% pull(gene_id))

plot_genes_by_group(cds_day7,
                    day7_cluster_list,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)
ggsave(filename = paste0(monocle_prefix, "UMAP_day7_cluster_markers.pdf"),
       width = 10, height = 20)

for(i in top_marker_cluster_day7$cell_group){
  GO = ensembldb::select(org.Hs.eg.db, keys = (top_marker_cluster_day7 %>% dplyr::filter(cell_group == i))$gene_id, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  go = goana(GO$ENTREZID, species = "Hs", convert = T)
  res = topGO(go, ontology = "BP", number = 30, truncate.term = 50)
  ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
    labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of BP of cluster ", i, " genes"))
  ggsave(filename = paste0(monocle_prefix, "GO_day7_cluster ", i, "_markers.pdf"),
         width = 8, height = 12)
}

generate_garnett_marker_file(top_marker_cluster_day7, file="/home/shared_folder/marker_file_day7.txt")

#select pseudotime root
cds_day4 <- order_cells(cds_day4)
p1 = plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "pseudotime")
p1
pseudotime = data.frame(p1$data$sample_name, libID = p1$data$sample_name, pseudotime = p1$data$cell_color, row.names = 1)
colData(cds_day4)$pseudotime = pseudotime$pseudotime

cds_day7 <- order_cells(cds_day7)
p1 = plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "pseudotime")
p1
pseudotime = data.frame(p1$data$sample_name, libID = p1$data$sample_name, pseudotime = p1$data$cell_color, row.names = 1)
colData(cds_day7)$pseudotime = pseudotime$pseudotime


##differential gene expression analysis
GSEA = c("GO:0061337", "GO:0086003", "GO:0055007", "GO:0010657", "GO:0051450", "GO:0030029", "GO:0000278")

#day4
p1 = plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster")
p1
ggsave(filename = paste0(monocle_prefix, "clusters_day4.pdf"),
       width = 20, height = 8)

clusters = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
colData(cds_day4)$cluster_day4 = clusters$cluster

for(i in unique(clusters$cluster)){
  clusters = clusters %>% dplyr::mutate(cfr = ifelse(cluster == i, i, 0))
  colData(cds_day4)$cfr = clusters$cfr
  gene_fits <- fit_models(cds_day4, expression_family = "poisson", model_formula_str = "~cfr")
  fit_coefs <- coefficient_table(gene_fits)
  fit_coefs = fit_coefs %>% dplyr::filter(term != "(Intercept)")
  diff = fit_coefs %>% dplyr::filter (q_value < 0.05) %>%
    dplyr::select(gene_short_name, term, q_value, estimate, test_val)
  
  diff_pos = diff %>% dplyr::filter(estimate > 0) %>% dplyr::arrange(abs(test_val))
  top_pos = diff_pos[c(1:12),]
  cds_subset <- cds_day4[
    cds_day4@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_pos$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_day4", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_day4_cluster_", i, "_pos.pdf"),
         width = 10, height = 24)
  plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = c(diff_pos$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_day4_cluster_", i, "_pos.pdf"),
         width = 20, height = 8)
  
  diff_neg = diff %>% dplyr::filter(estimate < 0) %>% dplyr::arrange(abs(test_val))
  top_neg = diff_neg[c(1:12),]
  cds_subset <- cds_day4[
    cds_day4@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_neg$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_day4", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_day4_cluster_", i, "_neg.pdf"),
         width = 10, height = 24)
  plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = c(diff_neg$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_day4_cluster_", i, "_neg.pdf"),
         width = 20, height = 8)
  
  diff = as.data.frame(diff) %>% dplyr::mutate(color = ifelse(estimate > 0, "up", ifelse(estimate < 0, "down", "other")))
  ggplot(diff, aes(x = estimate, y = -log10(q_value))) + geom_point(aes(color = color, size = 0.5, alpha = 0.5)) +
    scale_colour_manual(values = c("#FA0087", "#3283FE")) +
    ggrepel::geom_text_repel(aes(label = gene_short_name), max.overlaps = 20, size = 2, hjust = 1.2)
  ggsave(filename = paste0(monocle_prefix, "volcano_plot_cluster_", i, "_day4.pdf"),
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
  write.csv(res_gsea, file = paste0(monocle_prefix, "gsea_day4_cluster_", i, ".csv"))
  dotplot(gse, showCategory = 30, split = ".sign") + facet_grid(.~.sign)
  ggsave(filename = paste0(monocle_prefix, "dotplot_day4_cluster_", i, ".pdf"), width = 14, height = 20)
  
  for(j in GSEA){if(j %in% res_gsea$ID){
    res_j = res_gsea %>% dplyr::filter(ID == j);
    n = res_j$n;
    gseaplot(gse, by ="runningScore", color.line = "#009E73", title = paste0(gse$Description[n], " p = ", gse$p.adjust[n]), geneSetID = n) + theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7));
    ggsave(filename = paste0(monocle_prefix, "dotplot_day4_cluster_", i, "_", gse$Description[n], ".pdf"), width = 14, height = 8);
  }}
  
  write.csv(as.data.frame(fit_coefs) %>% dplyr::select(-c(model, model_summary, status)), paste0(monocle_prefix, "dge_day4_cluster_", i, ".csv"))
}

#day7
p1 = plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster")
p1
ggsave(filename = paste0(monocle_prefix, "clusters_day7.pdf"),
       width = 20, height = 8)

clusters = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
colData(cds_day7)$cluster_day7 = clusters$cluster

for(i in unique(clusters$cluster)){
  clusters = clusters %>% dplyr::mutate(cfr = ifelse(cluster == i, i, 0))
  colData(cds_day7)$cfr = clusters$cfr
  gene_fits <- fit_models(cds_day7, expression_family = "poisson", model_formula_str = "~cfr")
  fit_coefs <- coefficient_table(gene_fits)
  fit_coefs = fit_coefs %>% dplyr::filter(term != "(Intercept)")
  diff = fit_coefs %>% dplyr::filter (q_value < 0.05) %>%
    dplyr::select(gene_short_name, term, q_value, estimate, test_val)
  
  diff_pos = diff %>% dplyr::filter(estimate > 0) %>% dplyr::arrange(abs(test_val))
  top_pos = diff_pos[c(1:12),]
  cds_subset <- cds_day7[
    cds_day7@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_pos$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_day7", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_day7_cluster_", i, "_pos.pdf"),
         width = 10, height = 24)
  plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = c(diff_pos$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_day7_cluster_", i, "_pos.pdf"),
         width = 20, height = 8)
  
  diff_neg = diff %>% dplyr::filter(estimate < 0) %>% dplyr::arrange(abs(test_val))
  top_neg = diff_neg[c(1:12),]
  cds_subset <- cds_day7[
    cds_day7@rowRanges@elementMetadata@listData[["gene_short_name"]] %in% top_neg$gene_short_name
  ]
  plot_genes_violin(cds_subset, group_cells_by="cluster_day7", ncol=2, log_scale = TRUE)
  ggsave(filename = paste0(monocle_prefix, "violin_dge_day7_cluster_", i, "_neg.pdf"),
         width = 10, height = 24)
  plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = c(diff_neg$gene_short_name[1:4]), show_trajectory_graph = FALSE)
  ggsave(filename = paste0(monocle_prefix, "markers_dge_day7_cluster_", i, "_neg.pdf"),
         width = 20, height = 8)
  
  diff = as.data.frame(diff) %>% dplyr::mutate(color = ifelse(estimate > 0, "up", ifelse(estimate < 0, "down", "other")))
  ggplot(diff, aes(x = estimate, y = -log10(q_value))) + geom_point(aes(color = color, size = 0.5, alpha = 0.5)) +
    scale_colour_manual(values = c("#FA0087", "#3283FE")) +
    ggrepel::geom_text_repel(aes(label = gene_short_name), max.overlaps = 20, size = 2, hjust = 1.2)
  ggsave(filename = paste0(monocle_prefix, "volcano_plot_cluster_", i, "_day7.pdf"),
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
  write.csv(res_gsea, file = paste0(monocle_prefix, "gsea_day7_cluster_", i, ".csv"))
  dotplot(gse, showCategory = 30, split = ".sign") + facet_grid(.~.sign)
  ggsave(filename = paste0(monocle_prefix, "dotplot_day7_cluster_", i, ".pdf"), width = 14, height = 20)
  
  for(j in GSEA){if(j %in% res_gsea$ID){
    res_j = res_gsea %>% dplyr::filter(ID == j);
    n = res_j$n;
    gseaplot(gse, by ="runningScore", color.line = "#009E73", title = paste0(gse$Description[n], " p = ", gse$p.adjust[n]), geneSetID = n) + theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7));
    ggsave(filename = paste0(monocle_prefix, "dotplot_day7_cluster_", i, "_", gse$Description[n], ".pdf"), width = 14, height = 8);
  }}
  
  write.csv(as.data.frame(fit_coefs) %>% dplyr::select(-c(model, model_summary, status)), paste0(monocle_prefix, "dge_day7_cluster_", i, ".csv"))
}

##modules analysis
#subset for mature cardio
p1 = plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster")
p1
ggsave(filename = paste0(monocle_prefix, "clusters_day4.pdf"),
       width = 12, height = 8)

clusters = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
colData(cds_day4)$cluster_day4 = clusters$cluster

cardio_day4 = cds_day4[, colData(cds_day4) %>% subset(cluster_day4 %in% c("2", "4", "5", "6", "11")) %>% row.names]
meta_day4 = as.data.frame(cardio_day4@colData)
meta_day4 = meta_day4 %>% dplyr::mutate(type = paste0(shRNA, "_", knockdown))

p1 = plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster")
p1
ggsave(filename = paste0(monocle_prefix, "clusters_day7.pdf"),
       width = 12, height = 8)

clusters = data.frame(p1$data$sample_name, libID = p1$data$sample_name, cluster = p1$data$cell_color, row.names = 1)
colData(cds_day7)$cluster_day7 = clusters$cluster

cardio_day7 = cds_day7[, colData(cds_day7) %>% subset(cluster_day7 %in% c("1", "4")) %>% row.names]

#find modules
gene_module_day4 <- find_gene_modules(cardio_day4, resolution = 0.005)

ggplot(gene_module_day4) + geom_point(aes(x = dim_1, y = dim_2, color = module))
ggsave(filename = paste0(monocle_prefix, "genes_modules_umap_day4.pdf"),
       width = 8, height = 12)

for(i in unique(gene_module_day4$module)){
  genes_i = gene_module_day4 %>% dplyr::filter(module == i)
  GO = ensembldb::select(org.Hs.eg.db, keys = genes_i$id, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  go = goana(GO$ENTREZID, species = "Hs", convert = T)
  res = topGO(go, ontology = "BP", number = 30, truncate.term = 50)
  ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
    labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of BP of genes in ", as.character(i), " module"))
  ggsave(filename = paste0(monocle_prefix, "GO_module_", i, "_day4.pdf"), width = 14, height = 14)
}

cell_group_day4 <- data.frame(cell = row.names(colData(cardio_day4)), 
                              cluster = cardio_day4@colData@listData[["cluster_day4"]])
agg_mat_day4 <- aggregate_gene_expression(cardio_day4, gene_module_day4, cell_group_day4)

as.ggplot(pheatmap::pheatmap(agg_mat_day4, cluster_rows = TRUE, cluster_cols = TRUE, scale = "column", color = viridis::viridis(n = 25, option = "H")))
ggsave(filename = paste0(monocle_prefix, "heatmap_modules_clusters_day4.pdf"),
       width = 8, height = 12)


cell_group_KD <- data.frame(cell = row.names(colData(cardio_day4)), 
                            cluster = paste0(cardio_day4@colData@listData[["shRNA"]], "_", cardio_day4@colData@listData[["knockdown"]]))
agg_mat_KD <- aggregate_gene_expression(cardio_day4, gene_module_day4, cell_group_KD)

as.ggplot(pheatmap::pheatmap(agg_mat_KD, cluster_rows = TRUE, cluster_cols = TRUE, scale = "column", color = viridis::viridis(n = 25, option = "H")))
ggsave(filename = paste0(monocle_prefix, "heatmap_modules_KD_day4.pdf"),
       width = 8, height = 12)

gene_module_day7 <- find_gene_modules(cardio_day7, resolution = 0.005)

ggplot(gene_module_day7) + geom_point(aes(x = dim_1, y = dim_2, color = module))
ggsave(filename = paste0(monocle_prefix, "genes_modules_umap_day7.pdf"),
       width = 8, height = 12)

for(i in unique(gene_module_day7$module)){
  genes_i = gene_module_day7 %>% dplyr::filter(module == i)
  GO = ensembldb::select(org.Hs.eg.db, keys = genes_i$id, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  go = goana(GO$ENTREZID, species = "Hs", convert = T)
  res = topGO(go, ontology = "BP", number = 30, truncate.term = 50)
  ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
    labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of BP of genes in ", as.character(i), " module"))
  ggsave(filename = paste0(monocle_prefix, "GO_module_", i, "_day7.pdf"), width = 14, height = 14)
}

cell_group_day7 <- data.frame(cell = row.names(colData(cardio_day7)), 
                              cluster = cardio_day7@colData@listData[["cluster_day7"]])
agg_mat_day7 <- aggregate_gene_expression(cardio_day7, gene_module_day7, cell_group_day7)

as.ggplot(pheatmap::pheatmap(agg_mat_day7, cluster_rows = TRUE, cluster_cols = TRUE, scale = "column", color = viridis::viridis(n = 25, option = "H")))
ggsave(filename = paste0(monocle_prefix, "heatmap_modules_day7.pdf"),
       width = 8, height = 12)

cell_group_KD <- data.frame(cell = row.names(colData(cardio_day7)), 
                            cluster = paste0(cardio_day7@colData@listData[["shRNA"]], "_", cardio_day7@colData@listData[["knockdown"]]))
agg_mat_KD <- aggregate_gene_expression(cardio_day7, gene_module_day7, cell_group_KD)

as.ggplot(pheatmap::pheatmap(agg_mat_KD, cluster_rows = TRUE, cluster_cols = TRUE, scale = "column", color = viridis::viridis(n = 25, option = "H")))
ggsave(filename = paste0(monocle_prefix, "heatmap_modules_KD_day7.pdf"),
       width = 8, height = 12)

#plot modules in umap
module_dendro <- hclust(dist(agg_mat_day4))
gene_module_day4$module <- factor(gene_module_day4$module, 
                                  levels = row.names(agg_mat_day4)[module_dendro$order])
plot_cells(cardio_day4,
           genes = gene_module_day4,
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE)
ggsave(filename = paste0(monocle_prefix, "umap_modules_day4.pdf"),
       width = 8, height = 12)

module_dendro <- hclust(dist(agg_mat_day7))
gene_module_day7$module <- factor(gene_module_day7$module, 
                                  levels = row.names(agg_mat_day7)[module_dendro$order])
plot_cells(cardio_day7,
           genes = gene_module_day7,
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE)
ggsave(filename = paste0(monocle_prefix, "umap_modules_day7.pdf"),
       width = 8, height = 12)

##clusters enrichment

#fisher test
meta_day4 = as.data.frame(cds_day4@colData) %>% dplyr::mutate(type = paste0(shRNA, "_", knockdown)) %>% dplyr::mutate(n = 1) %>% dplyr::mutate(type = ifelse(knockdown == "CTR", "zCTR", type))
meta_day7 = as.data.frame(cds_day7@colData) %>% dplyr::mutate(type = paste0(shRNA, "_", knockdown)) %>% dplyr::mutate(n = 1) %>% dplyr::mutate(type = ifelse(knockdown == "CTR", "zCTR", type))

conditions = c("SCR_TET", "GATA4_TET", "CTCF_TET")

#day4
f_test_day4 = data.frame(condition = c(), cluster = c(), p.val = c())

for(j in conditions){
  
  for(i in unique(meta_day4$cluster_day4)){
    table = meta_day4 %>% dplyr::filter(type %in% c(j, "zCTR")) %>% dplyr::mutate(cluster = ifelse(cluster_day4 == i, i, "other")) %>%
      dplyr::select(c("cluster", "type"))
    table = table(table)
    test <- fisher.test(table, alternative = "two.sided")
    df = data.frame(condition = j, cluster = i, p.val = test$p.value)
    f_test_day4 = rbind(f_test_day4, df)
    
  }
}

f_test_day4$q.val=p.adjust(f_test_day4$p.val)

for(j in conditions){
  
  for(i in unique(meta_day4$cluster_day4)){
    table = meta_day4 %>% dplyr::filter(type %in% c(j, "zCTR")) %>% dplyr::mutate(cluster = ifelse(cluster_day4 == i, i, "other")) %>%
      dplyr::select(c("cluster", "type"))
    table = table(table)
    q = f_test_day4 %>% dplyr::filter(condition == j & cluster == i)
    
    pdf(paste0(monocle_prefix, "mosaic_cluster_day4_", i, "_", j, ".pdf"))
    mosaicplot(table, main = paste0("Mosaic plot day 4 cluster ", i, " ", j, "\np-value = ", q$q.val),
               color = okabe_pal)
    dev.off()
  }
}

df_day4 = data.frame(type = unique(meta_day4$type))
c_names = colnames(df_day4)
for(i in unique(meta_day4$cluster_day4)){
  meta_i = meta_day4 %>% dplyr::filter(cluster_day4 == i) %>% dplyr::select(type, n)
  meta_i = meta_i %>% dplyr::group_by(type) %>% dplyr::summarise(sum = sum(n))
  df_day4 = df_day4 %>% dplyr::left_join(meta_i, by = "type")
  c_names = append(c_names, i)
  colnames(df_day4) = c_names
}
rownames(df_day4) = df_day4$type
df_day4[is.na(df_day4)] = 0 
df_day4 = df_day4 %>% dplyr::select(-c("type"))
pdf(paste0(monocle_prefix, "mosaic_cluster_day4_all.pdf"))
mosaicplot(df_day4, main = "Mosaic plot day 4",
           color = okabe_pal)
dev.off()

meta_day4 = as.data.frame(cds_day4@colData) %>% dplyr::mutate(type = paste0(shRNA, "_", knockdown)) %>% dplyr::mutate(n = 1)
conditions = unique(meta_day4$type)
meta_day4 = meta_day4 %>% dplyr::select(c(type, cluster_day4, n))
perc_day4 = data.frame(cluster_day4 = unique(meta_day4$cluster_day4))

for(i in conditions){
  perc_i = meta_day4 %>% dplyr::filter(type == i)
  n = length(rownames(perc_i))
  perc_i = perc_i %>% dplyr::group_by(cluster_day4) %>% dplyr::summarise(count = sum(n)) %>% dplyr::mutate(count = count/n * 100)
  colnames(perc_i) = c("cluster_day4", paste0(i, "_%of_cells"))
  perc_day4 = perc_day4 %>% dplyr::left_join(perc_i, by = "cluster_day4")
  
  if(strsplit(i, "_")[[1]][2] != "CTR"){
    f_test = f_test_day4 %>% dplyr::filter(condition == i) %>% dplyr::select(cluster, p.val, q.val)
    colnames(f_test) = c("cluster_day4", paste0("p.val_", i), paste0("q.val", i))
    
    perc_day4 = perc_day4 %>% dplyr::left_join(f_test, by = "cluster_day4")
  }
}
write.csv(perc_day4, paste0("/home/shared_folder/tables/data_day4.csv"))

#day7
f_test_day7 = data.frame(condition = c(), cluster = c(), p.val = c())

for(j in conditions){
  
  for(i in unique(meta_day7$cluster_day7)){
    table = meta_day7 %>% dplyr::filter(type %in% c(j, "zCTR")) %>% dplyr::mutate(cluster = ifelse(cluster_day7 == i, i, "other")) %>%
      dplyr::select(c("cluster", "type"))
    table = table(table)
    test <- fisher.test(table, alternative = "two.sided")
    df = data.frame(condition = j, cluster = i, p.val = test$p.value)
    f_test_day7 = rbind(f_test_day7, df)
    
  }
}

f_test_day7$q.val=p.adjust(f_test_day7$p.val)

for(j in conditions){
  
  for(i in unique(meta_day7$cluster_day7)){
    table = meta_day7 %>% dplyr::filter(type %in% c(j, "zCTR")) %>% dplyr::mutate(cluster = ifelse(cluster_day7 == i, i, "other")) %>%
      dplyr::select(c("cluster", "type"))
    table = table(table)
    q = f_test_day7 %>% dplyr::filter(condition == j & cluster == i)
    
    pdf(paste0(monocle_prefix, "mosaic_cluster_day7_", i, "_", j, ".pdf"))
    mosaicplot(table, main = paste0("Mosaic plot day 7 cluster ", i, " ", j, "\np-value = ", q$q.val),
               color = okabe_pal)
    dev.off()
  }
}

df_day7 = data.frame(type = unique(meta_day7$type))
c_names = colnames(df_day7)
for(i in unique(meta_day7$cluster_day7)){
  meta_i = meta_day7 %>% dplyr::filter(cluster_day7 == i)
  meta_i = meta_i %>% dplyr::group_by(type) %>% dplyr::summarise(sum = sum(n))
  df_day7 = df_day7 %>% dplyr::left_join(meta_i, by = "type")
  c_names = append(c_names, i)
  colnames(df_day7) = c_names
}
rownames(df_day7) = df_day7$type
df_day7[is.na(df_day7)] = 0 
df_day7 = df_day7 %>% dplyr::select(-c("type"))
pdf(paste0(monocle_prefix, "mosaic_cluster_day7_all.pdf"))
mosaicplot(df_day7, main = "Mosaic plot day 7",
           color = okabe_pal)
dev.off()

meta_day7 = as.data.frame(cds_day7@colData) %>% dplyr::mutate(type = paste0(shRNA, "_", knockdown)) %>% dplyr::mutate(n = 1)
conditions = unique(meta_day7$type)
meta_day7 = meta_day7 %>% dplyr::select(c(type, cluster_day7, n))
perc_day7 = data.frame(cluster_day7 = unique(meta_day7$cluster_day7))

for(i in conditions){
  perc_i = meta_day7 %>% dplyr::filter(type == i)
  n = length(rownames(perc_i))
  perc_i = perc_i %>% dplyr::group_by(cluster_day7) %>% dplyr::summarise(count = sum(n)) %>% dplyr::mutate(count = count/n * 100)
  colnames(perc_i) = c("cluster_day7", paste0(i, "_%of_cells"))
  perc_day7 = perc_day7 %>% dplyr::left_join(perc_i, by = "cluster_day7")
  
  if(strsplit(i, "_")[[1]][2] != "CTR"){
    f_test = f_test_day7 %>% dplyr::filter(condition == i) %>% dplyr::select(cluster, p.val, q.val)
    colnames(f_test) = c("cluster_day7", paste0("p.val_", i), paste0("q.val", i))
    
    perc_day7 = perc_day7 %>% dplyr::left_join(f_test, by = "cluster_day7")
  }
}


write.csv(perc_day7, paste0("/home/shared_folder/tables/data_day7.csv"))

##plot markers
#pseudotime plots
plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "clusters_day4.pdf"),
       width = 10, height = 12)

plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "clusters_day7.pdf"),
       width = 10, height = 12)

plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "pseudotime")
ggsave(filename = paste0(monocle_prefix, "pseudotime_day4.pdf"),
       width = 10, height = 12)

plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "pseudotime")
ggsave(filename = paste0(monocle_prefix, "pseudotime_day7.pdf"),
       width = 10, height = 12)

##day4
#proliferation
gene_list = c("MKI67", "MCM6", "PCNA", "TOP2A")
plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day4_prol.pdf"),
       width = 10, height = 12)

#fibroblasts
gene_list = c("COL3A1", "COL1A1", "LUM", "FN1")
plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day4_fibro.pdf"),
       width = 10, height = 12)

#pluripotency
gene_list = c("POU5F1", "SOX2", "NANOG")
plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day4_plur.pdf"),
       width = 10, height = 12)

#mesoendoderm
gene_list = c("EOMES", "FOXA2", "SOX17", "APOB")
plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day4_mesoend.pdf"),
       width = 10, height = 12)

#angiogenesis
gene_list = c("KDR", "ERG", "TEK", "ESAM")
plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day4_angio.pdf"),
       width = 10, height = 12)

#cardio markers
gene_list = c("TTN", "RYR2", "SLC8A1")
plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day4_cardio.pdf"),
       width = 10, height = 12)

subset_day4 = cds_day4[, colData(cds_day4) %>% subset(cluster_day4 %in% c("2", "3", "4", "5", "6", "8")) %>% row.names]
plot_cells(subset_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "pseudotime")
ggsave(filename = paste0(monocle_prefix, "pseudotime_day4_cardio.pdf"),
       width = 10, height = 12)
cds_subset = subset_day4[rowData(subset_day4)$gene_short_name %in% gene_list]
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by = "cluster_day4",
                         min_expr=0.5)
ggsave(filename = paste0(monocle_prefix, "genes_pseudotime_day4_cardio.pdf"),
       width = 10, height = 18)

#plot by KD
subset_day4 = cds_day4[, colData(cds_day4) %>% subset(knockdown == "CTR") %>% row.names]
plot_cells(subset_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "subset_day4_CTR.pdf"),
       width = 10, height = 12)

subset_day4 = cds_day4[, colData(cds_day4) %>% subset(knockdown == "TET" & shRNA == "SCR") %>% row.names]
plot_cells(subset_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "subset_day4_SCR.pdf"),
       width = 10, height = 12)

subset_day4 = cds_day4[, colData(cds_day4) %>% subset(knockdown == "TET" & shRNA == "GATA4") %>% row.names]
plot_cells(subset_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "subset_day4_GATA4.pdf"),
       width = 10, height = 12)

subset_day4 = cds_day4[, colData(cds_day4) %>% subset(knockdown == "TET" & shRNA == "CTCF") %>% row.names]
plot_cells(subset_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "subset_day4_CTCF.pdf"),
       width = 10, height = 12)



#day7
#proliferation
gene_list = c("MKI67", "MCM6", "PCNA", "TOP2A")
plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day7_prol.pdf"),
       width = 10, height = 12)

#fibroblasts
gene_list = c("COL3A1", "COL1A1", "LUM", "FN1")
plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day7_fibro.pdf"),
       width = 10, height = 12)

#pluripotency
gene_list = c("POU5F1", "SOX2", "NANOG")
plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day7_plur.pdf"),
       width = 10, height = 12)

#mesoendoderm
gene_list = c("EOMES", "FOXA2", "SOX17", "APOB")
plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day7_mesoend.pdf"),
       width = 10, height = 12)

#angiogenesis
gene_list = c("KDR", "ERG", "TEK", "ESAM")
plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day7_angio.pdf"),
       width = 10, height = 12)

#cardio markers
gene_list = c("TTN", "RYR2", "SLC8A1")
plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", genes = gene_list, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
ggsave(filename = paste0(monocle_prefix, "markers_day7_cardio.pdf"),
       width = 10, height = 12)

subset_day7 = cds_day7[, colData(cds_day7) %>% subset(cluster_day7 %in% c("1", "4")) %>% row.names]
plot_cells(subset_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "pseudotime")
ggsave(filename = paste0(monocle_prefix, "pseudotime_day7_cardio.pdf"),
       width = 10, height = 12)
cds_subset = subset_day7[rowData(subset_day7)$gene_short_name %in% gene_list]
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by = "cluster_day7",
                         min_expr=0.5)
ggsave(filename = paste0(monocle_prefix, "genes_pseudotime_day7_cardio.pdf"),
       width = 10, height = 18)

#plot by KD
subset_day7 = cds_day7[, colData(cds_day7) %>% subset(knockdown == "CTR") %>% row.names]
plot_cells(subset_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "subset_day7_CTR.pdf"),
       width = 10, height = 12)

subset_day7 = cds_day7[, colData(cds_day7) %>% subset(knockdown == "TET" & shRNA == "SCR") %>% row.names]
plot_cells(subset_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "subset_day7_SCR.pdf"),
       width = 10, height = 12)

subset_day7 = cds_day7[, colData(cds_day7) %>% subset(knockdown == "TET" & shRNA == "GATA4") %>% row.names]
plot_cells(subset_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "subset_day7_GATA4.pdf"),
       width = 10, height = 12)

subset_day7 = cds_day7[, colData(cds_day7) %>% subset(knockdown == "TET" & shRNA == "CTCF") %>% row.names]
plot_cells(subset_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cluster", show_trajectory_graph=FALSE)
ggsave(filename = paste0(monocle_prefix, "subset_day7_CTCF.pdf"),
       width = 10, height = 12)


meta_day4 = as.data.frame(colData(cds_day4))
meta_day4 = meta_day4 %>% dplyr::mutate(cluster_name = ifelse(cluster_day4 == "1", "01_fibroblasts",
                                                              ifelse(cluster_day4 == "2", "02_mature_CMs1",
                                                                     ifelse(cluster_day4 == "3", "03_proliferating_cells",
                                                                            ifelse(cluster_day4 == "4", "04_immature_CMs",
                                                                                   ifelse(cluster_day4 == "5", "05_CTCF_only_CMs",
                                                                                          ifelse(cluster_day4 == "6", "06_mature_CMs2",
                                                                                                 ifelse(cluster_day4 == "7", "07_mesendoderm",
                                                                                                        ifelse(cluster_day4 == "8", "08_early_CM_progenitors",
                                                                                                               ifelse(cluster_day4 == "9", "09_undifferentiated_cells",
                                                                                                                      ifelse(cluster_day4 == "10", "10_endoderm",
                                                                                                                             ifelse(cluster_day4 == "11", "11_mature_CMs3",
                                                                                                                                    ifelse(cluster_day4 == "12", "12_endothelial_cells", "")))))))))))))
meta_day7 = as.data.frame(colData(cds_day7))
meta_day7 = meta_day7 %>% dplyr::mutate(cluster_name = ifelse(cluster_day7 == "1", "01_mature_CMs",
                                                              ifelse(cluster_day7 == "2", "02_undifferentiated_cells",
                                                                     ifelse(cluster_day7 == "3", "03_mesoderm",
                                                                            ifelse(cluster_day7 == "4", "04_immature_CMs",
                                                                                   ifelse(cluster_day7 == "5", "05_fibroblasts",
                                                                                          ifelse(cluster_day7 == "6", "06_endoderm",
                                                                                                 ifelse(cluster_day7 == "7", "07_proliferating_cells",
                                                                                                        ifelse(cluster_day7 == "8", "08_endothelial_cells",
                                                                                                               ifelse(cluster_day7 == "9", "09_mesendoderm",
                                                                                                                      ifelse(cluster_day7 == "10", "10_endothelial_cells2", "")))))))))))



colData(cds_day4)$cell_type = meta_day4$cluster_name
colData(cds_day7)$cell_type = meta_day7$cluster_name


#plots trajectory all
p1 = plot_cells(cds_day4, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cell_type", show_trajectory_graph=TRUE, label_cell_groups=FALSE)
p1
ggsave(paste0(monocle_prefix, "UMAP_traj_day4.pdf"), width = 12, height = 10)
meta_day4$UMAP1 =  p1$data$data_dim_1
meta_day4$UMAP2 =  p1$data$data_dim_2


p1 = plot_cells(cds_day7, norm_method = "log", reduction_method = "UMAP", color_cells_by = "cell_type", show_trajectory_graph=TRUE, label_cell_groups=FALSE)
p1
ggsave(paste0(monocle_prefix, "UMAP_traj_day7.pdf"), width = 12, height = 10)
meta_day7$UMAP1 =  p1$data$data_dim_1
meta_day7$UMAP2 =  p1$data$data_dim_2


#cardio_pseudotime
trajectory_day4 = meta_day4 %>% dplyr::filter(cluster_day4 %in% c("2", "3", "4", "5", "6", "8"))
p1 = ggplot(meta_day4) + geom_point(aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  scale_color_manual(values = pal)
p2 = ggplot() + geom_point(aes(x = meta_day4$UMAP1, y = meta_day4$UMAP2), color = "grey") +
  geom_point(aes(x = trajectory_day4$UMAP1, y = trajectory_day4$UMAP2, color = trajectory_day4$pseudotime)) + scale_color_viridis()
p1 + p2
ggsave(paste0(monocle_prefix, "clusters_pseudotime_day4.pdf"), width = 20, height = 12)

trajectory_day7 = meta_day7 %>% dplyr::filter(cluster_day7 %in% c("1", "4"))
p1 = ggplot(meta_day7) + geom_point(aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  scale_color_manual(values = pal)
p2 = ggplot() + geom_point(aes(x = meta_day7$UMAP1, y = meta_day7$UMAP2), color = "grey") +
  geom_point(aes(x = trajectory_day7$UMAP1, y = trajectory_day7$UMAP2, color = trajectory_day7$pseudotime)) + scale_color_viridis()
p1 + p2
ggsave(paste0(monocle_prefix, "clusters_pseudotime_day7.pdf"), width = 20, height = 12)

#downsampling
stat_day4 = meta_day4 %>% dplyr::group_by(sample) %>% dplyr::summarise(sum = sum(n))
down_day4 = data.frame(libID = c(), UMAP1 = c(), UMAP2 = c(), sample = c(), cell_type = c())
for(i in unique(meta_day4$sample)){
  meta_i = meta_day4 %>% dplyr::filter(sample == i)
  cells <- sample(x = meta_i$libID, size = 900, replace = FALSE)
  meta_i = meta_i %>% dplyr::filter(libID %in% cells) %>% dplyr::select(libID, UMAP1, UMAP2, sample, cell_type)
  down_day4 = rbind(down_day4, meta_i)
}

ggplot(down_day4) + geom_point(aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  facet_wrap(~forcats::fct_relevel(sample)) +
  scale_color_manual(values = pal)
ggsave(paste0(monocle_prefix, "UMAP_downsamples_day4.pdf"), width = 20, height = 18)

stat_day7 = meta_day7 %>% dplyr::group_by(sample) %>% dplyr::summarise(sum = sum(n))
down_day7 = data.frame(libID = c(), UMAP1 = c(), UMAP2 = c(), sample = c(), cell_type = c())
for(i in unique(meta_day7$sample)){
  meta_i = meta_day7 %>% dplyr::filter(sample == i)
  cells <- sample(x = meta_i$libID, size = 900, replace = FALSE)
  meta_i = meta_i %>% dplyr::filter(libID %in% cells) %>% dplyr::select(libID, UMAP1, UMAP2, sample, cell_type)
  down_day7 = rbind(down_day7, meta_i)
}

ggplot(down_day7) + geom_point(aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  facet_wrap(~forcats::fct_relevel(sample)) +
  scale_color_manual(values = pal)
ggsave(paste0(monocle_prefix, "UMAP_downsamples_day7.pdf"), width = 20, height = 18)

#save
save(cds_day4, cds_day7, meta_day4, meta_day7,
     file = paste0("/home/shared_folder/RData/cds_final.RData"))

##cell cycle analysis
meta_day4$n = 1
ggplot(meta_day4, aes(fill = Phase, y = n, x = sample_id)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("grey","#56B4E9","#ff6db6"))
ggsave(paste0(monocle_prefix, "percentage_Phase_day4.pdf"))
ggplot(meta_day4 %>% dplyr::filter(cluster_name %in% c("02_mature_CMs1","04_immature_CMs", "05_CTCF_only_CMs", "06_mature_CMs2", "11_mature_CMs3")), aes(fill = Phase, y = n, x = sample_id)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("grey","#56B4E9","#ff6db6"))
ggsave(paste0(monocle_prefix, "percentage_Phase_day4_cardio.pdf"))

meta_day7$n = 1
ggplot(meta_day7, aes(fill = Phase, y = n, x = sample_id)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("grey","#56B4E9","#ff6db6"))
ggsave(paste0(monocle_prefix, "percentage_Phase_day7.pdf"))
ggplot(meta_day7 %>% dplyr::filter(cluster_name %in% c("01_mature_CMs","04_immature_CMs")), aes(fill = Phase, y = n, x = sample_id)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("grey","#56B4E9","#ff6db6"))
ggsave(paste0(monocle_prefix, "percentage_Phase_day7_cardio.pdf"))

combs = list(c("SCR_CTR", "SCR_TET"), c("GATA4_CTR", "GATA4_TET"), c("CTCF_CTR", "CTCF_TET"))
ggplot(meta_day4, aes(x = type,  y = S.Score)) + geom_violin(aes(fill = type)) + geom_signif(comparisons = combs, test = "wilcox.test", step_increase = 0.05, size = 0.3, textsize = 3) +
  geom_boxplot(aes(alpha = 0.2)) + scale_fill_manual(values = samples)
ggsave(filename = paste0(monocle_prefix, "cc_all_S_day4_stats.pdf"))

ggplot(meta_day4, aes(x = type,  y = S.Score)) + geom_violin(aes(fill = type)) + geom_signif(comparisons = combs, test = "wilcox.test", step_increase = 0.05, size = 0.3, textsize = 3) +
  geom_boxplot(aes(alpha = 0.2)) + scale_fill_manual(values = samples)
ggsave(filename = paste0(monocle_prefix, "cc_all_G2M_day4_stats.pdf"))

ggplot(meta_day7, aes(x = type,  y = S.Score)) + geom_violin(aes(fill = type)) + geom_signif(comparisons = combs, test = "wilcox.test", step_increase = 0.05, size = 0.3, textsize = 3) +
  geom_boxplot(aes(alpha = 0.2)) + scale_fill_manual(values = samples)
ggsave(filename = paste0(monocle_prefix, "cc_all_S_day7_stats.pdf"))

ggplot(meta_day7, aes(x = type,  y = S.Score)) + geom_violin(aes(fill = type)) + geom_signif(comparisons = combs, test = "wilcox.test", step_increase = 0.05, size = 0.3, textsize = 3) +
  geom_boxplot(aes(alpha = 0.2)) + scale_fill_manual(values = samples)
ggsave(filename = paste0(monocle_prefix, "cc_all_G2M_day7_stats.pdf"))

write.csv(meta_day4 %>% dplyr::select(c(libID, sample_id, type, Phase, S.Score, G2M.Score, cluster_name)), paste0("/home/shared_folder/tables/cell_cycle_day4.csv"))
write.csv(meta_day7 %>% dplyr::select(c(libID, sample_id, type, Phase, S.Score, G2M.Score, cluster_name)), paste0("/home/shared_folder/tables/cell_cycle_day7.csv"))


##milo cluster enrichment analysis
counts <- assay(cds_day4, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(cds_day4) <- log2(t(t(counts)/size.factors) + 1)
seurat_day4 <- as.Seurat(cds_day4, data = NULL)

seurat_day4 <- FindVariableFeatures(seurat_day4, selection.method = "vst", nfeatures = 8000)
seurat_day4 = ScaleData(seurat_day4)
seurat_day4 = RunPCA(seurat_day4 , npcs = 30)
seurat_day4 = AddMetaData(seurat_day4, metadata = meta_day4)

seurat_SCR = seurat_day4[,seurat_day4$condition == "CTR"|seurat_day4$condition == "SCR_TET"]
seurat_GATA4 = seurat_day4[,seurat_day4$condition == "CTR"|seurat_day4$condition == "GATA4_TET"]
seurat_CTCF = seurat_day4[,seurat_day4$condition == "CTR"|seurat_day4$condition == "CTCF_TET"]

sce <- as.SingleCellExperiment(seurat_SCR)
milo <- Milo(sce)
set.seed(1)
milo <- buildGraph(milo, k = 30, d = 30)
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE)
plotNhoodSizeHist(milo)
ggsave(paste0(monocle_prefix, "SCR_milo_distribution_day4.pdf"))

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample = "sample_id")
design_df <- distinct(data.frame(milo@colData@listData)[,c("sample_id","condition")])
rownames(design_df) <- design_df$sample_id
design_df$condition <- factor(design_df$condition)

milo <- calcNhoodDistance(milo, d = 30)
milo_SCR = milo

da_results <- testNhoods(milo_SCR, design = ~ condition, design.df = design_df)
da_results <- annotateNhoods(milo_SCR, da_results, coldata_col = "cell_type")
da_results_SCR = da_results

ggplot(da_results_SCR, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
  labs(title = "Differential abundance per KD",
       x = "log2 Fold Change (mediums)",
       y = "-log10 SpatialFDR")
ggsave(paste0(monocle_prefix, "SCR_milo_volcano_day4.pdf"))

sce <- as.SingleCellExperiment(seurat_GATA4)
milo <- Milo(sce)
set.seed(1)
milo <- buildGraph(milo, k = 30, d = 30)
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE)
plotNhoodSizeHist(milo)
ggsave(paste0(monocle_prefix, "GATA4_milo_distribution_day4.pdf"))

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample = "sample_id")
design_df <- distinct(data.frame(milo@colData@listData)[,c("sample_id","condition")])
rownames(design_df) <- design_df$sample_id
design_df$condition <- factor(design_df$condition)

milo <- calcNhoodDistance(milo, d = 30)
milo_GATA4 = milo

da_results <- testNhoods(milo_GATA4, design = ~ condition, design.df = design_df)
da_results <- annotateNhoods(milo_GATA4, da_results, coldata_col = "cell_type")
da_results_GATA4 = da_results

ggplot(da_results_GATA4, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
  labs(title = "Differential abundance per KD",
       x = "log2 Fold Change (mediums)",
       y = "-log10 SpatialFDR")
ggsave(paste0(monocle_prefix, "GATA4_milo_volcano_day4.pdf"))

sce <- as.SingleCellExperiment(seurat_CTCF)
milo <- Milo(sce)
set.seed(1)
milo <- buildGraph(milo, k = 30, d = 30)
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE)
plotNhoodSizeHist(milo)
ggsave(paste0(monocle_prefix, "CTCF_milo_distribution_day4.pdf"))

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample = "sample_id")
design_df <- distinct(data.frame(milo@colData@listData)[,c("sample_id","condition")])
rownames(design_df) <- design_df$sample_id
design_df$condition <- factor(design_df$condition)

milo <- calcNhoodDistance(milo, d = 30)
milo_CTCF = milo

da_results <- testNhoods(milo_CTCF, design = ~ condition, design.df = design_df)
da_results <- annotateNhoods(milo_CTCF, da_results, coldata_col = "cell_type")
da_results_CTCF = da_results

ggplot(da_results_CTCF, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
  labs(title = "Differential abundance per KD",
       x = "log2 Fold Change (mediums)",
       y = "-log10 SpatialFDR")
ggsave(paste0(monocle_prefix, "CTCF_milo_volcano_day4.pdf"))

# beeswarm per cluster
beeswarm_SCR <-plotDAbeeswarm(da_results_SCR, group.by = "cell_type", alpha = 0.1)
beeswarm_df_SCR_day4 <- beeswarm_SCR$data
beeswarm_df_SCR_day4 <- beeswarm_df_SCR_day4 %>% dplyr::mutate(day = "day4", shRNA = "SCR")
beeswarm_GATA4 <-plotDAbeeswarm(da_results_GATA4, group.by = "cell_type", alpha = 0.1)
beeswarm_df_GATA4_day4 <- beeswarm_GATA4$data 
beeswarm_df_GATA4_day4 <- beeswarm_df_GATA4_day4 %>% dplyr::mutate(day = "day4", shRNA = "GATA4")
beeswarm_CTCF <-plotDAbeeswarm(da_results_CTCF, group.by = "cell_type", alpha = 0.1)
beeswarm_df_CTCF_day4 <- beeswarm_CTCF$data
beeswarm_df_CTCF_day4 <- beeswarm_df_CTCF_day4 %>% dplyr::mutate(day = "day4", shRNA = "CTCF")

df <- rbind(beeswarm_df_SCR_day4, beeswarm_df_GATA4_day4, beeswarm_df_CTCF_day4) %>%
  filter(!is.na(.data[["cell_type"]])) %>%
  mutate(
    cluster = .data[["cell_type"]],
    is_sig  = !is.na(SpatialFDR) & SpatialFDR < 0.05
  )

cluster_summary <- df %>%
  group_by(cell_type, day, shRNA) %>%
  summarise(
    n_nhoods         = n(),
    n_sig            = sum(is_sig, na.rm = TRUE),
    frac_sig         = ifelse(n_nhoods > 0, n_sig / n_nhoods, NA_real_),
    min_SpatialFDR   = suppressWarnings(min(SpatialFDR, na.rm = TRUE)),
    mean_logFC       = mean(logFC, na.rm = TRUE),
    median_logFC     = median(logFC, na.rm = TRUE),
    mean_logFC_sig   = if (any(is_sig)) mean(logFC[is_sig], na.rm = TRUE) else NA_real_,
    median_logFC_sig = if (any(is_sig)) median(logFC[is_sig], na.rm = TRUE) else NA_real_,
    prop_pos_sig     = if (any(is_sig)) mean(logFC[is_sig] > 0, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(desc(frac_sig), desc(abs(median_logFC_sig)))
cluster_summary$sig=ifelse(cluster_summary$median_logFC_sig > 0.5 & cluster_summary$frac_sig > 0.5 |
                                      cluster_summary$median_logFC_sig < 0.5*-1 &cluster_summary$frac_sig > 0.5,"*","")
cluster_summary_day4 = cluster_summary
write.csv(cluster_summary_day4, file = paste0("/home/shared_folder/tables/cluster_summary_Milo_day4.csv"))

ggplot(beeswarm_df_SCR_day4 %>% dplyr::left_join(cluster_summary_day4, by = c("cell_type", "shRNA")), aes(x = cell_type, y = -logFC*-1, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.7) +
  theme_classic(base_size = 14) +geom_hline(yintercept=2,linetype="dashed")+geom_hline(yintercept=-2,linetype="dashed") +
  labs(title = "Distribution of logFC per cluster SCR TET vs CTR",
       x = "", y = "log2 Fold Change comparing SCR TET vs CTR") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-13, 13)) +
  scale_fill_manual(values = pal) + geom_text(aes(x = cell_type, y = 12, label = sig))
ggsave(paste0(monocle_prefix, "SCR_milo_violin_day4.pdf"))

ggplot(beeswarm_df_GATA4_day4 %>% dplyr::left_join(cluster_summary_day4, by = c("cell_type", "shRNA")), aes(x = cell_type, y = -logFC*-1, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.7) +
  theme_classic(base_size = 14) +geom_hline(yintercept=2,linetype="dashed")+geom_hline(yintercept=-2,linetype="dashed") +
  labs(title = "Distribution of logFC per cluster GATA4 TET vs CTR",
       x = "", y = "log2 Fold Change comparing GATA4 TET vs CTR") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-13, 13)) +
  scale_fill_manual(values = pal) + geom_text(aes(x = cell_type, y = 12, label = sig))
ggsave(paste0(monocle_prefix, "GATA4_milo_violin_day4.pdf"))

ggplot(beeswarm_df_CTCF_day4 %>% dplyr::left_join(cluster_summary_day4, by = c("cell_type", "shRNA")), aes(x = cell_type, y = logFC*-1, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.7) +
  theme_classic(base_size = 14) +geom_hline(yintercept=2,linetype="dashed")+geom_hline(yintercept=-2,linetype="dashed") +
  labs(title = "Distribution of logFC per cluster CTCF TET vs CTR",
       x = "", y = "log2 Fold Change comparing CTCF TET vs CTR") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-13, 13)) +
  scale_fill_manual(values = pal) + geom_text(aes(x = cell_type, y = 12, label = sig))
ggsave(paste0(monocle_prefix, "CTCF_milo_violin_day4.pdf"))

meta_day7 = meta_day7 %>% dplyr::mutate(condition = ifelse(knockdown == "TET", paste0(shRNA, "_", knockdown), knockdown))
counts <- assay(cds_day7, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(cds_day7) <- log2(t(t(counts)/size.factors) + 1)
seurat_day7 <- as.Seurat(cds_day7, data = NULL)

seurat_day7 <- FindVariableFeatures(seurat_day7, selection.method = "vst", nfeatures = 8000)
seurat_day7 = ScaleData(seurat_day7)
seurat_day7 = RunPCA(seurat_day7 , npcs = 30)
seurat_day7 = AddMetaData(seurat_day7, metadata = meta_day7)

seurat_SCR = seurat_day7[,seurat_day7$condition == "CTR"|seurat_day7$condition == "SCR_TET"]
seurat_GATA4 = seurat_day7[,seurat_day7$condition == "CTR"|seurat_day7$condition == "GATA4_TET"]
seurat_CTCF = seurat_day7[,seurat_day7$condition == "CTR"|seurat_day7$condition == "CTCF_TET"]

sce <- as.SingleCellExperiment(seurat_SCR)
milo <- Milo(sce)
set.seed(1)
milo <- buildGraph(milo, k = 30, d = 30)
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE)
plotNhoodSizeHist(milo)
ggsave(paste0(monocle_prefix, "SCR_milo_distribution_day7.pdf"))

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample = "sample_id")
design_df <- distinct(data.frame(milo@colData@listData)[,c("sample_id","condition")])
rownames(design_df) <- design_df$sample_id
design_df$condition <- factor(design_df$condition)

milo <- calcNhoodDistance(milo, d = 30)
milo_SCR = milo

da_results <- testNhoods(milo_SCR, design = ~ condition, design.df = design_df)
da_results <- annotateNhoods(milo_SCR, da_results, coldata_col = "cell_type")
da_results_SCR = da_results

ggplot(da_results_SCR, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
  labs(title = "Differential abundance per KD",
       x = "log2 Fold Change (mediums)",
       y = "-log10 SpatialFDR")
ggsave(paste0(monocle_prefix, "SCR_milo_volcano_day7.pdf"))

sce <- as.SingleCellExperiment(seurat_GATA4)
milo <- Milo(sce)
set.seed(1)
milo <- buildGraph(milo, k = 30, d = 30)
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE)
plotNhoodSizeHist(milo)
ggsave(paste0(monocle_prefix, "GATA4_milo_distribution_day7.pdf"))

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample = "sample_id")
design_df <- distinct(data.frame(milo@colData@listData)[,c("sample_id","condition")])
rownames(design_df) <- design_df$sample_id
design_df$condition <- factor(design_df$condition)

milo <- calcNhoodDistance(milo, d = 30)
milo_GATA4 = milo

da_results <- testNhoods(milo_GATA4, design = ~ condition, design.df = design_df)
da_results <- annotateNhoods(milo_GATA4, da_results, coldata_col = "cell_type")
da_results_GATA4 = da_results

ggplot(da_results_GATA4, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
  labs(title = "Differential abundance per KD",
       x = "log2 Fold Change (mediums)",
       y = "-log10 SpatialFDR")
ggsave(paste0(monocle_prefix, "GATA4_milo_volcano_day7.pdf"))

sce <- as.SingleCellExperiment(seurat_CTCF)
milo <- Milo(sce)
set.seed(1)
milo <- buildGraph(milo, k = 30, d = 30)
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE)
plotNhoodSizeHist(milo)
ggsave(paste0(monocle_prefix, "CTCF_milo_distribution_day7.pdf"))

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample = "sample_id")
design_df <- distinct(data.frame(milo@colData@listData)[,c("sample_id","condition")])
rownames(design_df) <- design_df$sample_id
design_df$condition <- factor(design_df$condition)

milo <- calcNhoodDistance(milo, d = 30)
milo_CTCF = milo

da_results <- testNhoods(milo_CTCF, design = ~ condition, design.df = design_df)
da_results <- annotateNhoods(milo_CTCF, da_results, coldata_col = "cell_type")
da_results_CTCF = da_results

ggplot(da_results_CTCF, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
  labs(title = "Differential abundance per KD",
       x = "log2 Fold Change (mediums)",
       y = "-log10 SpatialFDR")
ggsave(paste0(monocle_prefix, "CTCF_milo_volcano_day7.pdf"))

# beeswarm per cluster
beeswarm_SCR <-plotDAbeeswarm(da_results_SCR, group.by = "cell_type", alpha = 0.1)
beeswarm_df_SCR_day7 <- beeswarm_SCR$data
beeswarm_df_SCR_day7 <- beeswarm_df_SCR_day7 %>% dplyr::mutate(day = "day7", shRNA = "SCR")
beeswarm_GATA4 <-plotDAbeeswarm(da_results_GATA4, group.by = "cell_type", alpha = 0.1)
beeswarm_df_GATA4_day7 <- beeswarm_GATA4$data 
beeswarm_df_GATA4_day7 <- beeswarm_df_GATA4_day7 %>% dplyr::mutate(day = "day7", shRNA = "GATA4")
beeswarm_CTCF <-plotDAbeeswarm(da_results_CTCF, group.by = "cell_type", alpha = 0.1)
beeswarm_df_CTCF_day7 <- beeswarm_CTCF$data
beeswarm_df_CTCF_day7 <- beeswarm_df_CTCF_day7 %>% dplyr::mutate(day = "day7", shRNA = "CTCF")

df <- rbind(beeswarm_df_SCR_day7, beeswarm_df_GATA4_day7, beeswarm_df_CTCF_day7) %>%
  filter(!is.na(.data[["cell_type"]])) %>%
  mutate(
    cluster = .data[["cell_type"]],
    is_sig  = !is.na(SpatialFDR) & SpatialFDR < 0.05
  )

cluster_summary <- df %>%
  group_by(cell_type, day, shRNA) %>%
  summarise(
    n_nhoods         = n(),
    n_sig            = sum(is_sig, na.rm = TRUE),
    frac_sig         = ifelse(n_nhoods > 0, n_sig / n_nhoods, NA_real_),
    min_SpatialFDR   = suppressWarnings(min(SpatialFDR, na.rm = TRUE)),
    mean_logFC       = mean(logFC, na.rm = TRUE),
    median_logFC     = median(logFC, na.rm = TRUE),
    mean_logFC_sig   = if (any(is_sig)) mean(logFC[is_sig], na.rm = TRUE) else NA_real_,
    median_logFC_sig = if (any(is_sig)) median(logFC[is_sig], na.rm = TRUE) else NA_real_,
    prop_pos_sig     = if (any(is_sig)) mean(logFC[is_sig] > 0, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(desc(frac_sig), desc(abs(median_logFC_sig)))
cluster_summary$sig=ifelse(cluster_summary$median_logFC_sig > 0.5 & cluster_summary$frac_sig > 0.5 |
                             cluster_summary$median_logFC_sig < 0.5*-1 &cluster_summary$frac_sig > 0.5,"*","")
cluster_summary_day7 = cluster_summary
write.csv(cluster_summary_day7, file = paste0("/home/shared_folder/tables/cluster_summary_Milo_day7.csv"))

ggplot(beeswarm_df_SCR_day7 %>% dplyr::left_join(cluster_summary_day7, by = c("cell_type", "shRNA")), aes(x = cell_type, y = -logFC*-1, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.7) +
  theme_classic(base_size = 14) +geom_hline(yintercept=2,linetype="dashed")+geom_hline(yintercept=-2,linetype="dashed") +
  labs(title = "Distribution of logFC per cluster SCR TET vs CTR",
       x = "", y = "log2 Fold Change comparing SCR TET vs CTR") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-13, 13)) +
  scale_fill_manual(values = pal) + geom_text(aes(x = cell_type, y = 12, label = sig))
ggsave(paste0(monocle_prefix, "SCR_milo_violin_day7.pdf"))

ggplot(beeswarm_df_GATA4_day7 %>% dplyr::left_join(cluster_summary_day7, by = c("cell_type", "shRNA")), aes(x = cell_type, y = -logFC*-1, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.7) +
  theme_classic(base_size = 14) +geom_hline(yintercept=2,linetype="dashed")+geom_hline(yintercept=-2,linetype="dashed") +
  labs(title = "Distribution of logFC per cluster GATA4 TET vs CTR",
       x = "", y = "log2 Fold Change comparing GATA4 TET vs CTR") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-13, 13)) +
  scale_fill_manual(values = pal) + geom_text(aes(x = cell_type, y = 12, label = sig))
ggsave(paste0(monocle_prefix, "GATA4_milo_violin_day7.pdf"))

ggplot(beeswarm_df_CTCF_day7 %>% dplyr::left_join(cluster_summary_day7, by = c("cell_type", "shRNA")), aes(x = cell_type, y = logFC*-1, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.7) +
  theme_classic(base_size = 14) +geom_hline(yintercept=2,linetype="dashed")+geom_hline(yintercept=-2,linetype="dashed") +
  labs(title = "Distribution of logFC per cluster CTCF TET vs CTR",
       x = "", y = "log2 Fold Change comparing CTCF TET vs CTR") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-13, 13)) +
  scale_fill_manual(values = pal) + geom_text(aes(x = cell_type, y = 12, label = sig))
ggsave(paste0(monocle_prefix, "CTCF_milo_violin_day7.pdf"))

save(beeswarm_df_SCR_day4, beeswarm_df_GATA4_day4, beeswarm_df_CTCF_day4, beeswarm_df_SCR_day7, beeswarm_df_GATA4_day7, beeswarm_df_CTCF_day7, cluster_summary_day4, cluster_summary_day7,
     file = paste0("/home/shared_folder/RData/milo_output.RData"))
