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
library(ggplotify)

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
dir.create(file.path("/home/shared_folder/Output", format(Sys.Date(), "%y%m%d"), "Seurat"))
seurat_prefix = paste0("/home/shared_folder/Output/", format(Sys.Date(), "%y%m%d"), "/Seurat/")

#palettes
okabe_pal=c("#E69F00","#56B4E9","#009E73","#f0E442","#0072B2","#D55E00","#CC79A7","#000000")
pal = c("#000000","#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff", "#920000","#924900","#db6d00","#24ff24","#ffff6d")


# save input files
matrix_dir ="/home/shared_folder/filtered_feature_bc_matrix"

barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "/features.tsv.gz")
matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")

feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

sample_info = read.csv("/home/shared_folder/aggregation.csv")
sample.names = sample_info$sample_id
sample_info$orig.ident = rownames(sample_info)

data = Read10X(
  data.dir = matrix_dir,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

#create Seurat object
data_seurat <- CreateSeuratObject(counts = data$`Gene Expression`, min.cells = 3, min.features = 200, names.delim = "-", names.field = 2)

#add meta data
OBJ_meta <- data_seurat@meta.data
OBJ_meta$libID = rownames(OBJ_meta)
OBJ_meta = OBJ_meta %>%
  dplyr::left_join(sample_info, by = "orig.ident") %>%
  dplyr::distinct(libID, .keep_all = T)

data_seurat <- AddMetaData(data_seurat, OBJ_meta, col.name = NULL)
pre_selection = data.frame(summary(as.factor(data_seurat@meta.data$sample_id)))

# Visualize QC metrics pre selection
data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
data_seurat[["percent.ribo"]] <- PercentageFeatureSet(data_seurat, pattern = "RPS")
QC=data_seurat[[c("percent.mt", "percent.ribo")]]
QC$libID = rownames(QC)
QC = QC %>% right_join(OBJ_meta, by = "libID")
QC = QC %>% dplyr::mutate(log10GenesPerUMI = log10(nFeature_RNA)/log10(nCount_RNA))
QC = QC %>% dplyr::mutate(selection = ifelse(percent.mt < 15, "kept", "discarded"))

VlnPlot(data_seurat, features = c("nCount_RNA"), group.by = "sample_id")
ggsave(filename = paste0(seurat_prefix, "violin_QC_RNA_nCounts.pdf"), width = 12, height = 7)
VlnPlot(data_seurat, features = c("nFeature_RNA"), group.by = "sample_id")
ggsave(filename = paste0(seurat_prefix, "violin_QC_RNA_features.pdf"), width = 12, height = 7)

QC %>%
  ggplot(aes(x = nCount_RNA, y = percent.mt, color = selection)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "MTPerCounts.pdf"),
       width = 12, height = 8)

QC %>%
  ggplot(aes(x = nCount_RNA, y = percent.ribo, color = selection)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "RiboPerCounts.pdf"),
       width = 12, height = 8)

QC %>%
  ggplot(aes(x = percent.mt, y = percent.ribo, color = selection)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "RiboPerMT.pdf"),
       width = 12, height = 8)

QC %>% 
  ggplot(aes(color = sample_id, x = nCount_RNA, fill = sample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  ylab("Cell density")
ggsave(filename = paste0(seurat_prefix, "density_CountsPerUMI.pdf"),
       width = 12, height = 8)

QC %>% 
  ggplot(aes(color = sample_id, x = nFeature_RNA, fill = sample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  ylab("Cell density")
ggsave(filename = paste0(seurat_prefix, "density_GenesPerUMI.pdf"),
       width = 12, height = 8)

QC %>% 
  ggplot(aes(color = sample_id, x = log10GenesPerUMI, fill = sample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10()
ggsave(filename = paste0(seurat_prefix, "density_log10ratio.pdf"),
       width = 12, height = 8)

QC %>% 
  ggplot(aes(x = sample_id, y = log10(nCount_RNA), fill = sample_id)) + 
  geom_boxplot()
ggsave(filename = paste0(seurat_prefix, "boxplot_CountsPerUMI.pdf"),
       width = 12, height = 8)

QC %>% 
  ggplot(aes(x = sample_id, y = log10(nFeature_RNA), fill = sample_id)) + 
  geom_boxplot()
ggsave(filename = paste0(seurat_prefix, "boxplot_GenesPerUMI.pdf"),
       width = 12, height = 8)

QC %>% 
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10()
ggsave(filename = paste0(seurat_prefix, "scatterplot_GenesAndCountsPerUMI.pdf"),
       width = 12, height = 8)


#filtering1
data_seurat <- subset(data_seurat, subset = percent.mt < 15)

VlnPlot(data_seurat, features = c("nCount_RNA"), group.by = "sample_id")
ggsave(filename = paste0(seurat_prefix, "violin_QC_RNA_nCounts_filtered1.pdf"), width = 12, height = 7)
VlnPlot(data_seurat, features = c("nFeature_RNA"), group.by = "sample_id")
ggsave(filename = paste0(seurat_prefix, "violin_QC_RNA_features_filtered1.pdf"), width = 12, height = 7)

OBJ_meta_2 <- data_seurat@meta.data
OBJ_meta_2$order = c(1:length(rownames(OBJ_meta_2)))
OBJ_meta_2 = OBJ_meta_2 %>% dplyr::arrange(nCount_RNA)
OBJ_meta_2$perc_counts = c(1:length(rownames(OBJ_meta_2)))/length(rownames(OBJ_meta_2)) * 100
threshold1_counts = (OBJ_meta_2 %>% dplyr::filter(perc_counts < 5) %>% dplyr::arrange(desc(perc_counts)))$nCount_RNA[1]
threshold2_counts = (OBJ_meta_2 %>% dplyr::filter(perc_counts > 99) %>% dplyr::arrange(perc_counts))$nCount_RNA[1]

OBJ_meta_2 = OBJ_meta_2 %>% dplyr::arrange(nFeature_RNA)
OBJ_meta_2$perc_genes = c(1:length(rownames(OBJ_meta_2)))/length(rownames(OBJ_meta_2)) * 100
threshold1_genes = (OBJ_meta_2 %>% dplyr::filter(perc_genes < 5) %>% dplyr::arrange(desc(perc_genes)))$nFeature_RNA[1]
threshold2_genes = (OBJ_meta_2 %>% dplyr::filter(perc_genes > 99) %>% dplyr::arrange(perc_genes))$nFeature_RNA[1]

OBJ_meta_2 = OBJ_meta_2 %>% dplyr::arrange(order)

QC_filtered1 = data_seurat[[c("percent.mt", "percent.ribo")]]
QC_filtered1$libID = rownames(QC_filtered1)
QC_filtered1 = QC_filtered1 %>% left_join(OBJ_meta_2, by = c("libID", "percent.mt", "percent.ribo"))
QC_filtered1 = QC_filtered1 %>% dplyr::mutate(log10GenesPerUMI = log10(nFeature_RNA)/log10(nCount_RNA))
QC_filtered1 = QC_filtered1 %>% dplyr::mutate(selection = ifelse(perc_counts < 1, "discarded", ifelse (perc_counts > 99, "discarded", ifelse(perc_genes < 1, "discarded", ifelse(perc_genes > 99, "discarded", "kept")))))

QC_filtered1 %>%
  ggplot(aes(x = nCount_RNA, y = percent.mt, color = selection)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "MTPerCounts_filtered1.pdf"),
       width = 12, height = 8)

QC_filtered1 %>%
  ggplot(aes(x = nCount_RNA, y = percent.ribo, color = selection)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "RiboPerCounts_filtered1.pdf"),
       width = 12, height = 8)

QC_filtered1 %>%
  ggplot(aes(x = percent.mt, y = percent.ribo, color = selection)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "RiboPerMT_filtered1.pdf"),
       width = 12, height = 8)

QC_filtered1 %>% 
  ggplot(aes(color = sample_id, x = nCount_RNA, fill = sample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  ylab("Cell density") +
  geom_vline(aes(xintercept = threshold1_counts), linetype = "dashed") +
  geom_vline(aes(xintercept = threshold2_counts), linetype = "dashed")
ggsave(filename = paste0(seurat_prefix, "density_CountsPerUMI_filtered1.pdf"),
       width = 12, height = 8)

QC_filtered1 %>% 
  ggplot(aes(color = sample_id, x = nFeature_RNA, fill = sample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  ylab("Cell density") +
  geom_vline(aes(xintercept = threshold1_genes), linetype = "dashed") +
  geom_vline(aes(xintercept = threshold2_genes), linetype = "dashed")
ggsave(filename = paste0(seurat_prefix, "density_GenesPerUMI_filtered1.pdf"),
       width = 12, height = 8)

QC_filtered1 %>% 
  ggplot(aes(color = sample_id, x = log10GenesPerUMI, fill = sample_id)) + 
  geom_density(alpha = 0.2) +
  geom_vline(aes(xintercept = log10(threshold1_genes)/log10(threshold1_counts)), linetype = "dashed") +
  geom_vline(aes(xintercept = log10(threshold2_genes)/log10(threshold2_counts)), linetype = "dashed")
  scale_x_log10()
ggsave(filename = paste0(seurat_prefix, "density_log10ratio_filtered1.pdf"),
       width = 12, height = 8)

QC_filtered1 %>% 
  ggplot(aes(x = sample_id, y = log10(nCount_RNA), fill = sample_id)) + 
  geom_hline(aes(yintercept = log10(threshold1_counts)), linetype = "dashed") +
  geom_hline(aes(yintercept = log10(threshold2_counts)), linetype = "dashed") +
  geom_boxplot()
ggsave(filename = paste0(seurat_prefix, "boxplot_CountsPerUMI_filtered1.pdf"),
       width = 12, height = 8)

QC_filtered1 %>% 
  ggplot(aes(x = sample_id, y = log10(nFeature_RNA), fill = sample_id)) + 
  geom_hline(aes(yintercept = log10(threshold1_genes)), linetype = "dashed") +
  geom_hline(aes(yintercept = log10(threshold2_genes)), linetype = "dashed") +
  geom_boxplot()
ggsave(filename = paste0(seurat_prefix, "boxplot_GenesPerUMI_filtered1.pdf"),
       width = 12, height = 8)

QC_filtered1 %>% 
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  geom_vline(aes(xintercept = threshold1_counts), linetype = "dashed") +
  geom_vline(aes(xintercept = threshold2_counts), linetype = "dashed") +
  geom_hline(aes(yintercept = threshold1_genes), linetype = "dashed") +
  geom_hline(aes(yintercept = threshold2_genes), linetype = "dashed") +
  scale_x_log10() + 
  scale_y_log10()
ggsave(filename = paste0(seurat_prefix, "scatterplot_GenesAndCountsPerUMI_filtered1.pdf"),
       width = 12, height = 8)

post_selection1 = data.frame(summary(as.factor(data_seurat@meta.data$sample_id)))

#filtering2
data_seurat <- AddMetaData(data_seurat, OBJ_meta_2, col.name = NULL)
data_seurat <- subset(data_seurat, subset = perc_counts < 99 & perc_counts > 5 & perc_genes < 99 & perc_genes > 5)

VlnPlot(data_seurat, features = c("nCount_RNA"), group.by = "sample_id")
ggsave(filename = paste0(seurat_prefix, "violin_QC_RNA_nCounts_filtered2.pdf"), width = 12, height = 7)
VlnPlot(data_seurat, features = c("nFeature_RNA"), group.by = "sample_id")
ggsave(filename = paste0(seurat_prefix, "violin_QC_RNA_features_filtered2.pdf"), width = 12, height = 7)

QC_filtered2 = data_seurat[[c("percent.mt", "percent.ribo")]]
QC_filtered2$libID = rownames(QC_filtered2)
QC_filtered2 = QC_filtered2 %>% left_join(OBJ_meta_2, by = c("libID", "percent.mt", "percent.ribo"))
QC_filtered2 = QC_filtered2 %>% dplyr::mutate(log10GenesPerUMI = log10(nFeature_RNA)/log10(nCount_RNA))

QC_filtered2 %>%
  ggplot(aes(x = nCount_RNA, y = percent.mt, color = sample_id)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "MTPerCounts_filtered2.pdf"),
       width = 12, height = 8)

QC_filtered2 %>%
  ggplot(aes(x = nCount_RNA, y = percent.ribo, color = sample_id)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "RiboPerCounts_filtered2.pdf"),
       width = 12, height = 8)

QC_filtered2 %>%
  ggplot(aes(x = percent.mt, y = percent.ribo, color = sample_id)) +
  geom_point()
ggsave(filename = paste0(seurat_prefix, "RiboPerMT_filtered2.pdf"),
       width = 12, height = 8)

QC_filtered2 %>% 
  ggplot(aes(color = sample_id, x = nCount_RNA, fill = sample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  ylab("Cell density")
ggsave(filename = paste0(seurat_prefix, "density_CountsPerUMI_filtered2.pdf"),
       width = 12, height = 8)

QC_filtered2 %>% 
  ggplot(aes(color = sample_id, x = nFeature_RNA, fill = sample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  ylab("Cell density")
ggsave(filename = paste0(seurat_prefix, "density_GenesPerUMI_filtered2.pdf"),
       width = 12, height = 8)

QC_filtered2 %>% 
  ggplot(aes(color = sample_id, x = log10GenesPerUMI, fill = sample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10()
ggsave(filename = paste0(seurat_prefix, "density_log10ratio_filtered2.pdf"),
       width = 12, height = 8)

QC_filtered2 %>% 
  ggplot(aes(x = sample_id, y = log10(nCount_RNA), fill = sample_id)) + 
  geom_boxplot()
ggsave(filename = paste0(seurat_prefix, "boxplot_CountsPerUMI_filtered2.pdf"),
       width = 12, height = 8)

QC_filtered2 %>% 
  ggplot(aes(x = sample_id, y = log10(nFeature_RNA), fill = sample_id)) + 
  geom_boxplot()
ggsave(filename = paste0(seurat_prefix, "boxplot_GenesPerUMI_filtered2.pdf"),
       width = 12, height = 8)

QC_filtered2 %>% 
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10()
ggsave(filename = paste0(seurat_prefix, "scatterplot_GenesAndCountsPerUMI_filtered2.pdf"),
       width = 12, height = 8)

post_selection2 = data.frame(summary(as.factor(data_seurat@meta.data$sample_id)))

#results
results = cbind(pre_selection, post_selection1, post_selection2)
colnames(results) = c("1unfiltered", "filtered1", "filtered2")
results$dataset = rownames(results)
melt = melt(results)
melt = melt %>% dplyr::mutate(sample = paste0(dataset, "_", variable))
ggplot(melt) + geom_boxplot(aes(x = sample, y = value, color = variable))
ggsave(filename = paste0(seurat_prefix, "filtering_results.pdf"),
       width = 12, height = 8)

#save
counts_post = GetAssayData(object = data_seurat, slot = "counts")
save(data_seurat, counts_post, QC_filtered2,
     file = paste0("/home/shared_folder/RData/seurat_filtered2.RData"))

rm(counts_post)

# normalize data
data_seurat_norm = NormalizeData(data_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# original expression distribution
raw_geneExp = as.vector(data_seurat_norm[['RNA']]@layers[["counts"]]) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
png(paste0(seurat_prefix, "histRaw.pdf"))
hist(raw_geneExp)
dev.off()

# expression distribution after normalization
logNorm_geneExp = as.vector(data_seurat_norm[['RNA']]@layers[["data"]]) %>% sample(10000)
logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
png(paste0(seurat_prefix, "histNorm.pdf"))
hist(logNorm_geneExp)
dev.off()

# find highly variable genes
data_seurat_norm = FindVariableFeatures(data_seurat_norm, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(data_seurat_norm), 20)

# plot variable features
plot1 = VariableFeaturePlot(data_seurat_norm) + 
  theme(legend.position="top")
LabelPoints(plot = plot1, points = top20, repel = TRUE) + 
  theme(legend.position="none")

ggsave(filename = paste0(seurat_prefix, "variableFeatures.pdf"),
       width = 12, height = 8)

#scaling
all.genes = rownames(data_seurat_norm)
data_seurat_norm = ScaleData(data_seurat_norm, features = all.genes)

logScaled_geneExp = as.vector(data_seurat_norm@assays[["RNA"]]@layers[["scale.data"]]) %>% sample(10000)
logScaled_geneExp = logScaled_geneExp[logScaled_geneExp != 0]
png(paste0(seurat_prefix, "histScaled.pdf"))
hist(logScaled_geneExp)
dev.off()

#PCA
data_seurat_norm = RunPCA(data_seurat_norm, features = VariableFeatures(object = data_seurat_norm))

as.ggplot(ElbowPlot(data_seurat_norm))
ggsave(filename = paste0(seurat_prefix, "ElbowPlot.pdf"),
       width = 12, height = 8)

as.ggplot(VizDimLoadings(data_seurat_norm, dims = 1:6, nfeatures = 30, col =,"#b66dff", reduction = "pca",
                         projected = FALSE, balanced = FALSE, ncol = NULL, combine = TRUE))

ggsave(filename = paste0(seurat_prefix, "VizDimLoadings.pdf"),
       width = 12, height = 20)

as.ggplot(DimPlot(data_seurat_norm, dims = c(1,2), reduction = "pca", group.by = c("sample_id")))
ggsave(filename = paste0(seurat_prefix, "DimPlot.pdf"),
       width = 12, height = 8)

#cell cycle
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

data_seurat_norm = CellCycleScoring(data_seurat_norm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

as.ggplot(RidgePlot(data_seurat_norm, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, log = TRUE, group.by = "sample_id"))
ggsave(filename = paste0(seurat_prefix, "CellCycle.pdf"),
       width = 12, height = 8)

#umap
data_seurat_norm <- FindNeighbors(data_seurat_norm, dims = 1:20)
data_seurat_norm <- FindClusters(data_seurat_norm, resolution = 0.4)
data_seurat_norm <- RunUMAP(data_seurat_norm, dims = 1:20)

as.ggplot(DimPlot(data_seurat_norm, reduction = "umap"))
ggsave(filename = paste0(seurat_prefix, "UMAP_all.pdf"),
       width = 20, height = 20)
as.ggplot(DimPlot(data_seurat_norm, reduction = "umap", group.by = "sample_id"))
ggsave(filename = paste0(seurat_prefix, "UMAP_samples.pdf"),
       width = 20, height = 20)
as.ggplot(DimPlot(data_seurat_norm, reduction = "umap", split.by = "sample_id"))
ggsave(filename = paste0(seurat_prefix, "UMAP_clusters_samples.pdf"),
       width = 12, height = 8)
as.ggplot(DimPlot(data_seurat_norm, reduction = "umap", group.by = c("Phase")))
ggsave(filename = paste0(seurat_prefix, "UMAP_cellCycle_all.pdf"),
       width = 20, height = 20)
as.ggplot(DimPlot(data_seurat_norm, reduction = "umap", split.by = "sample_id", group.by = c("Phase")))
ggsave(filename = paste0(seurat_prefix, "UMAP_cellCycle_samples.pdf"),
       width = 20, height = 20)

#save
meta_norm = as.data.frame(data_seurat_norm@meta.data)
save(meta_norm,
     file = paste0("/home/shared_folder/RData/seurat_meta_norm.RData"))