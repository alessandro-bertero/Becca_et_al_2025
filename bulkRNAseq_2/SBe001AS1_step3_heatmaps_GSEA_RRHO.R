#libraries
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggplotify)
library(readxl)
library(pheatmap)
devtools::install_github("antpiron/RedRibbon")
library(RedRibbon)
library(RNAseqQC)
library(DESeq2)
library(limma)

dir.create("plots")
theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 45, hjust = 1)))

samples = c("#FFB6DB", "#FA0087", "#B6DBFF", "#3283FE", "#808080", "#000000")

rawcounts_filtered <- read.csv("/home/shared_folder/Tables/Counts/rawcounts_filtered.csv")
rownames(rawcounts_filtered) = rawcounts_filtered$gene_ID
table <- read_excel("/home/shared_folder/table_DGEanalysis.xlsx")
get_genotype = function(str){
  return(strsplit(str, "_")[[1]][1])
}
table$genotype = as.character(lapply(table$sample_ID, get_genotype))
get_treatment = function(str){
  return(strsplit(str, "_")[[1]][2])
}
table$treatment = as.character(lapply(table$sample_ID, get_treatment))

meta_AT = table %>% dplyr::filter(cfr_RVvsAT == "AT")
counts_AT = rawcounts_filtered %>% dplyr::select(meta_AT$sample_ID)
dds_AT = make_dds(counts = counts_AT, metadata = meta_AT, ah_record = "AH89426")
vsd_AT = vst(dds_AT)
vsd_AT_corr = vsd_AT
assay(vsd_AT_corr) = limma::removeBatchEffect(assay(vsd_AT_corr), vsd_AT_corr$replicate)

meta_RV = table %>% dplyr::filter(cfr_RVvsAT == "RV")
counts_RV = rawcounts_filtered %>% dplyr::select(meta_RV$sample_ID)
dds_RV = make_dds(counts = counts_RV, metadata = meta_RV, ah_record = "AH89426")
vsd_RV = vst(dds_RV)
vsd_RV_corr = vsd_RV
vsd_RV_corr = vsd_RV
assay(vsd_RV_corr) = limma::removeBatchEffect(assay(vsd_RV_corr), vsd_RV_corr$replicate)

#PCA
pca_res_AT = plot_pca(vsd_AT_corr, show_plot = FALSE, n_feats = 1000)
pca_var_AT = as.data.frame(pca_res_AT[["var_exp"]])
colnames(pca_var_AT) = "pca_var"
pca_var_AT$pca = c(1:length(rownames(meta_AT)))
ggplot(pca_var_AT, aes(x = pca, y = pca_var)) + geom_point(size = 4, shape = 21)
ggsave(filename = paste0("/home/shared_folder/plots/scree_plot_AT.pdf"))

pca_AT = as.data.frame(pca_res_AT[["data"]])
ggplot(pca_AT, aes(x = PC1, y = PC2, shape = treatment)) + geom_point(aes(color = condition), size = 4) +
  scale_color_manual(values = samples) + ggtitle("PCA of the 1000 most variable genes before batch effect correction") +
  xlab(paste0(signif(pca_var_AT$pca_var[1], digits = 2), " %")) + ylab(paste0(signif(pca_var_AT$pca_var[2], digits = 2), " %"))
ggsave(filename = paste0("/home/shared_folder/plots/PCA_1_2_AT.pdf"))

loadings_AT = as.data.frame(pca_res_AT[["loadings"]])
write.csv(loadings_AT, file = paste0("/home/shared_folder/Tables/PCA/loadings_AT.csv"))

pca_res_RV = plot_pca(vsd_RV_corr, show_plot = FALSE, n_feats = 1000)
pca_var_RV = as.data.frame(pca_res_RV[["var_exp"]])
colnames(pca_var_RV) = "pca_var"
pca_var_RV$pca = c(1:length(rownames(meta_RV)))
ggplot(pca_var_RV, aes(x = pca, y = pca_var)) + geom_point(size = 4, shape = 21)
ggsave(filename = paste0("/home/shared_folder/plots/scree_plot_RV.pdf"))

pca_RV = as.data.frame(pca_res_RV[["data"]])
ggplot(pca_RV, aes(x = PC1, y = PC2, shape = treatment)) + geom_point(aes(color = condition), size = 4) +
  scale_color_manual(values = samples) + ggtitle("PCA of the 1000 most variable genes before batch effect correction") +
  xlab(paste0(signif(pca_var_RV$pca_var[1], digits = 2), " %")) + ylab(paste0(signif(pca_var_RV$pca_var[2], digits = 2), " %"))
ggsave(filename = paste0("/home/shared_folder/plots/PCA_1_2_RV.pdf"))

loadings_RV = as.data.frame(pca_res_RV[["loadings"]])
write.csv(loadings_RV, file = paste0("/home/shared_folder/Tables/PCA/loadings_RV.csv"))

#heatmap
TPMnormalized_counts_modules <- read_csv("/home/shared_folder/Table/Counts/TPMnormalized_counts_modules.csv") %>% dplyr::select(-"...1")
TPMcounts = melt(TPMnormalized_counts_modules)
colnames(TPMcounts) = c("gene_ID", "colors", "gene_name", "gene_type", "sample_ID", "TPM")
TPMcounts = TPMcounts %>% dplyr::left_join(table, by = "sample_ID")
markers = c("GATA4", "CTCF", "NEBL", "ACTN2", "TTN", "RYR2", "CTNNA3", "SLC8A1", "CCDC141", "LMO7", "MEF2A", "TGFB2", "CAMK2D", "CDH2", "ACTA2", "RGS5")
TPMmarkers = TPMcounts %>% dplyr::filter(gene_name %in% markers)
TPM_TET = TPMmarkers %>% dplyr::filter(treatment == "TET") %>% dplyr::select(c(gene_name, TPM, chamber, replicate, condition, shRNA))
colnames(TPM_TET) = c("gene_name", "TPM_TET", "chamber", "replicate", "condition_TET", "shRNA")
TPM_ctr = TPMmarkers %>% dplyr::filter(treatment == "ctr") %>% dplyr::select(c(gene_name, TPM, chamber, replicate, condition, shRNA))
colnames(TPM_ctr) = c("gene_name", "TPM_ctr", "chamber", "replicate", "condition_ctr", "shRNA")
ratio = TPM_TET %>% dplyr::left_join(TPM_ctr, by = c("gene_name", "chamber", "replicate", "shRNA")) %>% dplyr::mutate(FC = TPM_TET/TPM_ctr)

ratio_AT = ratio %>% dplyr::filter(chamber == "AT") %>% dplyr::group_by(shRNA, gene_name) %>% dplyr::summarise(ratio = mean(FC))
matrix_AT = dcast(ratio_AT, gene_name ~ shRNA)
rownames(matrix_AT)= matrix_AT$gene_name
matrix_AT = matrix_AT[markers,]
matrix_AT = matrix_AT %>% dplyr::select(-gene_name)

ratio_RV = ratio %>% dplyr::filter(chamber == "RV") %>% dplyr::group_by(shRNA, gene_name) %>% dplyr::summarise(ratio = mean(FC))
matrix_RV = dcast(ratio_RV, gene_name ~ shRNA)
rownames(matrix_RV)= matrix_RV$gene_name
matrix_RV = matrix_RV[markers,]
matrix_RV = matrix_RV %>% dplyr::select(-gene_name)

max = max(c(melt(matrix_AT)$value, melt(matrix_RV)$value))
min = min(c(melt(matrix_AT)$value, melt(matrix_RV)$value))
breaks = seq(min, max, length.out = 100)
colors = viridis::turbo(length(breaks) - 1)

as.ggplot(pheatmap(matrix_AT, breaks = breaks, color = colors, cluster_cols = FALSE, cluster_rows = FALSE))
ggsave("/home/shared_folder/plots/norm_heatmap_AT.pdf")

as.ggplot(pheatmap(matrix_RV, color = colors, breaks = breaks, cluster_cols = FALSE, cluster_rows = FALSE))
ggsave("/home/shared_folder/plots/norm_heatmap_RV.pdf")

##GSEA
ctrlvsCTCF.AT <- read.csv("/home/shared_folder/Tables/DGEs/ctrlvsCTCF.AT.csv", row.names=1) %>% dplyr::select(gene_ID, logFC)
colnames(ctrlvsCTCF.AT) = c("gene_ID", "logFC_CTCF_AT")
ctrlvsCTCF.RV <- read.csv("/home/shared_folder/Tables/DGEs/ctrlvsCTCF.RV.csv", row.names=1) %>% dplyr::select(gene_ID, logFC)
colnames(ctrlvsCTCF.RV) = c("gene_ID", "logFC_CTCF_RV")
ctrlvsGATA4.AT <- read.csv("/home/shared_folder/Tables/DGEs/ctrlvsGATA4.AT.csv", row.names=1) %>% dplyr::select(gene_ID, logFC)
colnames(ctrlvsGATA4.AT) = c("gene_ID", "logFC_GATA4_AT")
ctrlvsGATA4.RV <- read.csv("/home/shared_folder/Tables/DGEs/ctrlvsGATA4.RV.csv", row.names=1) %>% dplyr::select(gene_ID, logFC)
colnames(ctrlvsGATA4.RV) = c("gene_ID", "logFC_GATA4_RV")

#ctcf
library(clusterProfiler)
library(ensembldb)
library("org.Hs.eg.db", character.only = TRUE)
organism = "org.Hs.eg.db"

ctcfAT_FC = -ctrlvsCTCF.AT$logFC
names(ctcfAT_FC) = ctrlvsCTCF.AT$gene_ID
ctcfAT_FC = sort(ctcfAT_FC, decreasing = TRUE)
gse_ctcf_AT <- gseGO(geneList = ctcfAT_FC, 
                  ont ="BP", 
                  keyType = "ENSEMBL",
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 1, 
                  verbose = TRUE, 
                  OrgDb = organism, 
                  pAdjustMethod = "BH")
res_ctcfAT = gse_ctcf_AT@result
res_ctcfAT$n = c(1:length(rownames(res_ctcfAT)))

gseaplot(gse_ctcf_AT, by ="runningScore", title = paste0(gse_ctcf_AT$Description[2787], "\nNES = ", gse_ctcf_AT$NES[2787], "\np-value adjusted = ", gse_ctcf_AT$p.adjust[2787]), geneSetID = 3206, color.line = "blue") +
  theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7))
ggsave(filename = paste0("/home/shared_folder/plots/GSEA_ctcfAT.pdf"))

ctcfRV_FC = -ctrlvsCTCF.RV$logFC
names(ctcfRV_FC) = ctrlvsCTCF.RV$gene_ID
ctcfRV_FC = sort(ctcfRV_FC, decreasing = TRUE)
gse_ctcf_RV <- gseGO(geneList = ctcfRV_FC, 
                  ont ="BP", 
                  keyType = "ENSEMBL",
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 1, 
                  verbose = TRUE, 
                  OrgDb = organism, 
                  pAdjustMethod = "BH")
res_ctcfRV = gse_ctcf_RV@result
res_ctcfRV$n = c(1:length(rownames(res_ctcfRV)))

gseaplot(gse_ctcf_RV, by ="runningScore", title = paste0(gse_ctcf_RV$Description[32], "\nNES = ", gse_ctcf_RV$NES[32], "\np-value adjusted = ", gse_ctcf_RV$p.adjust[32]), geneSetID = 32, color.line = "red") +
  theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7))
ggsave(filename = paste0("/home/shared_folder/plots/GSEA_ctcfRV.pdf"))

#gata4
gata4AT_FC = -ctrlvsGATA4.AT$logFC
names(gata4AT_FC) = ctrlvsGATA4.AT$gene_ID
gata4AT_FC = sort(gata4AT_FC, decreasing = TRUE)
gse_gata4_AT <- gseGO(geneList = gata4AT_FC, 
                   ont ="BP", 
                   keyType = "ENSEMBL",
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = organism, 
                   pAdjustMethod = "BH")
res_gata4AT = gse_gata4_AT@result
res_gata4AT$n = c(1:length(rownames(res_gata4_AT)))

gseaplot(gse_gata4_AT, by ="runningScore", title = paste0(gse_gata4_AT$Description[27], "\nNES = ", gse_gata4_AT$NES[27], "\np-value adjusted = ", gse_gata4_AT$p.adjust[27]), geneSetID = 27, color.line = "blue") +
  theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7))
ggsave(filename = paste0("/home/shared_folder/plots/GSEA_gata4AT.pdf"))

gata4RV_FC = -ctrlvsGATA4.RV$logFC
names(gata4RV_FC) = ctrlvsGATA4.RV$gene_ID
gata4RV_FC = sort(gata4RV_FC, decreasing = TRUE)
gse_gata4_RV <- gseGO(geneList = gata4RV_FC, 
                   ont ="BP", 
                   keyType = "ENSEMBL",
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = organism, 
                   pAdjustMethod = "BH")
res_gata4RV = gse_gata4_RV@result
res_gata4RV$n = c(1:length(rownames(res_gata4RV)))

gseaplot(gse_gata4_RV, by ="runningScore", title = paste0(gse_gata4_RV$Description[9], "\nNES = ", gse_gata4_RV$NES[9], "\np-value adjusted = ", gse_gata4_RV$p.adjust[10]), geneSetID = 10, color.line = "red") +
  theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7))
ggsave(filename = paste0("/home/shared_folder/plots/GSEA_gata4RV.pdf"))


#RRHO
df = ctrlvsCTCF.AT %>% dplyr::left_join(ctrlvsGATA4.AT, by = "gene_ID")
colnames(df) = c("id", "a", "b")
df = df %>% dplyr::mutate(a = -a, b = -b)
rr <- RedRibbon(df, enrichment_mode="hyper-two-tailed")
quad <- quadrants(rr, algorithm="ea", permutation=TRUE, whole=FALSE)
gg <- ggRedRibbon(rr, quadrants=quad, repel.force = 250) + coord_fixed(ratio = 1, clip = "off")
gg
ggsave("/home/shared_folder/plots/AT_RRHO.pdf")

df = ctrlvsCTCF.RV %>% dplyr::left_join(ctrlvsGATA4.RV, by = "gene_ID")
colnames(df) = c("id", "a", "b")
df = df %>% dplyr::mutate(a = -a, b = -b)
rr <- RedRibbon(df, enrichment_mode="hyper-two-tailed")
quad <- quadrants(rr, algorithm="ea", permutation=TRUE, whole=FALSE)
gg <- ggRedRibbon(rr, quadrants=quad, repel.force = 250) + coord_fixed(ratio = 1, clip = "off")
gg
ggsave("/home/shared_folder/plots/RV_RRHO.pdf")

