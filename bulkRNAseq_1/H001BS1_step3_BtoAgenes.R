library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(ensembldb)
library(topGO)
library(limma)
library(org.Hs.eg.db)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(gtools)
library(readxl)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size=15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle=45, hjust = 1)))

day00.PC <- read.csv("/home/shared_folder/H9_data/day00.PC.bedGraph", header=FALSE, skip = 4, sep = "\t")
day00.PC = day00.PC %>% dplyr::mutate(V5 = paste0(V1, "_", (V3 - V2)/2 + V2)) %>%
  dplyr::mutate(V6 = ifelse(V4 > 0, "A", "B")) %>% dplyr::select(c(V5, V4, V6))
colnames(day00.PC) = c("bin", "1", "day00_comp")

day80.PC <- read.csv("/home/shared_folder/H9_data/day80.PC.bedGraph", header=FALSE, skip = 4, sep = "\t")
day80.PC = day80.PC %>% dplyr::mutate(V5 = paste0(V1, "_", (V3 - V2)/2 + V2)) %>%
  dplyr::mutate(V6 = ifelse(V4 > 0, "A", "B")) %>% dplyr::select(c(V5, V4, V6))
colnames(day80.PC) = c("bin", "6", "day80_comp")

PC_all = day00.PC %>% dplyr::left_join(day80.PC, by = "bin")

geneInfo <- read.table("/home/shared_folder/geneInfo.tab", quote="\"", skip = 1)
colnames(geneInfo) = c("geneID", "gene_name", "gene_type")

#tss annotation
tss <- read_csv("/home/shared_folder/H9_data/tss_hg19.csv")
tss = tss[, c(2:8)]
colnames(tss) = c("chr", "start", "end", "width", "strand", "gene_type", "geneID")
tss = tss %>% dplyr::mutate(tss = ifelse(strand == "+", start, end))
tss = tss %>% dplyr::select(c(geneID, chr, tss))
tss = tss %>% dplyr::left_join(geneInfo, by = "geneID") %>% dplyr::filter(gene_type == "protein_coding")

regions <- read.delim("/home/shared_folder/regions_abs.bed", header=FALSE)
colnames(regions) = c("chr", "start", "end", "n")
regions = regions %>% dplyr::mutate(bin = paste0(chr, "_", (end - start)/2 + start))
genes = c()
for(i in 1:length(rownames(regions))){
  tss_i = tss %>% dplyr::filter(chr == regions$chr[i])
  genes_i = ""
  for(j in 1:length(rownames(tss_i))){
    if(regions$start[i] < tss_i$tss[j] & regions$end[i] > tss_i$tss[j]){genes_i = paste0(genes_i, tss_i$gene_name[j], "_")}
  }
  genes = append(genes, genes_i)
}
regions$genes = genes
regions = regions %>% dplyr::select(c(bin, genes))
PC_all = PC_all %>% dplyr::left_join(regions, by = "bin")

#filtering of unidirectional changing B to A bins
melt = melt(PC_all)
dir.create("/home/shared_folder/plots")
melt$variable = as.numeric(melt$variable)
ggplot(melt, aes(x = variable, y = value, color = bin)) + geom_line() + geom_hline(aes(yintercept = 0), linetype = "dashed") + theme(legend.position="none")
ggsave("/home/shared_folder/plots/pre_selection.pdf")

PC_filtered = PC_all %>% dplyr::filter(day00_comp == "B" & day80_comp == "A")
melt = melt(PC_filtered)
melt$variable = as.numeric(melt$variable)
ggplot(melt, aes(x = variable, y = value, color = bin)) + geom_line() + geom_hline(aes(yintercept = 0), linetype = "dashed") + theme(legend.position="none")
ggsave("/home/shared_folder/plots/post_selection.pdf")

#significancy
dir.create("/home/shared_folder/plots")
sign = read.csv("/home/shared_folder/H9_data/differential_compartment.log10Padj.bedGraph", header = FALSE, skip = 4, sep = "\t")
sign = sign %>% dplyr::mutate(V5 = paste0(V1, "_", (V3 - V2)/2 + V2)) %>% dplyr::select(c(V5, V4))
colnames(sign) = c("bin", "pval")

PC_all = PC_all %>% dplyr::mutate(delta = PC_day80 - PC_day00) %>% dplyr::left_join(regions, by = "bin")
PC_all[is.na(PC_all)] = ""
PC_all = PC_all %>% dplyr::left_join(sign, by = "bin")
PC_all[is.na(PC_all)] = 0
PC_all = PC_all %>% dplyr::mutate(BtoA = ifelse(day00_comp == "B" & day80_comp == "A" & pval > -log10(0.05), "yes", "no")) %>% dplyr::mutate(BtoA_genes = ifelse(BtoA == "yes", genes, ""))

fun <- function(int) {
  return (2*(int)/(int-2))
}

PC_BtoA = PC_all %>% dplyr::filter(day00_comp == "B" & day80_comp == "A") %>% dplyr::mutate(BtoA = ifelse(pval > -log10(0.05) & delta > 2 & pval > fun(delta), "yes", "no")) %>%
  dplyr::mutate(BtoA_genes = ifelse(BtoA == "yes", genes, "")) %>% dplyr::filter(pval != 0)
ggplot(PC_BtoA, aes(x = delta, y = pval)) + geom_point(aes(fill = BtoA, color = BtoA), shape = 21, size = 4) +
  geom_vline(aes(xintercept = 2), linetype = "dashed") + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = BtoA_genes), max.overlaps = 60, size = 3, hjust = 1.2) +
  scale_fill_manual(values = c("grey", "#56B4E9")) + scale_color_manual(values = c("grey", "#006ddb")) +
  geom_function(fun = fun) + scale_y_continuous(limits = c(0, 210)) + scale_x_continuous(limits = c(1.5, 4))
ggsave("/home/shared_folder/plots/volcano.pdf")

#genes_overlapping
genes_BtoA = c()
for(i in PC_BtoA$BtoA_genes){
  x = strsplit(i, "_")
  for(j in x){genes_BtoA = append(genes_BtoA, j)}
}

BtoA_genes_table = geneInfo %>% dplyr::filter(gene_name %in% genes_BtoA)
write.csv(BtoA_genes_table, file = "/home/shared_folder/Tables/BtoA_genes.csv")

GO = ensembldb::select(org.Hs.eg.db, keys = genes_BtoA, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
go = goana(GO$ENTREZID, species = "Hs", convert = T)
res = topGO(go, ontology = "BP", number = 30, truncate.term = 50)
ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
  labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of BP of B to A genes"))
ggsave(filename = "/home/shared_folder/plots/BtoA_genes.pdf", width = 14, height = 20)

#H9 bulk
dir.create("/home/shared_folder/bulk")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3262nnn/GSM3262968/suppl/GSM3262968%5FRZY632%5FRNA%5FD0%5F1%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262969&format=file&file=GSM3262969%5FRZY638%5FRNA%5FD0%5F2%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262970&format=file&file=GSM3262970%5FRZY633%5FRNA%5FD2%5F1%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262971&format=file&file=GSM3262971%5FRZY639%5FRNA%5FD2%5F2%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262972&format=file&file=GSM3262972%5FRZY634%5FRNA%5FD5%5F1%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262973&format=file&file=GSM3262973%5FRZY640%5FRNA%5FD5%5F2%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262974&format=file&file=GSM3262974%5FRZY635%5FRNA%5FD7%5F1%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262975&format=file&file=GSM3262975%5FRZY641%5FRNA%5FD7%5F2%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262976&format=file&file=GSM3262976%5FRZY636%5FRNA%5FD15%5F1%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262977&format=file&file=GSM3262977%5FRZY642%5FRNA%5FD15%5F2%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262978&format=file&file=GSM3262978%5FRZY637%5FRNA%5FD80%5F1%2DchrM%2Erpkm%2Egz")
system2("wget", "-P /home/shared_folder/bulk ftp://ftp.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262979&format=file&file=GSM3262979%5FRZY643%5FRNA%5FD80%5F2%2DchrM%2Erpkm%2Egz")

system2("gzip", "-dk /home/shared_folder/bulk/GSM3262968_RZY632_RNA_D0_1-chrM.rpkm")
system2("gzip", "-dk /home/shared_folder/bulk/GSM3262969_RZY638_RNA_D0_2-chrM.rpkm")
bulk_day00_R1 <- read.delim("/home/shared_folder/bulk/GSM3262968_RZY632_RNA_D0_1-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY632_RNA_D0_1.nodup.bam"))
colnames(bulk_day00_R1) = c("gene_name", "length", "RPKM")
bulk_day00_R1 = bulk_day00_R1 %>% mutate(day00_R1_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day00_R1_TPM"))
bulk_day00_R2 <- read.delim("/home/shared_folder/bulk/GSM3262969_RZY638_RNA_D0_2-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY638_RNA_D0_2.nodup.bam"))
colnames(bulk_day00_R2) = c("gene_name", "length", "RPKM")
bulk_day00_R2 = bulk_day00_R2 %>% mutate(day00_R2_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day00_R2_TPM"))

system2("gzip", "-dk /home/shared_folder/bulk/GSM3262970_RZY633_RNA_D2_1-chrM.rpkm")
system2("gzip", "-dk /home/shared_folder/bulk/GSM3262971_RZY639_RNA_D2_2-chrM.rpkm")
bulk_day02_R1 <- read.delim("/home/shared_folder/bulk/GSM3262970_RZY633_RNA_D2_1-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY633_RNA_D2_1.nodup.bam"))
colnames(bulk_day02_R1) = c("gene_name", "length", "RPKM")
bulk_day02_R1 = bulk_day02_R1 %>% mutate(day02_R1_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day02_R1_TPM"))
bulk_day02_R2 <- read.delim("/home/shared_folder/bulk/GSM3262971_RZY639_RNA_D2_2-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY639_RNA_D2_2.nodup.bam"))
colnames(bulk_day02_R2) = c("gene_name", "length", "RPKM")
bulk_day02_R2 = bulk_day02_R2 %>% mutate(day02_R2_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day02_R2_TPM"))

system2("gzip", "-dk /home/shared_folder/bulk/GSM3262972_RZY634_RNA_D5_1-chrM.rpkm")
system2("gzip", "-dk /home/shared_folder/bulk/GSM3262973_RZY640_RNA_D5_2-chrM.rpkm")
bulk_day05_R1 <- read.delim("/home/shared_folder/bulk/GSM3262972_RZY634_RNA_D5_1-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY634_RNA_D5_1.nodup.bam"))
colnames(bulk_day05_R1) = c("gene_name", "length", "RPKM")
bulk_day05_R1 = bulk_day05_R1 %>% mutate(day05_R1_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day05_R1_TPM"))
bulk_day05_R2 <- read.delim("/home/shared_folder/bulk/GSM3262973_RZY640_RNA_D5_2-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY640_RNA_D5_2.nodup.bam"))
colnames(bulk_day05_R2) = c("gene_name", "length", "RPKM")
bulk_day05_R2 = bulk_day05_R2 %>% mutate(day05_R2_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day05_R2_TPM"))

system2("gzip", "-dk /home/shared_folder/bulk/GSM3262974_RZY635_RNA_D7_1-chrM.rpkm")
system2("gzip", "-dk /home/shared_folder/bulk/GSM3262975_RZY641_RNA_D7_2-chrM.rpkm")
bulk_day07_R1 <- read.delim("/home/shared_folder/bulk/GSM3262974_RZY635_RNA_D7_1-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY635_RNA_D7_1.nodup.bam"))
colnames(bulk_day07_R1) = c("gene_name", "length", "RPKM")
bulk_day07_R1 = bulk_day07_R1 %>% mutate(day07_R1_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day07_R1_TPM"))
bulk_day07_R2 <- read.delim("/home/shared_folder/bulk/GSM3262975_RZY641_RNA_D7_2-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY641_RNA_D7_2.nodup.bam"))
colnames(bulk_day07_R2) = c("gene_name", "length", "RPKM")
bulk_day07_R2 = bulk_day07_R2 %>% mutate(day07_R2_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day07_R2_TPM"))

system2("gzip", "-dk /home/shared_folder/bulk/GSM3262976_RZY636_RNA_D15_1-chrM.rpkm")
system2("gzip", "-dk /home/shared_folder/bulk/GSM3262977_RZY642_RNA_D15_2-chrM.rpkm")
bulk_day15_R1 <- read.delim("/home/shared_folder/bulk/GSM3262976_RZY636_RNA_D15_1-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY636_RNA_D15_1.nodup.bam"))
colnames(bulk_day15_R1) = c("gene_name", "length", "RPKM")
bulk_day15_R1 = bulk_day15_R1 %>% mutate(day15_R1_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day15_R1_TPM"))
bulk_day15_R2 <- read.delim("/home/shared_folder/bulk/GSM3262977_RZY642_RNA_D15_2-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY642_RNA_D15_2.nodup.bam"))
colnames(bulk_day15_R2) = c("gene_name", "length", "RPKM")
bulk_day15_R2 = bulk_day15_R2 %>% mutate(day15_R2_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day15_R2_TPM"))

system2("gzip", "-dk /home/shared_folder/bulk/GSM3262978_RZY637_RNA_D80_1-chrM.rpkm")
system2("gzip", "-dk /home/shared_folder/bulk/GSM3262979_RZY643_RNA_D80_2-chrM.rpkm")
bulk_day80_R1 <- read.delim("/home/shared_folder/bulk/GSM3262978_RZY637_RNA_D80_1-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY637_RNA_D80_1.nodup.bam"))
colnames(bulk_day80_R1) = c("gene_name", "length", "RPKM")
bulk_day80_R1 = bulk_day80_R1 %>% mutate(day80_R1_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day80_R1_TPM"))
bulk_day80_R2 <- read.delim("/home/shared_folder/bulk/GSM3262979_RZY643_RNA_D80_2-chrM.rpkm") %>% dplyr::select(c("Geneid", "Length", "bam.RZY643_RNA_D80_2.nodup.bam"))
colnames(bulk_day80_R2) = c("gene_name", "length", "RPKM")
bulk_day80_R2 = bulk_day80_R2 %>% mutate(day80_R2_TPM = (length * RPKM)/sum(RPKM) * 10^3) %>% dplyr::select(c("gene_name", "day80_R2_TPM"))

bulk = geneInfo %>% dplyr::left_join(bulk_day00_R1, by = "gene_name") %>% dplyr::left_join(bulk_day00_R2, by = "gene_name") %>% dplyr::left_join(bulk_day02_R1, by = "gene_name") %>% dplyr::left_join(bulk_day02_R2, by = "gene_name") %>%
  dplyr::left_join(bulk_day05_R1, by = "gene_name") %>% dplyr::left_join(bulk_day05_R2, by = "gene_name") %>% dplyr::left_join(bulk_day07_R1, by = "gene_name") %>% dplyr::left_join(bulk_day07_R2, by = "gene_name") %>%
  dplyr::left_join(bulk_day15_R1, by = "gene_name") %>% dplyr::left_join(bulk_day15_R2, by = "gene_name") %>% dplyr::left_join(bulk_day80_R1, by = "gene_name") %>% dplyr::left_join(bulk_day80_R2, by = "gene_name")
bulk = bulk %>% dplyr::filter(gene_type == "protein_coding") %>% dplyr::filter(gene_name %in% bulk_day00_R1$gene_name)

rep = c("R1", "R2", "R1", "R2", "R1", "R2", "R1", "R2", "R1", "R2", "R1", "R2")
rownames(bulk) = bulk$geneID
bulk = bulk %>% dplyr::select(-c(colnames(geneInfo)))
bulk = as.data.frame(limma::removeBatchEffect(bulk, rep))
bulk$geneID = rownames(bulk)
bulk = bulk %>% dplyr::left_join(geneInfo, by = "geneID")

bulk_sign = bulk %>% dplyr::filter(gene_name %in% genes_BtoA)
bulk_sign = bulk_sign %>% dplyr::mutate(sum = day00_R1_TPM + day00_R2_TPM + day02_R1_TPM + day02_R2_TPM + day05_R1_TPM + day05_R2_TPM + day07_R1_TPM + day07_R2_TPM + day15_R1_TPM + day15_R2_TPM + day80_R1_TPM + day80_R2_TPM) %>%
  dplyr::filter(sum != 0)
rownames(bulk_sign) = bulk_sign$gene_name
bulk_sign = bulk_sign %>% dplyr::select(-colnames(geneInfo)) %>% dplyr::select(-sum)
as.ggplot(pheatmap(bulk_sign, scale = "row", cluster_cols = FALSE, color = viridis::viridis(n = 25, option = "H")))
ggsave(filename = "/home/shared_folder/plots/heatmap.pdf", width = 14, height = 20)

#z_score
cal_z_score = function(x){(x - mean(x)) / sd(x)}
bulk_matrix = bulk
bulk_matrix = bulk_matrix %>% dplyr::mutate(sum = day00_R1_TPM + day00_R2_TPM + day02_R1_TPM + day02_R2_TPM + day05_R1_TPM + day05_R2_TPM + day07_R1_TPM + day07_R2_TPM + day15_R1_TPM + day15_R2_TPM + day80_R1_TPM + day80_R2_TPM) %>%
  dplyr::filter(sum != 0)
rownames(bulk_matrix) = bulk_matrix$geneID
bulk_matrix = bulk_matrix %>% dplyr::select(-colnames(geneInfo)) %>% dplyr::select(-sum)
bulk_matrix = as.data.frame(t(apply(as.data.frame(bulk_matrix), 1, cal_z_score)))
bulk_matrix = bulk_matrix %>% dplyr::mutate(day00 = (day00_R1_TPM + day00_R2_TPM)/2)
bulk_matrix = bulk_matrix %>% dplyr::mutate(day02 = (day02_R1_TPM + day02_R2_TPM)/2)
bulk_matrix = bulk_matrix %>% dplyr::mutate(day05 = (day05_R1_TPM + day05_R2_TPM)/2)
bulk_matrix = bulk_matrix %>% dplyr::mutate(day07 = (day07_R1_TPM + day07_R2_TPM)/2)
bulk_matrix = bulk_matrix %>% dplyr::mutate(day15 = (day15_R1_TPM + day15_R2_TPM)/2)
bulk_matrix = bulk_matrix %>% dplyr::mutate(day80 = (day80_R1_TPM + day80_R2_TPM)/2)
bulk_matrix$geneID = rownames(bulk_matrix)
bulk_matrix = bulk_matrix %>% dplyr::left_join(geneInfo, by = "geneID")

bulk_z = bulk_matrix
rownames(bulk_z) = bulk_z$geneID
bulk_z = bulk_z %>% dplyr::select(c(day00_R1_TPM, day00_R2_TPM, day02_R1_TPM, day02_R2_TPM, day05_R1_TPM, day05_R2_TPM, day07_R1_TPM, day07_R2_TPM, day15_R1_TPM, day15_R2_TPM, day80_R1_TPM, day80_R2_TPM))
bulk_z = as.data.frame(apply(as.data.frame(bulk_z), 2, cal_z_score))
bulk_z = bulk_z %>% dplyr::mutate(day00_z = (day00_R1_TPM + day00_R2_TPM)/2)
bulk_z = bulk_z %>% dplyr::mutate(day02_z = (day02_R1_TPM + day02_R2_TPM)/2)
bulk_z = bulk_z %>% dplyr::mutate(day05_z = (day05_R1_TPM + day05_R2_TPM)/2)
bulk_z = bulk_z %>% dplyr::mutate(day07_z = (day07_R1_TPM + day07_R2_TPM)/2)
bulk_z = bulk_z %>% dplyr::mutate(day15_z = (day15_R1_TPM + day15_R2_TPM)/2)
bulk_z = bulk_z %>% dplyr::mutate(day80_z = (day80_R1_TPM + day80_R2_TPM)/2)
bulk_z$geneID = rownames(bulk_z)

top_z = c()
for(i in bulk_z$geneID){
  bulk_i = bulk_z %>% dplyr::filter(geneID == i)
  list_i = c(bulk_i$day00_z, bulk_i$day02_z, bulk_i$day05_z, bulk_i$day07_z, bulk_i$day15_z, bulk_i$day80_z)
  names(list_i) = c("day00", "day02", "day05", "day07", "day15", "day80")
  list_i = sort(list_i)
  top = names(list_i)[6]
  top_z = append(top_z, top)
}

bulk_z$top = top_z
top_z = bulk_z %>% dplyr::select(c(geneID, top))

BtoA_genes = BtoA_genes %>% dplyr::left_join(top_z, by = "geneID")
ggplot(BtoA_genes) + geom_bar(aes(x = top, stat="count", fill = top))
ggsave("/home/shared_folder/plots/peaks_distribution.pdf", width = 7, height = 12)

annotation_row = BtoA_genes %>% dplyr::filter(gene_name %in% rownames(bulk_sign))
rownames(annotation_row) = annotation_row$gene_name
annotation_row = annotation_row %>% dplyr::select(top)

as.ggplot(pheatmap(bulk_sign, scale = "row", cluster_cols = FALSE, annotation_row = annotation_row, color = viridis::viridis(n = 25, option = "H")))
ggsave(filename = "/home/shared_folder/plots/heatmap_annotated.pdf", width = 14, height = 20)

#CTCF GATA4 analysis
table <- read_excel("/home/shared_folder/table_DGEanalysis.xlsx")
meta = table %>% dplyr::select(c("sample_ID", "condition", "replicate"))

counts_all = geneInfo
for(i in 1:length(rownames(table))){
  sample_i = read.table(table$counts[i], quote = "\"", comment.char = "", skip = 4) %>% dplyr::select("V1", "V4");
  colnames(sample_i) = c("geneID", table$sample_ID[i]);
  counts_all = counts_all %>% dplyr::left_join(sample_i, by = "geneID")
}

gl = read_csv("/home/shared_folder/gene_length.txt")
gl = as.data.frame(gl)
colnames(gl) = c("gene_length")
hgnc_symbol = c()
transcript_length = c()
for(i in gl$gene_length){row = strsplit(i, split = "\t"); hgnc_symbol = append(hgnc_symbol, row[[1]][1]); transcript_length = append(transcript_length, row[[1]][2])}
geneLength = data.frame(hgnc_symbol, transcript_length)
rownames(counts_all) = counts_all$geneID
counts_matrix = counts_all %>% dplyr::filter(gene_type == "protein_coding") %>% dplyr::select(-colnames(geneInfo))
counts_tpm = as.data.frame(ADImpute::NormalizeTPM(as.matrix(counts_matrix), sce = NULL, log = FALSE, tr_length = geneLength, scale = 1))
counts_tpm = as.data.frame(limma::removeBatchEffect(counts_tpm, meta$replicate))
counts_tpm$geneID = rownames(counts_tpm)
counts_tpm = counts_tpm %>% dplyr::left_join(geneInfo, by = "geneID")
counts_tpm = counts_tpm %>% dplyr::mutate(sum = SCR_ctr_1 + SCR_ctr_2 + SCR_tet_1 + SCR_tet_2 + GATA4_ctr_1 + GATA4_ctr_2 + GATA4_tet_1 + GATA4_tet_2 + CTCF_ctr_1 + CTCF_ctr_2 + CTCF_tet_1 + CTCF_tet_2) %>%
  dplyr::filter(sum != 0) %>% dplyr::select(-sum)

counts_BtoA = counts_tpm %>% dplyr::filter(geneID %in% BtoA_genes2$geneID) %>% dplyr::filter(gene_name %in% rownames(annotation_row))
annotation_row = BtoA_genes2 %>% dplyr::filter(gene_name %in% rownames(counts_BtoA))
rownames(annotation_row) = annotation_row$gene_name

counts_BtoA = counts_BtoA %>% dplyr::left_join(annotation_row, by = colnames(geneInfo)) %>% dplyr::arrange(top) %>% dplyr::select(-top)
annotation_row = annotation_row %>% dplyr::select(top)
rownames(counts_BtoA) = counts_BtoA$gene_name
counts_BtoA = counts_BtoA %>% dplyr::select(-colnames(geneInfo))

as.ggplot(pheatmap(counts_BtoA, scale = "row", cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = annotation_row, color = viridis::viridis(n = 25, option = "H")))
ggsave(filename = "/home/shared_folder/plots/CTCF_GATA4_all.pdf", width = 14, height = 20)

genes_day00 = (BtoA_genes2 %>% dplyr::filter(top == "day00"))$gene_name
genes_day02 = (BtoA_genes2 %>% dplyr::filter(top == "day02"))$gene_name
genes_day05 = (BtoA_genes2 %>% dplyr::filter(top == "day05"))$gene_name
genes_day07 = (BtoA_genes2 %>% dplyr::filter(top == "day07"))$gene_name
genes_day15 = (BtoA_genes2 %>% dplyr::filter(top == "day15"))$gene_name
genes_day80 = (BtoA_genes2 %>% dplyr::filter(top == "day80"))$gene_name

#gsea
CTCFtetvsCTCFctr <- read_csv("/home/shared_folder/Tables/DGEs/CTCFtetvsCTCFctr.csv")
GATA4tetvsGATA4ctr <- read_csv("/home/shared_folder/Tables/DGEs/GATA4tetvsGATA4ctr.csv")
SCRtetvsSCRctr <- read_csv("/home/shared_folder/Tables/DGEs/SCRtetvsSCRctr.csv")

#gsea
library(fgsea)
genes_late = c(genes_day15, genes_day80)
genes_mid = c(genes_day05, genes_day07)
genes_early = c(genes_day00, genes_day02)

pathways = list(BtoA_early = genes_early, 
                BtoA_mid = genes_mid,
                BtoA_late = genes_late)

#scr
ranks = SCRtetvsSCRctr$t
names(ranks) = SCRtetvsSCRctr$gene_name
fgseSCR <- fgsea(pathways = pathways, 
                 stats    = ranks,
                 scoreType = 'std',
                 minSize  = 1,
                 maxSize  = 500)

pdf(file = paste0( "/home/shared_folder/plots/gsea_SCR_grouped.pdf"), width = 20, height = 15)
plotGseaTable(pathways, stats = ranks, fgseaRes = fgseSCR, gseaParam = 0.5)
dev.off()

plotEnrichment(pathways[["BtoA_early"]], ranks) + labs(title = "SCR B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_SCR_early.pdf",width = 14, height = 8)

plotEnrichment(pathways[["BtoA_mid"]], ranks) + labs(title = "SCR B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_SCR_mid.pdf",width = 14, height = 8)

plotEnrichment(pathways[["BtoA_late"]], ranks) + labs(title = "SCR B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_SCR_late.pdf",width = 14, height = 8)

#gata4
ranks = -(GATA4tetvsGATA4ctr$t)
names(ranks) = GATA4tetvsGATA4ctr$gene_name
fgseGATA4 <- fgsea(pathways = pathways, 
                 stats    = ranks,
                 scoreType = 'std',
                 minSize  = 1,
                 maxSize  = 500)

pdf(file = paste0( "/home/shared_folder/plots/gsea_GATA4_grouped.pdf"), width = 20, height = 15)
plotGseaTable(pathways, stats = ranks, fgseaRes = fgseGATA4, gseaParam = 0.5)
dev.off()

plotEnrichment(pathways[["BtoA_early"]], ranks) + labs(title = "GATA4 B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_GATA4_early.pdf",width = 14, height = 8)

plotEnrichment(pathways[["BtoA_mid"]], ranks) + labs(title = "GATA4 B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_GATA4_mid.pdf",width = 14, height = 8)

plotEnrichment(pathways[["BtoA_late"]], ranks) + labs(title = "GATA4 B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_GATA4_late.pdf",width = 14, height = 8)

#ctcf
ranks = -(CTCFtetvsCTCFctr$t)
names(ranks) = CTCFtetvsCTCFctr$gene_name
fgseCTCF <- fgsea(pathways = pathways, 
                   stats    = ranks,
                   scoreType = 'std',
                   minSize  = 1,
                   maxSize  = 500)

pdf(file = paste0( "/home/shared_folder/plots/gsea_CTCF_grouped.pdf"), width = 20, height = 15)
plotGseaTable(pathways, stats = ranks, fgseaRes = fgseCTCF, gseaParam = 0.5)
dev.off()

plotEnrichment(pathways[["BtoA_early"]], ranks) + labs(title = "CTCF B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_CTCF_early.pdf",width = 14, height = 8)

plotEnrichment(pathways[["BtoA_mid"]], ranks) + labs(title = "CTCF B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_CTCF_mid.pdf",width = 14, height = 8)

plotEnrichment(pathways[["BtoA_late"]], ranks) + labs(title = "CTCF B to A")
ggsave(filename = "/home/shared_folder/plots/gsea_CTCF_late.pdf",width = 14, height = 8)

