library(dplyr)

#dataset1
dir.create("/home/shared_folder/dataset1")
system2("wget", "-P /home/shared_folder/dataset1 https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/929/E-MTAB-3929/Files/counts.txt")
system2("wget", "-P /home/shared_folder https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/929/E-MTAB-3929/Files/E-MTAB-3929.sdrf.txt")
counts1 <- read.delim("/home/shared_folder/dataset1/counts.txt")
E.MTAB.3929.sdrf <- read.delim("/home/shared_folder/dataset1/E-MTAB-3929.sdrf.txt")
meta1 = E.MTAB.3929.sdrf %>% dplyr::select(c("Source.Name", "Characteristics.developmental.stage.", "Characteristics.inferred.lineage.")) %>%
  dplyr::filter(Characteristics.inferred.lineage. %in% c("not applicable"))
rownames(counts1) = counts1$X
counts1 = counts1 %>% dplyr::select(meta1$Source.Name)
colnames(meta1) = c("cell", "timepoint", "phenotype")
meta1 = meta1 %>% dplyr::mutate(stage = "pre_implantation")

#dataset2
dir.create("/home/shared_folder/dataset2")
system2("git", "clone https://github.com/ScialdoneLab/human-gastrula-shiny.git /home/shared_folder/dataset2")
raw_reads <- readRDS("/home/shared_folder/dataset2/human-gastrula-shiny/raw_reads.rds")
umap <- readRDS("/home/shared_folder/dataset2/human-gastrula-shiny/umap.rds")
raw_reads$X = as.numeric(rownames(raw_reads)) - 1
rownames(raw_reads) = raw_reads$X
counts2 = as.data.frame(t(raw_reads %>% dplyr::select(-c("X"))))
meta2 = umap  %>% dplyr::mutate(timepoint = "day16") %>% dplyr::select(c("X", "timepoint", "cluster_id")) %>% dplyr::mutate(X = paste0("cell_", X)) %>%
  dplyr::filter(cluster_id %in% c("Primitive Streak", "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm"))
colnames(meta2) = c("cell", "timepoint", "phenotype")
meta2 = meta2 %>% dplyr::mutate(stage = "gastrula")
colnames(counts2) = paste0("cell_", colnames(counts2))
counts2 = counts2 %>% dplyr::select(meta2$cell)

#datset3
dir.create("/home/shared_folder/dataset3")
system2("wget", "-P /home/shared_folder/dataset3 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106118/suppl/GSE106118%5FUMI%5Fcount%5Fmerge%2Etxt%2Egz")
system2("wget", "-P /home/shared_folder/dataset3 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106118/suppl/GSE106118%5Fbarcode%5Finformation%2Etxt%2Egz")
GSE106118_UMI_count_merge <- read.delim("/home/shared_folder/dataset3/GSE106118_UMI_count_merge.txt")
GSE106118_barcode_information <- read.delim("/home/shared_folder/dataset3/GSE106118_barcode_information.txt")
get_timepoint = function(str){
  val = strsplit(str, "_")[[1]][1]
  val = strsplit(val, "")[[1]]
  if(length(val) == 4){time = paste0(val[3], val[4])}
  if(length(val) == 5){time = paste0(val[3], val[4], val[5])}
  return(time)
}
get_phenotype = function(str){
  val = strsplit(str, "_")[[1]][3]
  phenotype = strsplit(val, "[.]")[[1]][1]
  return(phenotype)
}
meta3 = GSE106118_barcode_information
counts3 = GSE106118_UMI_count_merge
meta3$phenotype = as.character(lapply(meta3$cell, get_phenotype))
meta3$timepoint = as.character(lapply(meta3$cell, get_timepoint))
meta3 = meta3 %>% dplyr::filter(phenotype %in% c("LA", "RA", "LV", "RV")) %>% dplyr::select(c("cell", "timepoint", "phenotype")) %>% dplyr::filter(cell %in% colnames(counts4))
meta3 = meta3 %>% dplyr::mutate(stage = "embrionic_heart")
rownames(counts3) = counts4$Gene
counts3 = counts3 %>% dplyr::select(-Gene) %>% dplyr::select(meta4$cell)

#intersection
common_genes = intersect(intersect(intersect(rownames(counts1), rownames(counts2)), rownames(counts3)), rownames(counts4))
counts1$genes = rownames(counts1)
counts1 = counts1 %>% dplyr::filter(genes %in% common_genes) %>% dplyr::select(-genes)
counts2$genes = rownames(counts2)
counts2 = counts2 %>% dplyr::filter(genes %in% common_genes) %>% dplyr::select(-genes)
counts3$genes = rownames(counts3)
counts3 = counts3 %>% dplyr::filter(genes %in% common_genes) %>% dplyr::select(-genes)

#save
save(meta1, counts1, meta2, counts2, meta3, counts3, file = "/home/shared_folder/preprocessed.RData")
