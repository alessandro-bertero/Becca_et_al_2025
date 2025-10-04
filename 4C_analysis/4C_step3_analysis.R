#libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(limma)
library(DescTools)
library(ggsignif)
library(Hmisc)
library(pheatmap)
library(ggplotify)

theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 45, hjust = 1)))

#upload bed files
bed_4c04_wt <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/WT_4C04.bedGraph", header=FALSE)
bed_4c04_c18 <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/C18_4C04.bedGraph", header=FALSE)
bed_4c04_c57 <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/C57_4C04.bedGraph", header=FALSE)
bed_4c05_wt <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/WT_4C05.bedGraph", header=FALSE)
bed_4c05_c18 <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/C18_4C05.bedGraph", header=FALSE)
bed_4c05_c57 <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/C57_4C05.bedGraph", header=FALSE)
bed_4c06_wt <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/WT_4C06.bedGraph", header=FALSE)
bed_4c06_c18 <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/C18_4C06.bedGraph", header=FALSE)
bed_4c06_c57 <- read.delim("/home/shared_folder/BedGraph_input_files/4C_tracks/C57_4C06.bedGraph", header=FALSE)
colnames(bed_4c04_wt) = c("chr", "start", "end", "wt_04")
colnames(bed_4c04_c18) = c("chr", "start", "end", "c18_04")
colnames(bed_4c04_c57) = c("chr", "start", "end", "c57_04")
colnames(bed_4c05_wt) = c("chr", "start", "end", "wt_05")
colnames(bed_4c05_c18) = c("chr", "start", "end", "c18_05")
colnames(bed_4c05_c57) = c("chr", "start", "end", "c57_05")
colnames(bed_4c06_wt) = c("chr", "start", "end", "wt_06")
colnames(bed_4c06_c18) = c("chr", "start", "end", "c18_06")
colnames(bed_4c06_c57) = c("chr", "start", "end", "c57_06")

tab_4c04 = bed_4c04_wt %>% dplyr::left_join(bed_4c04_c18, by = c("chr", "start", "end")) %>% dplyr::left_join(bed_4c04_c57, by = c("chr", "start", "end"))
tab_4c05 = bed_4c05_wt %>% dplyr::left_join(bed_4c05_c18, by = c("chr", "start", "end")) %>% dplyr::left_join(bed_4c05_c57, by = c("chr", "start", "end"))
tab_4c06 = bed_4c06_wt %>% dplyr::left_join(bed_4c06_c18, by = c("chr", "start", "end")) %>% dplyr::left_join(bed_4c06_c57, by = c("chr", "start", "end"))

tab_div0_4c04 = tab_4c04 %>% dplyr::mutate(sum = wt_04 + c18_04 + c57_04) %>% dplyr::filter(sum != 0) %>% dplyr::select(-sum) %>% dplyr::mutate(bin = paste0(chr, "_", start))
tab_div0_4c05 = tab_4c05 %>% dplyr::mutate(sum = wt_05 + c18_05 + c57_05) %>% dplyr::filter(sum != 0) %>% dplyr::select(-sum) %>% dplyr::mutate(bin = paste0(chr, "_", start))
tab_div0_4c06 = tab_4c06 %>% dplyr::mutate(sum = wt_06 + c18_06 + c57_06) %>% dplyr::filter(sum != 0) %>% dplyr::select(-sum) %>% dplyr::mutate(bin = paste0(chr, "_", start))

tab_tot = tab_div0_4c04 %>% dplyr::full_join(tab_div0_4c05, by = c("bin", "chr", "start", "end")) %>% dplyr::full_join(tab_div0_4c06, by = c("bin", "chr", "start", "end"))
tab_tot[is.na(tab_tot)] <- 0
rownames(tab_tot) = tab_tot$bin
pos_info = tab_tot %>% dplyr::select(bin, chr, start, end)
tab_matrix = tab_tot %>% dplyr::select(-colnames(pos_info))

tab_int = tab_div0_4c04 %>% dplyr::inner_join(tab_div0_4c05, by = c("bin", "chr", "start", "end")) %>% dplyr::inner_join(tab_div0_4c06, by = c("bin", "chr", "start", "end"))
rownames(tab_int) = tab_int$bin

#normalization
norm = tab_tot %>% dplyr::filter(chr == "chr2") %>% dplyr::select(-colnames(pos_info))
norm_value = colSums(norm)
tab_norm = as.data.frame(as.matrix(tab_matrix)/norm_value * length(rownames(tab_tot)))
tab_norm$bin = rownames(tab_norm)

#TTN plot
tab_TTN = tab_int %>% dplyr::filter(chr == "chr2") %>% dplyr::filter(start > 178522000) %>% dplyr::filter(end < 178830802)
tab_TTN = tab_TTN %>% dplyr::filter(start < 178802231)
tab_TTN = tab_TTN %>% dplyr::select(c("bin", "wt_04", "c18_04", "c57_04", "wt_05", "c18_05", "c57_05", "wt_06", "c18_06", "c57_06"))

tab_TTN_melt = melt(tab_TTN)
tab_TTN_melt = tab_TTN_melt %>% dplyr::left_join(pos_info, by = "bin")
ggplot(tab_TTN_melt) + 
  geom_area(aes(x = start, y = value), alpha = 0.5, color = 1, lwd = 0.5, linetype = 1) + 
  facet_wrap(~ variable) + geom_hline(aes(yintercept = 1000))
ggsave(filename = paste0("/home/shared_folder/analysis/areas_all.pdf"))

#batch correction
batch = c("4c04", "4c04", "4c04", "4c05", "4c05", "4c05", "4c06", "4c06", "4c06")
matrix_TTN = tab_TTN %>% dplyr::select(- bin)

TTN_corr = as.data.frame(limma::removeBatchEffect(matrix_TTN, batch))
TTN_corr$bin = rownames(TTN_corr)
tab_TTN_melt_corr = melt(TTN_corr)
tab_TTN_melt_corr = tab_TTN_melt_corr %>% dplyr::left_join(pos_info, by = "bin")
ggplot(tab_TTN_melt_corr) + 
  geom_area(aes(x = start, y = value), alpha = 0.5, color = 1, lwd = 0.5, linetype = 1) + 
  facet_wrap(~ variable) + geom_hline(aes(yintercept = 1000))
ggsave(filename = paste0("/home/shared_folder/analysis/areas_all_corr.pdf"))
TTN_corr = TTN_corr %>% dplyr::mutate(position = start + 0.5)

#save tracks
WT_track = TTN_corr %>% dplyr::mutate(WT = (wt_04 + wt_05 + wt_06)/3) %>% dplyr::select(c(chr, start, end, WT))
write.table(WT_track, file = "/home/shared_folder/TTN_plot_genome_tracks/TTN_wt.bedGraph", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
C18_track = TTN_corr %>% dplyr::mutate(C18 = (c18_04 + c18_05 + c18_06)/3) %>% dplyr::select(c(chr, start, end, WT))
write.table(C18_track, file = "/home/shared_folder/TTN_plot_genome_tracks/TTN_c18.bedGraph", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
C57_track = TTN_corr %>% dplyr::mutate(C57 = (c57_04 + c57_05 + c57_06)/3) %>% dplyr::select(c(chr, start, end, WT))
write.table(C57_track, file = "/home/shared_folder/TTN_plot_genome_tracks/TTN_c57.bedGraph", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#comparison with ChIP peaks
ENCFF692RPA_chip_ESC <- read.delim("/home/shared_folder/ChIP/ENCFF692RPA_chip_ESC.bed", header=FALSE)
chip = ENCFF692RPA_chip_ESC %>% dplyr::select(c("V1", "V2", "V3"))
colnames(chip) = c("chr", "start", "end")
chip = chip %>% dplyr::filter(chr == "chr2") %>% dplyr::filter(start >= TTN_corr$start[1] & end <= TTN_corr$end[1161])
chip = chip %>% dplyr::mutate(center = start + (end - start)/2)

#area under CTCF peaks
chip = chip %>% dplyr::arrange(desc(center))
combs = list(c("wt", "c18"), c("c18", "c57"), c("c57", "wt"))
areas_between = data.frame(wt = c(), c18 = c(), c57 = c(), batch = c(), region = c(), start = c(), end = c())

#1 2
TTN_peak_1_2 = TTN_corr %>% dplyr::filter(start >= chip$center[1])
TTN_area_wt_4c04 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$wt_04, method = "trapezoid")
TTN_area_c18_4c04 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$c18_04, method = "trapezoid")
TTN_area_c57_4c04 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$c57_04, method = "trapezoid")
TTN_area_wt_4c05 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$wt_05, method = "trapezoid")
TTN_area_c18_4c05 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$c18_05, method = "trapezoid")
TTN_area_c57_4c05 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$c57_05, method = "trapezoid")
TTN_area_wt_4c06 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$wt_06, method = "trapezoid")
TTN_area_c18_4c06 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$c18_06, method = "trapezoid")
TTN_area_c57_4c06 = AUC(x = TTN_peak_1_2$start, y = TTN_peak_1_2$c57_06, method = "trapezoid")

df = data.frame(wt = c(TTN_area_wt_4c04, TTN_area_wt_4c05, TTN_area_wt_4c06),
                c18 = c(TTN_area_c18_4c04, TTN_area_c18_4c05, TTN_area_c18_4c06),
                c57 = c(TTN_area_c57_4c04, TTN_area_c57_4c05, TTN_area_c57_4c06))
df$batch = c("4C_04", "4C_05", "4C_06")
df$region = rep(c("peaks_1_2"), times = 3)
df$start = rep(as.character(chip$center[1]), times = 3)
df$end = rep(c(178802230), times = 3)
areas_between = rbind(areas_between, df)

#2 3
TTN_peak_2_3 = TTN_corr %>% dplyr::filter(start >= chip$center[2] & start <= chip$center[1])
TTN_area_wt_4c04 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$wt_04, method = "trapezoid")
TTN_area_c18_4c04 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$c18_04, method = "trapezoid")
TTN_area_c57_4c04 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$c57_04, method = "trapezoid")
TTN_area_wt_4c05 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$wt_05, method = "trapezoid")
TTN_area_c18_4c05 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$c18_05, method = "trapezoid")
TTN_area_c57_4c05 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$c57_05, method = "trapezoid")
TTN_area_wt_4c06 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$wt_06, method = "trapezoid")
TTN_area_c18_4c06 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$c18_06, method = "trapezoid")
TTN_area_c57_4c06 = AUC(x = TTN_peak_2_3$start, y = TTN_peak_2_3$c57_06, method = "trapezoid")

df = data.frame(wt = c(TTN_area_wt_4c04, TTN_area_wt_4c05, TTN_area_wt_4c06),
                c18 = c(TTN_area_c18_4c04, TTN_area_c18_4c05, TTN_area_c18_4c06),
                c57 = c(TTN_area_c57_4c04, TTN_area_c57_4c05, TTN_area_c57_4c06))
df$batch = c("4C_04", "4C_05", "4C_06")
df$region = rep(c("peaks_2_3"), times = 3)
df$start = rep(c(chip$center[2]), times = 3)
df$end = rep(c(chip$center[1]), times = 3)
areas_between = rbind(areas_between, df)

#3 4
TTN_peak_3_4 = TTN_corr %>% dplyr::filter(start >= chip$center[3] & start <= chip$center[2])
TTN_area_wt_4c04 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$wt_04, method = "trapezoid")
TTN_area_c18_4c04 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$c18_04, method = "trapezoid")
TTN_area_c57_4c04 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$c57_04, method = "trapezoid")
TTN_area_wt_4c05 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$wt_05, method = "trapezoid")
TTN_area_c18_4c05 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$c18_05, method = "trapezoid")
TTN_area_c57_4c05 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$c57_05, method = "trapezoid")
TTN_area_wt_4c06 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$wt_06, method = "trapezoid")
TTN_area_c18_4c06 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$c18_06, method = "trapezoid")
TTN_area_c57_4c06 = AUC(x = TTN_peak_3_4$start, y = TTN_peak_3_4$c57_06, method = "trapezoid")

df = data.frame(wt = c(TTN_area_wt_4c04, TTN_area_wt_4c05, TTN_area_wt_4c06),
                c18 = c(TTN_area_c18_4c04, TTN_area_c18_4c05, TTN_area_c18_4c06),
                c57 = c(TTN_area_c57_4c04, TTN_area_c57_4c05, TTN_area_c57_4c06))
df$batch = c("4C_04", "4C_05", "4C_06")
df$region = rep(c("peaks_3_4"), times = 3)
df$start = rep(c(chip$center[3]), times = 3)
df$end = rep(c(chip$center[2]), times = 3)
areas_between = rbind(areas_between, df)

#4 5
TTN_peak_4_5 = TTN_corr %>% dplyr::filter(start >= chip$center[4] & start <= chip$center[3])
TTN_area_wt_4c04 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$wt_04, method = "trapezoid")
TTN_area_c18_4c04 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$c18_04, method = "trapezoid")
TTN_area_c57_4c04 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$c57_04, method = "trapezoid")
TTN_area_wt_4c05 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$wt_05, method = "trapezoid")
TTN_area_c18_4c05 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$c18_05, method = "trapezoid")
TTN_area_c57_4c05 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$c57_05, method = "trapezoid")
TTN_area_wt_4c06 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$wt_06, method = "trapezoid")
TTN_area_c18_4c06 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$c18_06, method = "trapezoid")
TTN_area_c57_4c06 = AUC(x = TTN_peak_4_5$start, y = TTN_peak_4_5$c57_06, method = "trapezoid")

df = data.frame(wt = c(TTN_area_wt_4c04, TTN_area_wt_4c05, TTN_area_wt_4c06),
                c18 = c(TTN_area_c18_4c04, TTN_area_c18_4c05, TTN_area_c18_4c06),
                c57 = c(TTN_area_c57_4c04, TTN_area_c57_4c05, TTN_area_c57_4c06))
df$batch = c("4C_04", "4C_05", "4C_06")
df$region = rep(c("peaks_4_5"), times = 3)
df$start = rep(c(chip$center[4]), times = 3)
df$end = rep(c(chip$center[3]), times = 3)
areas_between = rbind(areas_between, df)

#5 6
TTN_peak_5_6 = TTN_corr %>% dplyr::filter(start >= chip$center[5] & start <= chip$center[4])
TTN_area_wt_4c04 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$wt_04, method = "trapezoid")
TTN_area_c18_4c04 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$c18_04, method = "trapezoid")
TTN_area_c57_4c04 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$c57_04, method = "trapezoid")
TTN_area_wt_4c05 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$wt_05, method = "trapezoid")
TTN_area_c18_4c05 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$c18_05, method = "trapezoid")
TTN_area_c57_4c05 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$c57_05, method = "trapezoid")
TTN_area_wt_4c06 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$wt_06, method = "trapezoid")
TTN_area_c18_4c06 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$c18_06, method = "trapezoid")
TTN_area_c57_4c06 = AUC(x = TTN_peak_5_6$start, y = TTN_peak_5_6$c57_06, method = "trapezoid")

df = data.frame(wt = c(TTN_area_wt_4c04, TTN_area_wt_4c05, TTN_area_wt_4c06),
                c18 = c(TTN_area_c18_4c04, TTN_area_c18_4c05, TTN_area_c18_4c06),
                c57 = c(TTN_area_c57_4c04, TTN_area_c57_4c05, TTN_area_c57_4c06))
df$batch = c("4C_04", "4C_05", "4C_06")
df$region = rep(c("peaks_5_6"), times = 3)
df$start = rep(c(chip$center[5]), times = 3)
df$end = rep(c(chip$center[4]), times = 3)
areas_between = rbind(areas_between, df)

areas_between_peaks = areas_between %>% dplyr::mutate(len = as.numeric(end) - as.numeric(start))
write.csv(areas_between_peaks, "/home/shared_folder/tables_output/areas_between_peaks.csv")

#peaks interactions
table = data.frame(peak = c(), clone = c(), t = c(), p.val = c())
chip = chip %>% dplyr::arrange(start)
chip$id = c(6, 5, 4, 3, 2)
areas_peaks = data.frame(wt = c(), c18 = c(), c57 = c(), batch = c(), region = c(), start = c(), end = c())

for(i in chip$id){
  chip_i = chip %>% dplyr::filter(id == i)
  start_i = (TTN_corr %>% dplyr::filter(position <= chip_i$center[1]))$position[length(rownames(TTN_corr %>% dplyr::filter(position <= chip_i$start[1])))]
  end_i = (TTN_corr %>% dplyr::filter(position >= chip_i$center[1]))$position[1]
  TTN_i = TTN_corr %>% dplyr::filter(position >= start_i & position <= end_i)
  TTNmelt = melt(TTN_i %>% dplyr::select(-c(start, end, position, bin))) %>% dplyr::group_by(variable) %>% dplyr::summarise(value = mean(value))
  condition = c()
  for(j in TTNmelt$variable){condition = append(condition, strsplit(j, "_")[[1]][1])}
  TTNmelt$condition = condition
  replicate = c()
  for(j in TTNmelt$variable){replicate = append(replicate, strsplit(j, "_")[[1]][2])}
  TTNmelt$replicate = replicate
  
  ggplot(TTNmelt, aes(x = condition, y = value)) + geom_boxplot(aes(fill = condition, alpha = 0.5)) + geom_signif(comparisons = combs, test = "wilcox.test", step_increase = 0.1) +
    geom_point(aes(shape = replicate, size = 2)) + scale_fill_manual(values = c( "#E69F00", "#009292", "#000000"))
  ggsave(filename = paste0("/home/shared_folder/analysis/peak1_", i, ".pdf"))
  
  TTN_wt = (TTNmelt %>% dplyr::filter(condition == "wt") %>% dplyr::arrange(replicate))$value
  TTN_c18 = (TTNmelt %>% dplyr::filter(condition == "c18") %>% dplyr::arrange(replicate))$value
  TTN_c57 = (TTNmelt %>% dplyr::filter(condition == "c57") %>% dplyr::arrange(replicate))$value
  c18 = t.test(TTN_c18 - TTN_wt, mu = 0)
  c57 = t.test(TTN_c57 - TTN_wt, mu = 0)
  table = rbind(table, data.frame(peak = i, clone = "c18", t = c18[["statistic"]][["t"]], p.value = c18[["p.value"]]))
  table = rbind(table, data.frame(peak = i, clone = "c57", t = c57[["statistic"]][["t"]], p.value = c57[["p.value"]]))
  
  TTN_peak_i = TTN_i
  TTN_area_wt_4c04 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$wt_04, method = "trapezoid")
  TTN_area_c18_4c04 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$c18_04, method = "trapezoid")
  TTN_area_c57_4c04 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$c57_04, method = "trapezoid")
  TTN_area_wt_4c05 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$wt_05, method = "trapezoid")
  TTN_area_c18_4c05 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$c18_05, method = "trapezoid")
  TTN_area_c57_4c05 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$c57_05, method = "trapezoid")
  TTN_area_wt_4c06 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$wt_06, method = "trapezoid")
  TTN_area_c18_4c06 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$c18_06, method = "trapezoid")
  TTN_area_c57_4c06 = AUC(x = TTN_peak_i$start, y = TTN_peak_i$c57_06, method = "trapezoid")
  
  df = data.frame(wt = c(TTN_area_wt_4c04, TTN_area_wt_4c05, TTN_area_wt_4c06),
                  c18 = c(TTN_area_c18_4c04, TTN_area_c18_4c05, TTN_area_c18_4c06),
                  c57 = c(TTN_area_c57_4c04, TTN_area_c57_4c05, TTN_area_c57_4c06))
  df$batch = c("4C_04", "4C_05", "4C_06")
  df$region = rep(c(paste0("peak_", i)), times = 3)
  
  melt = melt(df)
  summ = melt %>% dplyr::group_by(variable) %>% dplyr::summarise(mean = mean(value))
  melt = melt %>% dplyr::left_join(summ, by = "variable")
  ggplot(melt, aes(x = variable, y = value)) + geom_boxplot(aes(fill = variable, alpha = 0.5)) + geom_point(aes(y = value, size = 1, shape = batch)) +
    scale_fill_manual(values = c("#000000", "#009292", "#E69F00")) + geom_signif(comparisons = combs, test = "wilcox.test", step_increase = 0.1)
  ggsave(filename = paste0("/home/shared_folder/analysis/AUC_boxplot_peak_", i, ".pdf"))
  ggplot(melt, aes(x = variable, y = value)) + geom_col(aes(y = mean/3, fill = variable, alpha = 0.5)) + geom_point(aes(y = value, size = 1, shape = batch)) +
    scale_fill_manual(values = c("#000000", "#009292", "#E69F00")) + geom_signif(comparisons = combs, test = "wilcox.test", step_increase = 0.1)
  ggsave(filename = paste0("/home/shared_folder/analysis/AUC_barchart_peak_", i, ".pdf"))
  
  df$start = rep(c(start_i), times = 3)
  df$end = rep(c(end_i), times = 3)
  areas_peaks = rbind(areas_peaks, df)
  
}
write.csv(areas_peaks, "/home/shared_folder/tables_output/areas_under_ChIP_peaks.csv")

areas_between_peaks <- read_csv("C:/Users/User/Downloads/areas_between_peaks.csv")
areas_between_peaks = areas_between_peaks %>% dplyr::mutate(diff_WT_C18 = wt - c18, diff_WT_C57 = wt - c57, diff_C18_C57 = c18 - c57)
for(i in unique(areas_between_peaks$region)){
  areas_i = areas_between_peaks %>% dplyr::filter(region == i) %>% dplyr::select(c(diff_WT_C18, diff_WT_C57, diff_C18_C57))
  melt = melt(areas_i)
  summ = melt %>% dplyr::group_by(variable) %>% dplyr::summarise(mean = mean(value), t = signif(t.test(value, mu = 0)[["p.value"]], digits = 3))
  melt = melt %>% dplyr::left_join(summ, by = "variable") %>% dplyr::mutate(t = ifelse(t < 0.05, paste0(t, " *"), t))
  ggplot(melt) + geom_col(aes(x = variable, y = mean/3, fill = variable)) + geom_point(aes(x = variable, y = value)) +
    scale_fill_manual(values = c("#009292", "#E69F00", "grey")) + geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(limits = c(-3300000, 3300000)) + xlab("") + ylab("AUC difference") +
    geom_label(aes(x = variable, y = mean, label = t), alpha = 0.5)
  ggsave(paste0("D:/Docker/4C/AUC_analysis/plots/diff_AUC_", i, ".pdf"))
}

areas_up_downstream = areas_between_peaks %>% dplyr::mutate(region = ifelse(region %in% c("peaks_1_2", "peaks_2_3", "peaks_3_4"), "before_peak4", "after_peak4")) %>%
  dplyr::group_by(region, batch) %>% dplyr::summarise(wt = sum(wt), c18 = sum(c18), c57 = sum(c57)) %>% dplyr::mutate(diff_WT_C18 = wt - c18, diff_WT_C57 = wt - c57, diff_C18_C57 = c18 - c57)
for(i in unique(areas_up_downstream$region)){
  areas_i = areas_up_downstream %>% dplyr::filter(region == i) %>% dplyr::select(c(diff_WT_C18, diff_WT_C57, diff_C18_C57))
  melt = melt(areas_i)
  summ = melt %>% dplyr::group_by(variable) %>% dplyr::summarise(mean = mean(value), t = signif(t.test(value, mu = 0)[["p.value"]], digits = 3))
  melt = melt %>% dplyr::left_join(summ, by = "variable") %>% dplyr::mutate(t = ifelse(t < 0.05, paste0(t, " *"), t))
  ggplot(melt) + geom_col(aes(x = variable, y = mean/3, fill = variable)) + geom_point(aes(x = variable, y = value)) +
    scale_fill_manual(values = c("#009292", "#E69F00", "grey")) + geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(limits = c(-4900000, 4900000)) + xlab("") + ylab("AUC difference") +
    geom_label(aes(x = variable, y = mean, label = t), alpha = 0.5)
  ggsave(paste0("D:/Docker/4C/AUC_analysis/plots/diff_AUC_", i, ".pdf"))
}

areas_peaks = areas_peaks %>% dplyr::mutate(diff_WT_C18 = wt - c18, diff_WT_C57 = wt - c57, diff_C18_C57 = c18 - c57)
for(i in unique(areas_peaks$region)){
  areas_i = areas_peaks %>% dplyr::filter(region == i) %>% dplyr::select(c(diff_WT_C18, diff_WT_C57, diff_C18_C57))
  melt = melt(areas_i)
  summ = melt %>% dplyr::group_by(variable) %>% dplyr::summarise(mean = mean(value), t = signif(t.test(value, mu = 0)[["p.value"]], digits = 3))
  melt = melt %>% dplyr::left_join(summ, by = "variable") %>% dplyr::mutate(t = ifelse(t < 0.05, paste0(t, " *"), t))
  ggplot(melt) + geom_col(aes(x = variable, y = mean/3, fill = variable)) + geom_point(aes(x = variable, y = value)) +
    scale_fill_manual(values = c("#009292", "#E69F00", "grey")) + geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("") + ylab("AUC difference") + geom_label(aes(x = variable, y = mean, label = t), alpha = 0.5)
  ggsave(paste0("D:/Docker/4C/AUC_analysis/plots/diff_AUC_", i, ".pdf"))
}

##virtual 4C correlation
#binning
bins = v4C %>% dplyr::select(c(chr, start, end)) %>% dplyr::mutate(bin = paste0(chr, "_", start))
rownames(tab_int) = tab_int$bin
matrix = tab_int %>% dplyr::select(wt_04, c18_04, c57_04, wt_05, c18_05, c57_05, wt_06, c18_06, c57_06)
batch = c("4c04", "4c04", "4c04", "4c05", "4c05", "4c05", "4c06", "4c06", "4c06")
matrix = as.data.frame(limma::removeBatchEffect(matrix, batch))
bins_info = tab_int %>% dplyr::select(chr, start, end, bin)
matrix$bin = rownames(matrix)
matrix = matrix %>% dplyr::left_join(bins_info, by = "bin")
bin_v4C = c()
for(i in 1:length(rownames(matrix))){
  v4C_i = bins %>% dplyr::filter(chr == matrix$chr[i] & start <= matrix$start[i] & end >= matrix$end[i])
  bin_v4C = append(bin_v4C, v4C_i$bin[1])
}
matrix$bin_v4C = bin_v4C

#virtual 4C correlation
undiff <- read.delim("/home/shared_folder/BedGraph_input_files/virtual4C/undiff_v4C.bedGraph", header=FALSE, skip = 1)
colnames(undiff) = c("chr", "start", "end", "undiff")
meso <- read.delim("/home/shared_folder/BedGraph_input_files/virtual4C/meso_v4C.bedGraph", header=FALSE, skip = 1)
colnames(meso) = c("chr", "start", "end", "meso")
prog <- read.delim("/home/shared_folder/BedGraph_input_files/virtual4C/prog_v4C.bedGraph", header=FALSE, skip = 1)
colnames(prog) = c("chr", "start", "end", "prog")
cardio <- read.delim("/home/shared_folder/BedGraph_input_files/virtual4C/cardio_v4C.bedGraph", header=FALSE, skip = 1)
colnames(cardio) = c("chr", "start", "end", "cardio")
v4C = undiff %>% dplyr::left_join(meso, by = c("chr", "start", "end")) %>% dplyr::left_join(prog, by = c("chr", "start", "end")) %>% dplyr::left_join(cardio, by = c("chr", "start", "end"))

new_matrix = matrix %>% dplyr::group_by(bin_v4C) %>% dplyr::summarise(wt_04 = mean(wt_04), wt_05 = mean(wt_05), wt_06 = mean(wt_06),
                                                                      c18_04 = mean(c18_04), c18_05 = mean(c18_05), c18_06 = mean(c18_06),
                                                                      c57_04 = mean(c57_04), c57_05 = mean(c57_05), c57_06 = mean(c57_06))
pos = c()
chr = c()
for(i in 1:length(rownames(new_matrix))){
  pos = append(pos, as.integer(strsplit(new_matrix$bin_v4C[i], "_")[[1]][2]) + 20000)
  chr = append(chr, strsplit(new_matrix$bin_v4C[i], "_")[[1]][1])
}
new_matrix$pos = pos
new_matrix$chr = chr

v4C = v4C %>% dplyr::mutate(bin_v4C = paste0(chr, "_", start))
v4C = v4C %>% dplyr::filter(bin_v4C %in% new_matrix$bin_v4C) %>% dplyr::select(-c(chr, start, end))

new_matrix = new_matrix %>% dplyr::filter(bin_v4C %in% v4C$bin_v4C)
new_matrix = new_matrix %>% dplyr::left_join(v4C, by = "bin_v4C")

new_matrix = new_matrix %>% dplyr::mutate(wt_04 = wt_04 / sum(wt_04)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(wt_05 = wt_05 / sum(wt_05)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(wt_06 = wt_06 / sum(wt_06)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(c18_04 = c18_04 / sum(c18_04)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(c18_05 = c18_05 / sum(c18_05)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(c18_06 = c18_06 / sum(c18_06)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(c57_04 = c57_04 / sum(c57_04)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(c57_05 = c57_05 / sum(c57_05)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(c57_06 = c57_06 / sum(c57_06)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(undiff = undiff / sum(undiff)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(meso = meso / sum(meso)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(prog = prog / sum(prog)* length(rownames(new_matrix)))
new_matrix = new_matrix %>% dplyr::mutate(cardio = cardio / sum(cardio)* length(rownames(new_matrix)))

new_matrix = new_matrix %>% dplyr::mutate(wt_04 = log2(ifelse(wt_04 < 1, 0, wt_04) + 1))
new_matrix = new_matrix %>% dplyr::mutate(wt_05 = log2(ifelse(wt_05 < 1, 0, wt_05) + 1))
new_matrix = new_matrix %>% dplyr::mutate(wt_06 = log2(ifelse(wt_06 < 1, 0, wt_06) + 1))
new_matrix = new_matrix %>% dplyr::mutate(c18_04 = log2(ifelse(c18_04 < 1, 0, c18_04) + 1))
new_matrix = new_matrix %>% dplyr::mutate(c18_05 = log2(ifelse(c18_05 < 1, 0, c18_05) + 1))
new_matrix = new_matrix %>% dplyr::mutate(c18_06 = log2(ifelse(c18_06 < 1, 0, c18_06) + 1))
new_matrix = new_matrix %>% dplyr::mutate(c57_04 = log2(ifelse(c57_04 < 1, 0, c57_04) + 1))
new_matrix = new_matrix %>% dplyr::mutate(c57_05 = log2(ifelse(c57_05 < 1, 0, c57_05) + 1))
new_matrix = new_matrix %>% dplyr::mutate(c57_06 = log2(ifelse(c57_06 < 1, 0, c57_06) + 1))
new_matrix = new_matrix %>% dplyr::mutate(undiff = log2(ifelse(undiff < 1, 0, undiff) + 1))
new_matrix = new_matrix %>% dplyr::mutate(meso = log2(ifelse(meso < 1, 0, meso) + 1))
new_matrix = new_matrix %>% dplyr::mutate(prog = log2(ifelse(prog < 1, 0, prog) + 1))
new_matrix = new_matrix %>% dplyr::mutate(cardio = log2(ifelse(cardio < 1, 0, cardio) + 1))

tab_coeff = data.frame(coeff = c(0:20))
rownames(tab_coeff) = tab_coeff$coeff
tab_coeff = tab_coeff %>% dplyr::select(-coeff)

wt04_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$wt_04 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(wt04_fun[["coefficients"]])$'wt04_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$wt_04 = coeff
wt05_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$wt_05 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(wt05_fun[["coefficients"]])$'wt05_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$wt_05 = coeff
wt06_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$wt_06 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(wt06_fun[["coefficients"]])$'wt06_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$wt_06 = coeff
c1804_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$c18_04 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(c1804_fun[["coefficients"]])$'c1804_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$c18_04 = coeff
c1805_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$c18_05 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(c1805_fun[["coefficients"]])$'c1805_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$c18_05 = coeff
c1806_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$c18_06 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(c1806_fun[["coefficients"]])$'c1806_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$c18_06 = coeff
c5704_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$c57_04 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(c5704_fun[["coefficients"]])$'c5704_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$c57_04 = coeff
c5705_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$c57_05 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(c5705_fun[["coefficients"]])$'c5705_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$c57_05 = coeff
c5706_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$c57_06 ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(c5706_fun[["coefficients"]])$'c5706_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$c57_06 = coeff

undiff_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$undiff ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(undiff_fun[["coefficients"]])$'undiff_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$undiff = coeff
meso_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$meso ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(meso_fun[["coefficients"]])$'meso_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$meso = coeff
prog_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$prog ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(prog_fun[["coefficients"]])$'prog_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$prog = coeff
cardio_fun = lm((new_matrix %>% dplyr::filter(chr == "chr2"))$cardio ~ poly((new_matrix %>% dplyr::filter(chr == "chr2"))$pos, degree = 20))
coeff = as.data.frame(cardio_fun[["coefficients"]])$'cardio_fun[["coefficients"]]'
coeff[is.na(coeff)] <- 0
tab_coeff$cardio = coeff

tab_coeff_d = tab_coeff
tab_coeff_d[1,] = tab_coeff[1,]*0
tab_coeff_d[2,] = tab_coeff[2,]*1
tab_coeff_d[3,] = tab_coeff[3,]*2
tab_coeff_d[4,] = tab_coeff[4,]*4
tab_coeff_d[5,] = tab_coeff[5,]*4
tab_coeff_d[6,] = tab_coeff[6,]*5
tab_coeff_d[7,] = tab_coeff[7,]*6
tab_coeff_d[8,] = tab_coeff[8,]*7
tab_coeff_d[9,] = tab_coeff[9,]*8
tab_coeff_d[10,] = tab_coeff[10,]*9
tab_coeff_d[11,] = tab_coeff[11,]*10
tab_coeff_d[12,] = tab_coeff[12,]*11
tab_coeff_d[13,] = tab_coeff[13,]*12
tab_coeff_d[14,] = tab_coeff[14,]*13
tab_coeff_d[15,] = tab_coeff[15,]*14
tab_coeff_d[16,] = tab_coeff[16,]*15
tab_coeff_d[17,] = tab_coeff[17,]*16
tab_coeff_d[18,] = tab_coeff[18,]*17
tab_coeff_d[19,] = tab_coeff[19,]*18
tab_coeff_d[20,] = tab_coeff[20,]*19
tab_coeff_d[21,] = tab_coeff[21,]*20

corr = rcorr(as.matrix(tab_coeff_d))
pearson = as.data.frame(corr[["r"]])
pearson = pearson[1:9, 10:13]
pearson$clone = rownames(pearson)
melt = melt(pearson)
condition = c()
for(i in 1:length(rownames(melt))){
  condition = append(condition, strsplit(melt$clone[i], "_")[[1]][1])
}
melt$condition = condition
ggplot(melt, aes(x = variable, y = value)) + geom_boxplot(aes(fill = condition)) + facet_wrap(~condition)
ggsave(filename = paste0("/home/shared_folder/analysis/correlation_rues2.pdf"))

