#------------------------------------------------------------------------------
#BULK RNA SEQ
#------------------------------------------------------------------------------

#packages
library(tidyr)
library(stringr)
library(grid)
library(viridis)
library(ggthemes)
library(dplyr)
library(tidyverse)
library(edgeR)
library(compareGroups)
library(dbscan)
library(reshape2)
library(scran)
library(fgsea)
library(RNAseqQC)
library(DESeq2)
library(ensembldb)
library(tibble)
library(ggplotify)
library(limma)
library(Glimma)
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggrepel)
library(ggsignif)
library(igraph)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(WGCNA)
library(readr)

#set up environment
theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 45, hjust = 1)))
options(stringsAsFactors = FALSE)
update_geom_defaults("point", aes(size = 4))
set.seed(1234)

pal = c("#3283FE", "#FA0087", "#009E73", "#FBE426", "#56B4E9", "#FEAF16", "#DEA0FD", "#1CBE4F", "#F6222E", "#1CFFCE", "#325A9B", "#AA0DFE","#D55E00", "#2ED9FF", "#f0E442", "#1C8356", "#0072B2", "#CC79A7")

current_time = format(Sys.time(), "%y%m%d")
dir.create(paste0("/home/shared_folder/Output_", current_time))


#upload settings

print("uploading settings file")

settings = as.data.frame(read_excel("/home/shared_folder/settings.xlsx", col_names = FALSE))
rownames(settings) = settings$...1
settings = settings %>% dplyr::select(-c("...1", "...2"))
settings = as.data.frame(t(settings))

if(length(unique(na.omit(settings$trimming))) != 1 | as.logical(settings$trimming[1]) == "NA"){
  print("invalid argument for trimming"); stop()
} else{trimming = as.logical(settings$trimming[1])}

if(trimming == TRUE){
  if(length(unique(na.omit(settings$fastqc))) != 1 | as.logical(settings$fastqc[1]) == "NA"){
    print("invalid argument for fastqc"); stop()
  } else{fastqc = as.logical(settings$fastqc[1]);
  if(fastqc == TRUE){fqc = " --fastqc"}
  else if(fastqc == FALSE){fqc = ""}}
  if(length(unique(na.omit(settings$mode))) != 1 | !as.character(settings$mode[1]) %in% c("single", "paired")){
    print("invalid argument for mode"); stop()
  } else{mode = as.character(settings$mode[1])};
  bases = strsplit(settings$adapter[1], "")
  invalid_bases = c()
  for(i in bases[[1]]){if(!i %in% c("A", "T", "G", "C")){invalid_bases = append(invalid_bases, i)}}
  if(length(unique(na.omit(settings$adapter))) != 1 | length(invalid_bases) != 0){
    print("invalid argument for adapter"); stop()
  } else{adapter = as.character(settings$adapter[1])};
  if(length(unique(na.omit(settings$clip_5_1))) != 1 | as.integer(settings$clip_5_1[1]) == "NA"){
    print("invalid argument for clip_5_1"); stop()
  } else{clip_5_1 = as.integer(settings$clip_5_1[1])};
  if(length(unique(na.omit(settings$clip_3_1))) != 1 | as.integer(settings$clip_3_1[1]) == "NA"){
    print("invalid argument for clip_3_1"); stop()
  } else{clip_3_1 = as.integer(settings$clip_3_1[1])};
  if(mode == "paired"){
    bases = strsplit(settings$adapter2[1], "")
    invalid_bases = c()
    for(i in bases[[1]]){if(!i %in% c("A", "T", "G", "C")){invalid_bases = append(invalid_bases, i)}}
    if(length(unique(na.omit(settings$adapter2))) != 1 | length(invalid_bases) != 0){
      print("invalid argument for adapter2"); stop()
    } else{adapter2 = as.character(settings$adapter2[1])};
    if(length(unique(na.omit(settings$clip_5_2))) != 1 | as.integer(settings$clip_5_2[1]) == "NA"){
      print("invalid argument for clip_5_2"); stop()
    } else{clip_5_2 = as.integer(settings$clip_5_2[1])};
    if(length(unique(na.omit(settings$clip_3_2))) != 1 | as.integer(settings$clip_3_2[1]) == "NA"){
      print("invalid argument for clip_3_2"); stop()
    } else{clip_3_2 = as.integer(settings$clip_3_2[1])};
  }
}

if(length(unique(na.omit(settings$indexing))) != 1 | as.logical(settings$indexing[1]) == "NA"){
  print("invalid argument for indexing"); stop()
} else{indexing = as.logical(settings$indexing[1])}

if(length(unique(na.omit(settings$alignment))) != 1 | as.logical(settings$alignment[1]) == "NA"){
  print("invalid argument for alignment"); stop()
} else{alignment = as.logical(settings$alignment[1])}

if(indexing == TRUE){
  if(length(unique(na.omit(settings$gtf))) != 1 | file.exists(settings$gtf[1]) == FALSE){
    print("invalid argument for gtf"); stop()
  } else{gtf = as.character(settings$gtf[1])};
  if(length(unique(na.omit(settings$fasta))) != 1 | file.exists(settings$fasta[1]) == FALSE){
    print("invalid argument for fasta"); stop()
  } else{fasta = as.character(settings$fasta[1])};
} else{if(alignment == TRUE){
  if(length(unique(na.omit(settings$genome_dir))) != 1 | dir.exists(settings$genome_dir[1]) == FALSE){
    print("invalid argument for genom_dir"); stop()
  } else{genome_dir = as.character(settings$genome_dir[1])}};
  if(length(unique(na.omit(settings$gene_length))) != 1 | file.exists(settings$gene_length[1]) == FALSE){
    print("invalid argument for gene_length"); stop()
  } else{gene_length = as.character(settings$gene_length[1])};
  if(length(unique(na.omit(settings$gene_info))) != 1 | file.exists(settings$gene_info[1]) == FALSE){
    print("invalid argument for gene_info"); stop()
  } else{geneInfo = as.character(settings$gene_info[1]); geneInfo = read.table(geneInfo, quote = "\"", comment.char = "", skip = 1)}
}

if(alignment == TRUE){if(length(unique(na.omit(settings$mode))) != 1 | !as.character(settings$mode[1]) %in% c("single", "paired")){
  print("invalid argument for mode"); stop()
} else{mode = as.character(settings$mode[1])}}

#table

print("uploading table file")

table = read_excel("/home/shared_folder/table.xlsx")
if(alignment == TRUE){
  if(mode == "single"){
    if(!"fastq" %in% colnames(table)){
      print("fastq column is missing"); stop()
    } else{
      for(i in table$fastq){
        if(file.exists(i) == FALSE){
          print("incorrect path in fastq column"); stop()
        }
      }
    }
  } else if(mode == "paired"){
    if(!"fastq1" %in% colnames(table) | !"fastq2" %in% colnames(table)){
      print("fastq1 or fastq2 columns are missing"); stop()
    } else{
      for(i in table$fastq1){if(file.exists(i) == FALSE){
        print("incorrect path in fastq1 column"); stop()}
      }
      for(i in table$fastq2){if(file.exists(i) == FALSE){
        print("incorrect path in fastq2 column"); stop()}
      }
    }
  }
} else if(alignment == FALSE){
  if(!"counts" %in% colnames(table)){
    print("counts column is missing"); stop()
  } else{
    for(i in table$counts){
      if(file.exists(i) == FALSE){
        print("incorrect path in counts column"); stop()
      }
    }
  }
}

#trimming
if(trimming == TRUE){
  print("starting trimming")
  
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/trimmed"));
  out_dir = paste0("/home/shared_folder/Output_", current_time, "/trimmed");
  if(mode == "single"){
    if(clip_3_1 == 0){
      clipping_R1_3 = ""
    } else {
      clipping_R1_3 = paste0(" --three_prime_clip_R1 ", as.character(clip_3_1))
    }
    if(clip_5_1 == 0){
      clipping_R1_5 = ""
    } else {
      clipping_R1_5 = paste0(" --clip_R1 ", as.character(clip_5_1))
    }
    
    for(i in 1:length(rownames(table))){fastq = table$fastq[i]; sample_ID = table$sample_ID[i]
      system2("/usr/bin/trim_galore", paste0("--adapter ", as.character(adapter),  clipping_R1_5, clipping_R1_3, " --gzip", fqc, " --basename ", sample_ID, " --output_dir ", as.character(out_dir), " ", as.character(fastq)))
    }
  } else if(mode == "paired"){
    if(clip_3_1 == 0){
      clipping_R1_3 = ""
    } else {
      clipping_R1_3 = paste0(" --three_prime_clip_R1 ", as.character(clip_3_1))
    }
    if(clip_5_1 == 0){
      clipping_R1_5 = ""
    } else {
      clipping_R1_5 = paste0(" --clip_R1 ", as.character(clip_5_1))
    }
    if(clip_3_2 == 0){
      clipping_R2_3 = ""
    } else {
      clipping_R2_3 = paste0(" --three_prime_clip_R2 ", as.character(clip_3_2))
    }
    if(clip_5_2 == 0){
      clipping_R2_5 = ""
    } else {
      clipping_R2_5 = paste0(" --clip_R2 ", as.character(clip_5_2))
    }
    for(i in 1:length(rownames(table))){fastq1 = table$fastq1[i]; fastq2 = table$fastq2[i]; sample_ID = table$sample_ID[i]
      system2("/usr/bin/trim_galore", paste0("--adapter ", as.character(adapter), " --adapter2 ", as.character(adapter2), clipping_R1_5, clipping_R1_3, clipping_R2_5, clipping_R2_3," --gzip --paired", fqc, " --basename ", sample_ID, " --output_dir ", as.character(out_dir), " ", as.character(fastq1), " ", as.character(fastq2)))
    }
  }
  print("finishing trimming")
}

#indexing
if(indexing == TRUE){
  print("starting indexing")
  
  system2("mkdir", "/home/genome")
  system2("chmod", "777 /home/genome")
  system2("/STAR", paste0("--runThreadN 4 --runMode genomeGenerate --limitGenomeGenerateRAM=150000000000 --genomeDir /home/genome --genomeFastaFiles ", fasta, " --sjdbGTFfile ", gtf));
  system2("rm", paste0("-r /home/shared_folder/Output_", current_time, "/genome"));
  system2("mv", paste0("/home/genome /home/shared_folder/Output_", current_time));
  genome_dir = paste0("/home/shared_folder/Output_", current_time, "/genome");
  system2("gtftools", paste0("-l /home/shared_folder/Output_", current_time, "/gene_length.txt ", gtf));
  gene_length = paste0("/home/shared_folder/Output_", current_time, "/gene_length.txt");
  geneInfo = read.table(paste0(genome_dir, "/geneInfo.tab"), quote = "\"", comment.char = "", skip = 1)
  
  print("finishing indexing")
} else{
  system2("wget", genome_dir)
}

colnames(geneInfo) = c("gene_ID", "gene_name", "gene_type")

#alignment
if(alignment == TRUE){
  print("starting alignment")
  
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/counts"));
  if(mode == "single"){
    for(i in 1:length(rownames(table))){
      sample_ID = table$sample_ID[i];
      counts_dir = paste0("/home/shared_folder/Output_", current_time, "/counts/", sample_ID);
      dir.create(counts_dir);
      system2("chmod", paste0("777 ", counts_dir));
      if(trimming == FALSE){fastq = table$fastq[i]} else if(trimming == TRUE){
        fastq = paste0(out_dir, "/", sample_ID, "_trimmed.fq.gz")};
      system2("/STAR", paste0("--runThreadN 4 --genomeDir ", genome_dir, " --readFilesIn ", fastq, " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix ", counts_dir, "/", sample_ID))}
  } else if(mode == "paired"){
    for(i in 1:length(rownames(table))){sample_ID = table$sample_ID[i];
    counts_dir = paste0("/home/shared_folder/Output_", current_time, "/counts/", sample_ID);
    if(trimming == FALSE){fastq1 = table$fastq1[i]; fastq2 = table$fastq2[i]} else if(trimming == TRUE){
      fastq1 = paste0(out_dir, "/", sample_ID, "_val_1.fq.gz");
      fastq2 = paste0(out_dir, "/", sample_ID, "_val_2.fq.gz")};
    system2("/STAR", paste0("--runThreadN 4 --genomeDir ", genome_dir, " --readFilesIn ", fastq1, " ", fastq2, " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix ", counts_dir, "/", sample_ID))}
  }
  print("finishing alignment")
}


