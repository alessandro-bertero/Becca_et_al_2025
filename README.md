# Becca_et_al_2025
Code repository for the bioinformatic analysis of Becca et al 2025 "Opposing CTCF and GATA4 activities set the pace of chromatin topology remodeling during cardiomyogenesis"

In order to reproduce our pipelines follow the steps below from downloading the supporting data files, pulling our docker containers and running the pipelines.
The pipelines in this repository are relative to the publication figures 1-6 and supplementary

### HiC_data_integration
HiC data from different developmental stages of cardiac differentiation in various cell lines were analyzed to investigate the compartment changes occurring on the TTN locus. The datasets were dowloaded from public repositories (4DN Data Portal and GEO). HiC data were reallinged (for RUES2 data only) and ICE-normalized with HiC-pro, and preprocessed with HiCExplorer to change the data format. Preprocessed files were used for differential compartment analysis with dcHiC, and for 4C analysis from the TTN promoter with the HiContacts R package. HiC tracks were visualized with pyGenomeTracks.
**Relative to Fig 1 and S1**

### scRNAseq_data_integration
Single cell data from different stages of human embryos development were integrated to analyze the changes in gene expression of CTCF and GATA4 across cardiac differentiation. With the first script (step1_preprocessing.R) three datasets are downloaded as scRNAseq counts with the respective meta data, and the they are preprocessed in order to be combined in a single object containing only the cells related to the mesodermal-cardiac lineage. With the second script (step2_analysis.R) the resulting sigle cell object is log2 normalized, corrected to account for the batch of the original dataset with Seurat, and analyzed with Monocle3 to identify the main clusters and to define the pseudotime trajectory related to cardiac differentiation. Gene expression patterns of the genes of interest and of some developmental markers were finally analyzed across the selected trajectory.
**Relative to Fig 1 and S1**

### 4C_analysis
4C data were analysed to compare the iteraction of the TTN promoter with the TTN gene body upon KO of 2 CTCF binding sites. Reads were alligned using o the hg38 human genome using the pipeline provided by Krijger P. H. L. et al., Methods (2020) with the 4C_step1_alignment.R, and the resulting bigWig files were converted to BedGraph for easy data integration on R with the 4C_step2_fileconversion.sh script. The 4C_step3_analysis.R script was used to normalize the tracks on the library size and to correct for the replicate batch. Counts realtive to the TTN gene were filtered and converted to BedGraph files for visualization with pyGenomeTracks. The resulting track was integrated with ChIP data of hESCs cells (ENCFF332TNJ) to obtain the area under the curbe (AUC) of the 4C profile under each peak and betwen each pair of consecutive peaks. The fifferences in the AUC of each regions for eaach pair of conditions were evaluated with a t.test. The 4C trackon chromosome2 was then binnd at 40kb, fitted into a polinomial function and the same fitting was performed for the RUES2 virtual 4C track (from the same viewpoint, retrieved from HiC data) of the different stages of cardiac development. The results were compared by computing the Pearson correlation of the second derivate of each pair of fucntions.
**Relative to Fig 3 and S3**

### bulk_RNAseq_2d_cardiomyocytes_with_KD (bulkRNAseq_1)
Bulk RNAseq data of 2D cardiomyocytes were analyzed to investigate the gene expression changes of B to A genes upon GATA4 and CTCF KD. With the first script (H001BS1_step1_preprocessing.R) raw reads are trimmed using trim-galore, and aligned on the GRCh38.112 genome using STAR. With the second script (H001BS1_step2_DGEanalysis) the counts are combined, filtered with DESeq2, CPM and TPM normalized, batch corrected and used for differential gene expression analysis with limma comapring each TET condition with its isogenic control. GSEA analysis was then performed on differentially expressed genes to analyze the main GO terms and Kegg pathways related to cardiac development. The thirs script (H001BS1_step3_BtoAgenes.R) H9 compartment and gene expression data were combined to identify the genes with the strongest B to A switch from day 0 to day 80 of cardiac development, and with the highest gene expression peak after day 15. These genes were used as a new gene set for GSEA analysis on the ranking resulting from the DGE analysis performed with the previous step.
**Relative to Fig 2 and S2**

### bulk_RNAseq_organoids_with_KD (bulkRNAseq_2)
Bulk RNAseq data of 3D right ventricle (RV) and atrial (AT) cardiac organoids were analyze to investigate the anticorrelation in gene expression changes upon GATA4 and CTCF KD. With the first script (SBe001AS1_step1_preprocessing.sh) raw reads are trimmed using trim-galore, and aligned on the GRCh38.111 genome using STAR. With the second script (SBe001AS1_step2_DGEanalysis.R) the counts are combined, filtered with DESeq2, CPM and TPM normalized, batch corrected and used for differential gene expression analysis with limma comapring each TET condition with its isogenic control. With the third script (SBe001AS1_step3_heatmaps_GSEA_RRHO.R) data were separated for the 2 cardiac chambers to perform a PCA analysis, and log2 TPM counts were analyzed for B to A genes in the 2 chambers. The ranking obtained from the DGE analysis performed in the previous step was used for GSEA analysis to investigate GO terms related to heart process fucntions, and for RRHO analysis to compute the global correlation of the gene expression perturbations in each chamber.
**Relative to Fig 4 and S4**

### scRNAseq_cardiac_organoids (scRNAseq_1)
Single cell RNAseq data of cardiac organoids were analysed to compare the specific phenotypic commitment of 3 differentiation protocols inducing the right ventricle (RV), left ventricle (LV) and atrium (AT) fates respectively. Raw reads were alligned (multi.sh) and combined (aggr.sh) with cellranger8, and scRNAseq counts were filtered with Seurat with the AB001_step1_Seurat_Filtering.R script to remove cells with few counts or expressed genes and with high expression of mitochondrial genes. With the AB001_step2_Monocle_SeparateChambersAnalysis.R, cells differentiated with each protol were separated in a different single cell dataset, and each dataset was analysed with Monocle3 to perform log2 normalization, dimensional reduction, clustering, and pseudotime analysis. Each cluster was further analysed with GO analysis on the top 10 expressed markers, with DGE analysis compared to the others and with gene modules detection to identify the cell type related to each cluster. Where needed new subclusters were created according to cell-type specific markers. With the last script (AB001_step3_Monocle_JointChambersAnalysis.R) the annotated cells were combined in a new unique single cell dataset. The dataset was again processed with dimensional reduction, and the main markers were visualized on the resulting UMAP. Finally the clusters classification was validated with Seurat by transfering the cell types label from the dataset retrieved form the GSE263193 accession number.
**Relative to Fig S4**

### scRNAseq_cardiac_organoids_with_KD (scRNAseq_2)
Single cell RNAseq data of left ventricle (LV) cardiac organoids were analysed to compare the cell types distribution upon GATA4 and CTCF KD. Raw reads were alligned and combined (aggr.sh) with cellranger8, and scRNAseq counts were filtered with the AB002_SeuratPrepocessing.R script to remove cells with few counts or expressed genes and with high expression of mitochondrial genes. The dataset was further processed with the Seurat pipeline in order to retrieve the metadata related to the cell cylce phase. With the AB002_step2_MonocleAnalysisAll.R, the single cell dataset was futher filtered to remove some artifacts, and it was preprocessed, and log2 normalized with Monocle3. With the AB002_step3_MonocleAnalysisSeparate.R the dataset was separated by timepoint (day4 and day7), and each resulting Monocle object was processed to obtain dimensional reduction, clustring, and pseudotime trajectories. Cluster idetity was defined upon combination of the GO of the top expressed markers, and GSEA (clusterProfiler) on differentially expressed genes of each cluster compared to the others. The resulting inferred cell types were validated with a gene module analysis and a visualization of the gene expression opatterns of validated markers. The enrichment of cells in each cluster in the TET compared to the control condition of each KD was evaluated with both a Fisher's exact test and a KNN graph analysis using MiloR. Differences in the cell cycle were evaluated with a Wilcoxon test based on Seurat's annotation.
**Relative to Fig 5, Fig 6, and Fig S5**

# Download the supporting datasets from the Zenodo repository:
Use wget or download directly on [Zenodo](https://doi.org/10.5281/zenodo.17256128).

    #scratch HiC data integration
    wget https://zenodo.org/records/17256128/files/scratch_HiC_H9.zip?download=1
    wget https://zenodo.org/records/17256128/files/scratch_HiC_RUES2.zip?download=1
    wget https://zenodo.org/records/17256128/files/scratch_HiC_WTC11.zip?download=1
    
    #scratch 4C 
    wget https://zenodo.org/records/17256128/files/scratch_4C.zip?download=1

    #scratch_bulk
    wget https://zenodo.org/records/17256128/files/scratch_bulk.zip?download=1

    #scratch_scRNAseq_1
    wget https://zenodo.org/records/17256128/files/scratch_scRNAseq1.zip?download=1
    
    #scratch_scRNAseq_2
    wget https://zenodo.org/records/17256128/files/scratch_scRNAseq2.zip?download=1

# Run the scripts within a Docker container:

The analysis can be performed within a Docker environment, which uses RStudio Server (2023.12.1 Build 402 Ocean storm release). The Docker image can be downloaded from our repository hedgelab/rstudio-hedgelab (rstudio password is rstudio when omitted).

### HiC_data_integration
Follow these steps to pull and run the Docker image for the HiC data analysis from the terminal:

    dockered
    docker pull hedgelab/hic_pipeline:image15
    
For H9 HiC data integration unzip the scratch_HiC_H9, add the H9_Rscript.R script to the unzipped folder, and in the parent folder proceed with the following steps:

    docker run -d -v ./scratch_HiC_H9:/home/shared_folder --name=NAME_CONTAINER hedgelab/hic_pipeline:image15
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/H9_Rscript.R
    docker exec -it NOME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/genomeTrack/track_TTN --region chr2:179350000-179700000 --outFileName /home/shared_folder/genomeTrack/TTN_H9_final.pdf

For WTC11 HiC data integration unzip the scratch_HiC_WTC11 and add the WTC11_Rscript.R script to the unzipped folder. Create a subfolder ./mcool and download the .mcool files from the 4DN Data Portal with the following Experiment Set Codes: 4DNFIUMP8ZZ6, 4DNFIK5SQFY2, 4DNFIUD4ECSX, 4DNFI4N2S9QW, 4DNFI4MXUCUV, 4DNFI1LYV4VR, 4DNFIX5V9GBA, 4DNFIJPZ4ASS, 4DNFIM3UJ1XI, 4DNFIRHM7URR, 4DNFIYBBPZ1V, 4DNFIQ57GW6C. Move to the parent folder and proceed with the following steps:

    docker run -d -v ./scratch_HiC_WTC11:/home/shared_folder --name=NAME_CONTAINER hedgelab/hic_pipeline:image15
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/WTC11_Rscript.R
    docker exec -it NOME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/genome_track/track_atrial --region chr2:178500000-178850000 --outFileName /home/shared_folder/genome_track/WTC11_TTN_atrial.pdf
    docker exec -it NOME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/genome_track/track_ventricular --region chr2:178500000-178850000 --outFileName /home/shared_folder/genome_track/WTC11_TTN_ventricular.pdf

For RUES2 HiC data integration unzip the scratch_HiC_RUES2 and add the RUES2_HiC_preprocessing.sh and RUES2_Rscript.R scripts to the unzipped folder. Move to the parent folder and proceed with the following steps:

    docker run -d -v ./scratch_HiC_RUES2:/home/shared_folder --name=NAME_CONTAINER hedgelab/hic_pipeline:image15
    docker exec -it NAME_CONTAINER /home/shared_folder/RUES2_HiC_preprocessing.sh
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/RUES2_Rscript.R
    docker exec -it NOME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/genomeTrack/track_TTN_RUES2 --region chr2:178500000-178850000 --outFileName /home/shared_folder/genomeTrack/TTN_RUES2_final.pdf

Refer to our GitHub [Hi-c pipeline](https://github.com/sara-bianchi/HiC_pipeline)

### 4C analysis
Follow these steps to pull and run the Docker image for the 4C data analysis from the terminal:

    dockered
    docker pull hedgelab/4c:image1

Unzip the scratch_4C.zip folder, and download the fastq files from ENA (ACCESSION NUMBER) in a fastq subfolder. Move the 4C scripts into the scratch folder and run the following commands to perform the analysis.

    docker run -d -v ./scratch_4C:/home/shared_folder --privileged=true --name=NAME_CONTAINER hedgelab/4c:image1
    docker exec -it NOME_CONTAINER wget -P /home/shared_folder https://zenodo.org/records/17256128/files/bowtie2_genome.zip?download=1
    docker exec -it NOME_CONTAINER Rscript 4C_step1_alignment.R --vpFile=./home/shared_fodler/VPinfo.txt --confFile /home/shared_folder/conf_plots.yml --fqFolder=./home/shared_folder/fastq --outFolder=./home/shared_folder/output --cores 8 --plot --wig
    docker exec -it NOME_CONTAINER /home/shared_folder/4C_step2_fileconversion.sh
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/4C_step3_analysis.R
    docker exec -it NAME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/TTN_plot_genomeTracks/TTN_track --region chr2:178500000-178850000 --outFileName /home/shared_folder/TTN_plot_genomeTracks/4C_clones.pdf


### scRNAseq analysis
Follow these steps to pull and run the Docker image for the scRNAseq data analysis from the terminal:

    dockered
    docker pull hedgelab/cellranger8:latest
    docker pull hedgelab/rstudio-hedgelab:media

For single_cell_data_integration move the single_cell_data_integration scripts in a ./scratch_single_cell folder and run the following commands:

    docker run -d -itv ./scratch_single_cell folder:/home/shared_folder --privileged=true --name=NAME_CONTAINER hedgelab/rstudio-hedgelab:media
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/step1_preprocessing.R
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/step2_analysis.R
    
For scRNAseq_1 analysis unzip the scratch_scRNAseq_1.zip folder, and download the fastq files from ArrayExpress (ACCESSION NUMBER) in a ./fastq subfolder. Move the scRNAseq_1 scripts into the scratch folder and run the following commands to perform the analysis.

    docker run -itv ./scratch_scRNAseq_1:/home/shared_folder hedgelab/cellranger8 /home/shared_folder/multi.sh
    docker run -itv ./scratch_scRNAseq_1:/home/shared_folder hedgelab/cellranger8 /home/shared_folder/aggr.sh
    docker run -d -itv ./scratch_scRNAseq_1:/home/shared_folder --privileged=true -p 8787:8787 \\
    -e PASSWORD=<your_password> -e USER=rstudio --name=NAME_CONTAINER hedgelab/rstudio-hedgelab:media
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/AB001_step1_Seurat_Filtering.R
    docker exec -it NAME_CONTAINER rstudio-server start

For scRNAseq_2 analysis unzip the scratch_scRNAseq_2.zip folder, and download the fastq files from ArrayExpress (ACCESSION NUMBER) in a ./fastq subfolder. Move the scRNAseq_1 scripts into the scratch folder and run the following commands to perform the analysis.

    docker run -itv ./scratch_scRNAseq_2:/home/shared_folder hedgelab/cellranger8 /cellranger-8.0.1/bin/cellranger multi --id=AB002_multi --output-dir=/home/shared_folder/Results --csv=/home/shared_folder/config_multi.csv
    docker run -itv ./scratch_scRNAseq_2:/home/shared_folder hedgelab/cellranger8 /home/shared_folder/aggr.sh
    docker run -d -itv ./scratch_scRNAseq_2:/home/shared_folder --privileged=true -p 8787:8787 \\
    -e PASSWORD=<your_password> -e USER=rstudio --name=NAME_CONTAINER hedgelab/rstudio-hedgelab:media
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/AB002_step1_SeuratPrerocessing.R
    docker exec -it NAME_CONTAINER rstudio-server start

Replace <your_password> with your desired password, if it is omitted the password will be "rstudio". This command maps port 8787 on your host machine to port 8787 in the container, allowing you to access RStudio via your browser at http://localhost:8787. The USER=rstudio part ensures you'll log in as the rstudio user, and the PASSWORD variable sets the password you'll use to log in. Then use Rstudio through a browser to execute the AB001_Monocle_SeparateChambersAnalysis.R and AB001_step3_JointChambersAnalysis.R for scRNAseq_1, and AB002_MoncleAnalysisAll and AB002_MonocleAnalysisSeparate.R for scRNAseq_2.

### bulkRNAseq analysis
Follow these steps to pull and run the Docker image for the scRNAseq data analysis from the terminal:

    dockered
    docker pull hedgelab/bulk_image:image3

For bulkRNAseq_1 analysis unzip the scratch_bulkRNAseq_1.zip folder, and download the fastq files from ArrayExpress (ACCESSION NUMBER) in a ./fastq subfolder. Move the bulkRNAseq_1 scripts into the scratch folder and run the following commands to perform the analysis.

    docker run -d -v /the/folder/you/want/to/share:/home/shared_folder --privileged=true --name=NAME_CONTAINER hedgelab/bulk_image:image3
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/H00BS1_step1_preprocessing.R
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/H00BS1_step2_DGEanalysis.R
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/H00BS1_step3_BtoAgenes.R

For bulkRNAseq_1 analysis unzip the scratch_bulkRNAseq_2.zip folder, and download the fastq files from ArrayExpress (ACCESSION NUMBER) in a ./fastq subfolder. Move the bulkRNAseq_2 scripts into the scratch folder and run the following commands to perform the analysis.

    docker run -d -v /the/folder/you/want/to/share:/home/shared_folder --privileged=true --name=NAME_CONTAINER hedgelab/bulk_image:image3
    docker exec -it NAME_CONTAINER /home/shared_folder/SBe001AS1_step1_preprocessing.sh
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/SBe001AS1_step2_DGEanalysis.R
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/SBe001AS1_step3_heatmaps_GSEA_RRHO.R

Refer to our GitHub [bulk RNAseq pipeline](https://github.com/sara-bianchi/Bulk_RNA_seq_pipeline)

