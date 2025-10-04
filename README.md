# Becca_et_al_2025
Code repository for the bioinformatic analysis of Becca et al 2025 "Opposing CTCF and GATA4 activities set the pace of chromatin topology remodeling during cardiomyogenesis"

In order to reproduce our pipelines follow the steps below from downloading the supporting data files, pulling our docker containers and running the pipelines.
The pipelines in this repository are relative to the publication figures 1-6 and supplementary

### HiC_data_integration
XXXadd a breef description here of the steps.
**Relative to Fig 1 and S1**

### scRNAseq_data_integration
XXXadd a breef description here of the steps.
**Relative to Fig 1 and S1**

### bulk_RNAseq_2d_cardiomyocytes (bulkRNAseq_1)
XXXadd a breef description here of the steps.
**Relative to Fig 2 and S2**

### 4C_analysis
XXXadd a breef description here of the steps.
**Relative to Fig 3 and S3**

### bulk_RNAseq_organoids_with_KD (bulkRNAseq_2)
XXXadd a breef description here of the steps.
**Relative to Fig 4 and S4**

### scRNAseq_cardiac_organoids (scRNAseq_1)
XXXadd a breef description here of the steps.
**Relative to Fig S4**

### scRNAseq_cardiac_organoids_with_KD (scRNAseq_2)
XXXadd a breef description here of the steps.
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

# HiC_data_integration
Follow these steps to pull and run the Docker image for the HiC data analysis from the terminal:

    dockered
    docker pull hedgelab/hic_pipeline:image15
    
For H9 HiC data integration unzip the scratch_HiC_H9, add the H9_Rscript.R script to the unzipped folder, and in the parent folder proceed with the following steps:

    docker run -d -itv ./scratch_HiC_H9:/home/shared_folder --name=NAME_CONTAINER hedgelab/hic_pipeline:image15
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/H9_Rscript.R
    docker exec -it NOME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/genomeTrack/track_TTN --region chr2:179350000-179700000 --outFileName /home/shared_folder/genomeTrack/TTN_H9_final.pdf

For WTC11 HiC data integration unzip the scratch_HiC_WTC11 and add the WTC11_Rscript.R script to the unzipped folder. Create a subfolder ./mcool and download the .mcool files from the 4DN Data Portal with the following Experiment Set Codes: 4DNFIUMP8ZZ6, 4DNFIK5SQFY2, 4DNFIUD4ECSX, 4DNFI4N2S9QW, 4DNFI4MXUCUV, 4DNFI1LYV4VR, 4DNFIX5V9GBA, 4DNFIJPZ4ASS, 4DNFIM3UJ1XI, 4DNFIRHM7URR, 4DNFIYBBPZ1V, 4DNFIQ57GW6C. Move to the parent folder and proceed with the following steps:

    docker run -d -itv ./scratch_HiC_WTC11:/home/shared_folder --name=NAME_CONTAINER hedgelab/hic_pipeline:image15
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/WTC11_Rscript.R
    docker exec -it NOME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/genome_track/track_atrial --region chr2:178500000-178850000 --outFileName /home/shared_folder/genome_track/WTC11_TTN_atrial.pdf
    docker exec -it NOME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/genome_track/track_ventricular --region chr2:178500000-178850000 --outFileName /home/shared_folder/genome_track/WTC11_TTN_ventricular.pdf

For RUES2 HiC data integration unzip the scratch_HiC_RUES2 and add the RUES2_HiC_preprocessing.sh and RUES2_Rscript.R scripts to the unzipped folder. Move to the parent folder and proceed with the following steps:

    docker run -d -itv ./scratch_HiC_RUES2:/home/shared_folder --name=NAME_CONTAINER hedgelab/hic_pipeline:image15
    docker exec -it NAME_CONTAINER ./home/shared_folder/RUES2_HiC_preprocessing.sh
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/RUES2_Rscript.R
    docker exec -it NOME_CONTAINER pyGenomeTracks --tracks /home/shared_folder/genomeTrack/track_TTN_RUES2 --region chr2:178500000-178850000 --outFileName /home/shared_folder/genomeTrack/TTN_RUES2_final.pdf

Refer to our GitHub [Hi-c pipeline](https://github.com/sara-bianchi/HiC_pipeline)

Follow these steps to pull and run the Docker image for the 4C data analysis from the terminal:

    dockered
    docker pull hedgelab/4c:image1

    docker run -d -itv /the/folder/you/want/to/share:/scratch --name=NAME_CONTAINER hedgelab/hic_pipeline:image15
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/Rscript.R

Refer to our GitHub [Hi-c pipeline](https://github.com/sara-bianchi/HiC_pipeline)

Follow these steps to pull and run the Docker image for the scRNAseq data analysis from the terminal:


    dockered
    docker pull hedgelab/rstudio-hedgelab:media

    docker run -d -itv /the/folder/you/want/to/share:/home/shared_folder --privileged=true -p 8787:8787 \\
    -e PASSWORD=<your_password> -e USER=rstudio --name=NAME_CONTAINER hedgelab/rstudio-hedgelab:media

    docker exec -it NAME_CONTAINER rstudio-server start


Replace <your_password> with your desired password, if it is omitted the password will be "rstudio". This command maps port 8787 on your host machine to port 8787 in the container, allowing you to access RStudio via your browser at http://localhost:8787. The USER=rstudio part ensures you'll log in as the rstudio user, and the PASSWORD variable sets the password you'll use to log in. Then use Rstudio through a browser.
Share the folder containing the scripts and the "scratch" downloaded as suggested below.

NB when analyzing the bulk RNAseq data pull image: hedgelab/bulk_image:image3 with 

    dockered
    docker pull hedgelab/bulk_image:image3

    docker run -d -v /the/folder/you/want/to/share:/home/shared_folder --privileged=true -p 8787:8787 \\
    -e PASSWORD=<your_password> -e USER=rstudio --name=NAME_CONTAINER hedgelab/bulk_image:image3

    docker exec -it NAME_CONTAINER rstudio-server start

Refer to our GitHub [bulk RNAseq pipeline](https://github.com/sara-bianchi/Bulk_RNA_seq_pipeline)

