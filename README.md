# Becca_et_al_2025
Code repository for the bioinformatic analysis of Becca et al 2025 "Opposing CTCF and GATA4 activities set the pace of chromatin topology remodeling during cardiomyogenesis"

In order to reproduce our pipelines follow the steps below from downloading the supporting data files, pulling our docker containers and running the pipelines.
The pipelines in this repository are relative to the publication figures 3-5 and supplementary

4C_analysis
XXXadd a breef description here of the steps.

Relative to Fig 3 and S3

bulk_RNAseq
XXXadd a breef description here of the steps.
Relative to Fig 4

scRNAseq_cardiac_organoids (scRNAseq_1)
XXXadd a breef description here of the steps.

Relative to Fig S4

scRNAseq_cardiac_organoids_with_KD (scRNAseq_2)
XXXadd a breef description here of the steps.

Relative to Fig 5 and Fig S5

# Run the scripts within a docker container:

The analysis can be performed within a Docker environment, which uses RStudio Server (2023.12.1 Build 402 ìOcean stormî release). The Docker image can be downloaded from our repository hedgelab/rstudio-hedgelab (rstudio password is rstudio when omitted).

Follow these steps to pull and run the Docker image for the scRNAseq data abìnalysis from the terminal:


    dockered
    docker pull hedgelab/rstudio-hedgelab:iPS2seq

    docker run -d -itv /the/folder/you/want/to/share:/scratch --privileged=true -p 8787:8787 \\
    -e PASSWORD=<your_password> -e USER=rstudio --name=NAME_CONTAINER hedgelab/rstudio-hedgelab:iPS2seq

    docker exec -idt NAME_CONTAINER rstudio-server start


Replace <your_password> with your desired password, if it is omitted the password will be "rstudio". This command maps port 8787 on your host machine to port 8787 in the container, allowing you to access RStudio via your browser at http://localhost:8787. The USER=rstudio part ensures you'll log in as the rstudio user, and the PASSWORD variable sets the password you'll use to log in. Then use Rstudio through a browser.
Share the folder containing the scripts and the "scratch" downloaded as suggested below.

NB when analyzing the bulk RNAseq data pull image: hedgelab/bulk_image:image3 with 

    dockered
    docker pull hedgelab/bulk_image:image3

Follow these steps to pull and run the Docker image for the 4C data analysis from the terminal:


    dockered
    docker pull hedgelab/hic_pipeline:image12

    docker run -d -itv /the/folder/you/want/to/share:/scratch --name=NAME_CONTAINER hedgelab/hic_pipeline:image12
    docker exec -it NAME_CONTAINER Rscript /home/shared_folder/Rscript.R

Refer to our github [Hi-c pipeline](https://github.com/sara-bianchi/HiC_pipeline)


# Download the supporting datasets from Zenodo repository:
Use wget or download directly on [Zenodo](https://doi.org/10.5281/zenodo.17256128).

    #scratch 4C 
    wget https://zenodo.org/records/17256128/files/scratch_4C.zip?download=1

    #scratch_bulk
    wget https://zenodo.org/records/17256128/files/scratch_bulk.zip?download=1

    #scratch_scRNAseq_1
    wget https://zenodo.org/records/17256128/files/scratch_scRNAseq1.zip?download=1
    
    #scratch_scRNAseq_2
    wget https://zenodo.org/records/17256128/files/scratch_scRNAseq2.zip?download=1

# Step pull the docker containers:

Alternatively, the repository can be download manually and loaded as follows: 
