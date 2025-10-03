# Becca_et_al_2025
Code repository for the bioinformatic analysis of Becca et al 2025 "Opposing CTCF and GATA4 activities set the pace of chromatin topology remodeling during cardiomyogenesis"

In order to reproduce our pipelines follow the steps below from downloading the supporting data files, pulling our docker containers and running the pipelines.
The pipelines in this repository are relative to the publication figures XXXXX

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


# Step 1:
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
