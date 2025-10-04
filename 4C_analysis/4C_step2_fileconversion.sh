#!/bin/bash

#docker conda
#conda create --name bigwig ucsc-bigwigtobedgraph -c bioconda
#conda init bash #then close shell and restart
#conda activate bigwig
#mkdir /home/shared_folder/analysis/BedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C04_WT_WIN21.bw /home/shared_folder/BedGraph/WT_4C04.bedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C04_C18_WIN21.bw /home/shared_folder/BedGraph/C18_4C04.bedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C04_C57_WIN21.bw /home/shared_folder/BedGraph/C57_4C04.bedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C05_WT_WIN21.bw /home/shared_folder/BedGraph/WT_4C05.bedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C05_C18_WIN21.bw /home/shared_folder/BedGraph/C18_4C05.bedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C05_C57_WIN21.bw /home/shared_folder/BedGraph/C57_4C05.bedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C06_WT_WIN21.bw /home/shared_folder/BedGraph/WT_4C06.bedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C06_C18_WIN21.bw /home/shared_folder/BedGraph/C18_4C06.bedGraph
bigWigToBedGraph /home/shared_folder/output/BWIG/4C06_C57_WIN21.bw /home/shared_folder/BedGraph/C57_4C06.bedGraph
