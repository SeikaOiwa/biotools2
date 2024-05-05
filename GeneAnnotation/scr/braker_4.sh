#!/bin/bash

# path設定
export PATH=/home/seika_oiwa_3590/notebooks/GeneAnnotation/BRAKER/scripts/:$PATH
export PATH=/home/seika_oiwa_3590/notebooks/GeneAnnotation/GeneMark-ETP/tools/:$PATH
export GENEMARK_PATH=/home/seika_oiwa_3590/notebooks/GeneAnnotation/GeneMark-ETP/bin
export AUGUSTUS_CONFIG_PATH=/home/seika_oiwa_3590/notebooks/GeneAnnotation/Augustus/config
export AUGUSTUS_BIN_PATH=/home/seika_oiwa_3590/notebooks/GeneAnnotation/Augustus/bin
export AUGUSTUS_SCRIPTS_PATH=/home/seika_oiwa_3590/notebooks/GeneAnnotation/Augustus/scripts
export BAMTOOLS_PATH=/home/seika_oiwa_3590/notebooks/GeneAnnotation/bamtools/bin
export TSEBRA_PATH=/home/seika_oiwa_3590/notebooks/GeneAnnotation/TSEBRA/bin

# 引数情報
# $1: strain name
# $2: path to genome
# $3: path to bam
# $4: path to protein_db

# braker予測
braker.pl --species=$1 --genome=$2 --bam=$3 --prot_seq=$4 --threads 8 --gff3

