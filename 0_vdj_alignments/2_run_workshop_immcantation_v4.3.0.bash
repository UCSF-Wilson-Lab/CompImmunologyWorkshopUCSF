#!/usr/bin/bash

# Goal: Generate results for 10X MS BCR & TCR data
#  - there are two BCR and TCR samples
#  - Immcantation v4.3.0 supports the analysis of these data (v4.5.0 had issues with BLASTP)

# Immcantation Documentation:
# https://immcantation.readthedocs.io/en/stable/

# Install container on your computer:
# docker pull immcantation/suite:4.3.0

# 1. Make sure data exists within container (use 'exit' command to exit container):
# docker run --name WORKSHOP -it immcantation/suite:4.3.0 bash 

# Arguments
# - Input formatted sample name
DATA_DIR=/data/input_immcantation
READS_BCR=$DATA_DIR\/input_bcr_filtered_contig.fasta
READS_TCR=$DATA_DIR\/input_tcr_filtered_contig.fasta
ANNOTATIONS_BCR=$DATA_DIR\/input_bcr_filtered_contig_annotations.csv
ANNOTATIONS_TCR=$DATA_DIR\/input_tcr_filtered_contig_annotations.csv
OUT_DIR_BCR=/data/results_immcantation/BCR_CSFPB
OUT_DIR_TCR=/data/results_immcantation/TCR_CSFPB
MODEL=aa
NPROC=4

# Run commands from within docker

# >>> BCR
changeo-10x -s $READS_BCR -a $ANNOTATIONS_BCR -x 0.15 -n BCR_CSFPB -o $OUT_DIR_BCR -p $NPROC -t ig -m $MODEL -f airr

# >>> TCR
changeo-10x -s $READS_TCR -a $ANNOTATIONS_TCR -x 0.00 -n TCR_CSFPB -o $OUT_DIR_TCR -p $NPROC -t tr -m $MODEL -f airr

# >>> process TCR Alpha chain results


cd $OUT_DIR_TCR
DefineClones.py -d TCR_CSFPB_light_productive-T.tsv --act first --model $MODEL --dist 0.00
CreateGermlines.py -d TCR_CSFPB_light_productive-T_clone-pass.tsv -o TCR_CSFPB_light_germ-pass.tsv -g dmask --cloned -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_TRA*

# >>> process BCR light chain results
cd $OUT_DIR_BCR
DefineClones.py -d BCR_CSFPB_light_productive-T.tsv --act first --model $MODEL --dist 0.15
CreateGermlines.py -d BCR_CSFPB_light_productive-T_clone-pass.tsv -o BCR_CSFPB_light_germ-pass.tsv -g dmask --cloned -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IG*


