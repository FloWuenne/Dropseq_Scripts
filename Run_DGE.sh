#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=12,pmem=4gb
#PBS -r n
#PBS -A nfk-724-aa
#PBS -N Dropseq_DGE


############################
scratch_dir=""
file_prefix=""
barcode_file=$scratch_dir/""
num_barcodes_to_discover=

############################
dropseq_tools="/gs/project/nfk-724-aa/Programs/Drop-seq_tools-1.12"

sample_BAM=$scratch_dir/$file_prefix".aligned.tagged.cleaned.bam"
sample_dge_discovered=$scratch_dir/$file_prefix".aligned.tagged.cleaned.discovered_BCs.DGE.txt"
sample_dge_bc_file=$scratch_dir/$file_prefix".aligned.tagged.cleaned.selected_STAMPS.DGE.txt"
discovered_sample_summary=$scratch_dir/$file_prefix".aligned.tagged.discovered_BCs.DGE.summary.txt"
selected_sample_summary=$scratch_dir/$file_prefix".aligned.tagged.selected_BCs.DGE.summary.txt"


## Perform DGE counting for species mixing experiment

## Use CORE Barcodes for STAMP selection
$dropseq_tools/DigitalExpression I=$sample_BAM O=$sample_dge_discovered SUMMARY=$sample_summary NUM_CORE_BARCODES=$num_barcodes_to_discover

## Use a list of barcodes for STAMP selection
$dropseq_tools/DigitalExpression I=$sample_BAM O=$sample_dge SUMMARY=$sample_summary CELL_BC_FILE=$barcode_file
