#!/bin/bash
#PBS -l walltime=15:00:00
#PBS -l nodes=1:ppn=12,pmem=4gb
#PBS -r n
#PBS -A 
#PBS -N Dropseq_pipeline

## Load modules on guillimin
## Main MUGQIC module
module use /cvmfs/soft.mugqic/CentOS6/modulefiles/mugqic
module load star
module load mugqic/fastqc

######################

## Read in variables from Dropseq_config.txt
## Change the name of the config file here if other than Dropseq_config.txt
source ./Dropseq_config.txt

#######################

#### File prefixes and output filenames

## Define the file prefix that will be the name of each subsequent file created in this pipeline 

## Tag Cell Barcodes
unaligned_bam=$file_prefix".unaligned.bam"
unaligned_bam_tagged_CB=$file_prefix".unaligned_tagged_Cell.bam"
unaligned_bam_tagged_CB_summary=$qc_stat_dir/$file_prefix".unaligned_tagged_Cellular.bam_summary.txt"

## Tag Molecular Barcodes
unaligned_bam_tagged_CB_MB=$file_prefix".unaligned_tagged_CellMolecular.bam"
unaligned_bam_tagged_CB_MB_summary=$qc_stat_dir/$file_prefix".unaligned_tagged_Molecular.bam_summary.txt"

## FilterBAM
unaligned_bam_tagged_CB_MB_FB=$file_prefix".unaligned_tagged_CellMolecular.filteredBAM.bam"

## Trim 5' primer sequence
unaligned_bam_tagged_CB_MB_trimmed_smart=$file_prefix".unaligned_tagged_CellMolecular_trimmed_smart.bam"
unaligned_bam_tagged_CB_MB_trimmed_smart_summary=$qc_stat_dir/$file_prefix".adapter_trimming_report.txt"

## Trim 3' polyA sequence
unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered=$file_prefix".unaligned_tagged_CellMolecular_trimmed_smart_poly_filtered.bam"
unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered_summary=$qc_stat_dir/$file_prefix".polyA_trimming_report.txt"

## BAM to FASTQ
unaligned_fastq_mc_tagged_polyA_filtered=$file_prefix".unaligned_tagged_CellMolecular_trimmed_smart_poly_filtered.fastq"

## Aligned Bam sorted by queryname
aligned_bam_sorted=$file_prefix".aligned.sorted.bam"

## Merged BAM
merged_bam=$file_prefix".aligned.merged.bam"

## Exon tagged BAM
gene_exon_tagged_bam=$file_prefix".aligned.merged.star_gene_exon_tagged.bam"

## Cleaned Bam (From Bead synthesis errors)
cleaned_bam=$file_prefix".aligned.tagged.cleaned.bam"
synthesis_stats=$qc_stat_dir/$file_prefix".synthesis.stats.txt"
synthesis_summary=$qc_stat_dir/$file_prefix".synthesis.summary.txt"

## BamTagHistogram
bamtaghistogram=$qc_stat_dir/$file_prefix".readcount.txt.gz"

## QC Metrics
UMI_by_gene=$qc_stat_dir/"UMI_by_gene_dist.tab"
single_cell_metrics=$qc_stat_dir/"Single_Cell_Quality_metrics.tab"

## Picard command 
## Will create a tmp folder in the outfile dir
java_cmd="java -jar -Xmx4g -Djava.io.tmpdir=$outfile_dir/tmp"

tbwrse_ex="TagBamWithReadSequenceExtended"
filter_bam="FilterBAM"
trim_start="TrimStartingSequence"
polyatrimmer="PolyATrimmer"
tagreadwithgeneexon="TagReadWithGeneExon"
detectbeadsynthesiserrors="DetectBeadSynthesisErrors"
gatherMolbarDistbyGene="GatherMolecularBarcodeDistributionByGene"
singlecellmetrics="SingleCellRnaSeqMetricsCollector"


################ Run the Dropseq pipeline on the data specified above

## If piping is activated, pipe as much as possible and delete intermedite files

if (($piping == 1))
then
	## Run Fastqc on original fastq files
        "fastqc" $fastq1 $fastq2 --outdir $outfile_dir

        ## Fastq --> Sam
        # Using picards FastqToSam command
        $java_cmd $PICARD_HOME FastqToSam F1=$fastq1 F2=$fastq2 SM=$file_prefix O=$outfile_dir/$unaligned_bam SO=queryname

        $dropseq_tools/$tbwrse_ex INPUT=$outfile_dir/$unaligned_bam OUTPUT=/dev/stdout SUMMARY=$outfile_dir/$unaligned_bam_tagged_CB_summary BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 COMPRESSION_LEVEL=0 |\
        $dropseq_tools/$tbwrse_ex INPUT=/dev/stdin OUTPUT=/dev/stdout SUMMARY=$outfile_dir/$unaligned_bam_tagged_CB_MB_summary BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 COMPRESSION_LEVEL=0  | \
	$dropseq_tools/$filter_bam TAG_REJECT=XQ INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \ 
        $dropseq_tools/$trim_start INPUT=/dev/stdin OUTPUT=/dev/stdout OUTPUT_SUMMARY=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_summary  SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5 COMPRESSION_LEVEL=0 | \
        $dropseq_tools/$polyatrimmer INPUT=/dev/stdin OUTPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered OUTPUT_SUMMARY=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered_summary MISMATCHES=0 NUM_BASES=6

	## Sam --> Fastq
        ## Using picards SamToFastq function
        $java_cmd $PICARD_HOME SamToFastq INPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered FASTQ=$outfile_dir/$unaligned_fastq_mc_tagged_polyA_filtered

        #### STAR alignment
        STAR --runMode alignReads --genomeDir $STAR_index --sjdbGTFfile $gtf_file --readFilesIn $outfile_dir/$unaligned_fastq_mc_tagged_polyA_filtered --outFileNamePrefix $outfile_dir/$file_prefix --runThreadN 12

        ## SortSam
        $java_cmd $PICARD_HOME SortSam I=$outfile_dir/$file_prefix"Aligned.out.sam" O=$outfile_dir/$aligned_bam_sorted SO=queryname

	# MergeBam
        $java_cmd $PICARD_HOME MergeBamAlignment REFERENCE_SEQUENCE=$reference_fasta UNMAPPED_BAM=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered ALIGNED_BAM=$outfile_dir/$aligned_bam_sorted OUTPUT=/dev/stdout INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false COMPRESSION_LEVEL=0| \
        dropseq_tools/$tagreadwithgeneexon I=/dev/stdin O=$outfile_dir/$gene_exon_tagged_bam ANNOTATIONS_FILE=$refflat TAG=GE

        # Detect Bead Synthesis Errors
        # DetectBeadSynthesisErrors
        $dropseq_tools/$detectbeadsynthesiserrors I=$outfile_dir/$gene_exon_tagged_bam O=$outfile_dir/$cleaned_bam OUTPUT_STATS=$outfile_dir/$synthesis_stats SUMMARY=$outfile_dir/$synthesis_summary NUM_BARCODES=$NUM_of_barcodes PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

	# Cleanup	
	# Delete files that are not longer needed
	rm $outfile_dir/$unaligned_bam
	rm $outfile_dir/$unaligned_fastq_mc_tagged_polyA_filtered
	rm $outfile_dir/$file_prefix"Aligned.out.sam"
	rm $outfile_dir/$aligned_bam_sorted
	rm $outfile_dir/$unaligned_fastq_mc_tagged_polyA_filtered
	rm $outfile_dir/$gene_exon_tagged_bam
	
	# Calculate reads per Cell Barcode or knee plot
        $dropseq_tools/BAMTagHistogram I=$outfile_dir/$cleaned_bam O=$outfile_dir/$bamtaghistogram TAG=XC

        # Run a couple of diagnostic scripts to make check what the UMI distribution etc looks like
        # Gather Molecular Barcode distributions by Gene
        $dropseq_tools/$gatherMolbarDistbyGene I=$outfile_dir/$cleaned_bam O=$outfile_dir/$UMI_by_gene NUM_CORE_BARCODES=$NUM_of_barcodes

        # Calculate single cell metrics
        $dropseq_tools/$singlecellmetrics I=$outfile_dir/$cleaned_bam O=$outfile_dir/$single_cell_metrics ANNOTATIONS_FILE=$gtf_file NUM_CORE_BARCODES=$NUM_of_barcodes


else
## Non piping 

	## Run Fastqc on original fastq files
	"fastqc" $fastq1 $fastq2 --outdir $outfile_dir

	## Fastq --> Sam
	# Using picards FastqToSam command
	$java_cmd $PICARD_HOME FastqToSam F1=$fastq1 F2=$fastq2 SM=$file_prefix O=$outfile_dir/$unaligned_bam SO=queryname

	# Dropseq-tools programs
	# TagBamWithReadSequenceExtended
	# 1) Extract Cell Barcode
	$dropseq_tools/$tbwrse_ex INPUT=$outfile_dir/$unaligned_bam OUTPUT=$outfile_dir/$unaligned_bam_tagged_CB SUMMARY=$outfile_dir/$unaligned_bam_tagged_CB_summary BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1
	
	# 2) Extract Molecular Barcode
	$dropseq_tools/$tbwrse_ex INPUT=$outfile_dir/$unaligned_bam_tagged_CB OUTPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB SUMMARY=$outfile_dir/$unaligned_bam_tagged_CB_MB_summary BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1	

	# FilterBam	
	$dropseq_tools/$filter_bam TAG_REJECT=XQ INPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB OUTPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB_FB COMPRESSION_LEVEL=0 

	# TrimStartingSequence 
	$dropseq_tools/$trim_start INPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB_FB OUTPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart OUTPUT_SUMMARY=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_summary SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5

	## PolyATrimmer
	$dropseq_tools/$polyatrimmer INPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart OUTPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered OUTPUT_SUMMARY=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered_summary MISMATCHES=0 NUM_BASES=6

	## Sam --> Fastq
	## Using picards SamToFastq function
	$java_cmd $PICARD_HOME SamToFastq INPUT=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered FASTQ=$outfile_dir/$unaligned_fastq_mc_tagged_polyA_filtered

	#### STAR alignment
	STAR --runMode alignReads --genomeDir $STAR_index --sjdbGTFfile $gtf_file --readFilesIn $outfile_dir/$unaligned_fastq_mc_tagged_polyA_filtered --outFileNamePrefix $outfile_dir/$file_prefix --runThreadN 12
	
	## SortSam
	$java_cmd $PICARD_HOME SortSam I=$outfile_dir/$file_prefix"Aligned.out.sam" O=$outfile_dir/$aligned_bam_sorted SO=queryname


	# MergeBam
	$java_cmd $PICARD_HOME MergeBamAlignment REFERENCE_SEQUENCE=$reference_fasta UNMAPPED_BAM=$outfile_dir/$unaligned_bam_tagged_CB_MB_trimmed_smart_poly_filtered ALIGNED_BAM=$outfile_dir/$aligned_bam_sorted OUTPUT=$outfile_dir/$merged_bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false


	# TagReadWithGeneExon
	$dropseq_tools/$tagreadwithgeneexon I=$outfile_dir/$merged_bam O=$outfile_dir/$gene_exon_tagged_bam ANNOTATIONS_FILE=$refflat TAG=GE

	# Detect Bead Synthesis Errors
	# DetectBeadSynthesisErrors
	$dropseq_tools/DetectBeadSynthesisErrors I=$outfile_dir/$gene_exon_tagged_bam O=$outfile_dir/$cleaned_bam OUTPUT_STATS=$outfile_dir/$synthesis_stats SUMMARY=$outfile_dir/$synthesis_summary NUM_BARCODES=$NUM_of_barcodes PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC


	# Calculate reads per Cell Barcode for knee plot
	$dropseq_tools/BAMTagHistogram I=$outfile_dir/$cleaned_bam O=$outfile_dir/$bamtaghistogram TAG=XC

	# Run a couple of diagnostic scripts to make check what the UMI distribution etc looks like
	# Gather Molecular Barcode distributions by Gene
	$dropseq_tools/$gatherMolbarDistbyGene I=$outfile_dir/$cleaned_bam O=$outfile_dir/$UMI_by_gene NUM_CORE_BARCODES=$NUM_of_barcodes

	# Calculate single cell metrics
	$dropseq_tools/$singlecellmetrics I=$outfile_dir/$cleaned_bam O=$outfile_dir/$single_cell_metrics ANNOTATIONS_FILE=$gtf_file NUM_CORE_BARCODES=$NUM_of_barcodes

fi
