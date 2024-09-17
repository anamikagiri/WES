#!/bin/bash
#SBATCH -A giria
#SBATCH -J step1
#SBATCH --mail-type=END
#SBATCH --mail-user=anamika.giri@dzne.de
#SBATCH --export=ALL
#SBATCH --partition=work
#SBATCH -n 1
#SBATCH --cpus-per-task 12
#SBATCH --exclude cluster-108



# To run, type sh run_step_1.sh in the command line. This will create a batch script for each sample. This script will create a gvcf file for each sample in the SAMPLEDIR directory.

# If necessary, arguments should be modified in this file.

# All (indexed) resources downloaded from the GATK bundle ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle

# SAMPLEDIR: Directory were the .fastq files are located. For each sample one  .fastq file per read must be present. Name of the files should be sampleID_R1/R2.fastq.gz

# ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. For more info about this argument please visit http://www.usadellab.org/cms/?page=trimmomatic.

	#- fastaWithAdaptersEtc: specifies the path to a fasta file containing all the adapters, PCR sequences etc. The naming of the various sequences within this file determines how they are used. See below.
	
	#- seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed.
	
	#- palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
	
	#- simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.

# LEADING: Integer. Specifies the minimum quality required to keep a base at the beginning of a read. For more info about this argument please visit http://www.usadellab.org/cms/?page=trimmomatic

# TRAILING: Integer.  Specifies the minimum quality required to keep a base at the end of a read. For more info about this argument please visit http://www.usadellab.org/cms/?page=trimmomatic

# SLIDING_WINDOW: windowSize:requiredQuality. windowSize specifies the number of bases to average across, while requiredQuality specifies the average quality required. For more info about this argument please visit http://www.usadellab.org/cms/?page=trimmomatic

# MINLEN: Integer. Specifies the minimum length of reads to be kept. For more info about this argument please visit http://www.usadellab.org/cms/?page=trimmomatic

# BUILD: Genome buld to be used in the pipeline. At this moment only hg19 is supported. Please type hg19.

# BAITFILE: Path to a .bed file with information about the genomic regions captured.

# PADDED: Integer. Number of bases upstream and downstream of the captured regions in which variants should be called. 100bp are recommended.

# MEM: Integer. Maximum memory use for java (in gigabites). Each node in the cluster has 190GB of memory.

# PROC: Maximum number of processor to be used. Usually 12.

# DEPTH_OF_COVERAGE: yes/no. Chose if you want to calculate depth of corage in capture regions per samples. Necessary if XHMM will be run downstream downstream.

# REMOVEFASTQ: Chose whether original .fastq file should be removed from SAMPLEDIR directory. This is recommended for computer memory reasons.

# OUTPUT: Name of output folder proceeded with "/".

#-----------------------------------------------------------#
# Argument variables passed in. User-modifiable parameters. #
#-----------------------------------------------------------#
SAMPLEDIR=Fastq_files/
ILLUMINACLIP=NGS_tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10
LEADING=15
TRAILING=15
SLIDING_WINDOW=4:20
MINLEN=50
BUILD=hg19
BAITFILE=Resources/Intervals/SureSelect_Human_All_Exon_V5_Covered.bed
PADDED=100
MEM=190
PROC=12
DEPTH_OF_COVERAGE=yes
OUTPUT=Pipeline_ouput/
REMOVEFASTQ=no

echo "Arguments:"
echo "SAMPLEDIR=$SAMPLEDIR"
echo "ILLUMINACLIP=$ILLUMINACLIP"
echo "LEADING=$LEADING"
echo "TRAILING=$TRAILING"
echo "SLIDING_WINDOW=$SLIDING_WINDOW"
echo "MINLEN=$MINLEN bp"
echo "BUILD=$BUILD. ucsc.hg19.fasta will be used as reference"
echo "BAITFILE=$BAITFILE"
echo "PADDED=$PADDED bp"
echo "MEM=$MEM GB"
echo "PROC=$PROC processors"
echo "DEPTH_OF_COVERAGE=$DEPTH_OF_COVERAGE"
echo "OUTPUT=$OUTPUT"
echo "REMOVEFASTQ=$REMOVEFASTQ"


#-----------------#
# Define software #
#-----------------#
GATK=NGS_tools/GenomeAnalysisTK-3.8/GenomeAnalysisTK.jar
TRIMMOMATIC=NGS_tools/Trimmomatic-0.36/trimmomatic-0.36.jar
FASTQC=NGS_tools/FastQC/fastqc
BWA=NGS_tools/bwa-0.7.15/bwa
PICARD=NGS_tools/picard-tools-2.6.0/picard.jar


#------------------------------------------#
# Check and setup of the genome build file #
#------------------------------------------#
if [ "$BUILD" = "hg19" ]
then
	REFERENCE=Resources/hg19/ucsc.hg19.fasta
	DBSNP=Resources/SNPs/dbsnp_138.hg19.excluding_sites_after_129.vcf
	INDEL_1KG=Resources/indels/1000G_phase1.indels.hg19.sites.vcf
	INDEL_MILLS=Resources/indels/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
	HAPMAP=Resources/SNPs/hapmap_3.3.hg19.sites.vcf
	OMNI=Resources/SNPs/1000G_omni2.5.hg19.sites.vcf
	SNP_1KG=Resources/SNPs/1000G_phase1.snps.high_confidence.hg19.sites.vcf
else
echo "`date`: At this moment, only only UCSC hg19 is implemented. Please select BUILD = hg19"	
fi


#---------------------------------------#
# Check if OUTPUT folder already exists #
#---------------------------------------#
if [ -d "$OUTPUT" ]
then
	echo "`date`: WARNING: OUPUT folder already exist. Please check there is no overlap between this sample and samples existing in this folder."
else
	echo "`date`: OUTPUT folder does not exist."
fi


#---------------------------------------------------# 
# Extract sample names and create ouput directories #
#---------------------------------------------------#
echo "`date`: Extracting sample name from FASTQ files"
base=$(basename "${1%_R1.fastq.gz}")
echo "`date`: processing ${base}"

echo "`date`: Creating output directories"
mkdir "$OUTPUT"
mkdir "$OUTPUT"QC_reports/
mkdir "$OUTPUT"BAM_files/
mkdir "$OUTPUT"gVCF_files/
mkdir "$OUTPUT"Per_sample_depth_of_coverage/
mkdir "$OUTPUT"tmp/
mkdir "$OUTPUT"tmp/"${base}"
echo "`date`: Output directories created successfully."


#----------------------------------------#
# Trim Fastq files based on base quality #
#----------------------------------------#
echo "`date`: Run Trimmomatic to trim .fastq files based on base quality"
srun java -Xmx"$MEM"g -jar "$TRIMMOMATIC" \
PE \
-threads "$PROC" \
"$SAMPLEDIR${base}_R1.fastq.gz" \
"$SAMPLEDIR${base}_R2.fastq.gz" \
"$SAMPLEDIR${base}_R1_trimmed.fastq.gz" \
"$SAMPLEDIR${base}_R1_unpaired.fastq.gz" \
"$SAMPLEDIR${base}_R2_trimmed.fastq.gz" \
"$SAMPLEDIR${base}_R2_unpaired.fastq.gz" \
ILLUMINACLIP:"$ILLUMINACLIP" \
LEADING:"$LEADING" \
TRAILING:"$TRAILING" \
SLIDINGWINDOW:"$SLIDING_WINDOW" \
MINLEN:"$MINLEN"

wait
echo "`date`: Trimmomatic -> Fastq files trimmed successfully"


#------------------------------------------------------------#
# Create quality control reports for each trimmed fastq file #
#------------------------------------------------------------#
echo "`date`: Create quality control reports for trimmed .fastq files"
srun "$FASTQC" $SAMPLEDIR"${base}"_R1_trimmed.fastq.gz \
-t "$PROC" \
-o "$OUTPUT"QC_reports/

wait

srun "$FASTQC" $SAMPLEDIR"${base}"_R2_trimmed.fastq.gz \
-t "$PROC" \
-o "$OUTPUT"QC_reports/

wait
echo "`date`: FastQC -> Quality control reports created for each _trimmed.fastq.gz file"


# -------------------------------------------------------#
# Run BWA mem on the .fastq files to generate .sam files #
# -------------------------------------------------------#
echo "`date`: Run BWA mem on the .fastq files to generate .sam files"
srun "$BWA" mem \
-t "$PROC" \
-M \
-R "@RG\tID:"${base}"\tPL:ILLUMINA\tLB:"${base}"\tSM:"${base}"\tCN:DZNE" \
   "$REFERENCE" "$SAMPLEDIR${base}_R1_trimmed.fastq.gz" "$SAMPLEDIR${base}_R2_trimmed.fastq.gz" > "$OUTPUT"BAM_files/"${base}.sam"
   
wait
echo "`date`: BWA mem -> Alignment to reference genome completed"


#--------------------------------------------#
# Clean up trimmed and unpaired .fastq files #
#--------------------------------------------#
if [ -f "$OUTPUT"BAM_files/"${base}".sam ]
then
	rm $SAMPLEDIR"${base}"_R1_trimmed.fastq.gz
	rm $SAMPLEDIR"${base}"_R2_trimmed.fastq.gz
	rm $SAMPLEDIR"${base}"_R1_unpaired.fastq.gz
	rm $SAMPLEDIR"${base}"_R2_unpaired.fastq.gz
	echo "`date`: Trimmed FastQ files removed"
else 
	echo ""$OUTPUT"BAM_files/"${base}".sam not found"
fi

wait


# -----------------------------------------------------------#
# Run Picard's SortSam to generate .bam files and sort reads #
# -----------------------------------------------------------#
echo "`date`:  Run Picard to convert SAM to BAM and sort reads" 

module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$PICARD" SortSam \
INPUT="$OUTPUT"BAM_files/"${base}".sam \
OUTPUT="$OUTPUT"BAM_files/"${base}"_sorted_reads.bam \
SORT_ORDER=coordinate \
TMP_DIR="$OUTPUT"tmp/"${base}"/ \
VALIDATION_STRINGENCY=SILENT \
VERBOSITY=ERROR \
MAX_RECORDS_IN_RAM=2500000

wait 
echo "`date`:  Picard -> .sam to .bam conversion and read sorting completed."


#----------------------------------------------------#
# Clean up sam files and temporary files within /tmp #
#----------------------------------------------------#
if [ -f "$OUTPUT"BAM_files/"${base}"_sorted_reads.bam ]
then
	rm "$OUTPUT"BAM_files/"${base}".sam
	echo "`date`: .sam file removed"
else 
	echo ""$OUTPUT"BAM_files/"${base}"_sorted_reads.bam not found"
fi

rm "$OUTPUT"tmp/"${base}"/*
wait


# -----------------------------------------------------------------#
# Run Picard's MarkDuplicates to flag PCR duplicates in .bam files #
# -----------------------------------------------------------------#
echo "`date`:  Run Picard's MarkDuplicates to flag PCR duplicates in .bam files"

module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$PICARD" MarkDuplicates \
INPUT="$OUTPUT"BAM_files/"${base}"_sorted_reads.bam \
OUTPUT="$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam \
METRICS_FILE="$OUTPUT"BAM_files/"${base}"_dedup_metrics.txt \
REMOVE_DUPLICATES=false \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR="$OUTPUT"tmp/"${base}"/ \
MAX_RECORDS_IN_RAM=2500000 \
VERBOSITY=ERROR

wait
echo "`date`:  Picard -> Duplicates marked"


#------------------------------------------------------------------#
# Clean up intermediate .bam files and temporary files within /tmp #
#------------------------------------------------------------------#
if [ -f "$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam ]
then 
	rm "$OUTPUT"BAM_files/"${base}"_sorted.bam
	echo "`date`: Intermediate .bam file removed"
else 
	echo ""$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam not found"
fi

rm "$OUTPUT"tmp/"${base}"/*
wait


# -------------------------------------------#
# Index sorted, duplicates-marked .bam files #
# -------------------------------------------#
echo "`date`:  Run Picard's BuildBamIndex for indexing of the sorted, duplicates-marked .bam files"

module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$PICARD" BuildBamIndex \
INPUT="$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam \
TMP_DIR="$OUTPUT"tmp/"${base}"/ \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=2500000 \
VERBOSITY=ERROR
	
wait
echo "`date`:  Picard -> Index created"


#----------------------------------#
# Generate Diagnose Target reports #
#----------------------------------#
echo "`date`:  Generate Diagnose Target reports"

module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T DiagnoseTargets \
-R "$REFERENCE" \
-I "$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam \
-L "$BAITFILE" \
-o "$OUTPUT"QC_reports/"${base}".DiagnoseTargets

wait
echo "`date`: GATK -> Diagnose Target reports created"


#---------------------------------------------------------------------------#
# Calculate per-sample depth of coverage (useful for CNV calling with XHMM) #
#---------------------------------------------------------------------------#
if [ "$DEPTH_OF_COVERAGE" = "yes" ]
then
	module load java/1.8.0
	srun java -Xmx"$MEM"g -jar "$GATK" \
	-T DepthOfCoverage \
	-I "$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam \
	-L "$BAITFILE" \
	-R "$REFERENCE" \
	-dt BY_SAMPLE \
	-dcov 5000 \
	-l INFO \
	--omitDepthOutputAtEachBase \
	--omitLocusTable \
	--minBaseQuality 0 \
	--minMappingQuality 20 \
	--start 1 \
	--stop 5000 \
	--nBins 200 \
	--includeRefNSites \
	--countType COUNT_FRAGMENTS \
	-o "$OUTPUT"Per_sample_depth_of_coverage/"${base}"

	echo "`date`: GATK -> Per-sample depth of coverage calculated"

else 
	echo "`date`: Per-sample depth of coverage not calculated"
fi

wait


#---------------------------------------------------------------#
# BSQR: Analyze patterns of covariation in the sequence dataset #
#---------------------------------------------------------------#
echo "`date`: Generate recalibration table based on user-specified covariates"

module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T BaseRecalibrator \
-nct "$PROC" \
-S SILENT \
-R "$REFERENCE" \
-L "$BAITFILE" \
-ip "$PADDED" \
-knownSites "$DBSNP" \
-knownSites "$INDEL_1KG" \
-knownSites "$INDEL_MILLS" \
-l ERROR \
-I "$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam \
-o "$OUTPUT"BAM_files/"${base}".recal_data.table
	
wait
echo "`date`: GATK -> Recalibration table created"


#---------------------------------------------------------#
# BSQR: Analyze covariation remaining after recalibration #
#---------------------------------------------------------#
echo "`date`: Analyze covariation remaining after recalibration"

module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T BaseRecalibrator \
-nct "$PROC" \
-S SILENT \
-R "$REFERENCE" \
-L "$BAITFILE" \
-ip "$PADDED" \
-knownSites "$DBSNP" \
-knownSites "$INDEL_1KG" \
-knownSites "$INDEL_MILLS" \
-BQSR "$OUTPUT"BAM_files/"${base}".recal_data.table \
-l ERROR \
-I "$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam \
-o "$OUTPUT"BAM_files/"${base}".post_recal_data.table
	
wait
echo "`date`: GATK -> Post-recalibration table created"


#-------------------------------------------------#
# BSQR: Generate before/after recalibration plots #
#-------------------------------------------------#
echo "`date`: Second pass to analyze covariation remaining after recalibration"

module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T AnalyzeCovariates \
-R "$REFERENCE" \
-L "$BAITFILE" \
-before "$OUTPUT"BAM_files/"${base}".recal_data.table \
-after "$OUTPUT"BAM_files/"${base}".post_recal_data.table \
-plots "$OUTPUT"BAM_files/"${base}".recalibration_plots.pdf

wait
echo "`date`: GATK -> before/after recalibration plots generated"


#-----------------------------------------------------#
# BSQR: Apply the recalibration to your sequence data #
#-----------------------------------------------------#
echo "`date`:  Apply the recalibration to your sequence data"

module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T PrintReads \
-nct "$PROC" \
-S SILENT \
-R "$REFERENCE" \
-L "$BAITFILE" \
-l ERROR \
-I "$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam \
-BQSR "$OUTPUT"BAM_files/"${base}".recal_data.table \
-o "$OUTPUT"BAM_files/"${base}"_sorted_dedup_recal_reads.bam
	
wait
echo "`date`:  GATK -> Sequence data recalibrated"


#------------------------------------------------------------------#
# Clean up intermediate .bam files and temporary files within /tmp #
#------------------------------------------------------------------#
if [ -f "$OUTPUT"BAM_files/"${base}"_sorted_dedup_recal_reads.bai ]
then 
	rm "$OUTPUT"BAM_files/"${base}"_sorted_reads.bam
	rm "$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bam
	rm "$OUTPUT"BAM_files/"${base}"_sorted_dedup_reads.bai
	rm "$OUTPUT"BAM_files/"${base}".recal_data.table
	rm "$OUTPUT"BAM_files/"${base}".post_recal_data.table
	rm "$OUTPUT"BAM_files/"${base}"_dedup_metrics.txt
	echo "`date`: Intermediate .bam files and duplication and recalibration tables removed"
else 
	echo ""$OUTPUT"BAM_files/BAM_files/"${base}"_sorted_dedup_recal_reads.bai not found"
fi

rm -r "$OUTPUT"tmp/"${base}"
wait

#--------------------#
# Collect HS metrics #
#--------------------#
echo "`date`:  Run Picard's CollectHsMetrics to calculate .bam file statistics"

if [ -f "${BAITFILE%.bed}".interval_list ]
then 
	echo "`date`:  Picard interval list detected for the BAITFILE used"
else
	echo "`date`:  Picard interval list not detected for the BAITFILE used. BedToIntervalList will be used to generate it."

	module load java/1.8.0
	srun java -Xmx"$MEM"g -jar "$PICARD" BedToIntervalList \
	I="$BAITFILE" \
	O="${BAITFILE%.bed}".interval_list \
	SD="${REFERENCE%.fasta}.dict"
fi
    
module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$PICARD" CollectHsMetrics \
INPUT="$OUTPUT"BAM_files/"${base}"_sorted_dedup_recal_reads.bam \
R="$REFERENCE" \
BAIT_INTERVALS="${BAITFILE%.bed}".interval_list \
TARGET_INTERVALS="${BAITFILE%.bed}".interval_list \
OUTPUT="$OUTPUT"QC_reports/"${base}"_HS_metrics.txt \
TMP_DIR="$OUTPUT"tmp/"${base}"/

wait
echo "`date`:  Picard -> HS metrics collected"


#--------------------------------------------#
# Run the Haplotype Caller on each BAM file. #
#--------------------------------------------#
echo "`date`:  Run GATK's Haplotype Caller to create a gVCF file for each sample"
module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T HaplotypeCaller \
-nct "$PROC" \
-R "$REFERENCE" \
-I "$OUTPUT"BAM_files/"${base}"_sorted_dedup_recal_reads.bam \
-L "$BAITFILE" \
-ip "$PADDED" \
--genotyping_mode DISCOVERY \
--emitRefConfidence GVCF \
--dbsnp "$DBSNP" \
-o "$OUTPUT"gVCF_files/"${base}".raw.g.vcf

wait
echo "`date`:  GATK -> Genotyping with HaplotypeCaller done. A gVCF file has been created for ${base}"


#---------------------#
# Remove .fastq files #
#---------------------#
if [ "$REMOVEFASTQ" = "yes" ]
then
	if [ -f "$OUTPUT"gVCF_files/"${base}".raw.g.vcf ]
		then
		rm $SAMPLEDIR"${base}"_R1.fastq.gz
		rm $SAMPLEDIR"${base}"_R2.fastq.gz
		echo "`date`: "${base}"_R1.fastq.gz file removed"
		echo "`date`: "${base}"_R2.fastq.gz file removed"
		else "`date`: "${base}".raw.snps.indels.g.vcf not found. Corresponding .fastq file will not be removed."
		fi
else 
	echo ".fastq files will not be removed from $SAMPLEDIR"
fi

wait
