#!/bin/bash
#SBATCH -A giria
#SBATCH -J step2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=anamika.giri@dzne.de
#SBATCH --export=ALL
#SBATCH --cpus-per-task 64
#SBATCH -N 1
#SBATCH --nodelist compute-03




# To run, type sbatch step2_joint_genotyping_and_recalibration.sh in the command line.

# If necessary, arguments should be modified in this file.

# All (indexed) resources downloaded from the GATK bundle ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle

# GVCF: Text file containing a list of locations of all gVCFs. Extension of file should be .list

# BUILD: Genome buld to be used in the pipeline. At this moment only hg19 is supported. Please type hg19.

# MAX_ALTERNATE_ALLELES: Maximum number of permited alternate alleles. Default is 6.

# VQSLOD: Tranche sensitivity threshold with ApplyRecalibration, expressed as a percentage (e.g. 99.9%).The program looks at what is the VQSLOD value above which 99.9% of the variants in the training callset are included. It then takes that value of VQSLOD and uses it as a threshold to filter your variants. Variants that are above the threshold pass the filter, so the FILTER field will contain PASS. Variants that are below the threshold will be filtered out; they will be written to the output file, but in the FILTER field they will have the name of the tranche they belonged to.

# TITV: The expected novel Ti/Tv ratio to use when calculating FDR tranches and for display on the optimization curve output figures. Default: 2.15.

# MEM: Integer. Maximum memory use for java (in gigabites). Each node in the cluster has 190GB of memory.

# PROC: Maximum number of processor to be used. Usually 64.

# HARD_FILTERING: yes/no. Perform genotype and variant QC.

# MISSINGNESS: yes/no. Filter variants with a missingness above the theshold specified with MAX_MISSING argument. A file containing a list of the male samples in the dataset must be specified.

# MAX_MISSING: Maximum ratio of missingness (0 - 1) of a variant across all samples.

# snpEff: yes/no. Perform default annotation of variants using snpEff. Please note that annotations will be included in the "classic" snpEFF format (ANN=) and that bgzip comression will be performed.

# MALE_SAMPLES: Text file containing a list of male samples in the dataset.

# OUTPUT: Folder where the output should be written.


#-----------------------------------------------------------#
# Argument variables passed in. User-modifiable parameters. #
#-----------------------------------------------------------#
GVCF=gVCF_list.list
BUILD=hg19
MAX_ALTERNATE_ALLELES=7
VQSLOD=99.9
TITV=2.8
MEM=1000
PROC=16
HARD_FILTERING=yes
MISSINGNESS=no
MAX_MISSING=0.1
MALE_SAMPLES=no
snpEff=yes
OUTPUT=VCF_files/WES_FTD_2017/

echo "Arguments:"
echo "GVCF=$GVCF"
echo "BUILD=$BUILD. ucsc.hg19.fasta will be used as reference"
echo "MAX_ALTERNATE_ALLELES=$MAX_ALTERNATE_ALLELES"
echo "VQSLOD=$VQSLOD"
echo "TITV=$TITV"
echo "VQSLOD=$VQSLOD"
echo "MEM=$MEM GB"
echo "PROC=$PROC processors"
echo "HARD_FILTERING=$HARD_FILTERING"
echo "MISSINGNESS=$MISSINGNESS"
echo "MAX_MISSING=$MAX_MISSING"
echo "MALE_SAMPLES=$MALE_SAMPLES"
echo "OUTPUT=$OUTPUT"

#-----------------#
# Define software #
#-----------------#
GATK=NGS_tools/GenomeAnalysisTK-3.8/GenomeAnalysisTK.jar
VCFTOOLS=./NGS_tools/vcftools/bin/vcftools
BCFTOOLS=./NGS_tools/bcftools/bcftools
SNPEFF=NGS_tools/snpEff_latest_core/snpEff/snpEff.jar
BGZIP=NGS_tools/tabix-0.2.6/bgzip
TABIX=NGS_tools/tabix-0.2.6/tabix

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
	echo "At this moment, only only UCSC hg19 is implemented. Please select BUILD = hg19"	
fi


#----------------------------------------------#
# Perform joint genotyping using GenotypeGVCFs #
#----------------------------------------------#
echo "`date`:  Run GATK's GenotypeGVCFs to perform the joint genotyping"
module load java/1.8.0
srun java -Xmx"$MEM"g -Djava.io.tmpdir=Temporal/ -jar "$GATK" \
-T GenotypeGVCFs \
-R "$REFERENCE" \
--variant "$GVCF" \
--dbsnp "$DBSNP" \
--disable_auto_index_creation_and_locking_when_reading_rods \
--max_alternate_alleles "$MAX_ALTERNATE_ALLELES" \
-nt "$PROC" \
-o "$OUTPUT"All_samples_raw.snps.indels.vcf

rm Temporal/*

echo "`date`:  GATK -> Joint genotyping done. A single .vcf file has been created for all samples."


#-------------------------------------------------------------------#	
# Variant recalibration for SNPs: Build the SNP recalibration model #
#-------------------------------------------------------------------#	
echo "`date`:  Run VSQR to perform variant filtering for SNPs"
module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T VariantRecalibrator \
-R "$REFERENCE" \
-input "$OUTPUT"All_samples_raw.snps.indels.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$HAPMAP" \
-resource:omni,known=false,training=true,truth=true,prior=12.0 "$OMNI" \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 "$SNP_1KG" \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
-recalFile "$OUTPUT"All_samples_recalibrate_SNP.recal \
-tranchesFile "$OUTPUT"All_samples_recalibrate_SNP.tranches \
-rscriptFile "$OUTPUT"All_samples_recalibrate_SNP.plots.R \
-an QD \
-an MQ \
-an FS \
-an SOR \
-an MQRankSum \
-an ReadPosRankSum \
-titv "$TITV" \
-mode SNP \
-nt "$PROC" \
--max_attempts 5 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0

echo "`date`:  GATK -> SNP recalibration modeling finished."


#--------------------------------------------------------------------------#
# Variant recalibration for SNPs: Apply the desired level of recalibration #
#--------------------------------------------------------------------------#
echo "`date`:  Apply recalibration to SNPs"
module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T ApplyRecalibration \
-R "$REFERENCE" \
-input "$OUTPUT"All_samples_raw.snps.indels.vcf \
-mode SNP \
--ts_filter_level "$VQSLOD" \
-recalFile "$OUTPUT"All_samples_recalibrate_SNP.recal \
-tranchesFile "$OUTPUT"All_samples_recalibrate_SNP.tranches \
-o "$OUTPUT"All_samples_recalibrated_snps_raw_indels.vcf

echo "`date`:  GATK -> SNPs recalibrated succesfully."


#-----------------------------------------------------------------------#
# Variant recalibration for indels: Build the indel recalibration model #
#-----------------------------------------------------------------------#
echo "`date`:  Run VSQR to perform variant filtering for indels"
module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T VariantRecalibrator \
-R "$REFERENCE" \
-input "$OUTPUT"All_samples_recalibrated_snps_raw_indels.vcf \
--maxGaussians 4 \
-resource:mills,known=false,training=true,truth=true,prior=12.0 "$INDEL_MILLS" \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
-an QD \
-an FS \
-an SOR \
-an MQRankSum \
-an ReadPosRankSum \
-mode INDEL \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile "$OUTPUT"All_samples_recalibrate_INDEL.recal \
-tranchesFile "$OUTPUT"All_samples_recalibrate_INDEL.tranches \
-nt "$PROC" \
-rscriptFile "$OUTPUT"All_samples_recalibrate_INDEL.plots.R 

echo "`date`:  GATK -> Indel recalibration modeling finished."


#----------------------------------------------------------------------------#
# Variant recalibration for indels: Apply the desired level of recalibration #
#----------------------------------------------------------------------------#
echo "`date`:  Apply recalibration to indels"
module load java/1.8.0
srun java -Xmx"$MEM"g -jar "$GATK" \
-T ApplyRecalibration \
-R "$REFERENCE" \
-input "$OUTPUT"All_samples_recalibrated_snps_raw_indels.vcf \
-mode INDEL \
--ts_filter_level "$VQSLOD" \
-recalFile "$OUTPUT"All_samples_recalibrate_INDEL.recal \
-tranchesFile "$OUTPUT"All_samples_recalibrate_INDEL.tranches \
-o "$OUTPUT"All_samples_recalibrated_variants.vcf

echo "`date`:  GATK -> Indels recalibrated succesfully."


#-------------------------------------------#
# Normalize VCF and split multiallelic loci #
#-------------------------------------------#
echo "`date`:  Normalize indels and split multiallelic loci. All variant IDs will be changed to CHROM:POS:REF:ALT"
srun "$BCFTOOLS" norm -Ou -m -any "$OUTPUT"All_samples_recalibrated_variants.vcf |
srun "$BCFTOOLS" norm -Ou -f "$REFERENCE" |
srun "$BCFTOOLS" annotate -Ov -I '%CHROM:%POS:%REF:%ALT' > "$OUTPUT"All_samples_recalibrated_and_normalized_variants.vcf 

echo "`date`:  bcftools -> Indels normalized and  multiallelic loci splitted. Variant IDs changed to CHROM:POS:REF:ALT"


#---------------------------#
# Delete intermediate files #
#---------------------------#
if [ -f "$OUTPUT"All_samples_recalibrated_and_normalized_variants.vcf   ]
then 
	rm "$OUTPUT"All_samples_recalibrate_SNP.recal
	rm "$OUTPUT"All_samples_recalibrate_SNP.recal.idx
	rm "$OUTPUT"All_samples_recalibrate_INDEL.recal
	rm "$OUTPUT"All_samples_recalibrate_INDEL.recal.idx
	rm "$OUTPUT"All_samples_raw.snps.indels.vcf
	rm "$OUTPUT"All_samples_raw.snps.indels.vcf.idx
	rm "$OUTPUT"All_samples_recalibrated_snps_raw_indels.vcf
	rm "$OUTPUT"All_samples_recalibrated_snps_raw_indels.vcf.idx
	
	echo "`date`:  File All_samples_recalibrated_and_normalized_variants.vcf found. All intermediate files deleted"
else 
	echo "All_samples_recalibrated_variants.vcf not found"
fi


#-------------------------------------------------#
# Annotate variant calls with context information #
#-------------------------------------------------#
if [ "$HARD_FILTERING" = "yes" ]
then
	echo "`date`:  Annotate variant calls with context information"
	module load java/1.8.0
	srun java -Xmx"$MEM"g -jar $GATK \
	-R "$REFERENCE" \
	-T VariantAnnotator \
	--variant "$OUTPUT"All_samples_recalibrated_and_normalized_variants.vcf \
	--useAllAnnotations \
	--out "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.vcf

	echo "`date`:  GATK -> Variant fully annotated."


	#------------------#
	# Variant QC: SNPs #
	#------------------#
	echo "`date`:  Apply hard filters to SNPs"
	module load java/1.8.0
	srun java -Xmx"$MEM"g -jar $GATK \
	-R "$REFERENCE" \
	-T SelectVariants \
	--variant "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.vcf \
	-selectType SNP \
	-o "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.SNP.vcf \
                
	module load java/1.8.0
	srun java -Xmx"$MEM"g -jar $GATK \
	-R "$REFERENCE" \
	-T VariantFiltration \
	--variant "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.SNP.vcf \
	--filterExpression 'QD < 2.0' \
	--filterName 'QD_filter' \
	--filterExpression 'MQ < 40.0' \
	--filterName 'MQ_filter' \
	--filterExpression 'FS > 60.0' \
	--filterName 'FS_filter' \
	--filterExpression 'SOR > 3.0' \
	--filterName 'SOR_filter' \
	--filterExpression 'MQRankSum < -12.5' \
	--filterName 'MQRankSum_filter' \
	--filterExpression 'ReadPosRankSum < -8.0' \
	--filterName 'ReadPosRankSum_filter' \
	-o "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.SNP.variant.QC.vcf \

	echo "`date`:  GATK -> Hard filtering applied to SNPs."


	#--------------------#
	# Variant QC: indels #
	#--------------------#
	echo "`date`:  Apply hard filters to indels"
	module load java/1.8.0
	srun java -Xmx"$MEM"g -jar $GATK \
	-R "$REFERENCE" \
	-T SelectVariants \
	--variant "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.vcf \
	-selectType INDEL \
	-o "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.INDEL.vcf

	module load java/1.8.0
	srun java -Xmx"$MEM"g -jar $GATK \
	-R "$REFERENCE" \
	-T VariantFiltration \
	--variant "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.INDEL.vcf \
	--filterExpression "QD < 2.0" \
	--filterName "QD_filter" \
	--filterExpression "ReadPosRankSum < -20.0" \
	--filterName "ReadPosRankSum_filter" \
	--filterExpression "FS > 200.0" \
	--filterName "FS_filter" \
	--filterExpression 'SOR > 10.0' \
	--filterName 'SOR_filter' \
	-o "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.INDEL.variant.QC.vcf

	echo "`date`:  GATK -> Hard filtering applied to indels."


	#----------------------------------#
	# Combine filtered SNPs and indels #
	#----------------------------------#
	module load java/1.8.0
	srun java -Xmx"$MEM"g -jar $GATK \
	-R "$REFERENCE" \
	-T CombineVariants \
	--variant "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.SNP.variant.QC.vcf \
	--variant "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.INDEL.variant.QC.vcf \
	--assumeIdenticalSamples \
	-o "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.variant.QC.vcf

	echo "`date`:  GATK -> Variant QC applied to both SNPs and indels."


	#-------------#
	# Genotype QC #
	#-------------#
	srun "$VCFTOOLS" \
	--vcf "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.variant.QC.vcf \
	--mac 1 \
	--remove-filtered-all \
	--minGQ 20 \
	--minDP 10 \
	--recode \
	--recode-INFO-all \
	--stdout | gzip -c > "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.vcf.gz
	
	echo "`date`:  VCF tools -> Genotype QC applied. A file named All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.vcf.gz has been created."


	#---------------------------#
	# Delete intermediate files #
	#---------------------------#
	if [ -f  "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.vcf.gz ]
	then 
		rm "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.vcf*
		rm "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.SNP.vcf*
		rm "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.SNP.variant.QC.vcf*
		rm "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.INDEL.vcf*
		rm "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.INDEL.variant.QC.vcf*
		rm "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.variant.QC.vcf
		rm "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.variant.QC.vcf.idx

		echo "`date`:  File All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.vcf.gz found. All intermediate files deleted"
	else 
		echo "All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.vcf.gz not found"
	fi

else
	echo "`date`:  Variant and genotype QC has not been performed."
fi


#-----------------#
# Missingness QC  #
#-----------------#
if [ "$MISSINGNESS" = "yes" ]
then
	srun "$VCFTOOLS" \
	--gzvcf "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.vcf.gz \
	--keep "$MALE_SAMPLES" \
	--chr chrY \
	--missing-site --out "$OUTPUT"MALE_SAMPLES

	srun "$VCFTOOLS" \
	--gzvcf "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.vcf.gz \
	--missing-site \
	--out "$OUTPUT"ALL

	cat "$OUTPUT"ALL.lmiss | grep -v chrY | grep -v CHR | awk '$6>'$MAX_MISSING'{print$1,$2}' > "$OUTPUT"Sites_with_missingness_over_"$MAX_MISSING".txt

	cat "$OUTPUT"MALE_SAMPLES.lmiss | grep -v CHR | awk '$6>'$MAX_MISSING'{print$1,$2}' >> "$OUTPUT"Sites_with_missingness_over_"$MAX_MISSING".txt

	srun "$VCFTOOLS" \
	--gzvcf "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.vcf.gz \
	--mac 1 \
	--exclude-positions "$OUTPUT"Sites_with_missingness_over_"$MAX_MISSING".txt \
	--recode \
	--recode-INFO-all \
	--stdout | gzip -c > "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Max_missing_"$MAX_MISSING".Mac1.vcf.gz

	echo "`date`:  VCF tools -> Missingness QC finished. File All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.vcf.gz.Max_missing_${MAX_MISSING}.Mac1.vcf.gz has been created."

else
	echo "`date`: Missingness QC not performed"
fi


#-----------------------------------------#
# snpEff annotation and bgzip compression #
#-----------------------------------------#
if [ "$snpEff" = "yes" ]
then
	module load java/1.8.0
	echo "`date`:  Basic annotation with snpEff"
	
	if [ -f "$OUTPUT"All_samples_recalibrated_and_normalized_variants.vcf   ]
 		srun java -Xmx"$MEM"g -jar "$SNPEFF" -i vcf -o vcf hg19 -classic -formatEff "$OUTPUT"All_samples_recalibrated_and_normalized_variants.vcf \
   		| "$BGZIP" -c > "$OUTPUT"All_samples_recalibrated_and_normalized_variants.snpEff_annotated.vcf 
   		"$TABIX" -p vcf "$OUTPUT"All_samples_recalibrated_and_normalized_variants.snpEff_annotated.vcf 
	else
		echo "`date`: All_samples_recalibrated_and_normalized_variants.vcf does not exist"
	fi
	
			
	if [ -f "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Mac1.vcf.gz   ]
	then
   		srun java -Xmx"$MEM"g -jar "$SNPEFF" -i vcf -o vcf hg19 -classic -formatEff "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Mac1.vcf.gz \
   		| "$BGZIP" -c > "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.snpEff_annotated.vcf.gz
		"$TABIX" -p vcf "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.mac1.snpEff_annotated.vcf.gz
	else
		echo "`date`: All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Mac1.snpEff_annotated.vcf.gz does not exist"
	fi
	
	if [ -f "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Max_missing_"$MAX_MISSING".Mac1.vcf.gz   ]
	then
		srun java -Xmx"$MEM"g -jar "$SNPEFF" -i vcf -o vcf hg19 -classic -formatEff "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Max_missing_"$MAX_MISSING".Mac1.vcf.gz \
   		| "$BGZIP" -c > "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Max_missing_"$MAX_MISSING".Mac1.snpEff_annotated.vcf.gz
		"$TABIX" -p vcf "$OUTPUT"All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Max_missing_"$MAX_MISSING".Mac1.snpEff_annotated.vcf.gz
	else
		echo "`date`: All_samples_recalibrated_and_normalized_variants.annot.GT.variant.QC.Max_missing_"$MAX_MISSING".Mac1.snpEff_annotated.vcf.gz does not exist"
	fi	
	
	
else
	echo "`date`: snpEff annotationn and bgzip compression will not be performed"
fi
