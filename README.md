# Whole Exome Sequencing
Pipeline to analyze Whole Exome Sequencing data

### Pre-processing and Quality Control of Raw Data ####
Input: Raw FASTQ files.
Tools: FastQC and Trimmomatic (for trimming).

Steps:
Quality Control: Run FastQC on raw FASTQ files to assess quality (per base sequence quality, GC content, adapter contamination, etc.).

Trimming: Use Trimmomatic to remove low-quality bases and adapter sequences from reads.

### Alignment to Reference Genome ###
Input: Trimmed FASTQ files.
Tools: BWA (Burrows-Wheeler Aligner), Samtools, Picard.

Steps:
Alignment: Align reads to the reference genome using BWA-MEM.
Sort BAM: Use Samtools or Picard to sort the resulting BAM file.

Mark Duplicates: Use Picard to mark duplicate reads.

### Post-alignment Processing ###
Input: Deduplicated BAM file.
Tools: GATK (Genome Analysis Toolkit), Picard.

Steps:
Base Quality Score Recalibration (BQSR): Perform BQSR using GATK's BaseRecalibrator and ApplyBQSR.
Index BAM: Index the BAM file for faster access during variant calling.

### Variant Calling ###
Input: Recalibrated BAM file.
Tools: GATK (HaplotypeCaller).

Steps:
Variant Calling: Use GATK's HaplotypeCaller in gVCF mode to call variants for each sample.
Joint Genotyping: Combine individual gVCF files using GATK’s GenotypeGVCFs.

## Variant Recalibration, Filtering and Annotation ###
Input: VCF file from variant calling.
Tools: GATK (VariantRecalibrator, VariantFiltration), SnpEff.
