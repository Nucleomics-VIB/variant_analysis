#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2019-12-12
# from: https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020

## Where are we?
workdir="/data/NC_projects/4120_ARuiz_chicken_variants"

## remap reads from SRR17149558 to gallus_gallus GRCg6a.105
cd ${workdir}

# create tmp folder
mkdir -p tmpfiles

# multithreading, adapt to your cpu count
samtoolsthr=24

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

reference_fa="/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa"

# base quality score recalibration
dbsnp="/data/biodata/references/galGal6a/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz"

infolder=gatk4_BP
outfolder=gatk_BP_recal
mkdir -p ${outfolder}

# more records in RAM speeds up when enough RAM is present
recinram=10000000

#################
# function
#################

runall () {

# get bam from the loop
bam=$1
pfx=$(basename ${bam%_mrkdup_srt-tags.bam})
recalbamfile=${pfx}_recal.bam

#############################################
# GATK4 BAM RECALIBRATION
#############################################

if [ ! -f gatk_preprocessing/recalibration_done ]; then

# compute table before
java ${javaopts} -jar $GATK/gatk.jar \
	BaseRecalibrator \
	-I ${bam} \
	-R ${reference_fa} \
	--known-sites ${dbsnp} \
	-O ${outfolder}/${pfx}_recal_data.table \
	--tmp-dir tmpfiles

# apply recalibration table
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyBQSR \
	-I ${bam} \
	-R ${reference_fa} \
	-bqsr ${outfolder}/${pfx}_recal_data.table \
	-O ${outfolder}/${recalbamfile} \
	--interval-padding 100 \
	--add-output-sam-program-record \
	--use-original-qualities \
	--static-quantized-quals 10 \
	--static-quantized-quals 20 \
	--static-quantized-quals 30 \
	--tmp-dir tmpfiles

# compute table after
java ${javaopts} -jar $GATK/gatk.jar \
	BaseRecalibrator \
	-I ${outfolder}/${recalbamfile} \
	-R ${reference_fa} \
	--known-sites ${dbsnp} \
	-O ${outfolder}/${pfx}_recal_data_after.table \
	--tmp-dir tmpfiles

# create plots from both tables
java ${javaopts} -jar $GATK/gatk.jar \
	AnalyzeCovariates \
	-before ${outfolder}/${pfx}_recal_data.table \
	-after ${outfolder}/${pfx}_recal_data_after.table \
	-plots ${outfolder}/${pfx}_BQSR_report.pdf \
	-csv ${outfolder}/${pfx}_BQSR-report.csv \
	--tmp-dir tmpfiles

# Picard CollectMultipleMetrics on final BAM
java ${javaopts} -jar $PICARD/picard.jar \
	CollectMultipleMetrics \
	I=${outfolder}/${recalbamfile} \
	O=${outfolder}/${pfx}_multiple_metrics \
	R=${reference_fa} \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles

fi

}

for bam in ${infolder}/*_mrkdup_srt-tags.bam; do
runall ${bam} &
done

wait
