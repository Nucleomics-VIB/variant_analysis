#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2022-03-18
# from: https://github.com/BITS-VIB/NGS-Variant-Analysis-training-2020

# process data after elprep was ran

## Where are we?
workdir="/data/NC_projects/4120_ARuiz_chicken_variants"
cd ${workdir}

# create tmp folder
mkdir -p ${workdir}tmpfiles

# create date-tagged logs 
mkdir -p ${workdir}/logs
logfile=${workdir}logs/gatk_runlog.txt
cat /dev/null > ${logfile}

# record all script actions
#set -x
#exec > ${log} 2>&1
exec &> >(tee -i ${logfile})

# multithreading, adapt to your cpu count
bwathr=84
samtoolsthr=24

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

# references (indexed)
reference_fa=/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa
dbsnp=/data/biodata/references/galGal6a/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz

infolder=elprep

# add index to g.VCF files generated by elprep
for vcf in ${infolder}/[5W]*.g.vcf.gz; do
java ${javaopts} -jar $GATK/gatk.jar \
  IndexFeatureFile \
  -I ${vcf}
done

outfolder=gatk_variants
mkdir -p ${outfolder}

#############################################
# GATK4 MERGE GVCF VARIANT FILES
#############################################

# convert to multi-VCF (merge samples)
java ${javaopts} -jar $GATK/gatk.jar \
	CombineGVCFs \
	-R ${reference_fa} \
	-V ${infolder}/524.g.vcf.gz \
	-V ${infolder}/526.g.vcf.gz \
	-V ${infolder}/528.g.vcf.gz \
	-V ${infolder}/530.g.vcf.gz \
	-V ${infolder}/W201120E14-EMB.g.vcf.gz \
	-V ${infolder}/W291020E14-EMB.g.vcf.gz \
	-V ${infolder}/W291020E15-EMB.g.vcf.gz \
	-O ${outfolder}/combined_7.g.vcf.gz \
	--tmp-dir tmpfiles/

# Genotype (multi-)gVCF to (multi-)VCF and add dbSNP IDs
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
java ${javaopts} -jar $GATK/gatk.jar \
	GenotypeGVCFs \
	--reference ${reference_fa} \
	--variant ${outfolder}/combined_7.g.vcf.gz \
	--output ${outfolder}/combined_7.vcf.gz \
	--dbsnp ${dbsnp} \
	--use-new-qual-calculator \
	--tmp-dir tmpfiles/

#############################################
# GATK4 VARIANT HARD FILTERING
#############################################

# from https://www.nature.com/articles/s41422-020-0349-y
# - MQ < 25.0
# - QUAL < 40.0
# - MQ0 ≥ 4 && ((MQ0/(1.0*DP)) > 0.1)
# --cluster-size 3
# --cluster-window-size 10 (flag more than 3 clustered variants within 10bps +++)

# instructions and filters from:
# https://gatkforums.broadinstitute.org/gatk/discussion/23216/how-to-filter-variants-either-with-vqsr-or-by-hard-filtering
# Genepattern: QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30

# This produces a VCF with the same variant records now annotated with filter status. Specifically, if a record passes all the filters, it receives a PASS label in the FILTER column.
# A record that fails a filter #receives the filter name in the FILTER column, e.g. SOR3.
# If a record fails multiple filters, then each failing filter name appears in the FILTER column separated by semi-colons ; e.g. "MQRankSum-12.5;ReadPosRankSum-8".

infolder=${outfolder}
outfolder=gatk_varianthardfiltering
mkdir -p ${outfolder}

###########################################
# 1) hard-Filter SNPs on multiple metrics
###########################################

# produces a VCF with records with SNP-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
	SelectVariants \
	-V ${infolder}/combined_7.vcf.gz \
	-O ${outfolder}/combined_7_snp.vcf.gz \
	--select-type-to-include SNP

threshold=54.69

java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	-V ${outfolder}/combined_7_snp.vcf.gz \
	-O ${outfolder}/combined_7_snp_filtered.vcf.gz \
	--filter "QD < 2.0" --filter-name "QD2" \
	--filter "QUAL < 30.0" --filter-name "QUAL30" \
	--filter "SOR > 3.0" --filter-name "SOR3" \
	--filter "FS > 60.0" --filter-name "FS60" \
	--filter "MQ < 40.0" --filter-name "MQ40" \
	--filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  	--filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	--filter-expression "ExcessHet > ${threshold}" --filter-name "ExcessHet" \
    --cluster-size 3 \
    --cluster-window-size 10

# not applicable to chicken data
#	--filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
# 	--filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \

######################################################
# 2) hard-Filter INDELs and MIXED on multiple metrics
######################################################

# produces a VCF with records with INDEL-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
	SelectVariants \
	-V ${infolder}/combined_7.vcf.gz \
	-O ${outfolder}/combined_7_indel.vcf.gz
	--select-type-to-include INDEL \
	--select-type-to-include MIXED

threshold=54.69

java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	-V ${outfolder}/combined_7_indel.vcf.gz \
	-O ${outfolder}/combined_7_indel_filtered.vcf.gz \
	--filter "QD < 2.0" --filter-name "QD2" \
	--filter "QUAL < 30.0" --filter-name "QUAL30" \
	--filter "FS > 200.0" --filter-name "FS200" \
	--filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	--filter-expression "ExcessHet > ${threshold}" --filter-name "ExcessHet" \
    --cluster-size 3 \
    --cluster-window-size 10
	
# not applicable to chicken data
#	--filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

###########################################
## 3) merge SNP and Indel filtered calls	
###########################################

# combine filtered into one file using Picard
java ${javaopts} -jar $GATK/gatk.jar \
	MergeVcfs \
	-I ${outfolder}/combined_7_snp_filtered.vcf.gz \
	-I ${outfolder}/combined_7_indel_filtered.vcf.gz \
	-R ${reference_fa} \
	-O ${outfolder}/combined_7_snp_indel_filtered.vcf.gz

#############################################
# SnpEff ANNOTATION
#############################################

infolder=${outfolder}
outfolder="snpeff"
mkdir -p ${outfolder}

build="GRCg6a.105"

java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfolder}/gatk_hardfiltering_snpEff_summary.html \
	-nodownload \
	${build} \
	${infolder}/combined_7_snp_indel_filtered.vcf.gz | \
	bgzip -c > ${outfolder}/combined_7_snpeff.vcf.gz && \
	tabix -p vcf ${outfolder}/combined_7_snpeff.vcf.gz

#############################################
# FILTER AND CREATE CANDIDATE LIST
#############################################

# add hom/het/tot calls for each sample
# count variant in sample1 1:4 only, 5-7 are neutral
java ${javaopts} -jar $SNPEFF/SnpSift.jar caseControl \
	"++++000" ${outfolder}/combined_7_snpeff.vcf.gz | \
	bgzip -c > ${outfolder}/combined_7_snpeff_cnt.vcf.gz

# filter calls with the newly created INFO field
zcat ${outfolder}/combined_7_snpeff_cnt.vcf.gz | \
java ${javaopts} -jar $SNPEFF/SnpSift.jar filter \
	"(Cases[0] + Cases[1] == 1)  & \
	((exists ID) & ( ID =~ 'rs' ))" \
	> private_candidates.vcf

# attractive filters
#"( GEN[*].GQ > 60 )"
#" ( QUAL[*] >= 30 )"
#(DP >= 25)
#(countVariant() = 1) & (countRef() == 3)
#isVariant( GEN[0] )
#"(exists ID) & ( ID =~ 'rs' )"
#(( exists INDEL ) & (QUAL >= 20))
#ANN[*].EFFECT has 'missense_variant'
# & ( GEN[*].GQ > 60 )
#	( ANN[*].EFFECT has 'missense_variant' | ANN[*].EFFECT has 'synonymous_variant' )" \

# filter candidates from full VCF to be called in max 1 sample
#zcat ${outfolder}/combined_7_snpeff.vcf.gz | java ${javaopts} -jar $SNPEFF/SnpSift.jar filter \
#	"(countRef() == 3) & (countVariant() == 1)  & \
#	((exists ID) & ( ID =~ 'rs' ))" \
#	> private_candidates.vcf

# create excel table
java ${javaopts} -jar $SNPEFF/SnpSift.jar extractFields -s "," -e "." private_candidates.vcf \
  "CHROM" "POS" "ID" "REF" "ALT" "FILTER" \
  "DP" "MQ" "GEN[0].GT" "GEN[1].GT" "GEN[2].GT" "GEN[3].GT"  "GEN[4].GT" "GEN[5].GT" "GEN[6].GT" \
  "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
  > private_candidates.tsv

# cleanup leftovers
rm tmpfiles/*
