#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2022-03-28

# http://samtools.github.io/bcftools/howtos/variant-calling.html

workdir=$PWD
cd ${workdir}

outfolder="bcftools_results"
mkdir -p ${outfolder}

###########
# variables
###########

reference_fa="/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa"
dbsnp="/data/biodata/references/galGal6a/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz"
chrlist="/data/biodata/references/galGal6a/main_chr.txt"
echo $(seq 1 28; seq 30 33) MT W Z | tr " " "\n" > ${chrlist}
snpeff_build="GRCg6a.105"

javaopts="-Xms24g -Xmx24g"
bcft=4

#######################################
# call all main chromosomes in parallel
#######################################

while read chr; do

vcfout=bcftools_${chr}_multi.vcf.gz

(bcftools mpileup -Ou \
  -r ${chr} \
  -f ${reference_fa} \
  $(ls bwa_mappings/*.bam) \
  | bcftools call -Ov -mv \
  | bcftools norm -f ${reference_fa} - \
  | java ${javaopts} -jar $SNPEFF/SnpSift.jar filter "countVariant()>0" \
  | bcftools filter -sLowQual -e'%QUAL<10' -g3 -G10 -O z \
    -o ${outfolder}/${vcfout%.vcf.gz}_filt.vcf.gz - \
    && tabix -p vcf ${outfolder}/${vcfout%.vcf.gz}_filt.vcf.gz) &
  
done < ${chrlist}

wait

echo "# Please wait until all chromosomes are processed before running the next script"
