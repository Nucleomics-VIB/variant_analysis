# http://samtools.github.io/bcftools/howtos/variant-calling.html

reference_fa="/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa"
dbsnp="/data/biodata/references/galGal6a/ens_dbSNP_150/GRCg6a_dbSNP_150.vcf.gz"

# call from multiple bam and 1 chr
chr=1
bcftools mpileup -Ou \
  -r ${chr} \
  -f ${reference_fa} \
  bam_recalibration/524_mrkdup_srt_recal.bam \
  bam_recalibration/526_mrkdup_srt_recal.bam \
  bam_recalibration/528_mrkdup_srt_recal.bam \
  bam_recalibration/530_mrkdup_srt_recal.bam \
  | bcftools call -Ov -mv -o bcftools/bcftools_${chr}_multi_4.vcf

reference_fa="/data/biodata/references/galGal6a/Gallus_gallus.GRCg6a.dna.toplevel.fa"
while read chr; do
(bcftools mpileup -Ou \
  -r ${chr} \
  -f ${reference_fa} \
  bam_recalibration/524_mrkdup_srt_recal.bam \
  bam_recalibration/526_mrkdup_srt_recal.bam \
  bam_recalibration/528_mrkdup_srt_recal.bam \
  bam_recalibration/530_mrkdup_srt_recal.bam \
  bam_recalibration/W201120E14_mrkdup_srt_recal.bam \
  bam_recalibration/W291020E14_mrkdup_srt_recal.bam \
  bam_recalibration/W291020E15_mrkdup_srt_recal.bam \
  | bcftools call -Ov -mv -o bcftools/bcftools_${chr}_multi_7.vcf)&
done < bcftools/main_chr.txt

# normalize calls
bcft=4
vcf=bcftools/bcftools_${chr}_multi_4.vcf
vcfout=bcftools/$(basename ${vcf%.vcf}_norm.vcf.gz)
bcftools norm \
  --threads ${bcft} \
  -f ${reference_fa} \
  ${vcf} | \
  java ${javaopts} -jar $SNPEFF/SnpSift.jar filter "countVariant()>0" \
  | bgzip -c > ${vcfout} && tabix -p vcf ${vcfout}

# filter calls
bcftools filter -sLowQual -e'%QUAL<10' -g3 -G10 -O z -o ${vcfout%.vcf.gz}_filt.vcf.gz ${vcfout} \
  && tabix -p vcf ${vcfout%.vcf.gz}_filt.vcf.gz

# annotate with SnpEff
infile=${vcfout%.vcf.gz}_filt.vcf.gz
outfile=${infile%.vcf.gz}_snpeff.vcf.gz
build="GRCg6a.105"

java ${javaopts} -jar $SNPEFF/SnpSift.jar annotate \
  -id ${dbsnp} \
  ${infile} | \
  java ${javaopts} -jar $SNPEFF/SnpSift.jar caseControl "++++" | \
  java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfile%.vcf.gz}_summary.html \
	-nodownload \
	${build} - | \
	bgzip -c > ${outfile} && \
	tabix -p vcf ${outfile}

# candidates 524
zcat ${outfile} | \
java ${javaopts} -jar $SNPEFF/SnpSift.jar filter \
  "(( na FILTER ) | (FILTER = 'PASS')) & \
   ((exists ID) & ( ID =~ 'rs' )) & \
   (Cases[0] + Cases[1] == 1) & \
   isVariant( GEN[0] )" \
  | bgzip -c > ${outfile%.vcf.gz}_524.vcf.gz \
  && tabix -p vcf ${outfile%.vcf.gz}_524.vcf.gz
