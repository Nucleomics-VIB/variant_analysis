#!/bin/bash

#SBATCH --job-name="parabricks_deepvariant_germline"
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=4
#SBATCH --time=0-04:00:00
#SBATCH --mem=512G
#SBATCH --account=s04
#SBATCH --gres=gpu:4
#SBATCH --partition=gpu_l40s_64C_128T_1TB

# Aim: run the NVidia clara parabricks deepvariant_germline (v4.4.0) on a single read pair
# script: slurm_deepvariant_germline.sh
#   run with 'sbatch scripts/slurm_deepvariant_germline.sh reads/SRR29676022_1.fq.gz'
# author: SP@NC; 2025-01-17 v1.0

# Check if an argument is provided
if [ $# -eq 0 ]; then
    echo "Error: Please provide the path to the reads_1.fq.gz file. \
          Usage: sbatch scripts/slurm_deepvariant_germline.sh reads/<read1.fq.gz>"
    exit 1
fi

# set the argument to the path to read_1.fq.gz provided by the user
fq=${1}

# set WORKDIR
# obtained from PWD inside the interactive session WORKDIR="/home/VIB.LOCAL/stephane.plaisance/project_data"
export WORKDIR="/data/projects/s04/wgs_variant_analysis"

# set pfx
export pfx=$(basename "${fq%_1.fq.gz}")

# output folder for mappings
outbam="mappings"
mkdir -p "${outbam}"

# output folder for variants
outvcf="variants_XYhap"
mkdir -p "${outvcf}"

# output folder for logs
outlogs="logs_XYhap"
mkdir -p "${outlogs}"

#######################
# additional arguments
#######################

singimg="${WORKDIR}/bin/clara-parabricks_4.4.0-1.sif"

refidx="mRatBN7.2.fa"
knownsites="rattus_norvegicus.vcf.gz"
optdist=2500
platform="Illumina"

# deduce fq1
fq1=$(basename "${fq}")

# deduce fq2
fq2=${fq1/_1.fq.gz/_2.fq.gz}

# define read-group tag
RGTAG="@RG\tID:${pfx}\tLB:lib_${pfx}\tPL:Aviti\tSM:${pfx}\tPU:${pfx}"

# run with singularity

# use gpu
numgpu=4

# Your job commands go here
echo "Starting slurm_deepvariant_germline.sh job execution..." > ${WORKDIR}/${pfx}_slurm_resources.txt
echo $(date) >> ${WORKDIR}/${pfx}_slurm_resources.txt

singularity run --nv \
    --bind "${WORKDIR}:/workdir" \
    --bind "${TMPDIR}:/tmp" \
    --pwd "/workdir" \
    ${singimg} \
    pbrun deepvariant_germline \
    --gpusort \
    --gpuwrite \
    --num-gpus ${numgpu} \
    --bwa-options '-M -K 10000000' \
    --optical-duplicate-pixel-distance "${optdist}" \
    --read-group-sm "${pfx}" \
    --read-group-lb "lib_${pfx}" \
    --read-group-pl "${platform}" \
    --read-group-id-prefix "${pfx}" \
    --ref "/workdir/bwaidx/${refidx}" \
    --in-fq "/workdir/reads/${fq1}" "/workdir/reads/${fq2}" \
    --knownSites "/workdir/reference/${knownsites}" \
    --haploid-contigs X,Y \
    --out-bam "/workdir/${outbam}/${pfx}_mrkdup.bam" \
    --out-recal-file "/workdir/${outbam}/${pfx}_recal.txt" \
    --out-duplicate-metrics "/workdir/${outbam}/${pfx}_duplicate_metrics" \
    --out-variants "/workdir/${outvcf}/${pfx}_XYhap.g.vcf.gz" \
    --gvcf \
    --tmp-dir "/tmp"

# cleanup

# Rename output and error files to include ${pfx}
mv ${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out ${pfx}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out
mv ${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err ${pfx}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err

$(date) >> ${WORKDIR}/${pfx}_slurm_resources.txt
# Print resource usage for logging purposes
sacct -j $SLURM_JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,CPUTime,TotalCPU,ReqCPUS,AllocCPUS,AveRSS,MaxRSS,ReqMem,AveVMSize,MaxVMSize \
      >> ${WORKDIR}/${pfx}_slurm_resources.txt
echo "Cleanup complete. Exiting."

# terminate session
exit 0
