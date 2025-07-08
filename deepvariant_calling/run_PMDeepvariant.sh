#!/bin/bash

# call variants from long reads and a reference
# using PEPPER Margin Deepvariant
#
# Requirements:
# you will be asked your sudo password (no sudo no run!)
# run on a unix computer installed with
#   docker image: kishwars/pepper_deepvariant:r0.8
#   # get with: sudo docker pull kishwars/pepper_deepvariant:r0.8
#   samtools, minimap2,
#   reference genome fasta present and indexed
#   a minimum of 4 threads
#
# this script will write results in the current folder
#
# Stephane Plaisance (VIB-NC) 2022/06/01; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# version number
version="1.0, 2022_06_01"

usage='# Usage: run_PMDeepvariant.sh -f <fastq reads> -r <fasta reference>
# script version '${version}'
# [optional: -n <sample name|sample>]
# [optional: -t <threads|4>]
# [optional: -X <extra parameter string between quotes for the docker command>]'

while getopts "f:r:n:t:X:h" opt; do
  case $opt in
    f) reads=${OPTARG} ;;
    r) ref=${OPTARG} ;;
    n) opt_name=${OPTARG} ;;
    t) opt_thr=${OPTARG} ;;
    X) opt_extra=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2;
       exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; 
       exit 1 ;;
  esac
done

# set defaults
nthr=${opt_thr:-4}
stthr=4
extra=${opt_extra:-""}
pfx=${opt_name:-"sample"}

# check if all dependencies are present
declare -a arr=("minimap2" "samtools")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || \
  ( echo "# required ${prog} not found in PATH"
exit 1 )
done

# check for the docker image
res=$(sudo docker images | grep kishwars/pepper_deepvariant | cut -d " " -f 1)
if [ ! ${res} == "kishwars/pepper_deepvariant" ]; then
  (echo "docker image not found"; exit 1)
fi

# test if minimal arguments were provided
if [ -z "${reads}" ]; then
   echo "# no long reads provided!"
   echo "${usage}"
   exit 1
fi

if [ ! -f "${reads}" ]; then
    echo "${reads} file not found!";
    exit 1
fi

if [ -z "${ref}" ]
then
	echo "# no fasta reference provided"
	echo "${usage}"
	exit 1
fi

if [ ! -f "${ref}" ]; then
    echo "${ref} file not found!";
    exit 1
fi

# other parameters or defaults
workdir=$PWD

# create folders for results
mkdir -p input output/intermediate_files output/logs

#############################################
# map reads to the reference using minimap2 #
#############################################

reffile=$(basename ${ref})
cp ${ref} input/${reffile} && \
  samtools faidx input/${reffile}

# build the command
cmd="minimap2 \
  -a \
  -t ${nthr} \
  -x map-hifi \
  -R '@RG\tID:'${pfx}'\tSM:'${pfx} \
  --MD \
  input/${reffile} \
  ${reads} \
  > output/${pfx}_mappings.sam && \
    samtools sort -@${stthr} \
      output/${pfx}_mappings.sam \
      -O BAM \
      -o output/${pfx}_sorted.bam && \
        samtools index output/${pfx}_sorted.bam && \
          rm output/${pfx}_mappings.sam"

# show and execute
echo "# mapping reads to reference"
echo "# ${cmd}"
time eval ${cmd}

#################################
# run pepper margin deepvariant #
#################################

input_dir="${workdir}/input"
output_dir="${workdir}/output"
out_pfx="${pfx}_PEPPER_Margin_DeepVariant"

param_haplotag="allParams.haplotag.pb-hifi.hapDup.json"
param_phase="allParams.phase_vcf.pb-hifi.json"

# Run PEPPER-Margin-DeepVariant
echo "# running pepper_margin_deepvariant docker"
time sudo docker run \
  -u "$(id -u):$(id -g)" \
  -v "${input_dir}":"/input" \
  -v "${output_dir}":"/output" \
  kishwars/pepper_deepvariant:r0.8 \
    run_pepper_margin_deepvariant call_variant \
    -b "/output/${pfx}_sorted.bam" \
    -f "/input/${reffile}" \
    -o "/output" \
    -p "${out_pfx}"_q \
    -t "${nthr}" \
    -s ${pfx} \
    --keep_intermediate_bam_files \
    --sample_name "${pfx}" \
    --gvcf \
    --phased_output \
    --margin_haplotag_model "/opt/margin_dir/params/phase/${param_haplotag}" \
    --margin_phase_model "/opt/margin_dir/params/phase/${param_phase}" \
    --hifi \
    ${extra}

exit 0
