#!/bin/bash

#SBATCH --job-name=PARABRICKS-HELP
#SBATCH --partition=gpu_l40s_64C_128T_1TB
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH -–ntasks=1
#SBATCH -–cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=s04
#SBATCH --time=00:03:00

# singularity image
img='/data/projects/s04/wgs_variant_analysis/bin/clara-parabricks_4.4.0-1.sif'

singularity exec --nv ${img}

exit 0

