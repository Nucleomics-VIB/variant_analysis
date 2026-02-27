#!/bin/bash
#SBATCH --job-name=SMI-INFO
#SBATCH --partition=gpu_l40s_64C_128T_1TB
#SBATCH --gres=gpu:4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=s04
#SBATCH --time=0-1:00
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

nvidia-smi
