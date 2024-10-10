#!/bin/bash

#SBATCH -J SingleALLpatho

#SBATCH --account=investor

#SBATCH --ntasks-per-node=1

#SBATCH --qos=mid_investor

#SBATCH --mail-type=ALL

#SBATCH --partition=mid_investor

#SBATCH --mail-user=ilaydakaytaran@sabanciuniv.edu

#SBATCH --output=/cta/users/ilaydakaytaran/maf_analysis_2/outputs_and_errors/slurm-%j.out 

#SBATCH --error=/cta/users/ilaydakaytaran/maf_analysis_2/outputs_and_errors/slurm-%j.err

#SBATCH --array=1-2

module load python-3.8.0-gcc-9.2.0-gjpyy5a 

export TMPDIR=/cta/users/ilaydakaytaran/temp

study_ids=("acc_tcga_pan_can_atlas_2018" "ov_tcga_pan_can_atlas_2018")

current_id=${study_ids[$SLURM_ARRAY_TASK_ID - 1]}

python /cta/users/ilaydakaytaran/maf_analysis_2/singlegene_analysis/singlegene_pathogenic_analysis_all.py $current_id 