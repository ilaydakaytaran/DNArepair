#!/bin/bash

#SBATCH -J BackNoiseXPC

#SBATCH --account=investor

#SBATCH --ntasks-per-node=1

#SBATCH --qos=mid_investor

#SBATCH --mail-type=ALL

#SBATCH --partition=mid_investor

#SBATCH --mail-user=ilaydakaytaran@sabanciuniv.edu

#SBATCH --output=/cta/users/ilaydakaytaran/backgrounds/outputs_and_errors/slurm-%j.out 

#SBATCH --error=/cta/users/ilaydakaytaran/backgrounds/outputs_and_errors/slurm-%j.err

#SBATCH --array=1-32

module load python-3.8.0-gcc-9.2.0-gjpyy5a 


export TMPDIR=/cta/users/ilaydakaytaran/temp

study_ids=("laml_tcga_pan_can_atlas_2018" "acc_tcga_pan_can_atlas_2018" "blca_tcga_pan_can_atlas_2018" "brca_tcga_pan_can_atlas_2018"
           "cesc_tcga_pan_can_atlas_2018" "coadread_tcga_pan_can_atlas_2018" "chol_tcga_pan_can_atlas_2018" "dlbc_tcga_pan_can_atlas_2018"
		   "esca_tcga_pan_can_atlas_2018" "gbm_tcga_pan_can_atlas_2018" "hnsc_tcga_pan_can_atlas_2018" "kich_tcga_pan_can_atlas_2018" 
		   "kirc_tcga_pan_can_atlas_2018" "kirp_tcga_pan_can_atlas_2018" "lgg_tcga_pan_can_atlas_2018" "lihc_tcga_pan_can_atlas_2018"
		   "luad_tcga_pan_can_atlas_2018" "lusc_tcga_pan_can_atlas_2018" "meso_tcga_pan_can_atlas_2018" "ov_tcga_pan_can_atlas_2018"
		   "paad_tcga_pan_can_atlas_2018" "pcpg_tcga_pan_can_atlas_2018" "prad_tcga_pan_can_atlas_2018" "sarc_tcga_pan_can_atlas_2018"
		   "skcm_tcga_pan_can_atlas_2018" "stad_tcga_pan_can_atlas_2018" "tgct_tcga_pan_can_atlas_2018" "thca_tcga_pan_can_atlas_2018"
		   "thym_tcga_pan_can_atlas_2018" "ucec_tcga_pan_can_atlas_2018" "ucs_tcga_pan_can_atlas_2018" "uvm_tcga_pan_can_atlas_2018")
knees_indices=("5" "6" "11" "1" "8" "12" "1" "5" "7" "11" "9" "4" "10" "4" "16" "8" "19" "14" "5" "3" "8" "4" "5" "3" "20" "10" "4" "7" "5" "14" "5" "6")

current_id=${study_ids[$SLURM_ARRAY_TASK_ID - 1]}
current_knee=${knees_indices[$SLURM_ARRAY_TASK_ID - 1]} 
python /cta/users/ilaydakaytaran/backgrounds/background_check_xpc.py $current_id $current_knee





