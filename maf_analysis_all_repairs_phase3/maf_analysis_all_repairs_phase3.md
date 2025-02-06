# MAF_ANALYSIS_ALL_REPAIRS_PHASE3
The scripts under this folder are used for 
- checking the mutation count distributions of 33 TCGA studies
- creating a list of all background genes for each gene in of interest
- creating pseudo-repair gene lists with found background genes
- calculating the mutation percentage of the random lists and all DNA repair genes from mutation data of MC3 MAF file with all mutations (indicated with the word "background")
- calculating the mutation percentage of the random lists and all DNA repair genes from mutation data of MC3 MAF file with pathogenic mutations (indicated with the word "pathogenic")
- comparing the null case (random lists) with observed value (genes of interest) using right-tailed z-score test for significance
- visualising the findings

These analyses on MC3 MAF file is run for a gene list containing all DNA repair genes: Global Genome NER genes (NER-GG), Transciption-Coupled NER (NER-TCR), NER genes common for NER-GG and NER-TCR processes (NER-Common), Base Excision Repair (BER), Crosslink Repair (CLR), Double Strand Break Repair (DSB) and Mismatch Repair (MMR).

The gene list used for preparing all background genes consists of other human protein coding genes taken from NCBI.

The analysis results of significantly mutated studies are found in [significants](https://github.com/ilaydakaytaran/DNArepair/tree/main/maf_analysis_all_repairs_phase3/significants) directory.

# FILE CONTENTS

## [input_files](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_all_repairs_phase3/input_files/)
- **preparing_all_genes_dfs.ipynb**: Python notebook for obtaining the SwissProt IDs of all human protein coding genes from NCBI and recording their protein length in UniProt. Genes without a corresponding SwissProt ID are noted as “NotFound”, and 55 genes with multiple SwissProt IDs are assigned to a single ID manually. 

- **df_all_most_23032024.txt**: dataframe of all human protein coding genes in [ncbi_all_gene_result_27122023.txt](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/input_files/ncbi_all_gene_result_27122023.txt) with two additional columns “swissprot_id” and “protein_length”. 1306 genes out of 20627 have “NotFound” in swissprot_id column.

- **genes_without_swissprot_id_found.txt**: contains 1306 genes that did not have any assigned SwissProt ID. These genes were excluded from the downstream analysis since no known protein length was found. 
 
## [analysis](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_all_repairs_phase3/analysis/)
- **background_analysis_all.py**: Python script for finding background genes for all DNA repair genes (NER-GG, NER-TCR, NER-COMMON, BER, CLR, DSB, MMR) from non-repair human protein coding genes with known protein length, creating random lists and calculating observed mutation rate in all repair genes and expected mutation rate in random lists for a single study. No filtering is done on MC3 MAF files to select pathogenic mutations. 

- **background_analysis.slurm**: Slurm file used for submitting parallel HPC jobs for 33 cancer studies where the MC3 MAF file is not filtered to select pathological mutations.

- **pathogenic_analysis_all.py**: Python script for finding background genes for all DNA repair genes (NER-GG, NER-TCR, NER-COMMON, BER, CLR, DSB, MMR) non-repair human protein coding genes with known protein length, creating random lists and calculating observed mutation rate in all repair genes and expected mutation rate in random lists for a single study. Mutation filtering is done on MC3 MAF files to select pathogenic mutations. 

- **pathogenic_analysis.slurm**: Slurm file used for submitting parallel HPC jobs for 33 cancer studies where the MC3 MAF file is filtered to select pathological mutations.

## [visualisation](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_all_repairs_phase3/visualisation/)
- **background_visualisation.py** / **pathogenic_visualisation.py**: Python scripts for visualising the graphs of observed vs expected mutation rates as a line graph and calculating the P-values & Z-scores for a single study for the analyses. 

- **background_visualisation.slurm** / **pathogenic_visualisation.slurm**: Slurm files used for submitting parallel HPC jobs for 33 cancer studies.

## [z-score](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_all_repairs_phase3/z_score/)
- **z_score_plotting.py**: Python script for plotting the Z-scores of the cancer studies with a significant mutation rate in all DNA repair genes as a bar plot.
  
## [significants](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/significants/)
- **\*** **significant_scores_all_background.txt**: Text files containing the study names that were found to have significant mutations in all repair genes from the analyses run on MC3 MAF file without pathological filtering. The file is in tab-separated format, where the first column is the study name, second column is p-value and the third column is Z-score.

- **\*** **significant_scores_all.txt**: Text files containing the cancer study names that were found to have significant mutations in all repair genes from the analyses run on MC3 MAF file with pathological filtering. The file is in tab-separated format, where the first column is the study name, second column is p-value and the third column is Z-score. 

- **z_scores_all_background.pdf**: The plots of significant study names vs. Z-scores from all DNA repair genes analyses run on MC3 MAF file without pathological filtering.

- **z_scores_all.pdf**: The plots of significant studies vs. Z-scores from all DNA repair genes analyses run on MC3 MAF file with pathological filtering.


**\*** An empty file indicates an analysis that found no studies with significant mutations in all DNA repair genes.
