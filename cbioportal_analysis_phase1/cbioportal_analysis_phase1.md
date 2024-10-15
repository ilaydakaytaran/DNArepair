# CBIOPORTAL_ANALYSIS_PHASE1
The scripts under this folder are used for 
- creating a list of all background genes for each gene in of interest
- creating pseudo-NER gene lists with found background genes
- calculating the mutation percentage of the random lists
- calculating the mutation percentage of the genes of interest
- comparing the null case (random lists) with observed value (genes of interest) using right-tailed Z-score test for significance
- visualising the findings 

These analyses on mutation data of [TCGA PanCancer Atlas studies from cBioPortal](https://www.cbioportal.org/study/summary?id=laml_tcga_pan_can_atlas_2018%2Cacc_tcga_pan_can_atlas_2018%2Cblca_tcga_pan_can_atlas_2018%2Clgg_tcga_pan_can_atlas_2018%2Cbrca_tcga_pan_can_atlas_2018%2Ccesc_tcga_pan_can_atlas_2018%2Cchol_tcga_pan_can_atlas_2018%2Ccoadread_tcga_pan_can_atlas_2018%2Cdlbc_tcga_pan_can_atlas_2018%2Cesca_tcga_pan_can_atlas_2018%2Cgbm_tcga_pan_can_atlas_2018%2Chnsc_tcga_pan_can_atlas_2018%2Ckich_tcga_pan_can_atlas_2018%2Ckirc_tcga_pan_can_atlas_2018%2Ckirp_tcga_pan_can_atlas_2018%2Clihc_tcga_pan_can_atlas_2018%2Cluad_tcga_pan_can_atlas_2018%2Clusc_tcga_pan_can_atlas_2018%2Cmeso_tcga_pan_can_atlas_2018%2Cov_tcga_pan_can_atlas_2018%2Cpaad_tcga_pan_can_atlas_2018%2Cpcpg_tcga_pan_can_atlas_2018%2Cprad_tcga_pan_can_atlas_2018%2Csarc_tcga_pan_can_atlas_2018%2Cskcm_tcga_pan_can_atlas_2018%2Cstad_tcga_pan_can_atlas_2018%2Ctgct_tcga_pan_can_atlas_2018%2Cthym_tcga_pan_can_atlas_2018%2Cthca_tcga_pan_can_atlas_2018%2Cucs_tcga_pan_can_atlas_2018%2Cucec_tcga_pan_can_atlas_2018%2Cuvm_tcga_pan_can_atlas_2018) is run for 3 gene lists related to Nucleotide Excision Repair (NER): Global Genome NER genes (NER-GG), Transciption-Coupled NER (NER-TCR) and NER genes common for NER-GG and NER-TCR processes (NER-Common). In file names below, the respective identifiers for these lists are _gg_, _tcr_ and _common_ respectively (in place of the keyword _genelistname_).

The gene list used for preparing all background genes consists of other DNA damage repair genes, specifically Base Excision Repair (BER), Crosslink Repair (CLR), Double Strand Break Repair (DSB) and Mismatch Repair (MMR). Additional genes are selected from Gene Ontology term “positive regulation of multicellular organismal process” genes in the cases that no suitable genes were found from non-NER repair genes.

The analysis results of significantly mutated studies are found in [significants](https://github.com/ilaydakaytaran/DNArepair/tree/main/cbioportal_analysis_phase1/significants) directory.

# FILE CONTENTS

## [input_files](https://github.com/ilaydakaytaran/DNArepair/blob/main/cbioportal_analysis_phase1/input_files/)
- **df_all_repair_ner.txt**: Dataframe of all repair genes (Nucleotide Excision Repair, Base Excision Repair, Crosslink Repair, Mismatch Repair, Double Strand Break Repair) in tab-separated format, taken from corresponding Gene Ontology terms.
- **df_rna_ner.txt**: Dataframe of additional background genes in tab-separated format, taken from the Gene Ontology term "positive regulation of multicellular organismal process".
- **ncbi_all_gene_result.txt**: Dataframe of all human protein-coding genes in tab-separated format, taken from NCBI Genes database with the search *' "9606"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties] '*


## [check](https://github.com/ilaydakaytaran/DNArepair/blob/main/cbioportal_analysis_phase1/check/)
- **background_check.slurm**: Slurm file used for submitting parallel HPC jobs for 32 cancer studies.
- **background_check_*genelistname*.py**: Python script for finding background genes for the respective gene list genes from non-NER repair genes, creating random lists and calculating observed mutation rate in NER genes and expected mutation rate in random lists for a single study. For gene lists NER-TCR and NER-Common, additional Gene Ontology term “positive regulation of multicellular organismal process” genes are added to the pool of background genes to be able to represent some shorter genes in these lists.
- **background_check_ercc6.py**: Python script for finding background genes for a single gene of interest (ERCC6 in the script) from non-NER repair genes, creating random lists and calculating observed mutation rate in the gene of interest and expected mutation rate in random lists for a single study. Note that unlike all the other background gene finding scripts with a range criteria of +/- 5% of the NER protein length, this one uses a range criteria of +/- 10% during the selection of potential background genes for an NER gene. This alteration was done to obtain enough background genes during the analysis of a single gene (ERCC6).


## [visualisation](https://github.com/ilaydakaytaran/DNArepair/blob/main/cbioportal_analysis_phase1/visualisation/)
- **background_visualisation.slurm**: Slurm file used for submitting parallel HPC jobs for 32 cancer studies.
- **\*** **background_visualisation.py**: Python script for visualising the graphs of observed vs expected mutation rates as a line graph for a single study.


## [z_score](https://github.com/ilaydakaytaran/DNArepair/blob/main/cbioportal_analysis_phase1/z_score/)
- **\*** **z_score_plotting.py**: Python script for plotting the Z-scores of the cancer studies with a significant mutation rate in NER genes as a bar plot. 

## [significants](https://github.com/ilaydakaytaran/DNArepair/blob/main/cbioportal_analysis_phase1/significants/)
- **\*\*** **significant_z_scores_*genelistname*.txt**: Text files containing the study names that were found to have significant mutations in the respective gene list genes from the analyses run on mutation data from cBioPortal. The file is in tab-separated format, where the first column is the study name and second column is Z-score.

- **\*\*** **significant_z_scores_ercc6.txt**: Text file containing the study names that were found to have significant mutations in ERCC6 from the analyses run on mutation data from cBioPortal. The file is in tab-separated format, where the first column is the study name and second column is Z-score.

- **z_scores_*genelistname*.pdf**: The plots of significant study names vs. Z-scores from the respective gene list analyses run on mutation data from cBioPortal.

- **z_scores_ercc6.pdf**: The plot of significant studies vs. Z-scores from the ERCC6 analyses run on mutation data from cBioPortal.

**\*** The file names and the paths in the script need to be changed according to the gene list name used (common, gg, tcr or a single gene symbol like ercc6).

**\*\*** An empty file indicates an analysis that found no studies with significant mutations in the respective gene list.
