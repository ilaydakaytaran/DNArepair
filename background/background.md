# BACKGROUND.md
The scripts under this folder are used for 
- creating a list of all background genes for each gene in of interest
- creating pseudo-NER gene lists with found background genes
- calculating the mutation percentage of the random lists
- calculating the mutation percentage of the genes of interest
- comparing the null case (random lists) with observed value (genes of interest) using right-tailed z-score test for significance
- visualising the findings 

# FILE CONTENTS

## input_files
- [df_all_repair_ner.txt](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/input_files/df_all_repair_ner.txt): dataframe of all repair genes (Nucleotide Excision Repair, Base Excision Repair, Crosslink Repair, Mismatch Repair, Double Strand Break Repair) in tab-separated format, taken from corresponding Gene Ontology terms.
- [df_rna_ner.txt](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/input_files/df_rna_ner.txt): dataframe of additional background genes in tab-separated format, taken from the Gene Ontology term "positive regulation of multicellular organismal process".
- [ncbi_all_gene_result.txt](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/input_files/ncbi_all_gene_result.txt): dataframe of all human protein-coding genes in tab-separated format, taken from NCBI Genes database with the search *' "9606"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties] '*


## check
- [background_check.slurm](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/check/background_check.slurm): slurm file used for submitting parallel HPC jobs for 32 cancer studies.
- [background_check_gg.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/check/background_check_gg.py): python script for finding background genes for Nucleotide Excision Repair-GG genes, creating random lists and calculating observed mutation rate in NER genes and expected mutation rate in random lists for a single study.
- [background_check_common.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/check/background_check_common.py) / [background_check_tcr.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/check/background_check_tcr.py): python script for finding background genes for Nucleotide Excision Repair-Common/Nucleotide Excision Repair-TCR genes (with additional Gene Ontology term “positive regulation of multicellular organismal process” genes), eliminating short NER-common genes, creating random lists and calculating observed mutation rate in NER genes and expected mutation rate in random lists for a single study.
- [background_check_xpc.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/check/background_check_xpc.py): python script for finding background genes for a single gene of interest, creating random lists and calculating observed mutation rate in the gene of interest and expected mutation rate in random lists for a single study.


## visualisation
- [background_visualisation.slurm](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/visualisation/background_visualisation.slurm): slurm file used for submitting parallel HPC jobs for 32 cancer studies.
- [background_visualisation.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/visualisation/background_visualisation.py): python script for visualising the graphs of observed vs expected mutation rates as a line graph for a single study.


## z_score
- [z_score_plotting.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/z_score/z_score_plotting.py): python script for plotting the z-scores of the cancer studies with a significant mutation rate in NER genes as a bar plot. 
