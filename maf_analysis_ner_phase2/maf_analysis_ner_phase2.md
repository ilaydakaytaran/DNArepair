# MAF_ANALYSIS_NER_PHASE2
The scripts under this folder are used for 
- checking the mutation count distributions of 33 TCGA studies
- creating a list of all background genes for each gene in of interest
- creating pseudo-NER gene lists with found background genes
- calculating the mutation percentage of the random lists and the genes of interest from mutation data of MC3 MAF file with all mutations (indicated with the word "background")
- calculating the mutation percentage of the random lists and the genes of interest from mutation data of MC3 MAF file with pathogenic mutations (indicated with the word "pathogenic")
- comparing the null case (random lists) with observed value (genes of interest) using right-tailed z-score test for significance
- visualising the findings
- running the pathogenic analysis for each gene separately, for the studies which were found to have significant pathogenic mutations in analysed gene lists

These analyses on MC3 MAF file is run for 4 gene lists related to Nucleotide Excision Repair (NER): All NER genes, Global Genome NER genes (NER-GG), Transciption-Coupled NER (NER-TCR) and NER genes common for NER-GG and NER-TCR processes (NER-Common). In file names below, the respective identifiers for these lists are _all_, _gg_, _tcr_ and _common_ respectively (in place of the keyword _genelistname_).

The gene list used for preparing all background genes consists of other DNA damage repair genes, specifically Base Excision Repair (BER), Crosslink Repair (CLR), Double Strand Break Repair (DSB) and Mismatch Repair (MMR). Additional genes are selected from Gene Ontology term “positive regulation of multicellular organismal process” genes in the cases that no suitable genes were found from non-NER genes.

The analysis results of significantly mutated studies are found in [significants](https://github.com/ilaydakaytaran/DNArepair/tree/main/maf_analysis_ner_phase2/significants) directory.

# FILE CONTENTS

## [input_files](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/input_files/)
- **mc3_manifest.txt**: Manifest file to get the MC3 MAF file with GDC client. Using this file with GDC Client, the public MC3 MAF file was downloaded, containing mutation data for the samples of 33 TCGA studies. Three new columns were created to store Patient IDs, Sample IDs and Study IDs. Additionally, this public data has Entrez Gene ID column empty for all lines which is the main information we are using to identify genes. This column was filled by combining with Gene Symbol – Entrez ID pairs taken from NCBI.

- **genes_without_id.txt**: During the process of filling up the Entrez Gene ID column, 2633 symbols from the MC3 MAF file were unable to be matched to any ID. This file contains those genes with no corresponding gene ID which were eliminated from the analyses.

- **ncbi_all_gene_result_27122023.txt**: Dataframe of all human protein-coding genes in tab-separated format (redownloaded on 27 December 2023) taken from NCBI Genes database with the search *' "9606"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties] '*

## [mutation_counts](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/mutation_counts/)
- **mutation_count_analysis.py**: Python script for plotting the box plots of the mutation counts of each study from the MC3 MAF file.

- **mutation_counts_ordered.pdf** / **mutation_counts_ordered_log2.pdf**: The boxplots showing the mutation count distributions of 33 TCGA studies in increasing mean mutation counts. 

## [analysis](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/analysis/)
- **background_analysis_*genelistname*.py**: Python script for finding background genes for the respective gene list genes from non-NER repair genes, creating random lists and calculating observed mutation rate in the respective gene list and expected mutation rate in random lists for a single study. No filtering is done on MC3 MAF files to select pathogenic mutations. For gene lists All NER, NER-TCR and NER-Common, additional Gene Ontology term “positive regulation of multicellular organismal process” genes are added to the pool of background genes to be able to represent some shorter genes in these lists.

- **background_analysis.slurm**: Slurm file used for submitting parallel HPC jobs for 33 cancer studies where the MC3 MAF file is not filtered to select pathological mutations.

- **pathogenic_analysis_*genelistname*.py**: Python script for finding background genes for the respective gene list genes from non-NER repair genes, creating random lists and calculating observed mutation rate in the respective gene list and expected mutation rate in random lists for a single study. Mutation filtering is done on MC3 MAF files to select pathogenic mutations. For gene lists All NER, NER-TCR and NER-Common, additional Gene Ontology term “positive regulation of multicellular organismal process” genes are added to the pool of background genes to be able to represent some shorter genes in these lists.

- **pathogenic_analysis.slurm**: Slurm file used for submitting parallel HPC jobs for 33 cancer studies where the MC3 MAF file is filtered to select pathological mutations.

## [visualisation](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/visualisation/)
- **\*** **background_visualisation.py** / **\*** **pathogenic_visualisation.py**: Python scripts for visualising the graphs of observed vs expected mutation rates as a line graph and calculating the P-values & Z-scores for a single study for the analyses. 

- **background_visualisation.slurm** / **pathogenic_visualisation.slurm**: Slurm files used for submitting parallel HPC jobs for 33 cancer studies.

## [z-score](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/z_score/)
- **\*** **z_score_plotting.py**: Python script for plotting the Z-scores of the cancer studies with a significant mutation rate in a gene list as a bar plot.

## [significants](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/significants/)
- **\*\*** **significant_scores_*genelistname*_background.txt**: Text files containing the study names that were found to have significant mutations in the respective gene list genes from the analyses run on MC3 MAF file without pathological filtering. The file is in tab-separated format, where the first column is the study name, second column is p-value and the third column is Z-score.

- **\*\*** **significant_scores_*genelistname*.txt**: Text files containing the cancer study names that were found to have significant mutations in the respective gene list genes from the analyses run on MC3 MAF file with pathological filtering. The file is in tab-separated format, where the first column is the study name, second column is p-value and the third column is Z-score. 

- **z_scores_*genelistname*_background.pdf**: The plots of significant study names vs. Z-scores from the respective gene list analyses run on MC3 MAF file without pathological filtering.

- **z_scores_*genelistname*.pdf**: The plots of significant studies vs. Z-scores from the respective gene list analyses run on MC3 MAF file with pathological filtering.

## [singlegene_analysis](https://github.com/ilaydakaytaran/DNArepair/blob/main/maf_analysis_ner_phase2/singlegene_analysis/)
- **singlegene_pathogenic_analysis_*genelistname*.py**: Python script for finding background genes for each gene from the respective gene list (separately) from non-NER repair genes, creating random lists (consisting of a single gene each) and calculating observed mutation rate in the respective NER gene and expected mutation rate in random lists for a single study. Mutation filtering is done on MC3 MAF files to select pathogenic mutations. For gene lists All NER, NER-TCR and NER-Common, additional Gene Ontology term “positive regulation of multicellular organismal process” genes are added to the pool of background genes to be able to represent some shorter genes in these lists. This analysis is done only for studies that were found to have significant pathogenic mutations in NER genes from the pathogenic analysis of NER-GG genes.
  
  - Studies in the ALL NER single gene analysis: Bladder Urothelial Carcinoma, Cholangiocarcinoma, Esophageal Adenocarcinoma, Liver Hepatocellular Carcinoma and Ovarian Serous Cystadenocarcinoma.
  - Studies in the NER-Common single gene analysis: Adrenocortical Carcinoma and Ovarian Serous Cystadenocarcinoma.
  - Studies in the NER-GG single gene analysis: Bladder Urothelial Carcinoma, Cholangiocarcinoma, Esophageal Adenocarcinoma, Liver Hepatocellular Carcinoma and Ovarian Serous Cystadenocarcinoma).
  - Studies in the NER-TCR single gene analysis: Acute Myeloid Leukemia, Adrenocortical Carcinoma, Mesothelioma and Uterine Carcinosarcoma.
 
- **singlegene_pathogenic_analysis_*genelistname*.slurm**: Slurm files for submitting parallel HPC jobs with the respective cancer studies for the respective single gene analysis.

- **singlegene_pathogenic_visualisation_*genelistname*.py**: Python scripts for calculating the P-values & Z-scores for respective cancer studies of the respective single gene analysis and plotting scatter plots of genes with significant pathogenic mutations. The calculated values and the obtained plots are found in respective directories ([all](https://github.com/ilaydakaytaran/DNArepair/tree/main/maf_analysis_ner_phase2/singlegene_analysis/all) / [common](https://github.com/ilaydakaytaran/DNArepair/tree/main/maf_analysis_ner_phase2/singlegene_analysis/common) / [gg](https://github.com/ilaydakaytaran/DNArepair/tree/main/maf_analysis_ner_phase2/singlegene_analysis/gg) / [tcr](https://github.com/ilaydakaytaran/DNArepair/tree/main/maf_analysis_ner_phase2/singlegene_analysis/tcr)) under this directory. The file formats under these directories are as follows:
  
	- **p_values_*genelistname*.txt** / **z_scores_*genelistname*.txt**: Text files containing the calculated p-values/Z_scores from the respective single gene analysis. This file is in tab-separated format, where the first column contains the study names the analysis was run on, and the rest of the columns have the gene symbols of the gene list used in the analysis. Having a value of “NaN” indicates a gene whose standard deviation was calculated to be below 1e-15 (a standard deviation of 0 in practice), therefore no p-value/Z-score was calculated. 
	- **significant_genes_per_study_*genelistname*.pdf**: PDF files containing the scatter plots obtained from respective single gene analyses. x-axis has columns for study names the analysis was run on, while y-axis has -log10 transformation of the p-values calculated for genes from the gene list used in the analysis. Each scatter dot represents a gene from the gene list, where a redder dot indicates a gene with a lower p-value, therefore higher significance. Genes with non-significant pathogenic mutations are represented in grey. 

**\*** The file names and the paths in the script need to be changed according to the gene list name used (all, common, gg or tcr).

**\*\*** An empty file indicates an analysis that found no studies with significant mutations in the respective gene list.
