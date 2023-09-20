The scripts under this folder is used for 
- creating a list of all background genes for each gene in of interest
- creating pseudo-NER gene lists with found background genes
- calculating the mutation percentage of the random lists
- calculating the mutation percentage of the genes of interest
- comparing the null case (random lists) with observed value (genes of interest) using right-tailed z-score test for significance
- visualising the findings 

### FILE CONTENTS

- [background_check.slurm](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/background_check.slurm): slurm file used for submitting parallel HPC jobs for 32 cancer studies.
- [background_check_gg.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/background_check_gg.py): python script for finding background genes for Nucleotide Excision Repair-GG genes, creating random lists and calculating observed mutation rate in NER genes and expected mutation rate in random lists.
- [background_check_common.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/background_check_common.py) / [background_check_tcr.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/background_check_tcr.py): python script for finding background genes for Nucleotide Excision Repair-Common/Nucleotide Excision Repair-TCR genes (with additional Gene Ontology term “positive regulation of multicellular organismal process” genes), eliminating short NER-common genes, creating random lists and calculating observed mutation rate in NER genes and expected mutation rate in random lists.
- [background_check_xpc.py](https://github.com/ilaydakaytaran/DNArepair/blob/main/background/background_check_xpc.py): python script for finding background genes for a single gene of interest, creating random lists and calculating observed mutation rate in the gene of interest and expected mutation rate in random lists.
