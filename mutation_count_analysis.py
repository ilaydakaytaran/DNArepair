import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

study_ids = {
    "laml_tcga_pan_can_atlas_2018": "Acute Myeloid Leukemia", 
    "acc_tcga_pan_can_atlas_2018": "Adrenocortical Carcinoma", 
    "blca_tcga_pan_can_atlas_2018": "Bladder Urothelial Carcinoma",
    "brca_tcga_pan_can_atlas_2018": "Breast Invasive Carcinoma", 
    "cesc_tcga_pan_can_atlas_2018": "Cervical Squamous Cell Carcinoma", 
    "coadread_tcga_pan_can_atlas_2018": "Colorectal Adenocarcinoma",
    "chol_tcga_pan_can_atlas_2018": "Cholangiocarcinoma", 
    "dlbc_tcga_pan_can_atlas_2018": "Diffuse Large B-Cell Lymphoma", 
    "esca_tcga_pan_can_atlas_2018": "Esophageal Adenocarcinoma",
    "gbm_tcga_pan_can_atlas_2018": "Glioblastoma Multiforme", 
    "hnsc_tcga_pan_can_atlas_2018": "Head and Neck Squamous Cell Carcinoma", 
    "kich_tcga_pan_can_atlas_2018": "Kidney Chromophobe",
    "kirc_tcga_pan_can_atlas_2018": "Kidney Renal Clear Cell Carcinoma", 
    "kirp_tcga_pan_can_atlas_2018": "Kidney Renal Papillary Cell Carcinoma", 
    "lgg_tcga_pan_can_atlas_2018": "Brain Lower Grade Glioma",
    "lihc_tcga_pan_can_atlas_2018": "Liver Hepatocellular Carcinoma", 
    "luad_tcga_pan_can_atlas_2018": "Lung Adenocarcinoma", 
    "lusc_tcga_pan_can_atlas_2018": "Lung Squamous Cell Carcinoma",
    "meso_tcga_pan_can_atlas_2018": "Mesothelioma", 
    "ov_tcga_pan_can_atlas_2018": "Ovarian Serous Cystadenocarcinoma", 
    "paad_tcga_pan_can_atlas_2018": "Pancreatic Adenocarcinoma",
    "pcpg_tcga_pan_can_atlas_2018": "Pheochromocytoma and Paraganglioma", 
    "prad_tcga_pan_can_atlas_2018": "Prostate Adenocarcinoma", 
    "sarc_tcga_pan_can_atlas_2018": "Sarcoma",
    "skcm_tcga_pan_can_atlas_2018": "Skin Cutaneous Melanoma", 
    "stad_tcga_pan_can_atlas_2018": "Stomach Adenocarcinoma", 
    "tgct_tcga_pan_can_atlas_2018": "Testicular Germ Cell Tumors",
    "thca_tcga_pan_can_atlas_2018": "Thyroid Carcinoma", 
    "thym_tcga_pan_can_atlas_2018": "Thymoma", 
    "ucec_tcga_pan_can_atlas_2018": "Uterine Corpus Endometrial Carcinoma",
    "ucs_tcga_pan_can_atlas_2018": "Uterine Carcinosarcoma", 
    "uvm_tcga_pan_can_atlas_2018": "Uveal Melanoma",
    "tcga_read": "Rectum Adenocarcinoma"
}

# calculate mutation counts for each sample in each study
mutation_counts = {study_id: {} for study_id in study_ids.keys()}  # {study_id1 : {sample_id1: mut_count1, sample_id2: mutcount2, ...}, ... , study_id33 : {...}} 

for study_id in mutation_counts.keys():
  # get df of maf file
  with open("/cta/users/ilaydakaytaran/maf_analysis_2/cancer_studies/" + study_id + "/mc3_maf_" + study_id + ".txt") as f:
    maf_df = pd.read_csv(f, sep="\t", low_memory=False)
    maf_df.dropna(how='all', inplace=True)
    maf_df = maf_df.reset_index(drop=True)

  # get all sample ids in the study
  sample_ids = maf_df["sample_id"].unique()

  # find the mutation count of each sample in study
  mutations_in_samples = {}
  for sample_id in sample_ids:
    sample_df = maf_df[maf_df['sample_id'] == sample_id]
    mut_count = sample_df.shape[0]
    mutations_in_samples[sample_id] = mut_count
  mutation_counts[study_id] = mutations_in_samples

  # save mutations_in_samples
  df = pd.DataFrame(mutations_in_samples.values(), index=mutations_in_samples.keys())
  df.to_csv("/cta/users/ilaydakaytaran/maf_analysis_2/cancer_studies_background/" + study_id + "/mutations_in_samples_" + study_id + ".txt", sep="\t")

"""
# box plot visualisation - unordered

# extract study IDs and mutation counts
study_labels = []
mutation_counts_studies = []
for study_id in study_ids:
  study_labels.append(study_ids[study_id])
  sample_counts = mutation_counts[study_id]
  mutation_counts_studies.append(list(sample_counts.values()))
  #mutation_counts_studies.append(list(np.log2(i) for i in sample_counts.values()))

plt.figure(figsize=(30, 10))

# plot each boxplot
plt.boxplot(mutation_counts_studies, labels=study_labels, patch_artist=True)

plt.xlabel('Cancer Studies')
plt.ylabel('Mutation Counts')
#plt.ylabel('log2(Mutation Counts)')
plt.title('Box Plot of Mutation Counts for Different Cancer Studies')
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
plt.grid(True)
plt.tight_layout()  # Adjust layout to prevent clipping of labels
plt.savefig('/cta/users/ilaydakaytaran/maf_analysis_2/mutation_counts.pdf', format='pdf')
#plt.savefig('/cta/users/ilaydakaytaran/maf_analysis_2/mutation_counts_log2.pdf', format='pdf')

"""

# box plot visualisation - ordered

# calculate means of mutation counts for each study
#means = {study_id: np.mean(list(counts.values())) for study_id, counts in mutation_counts.items()}
means = {study_id: np.mean(np.log2(list(counts.values()))) for study_id, counts in mutation_counts.items()}

# sort study labels based on means
sorted_study_labels = [study_ids[study_id] for study_id, _ in sorted(means.items(), key=lambda x: x[1])]

# sort mutation counts accordingly and extract values
#sorted_mutation_counts = [list(mutation_counts[study_id].values()) for study_id in sorted(means, key=means.get)]
sorted_mutation_counts = [np.log2(list(mutation_counts[study_id].values())) for study_id in sorted(means, key=means.get)]

plt.figure(figsize=(30, 10))

# Plot each boxplot
plt.boxplot(sorted_mutation_counts, labels=sorted_study_labels, patch_artist=True)

plt.xlabel('Cancer Studies')
#plt.ylabel('Mutation Counts')
plt.ylabel('log2(Mutation Counts)')
plt.title('Box Plot of Mutation Counts for Different Cancer Studies (Ordered by Mean)')
plt.xticks(rotation=45, ha='right')
plt.grid(True)
plt.tight_layout()
#plt.savefig('/cta/users/ilaydakaytaran/maf_analysis_2/mutation_counts_ordered.pdf', format='pdf')
plt.savefig('/cta/users/ilaydakaytaran/maf_analysis_2/mutation_counts_ordered_log2.pdf', format='pdf')