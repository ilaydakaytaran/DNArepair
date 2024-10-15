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
    "uvm_tcga_pan_can_atlas_2018": "Uveal Melanoma"
}

# get significant z-scores
z_scores_gg = {}
for study_id in study_ids.keys():
    with open(study_id + "/p_value_ercc6.txt", "r") as f:
        line = f.readlines()[0]
        entries = line.replace(" ", "\t").strip().split("\t")
        p_value = float(entries[2])
        if p_value <= 0.05:
            z_score = float(entries[4])
            if z_score >= 0:
                z_scores_gg[study_ids[study_id]] = z_score
# save into a file
with open("significants/significant_z_scores_ercc6.txt", "w") as f:
    for study_id, z_score in z_scores_gg.items():
        f.write(study_id + "\t" + str(z_score) + "\n")


# visualise 
from matplotlib import pyplot as plt

plt.figure(figsize=(12,8))

study_names_gg = list(z_scores_gg.keys())
z_scores_gg_values = list(z_scores_gg.values())
x_indices = range(len(study_names_gg))
num_bars = len(study_names_gg)
bar_width = 2/len(study_names_gg)
x_indices_shifted = [i + bar_width * 0.5 for i in x_indices]  # Shifted x-coordinates
tick_positions = [i + bar_width / 2 for i in x_indices_shifted]

plt.bar(x_indices_shifted, z_scores_gg_values, color="#2a9d8f", width = 2/len(study_names_gg))
plt.xlabel('Study Name')
plt.ylabel('Z-score')
plt.title('Mutations - ERCC6 - Significant Z-Scores', fontweight='bold', fontsize = 14)
plt.xticks(tick_positions, study_names_gg, rotation=0, fontsize=9, ha='center')
plt.yticks(fontsize=10)
plt.grid(axis='y', linestyle='--', alpha=0.2, color="#2a9d8f")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('significants/z_scores_ercc6.pdf', format='pdf')

