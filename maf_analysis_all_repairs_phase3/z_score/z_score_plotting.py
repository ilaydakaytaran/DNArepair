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



# get significant z-scores
z_scores_all = {}
p_values_all = {}
for study_id in study_ids.keys():
    with open("./cancer_studies_background/" + study_id + "/p_value_all.txt", "r") as f:
        line = f.readlines()[0]
        entries = line.replace(" ", "\t").strip().split("\t")
        p_value = float(entries[2])
        if p_value <= 0.05:
            z_score = float(entries[4])
            if z_score >= 0:
                p_values_all[study_ids[study_id]] = p_value
                z_scores_all[study_ids[study_id]] = z_score
# save into a file
with open("significants/significant_scores_all_background.txt", "w") as f:
    for study in z_scores_all.keys():
        f.write(study + "\t" + str(p_values_all[study]) + "\t" + str(z_scores_all[study]) + "\n")


# visualise 
from matplotlib import pyplot as plt

if len(list(z_scores_all.keys())) > 0:
    plt.figure(figsize=(12,8))

    study_names_all = list(z_scores_all.keys())
    z_scores_all_values = list(z_scores_all.values())
    x_indices = range(len(study_names_all))

    plt.bar(x_indices, z_scores_all_values, color="#2a9d8f", width = 0.8)
    plt.xlabel('Study Name')
    plt.ylabel('Z-score')
    plt.title('Mutations - ALL REPAIR GENES - Significant Z-Scores', fontweight='bold', fontsize = 14)
    plt.xticks(x_indices, study_names_all, rotation=30, fontsize=9, ha='right')
    plt.yticks(fontsize=10)
    plt.grid(axis='y', linestyle='--', alpha=0.2, color="#2a9d8f")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig('significants/z_scores_all_background.pdf', format='pdf')

"""
# SAME CODE FOR PATHOGENIC ANALYSIS
# get significant z-scores
z_scores_all = {}
p_values_all = {}
for study_id in study_ids.keys():
    with open("./cancer_studies/" + study_id + "/p_value_all.txt", "r") as f:
        line = f.readlines()[0]
        entries = line.replace(" ", "\t").strip().split("\t")
        p_value = float(entries[2])
        if p_value <= 0.05:
            z_score = float(entries[4])
            if z_score >= 0:
                p_values_all[study_ids[study_id]] = p_value
                z_scores_all[study_ids[study_id]] = z_score
# save into a file
with open("significants/significant_scores_all.txt", "w") as f:
    for study in z_scores_all.keys():
        f.write(study + "\t" + str(p_values_all[study]) + "\t" + str(z_scores_all[study]) + "\n")


# visualise 
from matplotlib import pyplot as plt

if len(list(z_scores_all.keys())) > 0:
    plt.figure(figsize=(12,8))

    study_names_all = list(z_scores_all.keys())
    z_scores_all_values = list(z_scores_all.values())
    x_indices = range(len(study_names_all))

    plt.bar(x_indices, z_scores_all_values, color="#2a9d8f", width = 0.8)
    plt.xlabel('Study Name')
    plt.ylabel('Z-score')
    plt.title('Pathogenic Mutations - ALL REPAIR GENES - Significant Z-Scores', fontweight='bold', fontsize = 14)
    plt.xticks(x_indices, study_names_all, rotation=30, fontsize=9, ha='right')
    plt.yticks(fontsize=10)
    plt.grid(axis='y', linestyle='--', alpha=0.2, color="#2a9d8f")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig('significants/z_scores_all.pdf', format='pdf')
"""