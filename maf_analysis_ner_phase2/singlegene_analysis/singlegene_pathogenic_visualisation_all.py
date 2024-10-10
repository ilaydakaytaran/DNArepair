from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize


all_genes = ['CCNH', 'CDK7', 'CETN2', 'ERCC8', 'DDB1', 'DDB2', 'ERCC1', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERCC6', 
             'GTF2H1', 'GTF2H2', 'GTF2H3', 'GTF2H4', 'LIG1', 'MNAT1', 'POLD1', 'POLE', 'POLR2A', 'POLR2B', 'POLR2C', 
             'POLR2D', 'POLR2E', 'POLR2F', 'POLR2G', 'POLR2H', 'POLR2I', 'POLR2J', 'RAD23A', 'RAD23B', 'RPA1', 'RPA2', 
             'RPA3', 'XPA', 'XPC', 'USP7', 'CUL4A', 'RBX1', 'POLD3', 'RNF111', 'XAB2', 'UVSSA', 'POLR2M', 'POLR2J2', 'POLR2J3']


study_ids = ["acc_tcga_pan_can_atlas_2018", "ov_tcga_pan_can_atlas_2018"]
study_names = ["Adrenocortical Carcinoma", "Ovarian Serous Cystadenocarcinoma"] 

p_values_all = [] # store p_values list of each study in the same order as study_ids
z_scores_all = [] # store z_scores list of each study in the same order as study_ids
for study_id in study_ids:
  p_values = []  # store p values of the genes in the same order as all_genes
  z_scores = []  # store z scores of the genes in the same order as all_genes
  for all_gene in all_genes:
    df = pd.read_csv("./all/" + study_id + "/cancer_type_summary_sample_percentage_all_" + all_gene + "_" + study_id + ".txt", sep="\t")
    df.set_index("Unnamed: 0", inplace=True)
    # calculate sample percentages of the singlegene's random lists
    cancer_type_summary_all_singlegene_sample_percentage = {}
    cols = list(df.columns)
    for randomlist, row in df.iterrows():
      gene_dict = {}
      for col in cols:
        if not pd.isna(row[col]):
          gene_dict[col] = float(row[col])
      cancer_type_summary_all_singlegene_sample_percentage[randomlist] = gene_dict

    df = pd.read_csv("./all/" + study_id + "/ner_sample_percentages_all_" + all_gene + "_" + study_id + ".txt", sep="\t")
    df.set_index("Unnamed: 0", inplace=True)
    # calculate sample percentage of the singlegene
    ner_singlegene_sample_percentages = {}
    for index, row in df.iterrows():
      ner_singlegene_sample_percentages[index] = float(row["0"])

    # form zipped version of data for plotting
    random_lists = [] 
    random_list_percentages = []

    for random_list, percentage_dict in cancer_type_summary_all_singlegene_sample_percentage.items():
      random_lists.append(random_list)
      random_list_percentages.append(np.double(percentage_dict["cumulative"]))
    sorted_data = sorted(zip(random_list_percentages, random_lists), reverse=True)
    sorted_values, sorted_categories = zip(*sorted_data)

    # z-test for significance
    mean = np.mean(random_list_percentages)
    std_dev = np.std(random_list_percentages)
    if std_dev > 1e-15:
      z_score = (ner_singlegene_sample_percentages["cumulative"] - mean) / (std_dev)
      p_value = stats.norm.sf(abs(z_score))
    else:
      z_score = np.nan
      p_value = np.nan
    
    # store gene's values separately for study
    p_values.append(p_value)
    z_scores.append(z_score)

  p_values_all.append(p_values)
  z_scores_all.append(z_scores)

# save p_values_all
with open("./all/p_values_all.txt", "w") as f:
  f.write("ALL_genes\t")
  for gene in all_genes:
    f.write("\t" + gene)
  f.write("\n")

  for idx in range(len(study_ids)):
    study_id = study_ids[idx]
    f.write(study_id + "\t")
    for p_value in p_values_all[idx]:
      f.write("\t" + str(p_value))
    f.write("\n")

# save z_scores_all
with open("./all/z_scores_all.txt", "w") as f:
  f.write("ALL_genes\t")
  for gene in all_genes:
    f.write("\t" + gene)
  f.write("\n")

  for idx in range(len(study_ids)):
    study_id = study_ids[idx]
    f.write(study_id + "\t")
    for z_score in z_scores_all[idx]:
      f.write("\t" + str(z_score))
    f.write("\n")

# create a redness colormap
redness_cmap = LinearSegmentedColormap.from_list('custom_redness', [(1, 0.8, 0.8), (1, 0, 0)], N=256)
# define the normalization for your colormap based on the range of p-values
norm = Normalize(vmin=-np.log10(0.05), vmax=15)

# create a figure and axis
fig, ax = plt.subplots(figsize=(16, 28))

# plot the scatter plot for each gene in each study
for i, gene in enumerate(all_genes):
  # extract p-values and z-scores for the current gene
  gene_p_values = [study[i] for study in p_values_all]
  gene_z_scores = [study[i] for study in z_scores_all]
  
  # create a mask for significant genes
  significant_mask = [(p <= 0.05) and (z >= 0) for p, z in zip(gene_p_values, gene_z_scores)]
  significant_x_values = [j for j, is_significant in enumerate(significant_mask) if is_significant]
  significant_y_values = [-np.log10(p) for j, p in enumerate(gene_p_values) if j in significant_x_values]
  
  # scatter plot for significant genes
  scatter_sig = ax.scatter(significant_x_values, significant_y_values, c=significant_y_values, cmap=redness_cmap, norm=norm)

  # scatter plot for insignificant genes
  ax.scatter([j for j, is_significant in enumerate(significant_mask) if not is_significant], 
             [-np.log10(gene_p_values[j]) for j, is_significant in enumerate(significant_mask) if not is_significant], 
             color='grey')
  
  # annotate each significant point with the gene name
  for j, study in enumerate(study_ids):
    if significant_mask[j]:
      point_color = redness_cmap(norm(-np.log10(gene_p_values[j])))
      ax.annotate(gene, (j, -np.log10(gene_p_values[j])), textcoords="offset points", xytext=(30, -2), ha='center', color=point_color, fontsize=12)
    else:
      ax.annotate(gene, (j, -np.log10(gene_p_values[j])), textcoords="offset points", xytext=(30, -2), ha='center', color='grey', fontsize=12)

# customize the plot
for i in range(len(study_ids) - 1): 
  ax.vlines(i + 0.5, ymin=0, ymax=15, linestyle='dotted', color='grey', alpha=0.5)
ax.axhline(y=-np.log10(0.05), color='r', linestyle='--', label='Significance Limit (p-value <= 0.05)')
ax.set_xticks(range(len(study_names)))
ax.set_xticklabels(study_names, rotation=45, ha='right', fontsize=18)
ax.set_ylabel('-log10(P-values)', fontsize=20)
ax.set_title('Significant Genes in Different Cancer Studies - ALL', fontsize=25)
cbar = plt.colorbar(scatter_sig, ax=ax, label='P-Values')
# Calculate the ticks
start = -np.log10(0.05)
end = 15
num_ticks = 6
tick_interval = (end - start) / (num_ticks - 1)
tick_positions = [start + i * tick_interval for i in range(num_ticks)]
cbar.set_ticks(tick_positions)
cbar.set_ticklabels([str(val) for val in [0.05, 0.04, 0.03, 0.02, 0.01, 0.0]])
cbar.set_label('P-Values', rotation=270, labelpad=15, fontsize=20)
ax.legend()

# save the plot
plt.tight_layout()  # add this line for better layout adjustment
plt.savefig('all/significant_genes_per_study_all.pdf', format='pdf')
