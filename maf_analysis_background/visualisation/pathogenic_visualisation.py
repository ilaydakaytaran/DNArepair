import sys

def main(study_id):
  import pandas as pd
  import numpy as np

  df = pd.read_csv("./cancer_studies/" + study_id + "/cancer_type_summary_sample_percentage_gg_1000_" + study_id + ".txt", sep="\t")
  df.set_index("Unnamed: 0", inplace=True)

  # calculate sample percentages
  cancer_type_summary_gg_1000_sample_percentage = {}
  cols = list(df.columns)
  for randomlist, row in df.iterrows():
    gene_dict = {}
    for col in cols:
      if not pd.isna(row[col]):
        gene_dict[col] = float(row[col])
    cancer_type_summary_gg_1000_sample_percentage[randomlist] = gene_dict

  df = pd.read_csv("./cancer_studies/" + study_id + "/ner_sample_percentages_gg_" + study_id + ".txt", sep="\t")
  df.set_index("Unnamed: 0", inplace=True)

  ner_sample_percentages = {}
  for index, row in df.iterrows():
    ner_sample_percentages[index] = float(row["0"])

  # form zipped version of data for plotting
  random_lists = [] 
  random_list_percentages = []

  for random_list, percentage_dict in cancer_type_summary_gg_1000_sample_percentage.items():
    random_lists.append(random_list)
    random_list_percentages.append(np.double(percentage_dict["cumulative"]))
  sorted_data = sorted(zip(random_list_percentages, random_lists), reverse=True)
  sorted_values, sorted_categories = zip(*sorted_data)

  # visualise
  from matplotlib import pyplot as plt

  plt.figure(figsize=(17, 8))
  plt.bar(sorted_categories, sorted_values, color = "purple")

  plt.xlabel('Random Lists')
  plt.ylabel('Frequency (%)')
  plt.title(study_id + ' Pathogenic Mutations - GG - 1000 Lists', fontweight='bold', fontsize = 14)

  plt.xticks(rotation=-45, fontsize=10, ha='left')
  plt.yticks(fontsize=10)
  plt.grid(axis='y', linestyle='--', alpha=0.4, color="purple")

  plt.gca().spines['top'].set_visible(False)
  plt.gca().spines['right'].set_visible(False)

  plt.savefig("./cancer_studies/" + study_id + '/mut_cancer_type_summary_gg_1000_' + study_id + '.pdf', format='pdf')

  # create normal distribution
  from scipy import stats
  from scipy.stats import norm


  # z-test for significance
  mean = np.mean(random_list_percentages)
  std_dev = np.std(random_list_percentages)
  if std_dev > 1e-15:
    z_score = (ner_sample_percentages["cumulative"] - mean) / (std_dev)
    p_value = stats.norm.sf(abs(z_score))
  else:
    z_score = np.nan
    p_value = np.nan
  
  with open("./cancer_studies/" + study_id + "/p_value_gg.txt", "w") as f:
    f.write(study_id + "\tp_value: " + str(p_value) + "\tz_score: " + str(z_score) + "\n")

  # visualise
  plt.figure(figsize=(15,8))

  x_range = np.linspace(min(random_list_percentages), max(random_list_percentages), 10000)
  plt.plot(x_range, norm.pdf(x_range, mean, std_dev), 'r')

  plt.axvline(x=float(ner_sample_percentages["cumulative"]), color='g', linestyle='--')

  plt.xlabel('Pathogenic Mutation Percentage of Randomised Repair Genes')
  plt.ylabel('Frequency')
  plt.title('Normal Distribution - 1000 Lists - ' + study_id)
  legend_labels = ['Null Case Pathogenic Mutation Percentage', f'Observed Value\nZ-score: {z_score:.2f}\np-value: {p_value:.4f}']
  plt.legend(legend_labels)
  plt.savefig("./cancer_studies/" + study_id + '/mut_norm_dist_gg_1000_' + study_id + '.pdf', format='pdf')


if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("Usage: python script.py <study_id>")
    sys.exit(1)
  study_id = sys.argv[1]
  main(study_id)
