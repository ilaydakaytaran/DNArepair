import sys

def main(study_id, knee_index):
  # import modules
  import pandas as pd
  import numpy as np

  # get df of all protein coding genes
  with open("ncbi_all_gene_result.txt", "r") as f_in:
    df_ncbi = pd.read_csv(f_in, sep='\t')           # initial gene data coming from NCBI
  df_ncbi = df_ncbi[['GeneID', 'Symbol']]

  # get df of all repairs
  with open("df_all_repair_ner.txt", "r") as f_in:
    df_all_repair = pd.read_csv(f_in, sep="\t")  # df of NER, BER, DSB, MMR and CLR, with swissprotID and protein_length
  df_all_repair.drop("Unnamed: 0", axis=1, inplace=True)

  # get df of ercc6
  df_ercc6 = df_all_repair[df_all_repair["Symbol"] == "ERCC6"]
  df_ercc6.reset_index(drop=True, inplace=True)

  # get df of other repairs
  df_other_repair = df_all_repair[~(df_all_repair['Pathway'].str.contains("NER"))]
  df_other_repair.reset_index(drop=True, inplace=True)

  from cbio_py import cbio_mod as cb

  # obtain the mutated genes in study
  mutations_in_study = cb.getMutationsInMolecularProfile(study_id + "_mutations", study_id + "_all",
                                                        projection='DETAILED', return_type = 'dict', append='no')
                                                        # list of all mutations in all genes in all samples in this study

  # save mutations_in_study
  df = pd.DataFrame(mutations_in_study)
  df.to_csv(study_id + "/mutations_in_study_" + study_id + ".txt", sep="\t")


  # get all ncbi genes
  mutations_ncbi = {int(gene_id): [] for gene_id in df_ncbi["GeneID"]}   # {gene1: [mut_sample1, mut_sample2...], gene2: [mut_sample1, mut_sample2...]...
                        #  cumulative: [mut_sample1, mut_sample2...]}
  mutations_ncbi["cumulative"] = []
  for mut in mutations_in_study:
    mut_gene_id =  mut["gene"].entrezGeneId
    if mut_gene_id in mutations_ncbi:
      sample_id = mut["sampleId"]
      if sample_id not in mutations_ncbi[mut_gene_id]:
        mutations_ncbi[mut_gene_id].append(sample_id)
      if sample_id not in mutations_ncbi["cumulative"]:
        mutations_ncbi["cumulative"].append(sample_id)

  # get percentages of ncbi genes
  total_sample_num = len(mutations_ncbi['cumulative'])
  mutations_ncbi_percentages = {}
  for gene_id, samples in mutations_ncbi.items():
    mutations_ncbi_percentages[gene_id] = float(format((len(samples) / total_sample_num), ".10f"))

  df_ncbi["Percentage"] = float("NaN")
  for gene_id, percentage in mutations_ncbi_percentages.items():
    df_ncbi.loc[df_ncbi["GeneID"] == gene_id, "Percentage"] = percentage
  df_ncbi = df_ncbi.sort_values(by="Percentage", ascending=False)
  df_ncbi = df_ncbi.reset_index()

  # save df_ncbi
  df_ncbi.to_csv(study_id + "/df_ncbi_" + study_id + ".txt", sep="\t")

  # get the genes above the knee
  df_ncbi = df_ncbi.head(knee_index-1)


  # get all potential background genes
  all_background_genes = {}  # {gene: [backgene1, backgene2, ...] for gene in df_ercc6["Symbol"]} for study
  for gene in df_ercc6["GeneID"]:        # iterate genes in each study
    gene = int(gene)
    background_genes = []
    ner_length = int(df_ercc6[df_ercc6["GeneID"] == gene]["protein_length"].iloc[0])

    # set range
    range_min = ner_length * (1-0.1)
    range_max = ner_length * (1+0.1)

    for row in df_other_repair.itertuples(index=False):    # iterate the dataframe for suitable genes
      if int(row.GeneID) not in df_ncbi['GeneID'].values:
        other_length = int(row.protein_length)
        if range_min <= other_length and other_length <= range_max:
          background_genes.append(row.Symbol)
    all_background_genes[gene] = background_genes

  # save all_background_genes
  df = pd.DataFrame(all_background_genes.values(), index=all_background_genes.keys())
  df.to_csv(study_id + "/all_background_genes_ercc6_" + study_id + ".txt", sep="\t")

  # prepare randomised lists
  import random as r
  r.seed(1234)
  random_backgrounds = {} # {random_list1: [random_gene1...], random_list2: [random_gene1...], ... , random_list30 : [random_gene1...], ner_list: [ner_gene1...]}
  ner_list = ["ERCC6"]
  for i in range(len(all_background_genes[2074])):    # prepare random lists from repair genes
    random_list = [all_background_genes[2074][i]] 
    random_backgrounds["random_list"+str(i+1)] = random_list
  random_backgrounds["ner_list"] = ner_list

  # save random_backgrounds
  df = pd.DataFrame(random_backgrounds.values(), index = random_backgrounds.keys())
  df.to_csv(study_id + "/random_backgrounds_ercc6_" + study_id + ".txt", sep="\t")

  # obtain mutated samples for background genes
  cancer_type_summary_ercc6 = {} # {gene: [] for gene in df_ercc6["Gene"]} for study_id
  df = pd.DataFrame(mutations_in_study)
  for randomlist, genes in random_backgrounds.items():
    gene_dict = {gene: [] for gene in genes}
    gene_dict["samples_with_data"] = df['sampleId'].unique()     # get the unique IDs in study with mutation data

    gene_dict["cumulative"] = []
    for mut in mutations_in_study:
      gene = mut["gene"].hugoGeneSymbol
      if gene in gene_dict:
        sample_id = mut["sampleId"]
        if sample_id not in gene_dict[gene]:
          gene_dict[gene].append(sample_id)
        if sample_id not in gene_dict["cumulative"]:
          gene_dict["cumulative"].append(sample_id)
    if randomlist != "ner_list":
      cancer_type_summary_ercc6[randomlist] = gene_dict

  # save cancer_type_summary_ercc6
  df = pd.DataFrame(cancer_type_summary_ercc6.values(), index = cancer_type_summary_ercc6.keys())
  df.to_csv(study_id + "/cancer_type_summary_ercc6_" + study_id + ".txt", sep="\t")

  # obtain mutated samples for ner genes
  ner_gene_dict = {gene: [] for gene in df_ercc6["Symbol"]}
  df = pd.DataFrame(mutations_in_study)
  ner_gene_dict["samples_with_data"] = df['sampleId'].unique()     # get the unique IDs in study with mutation data

  ner_gene_dict["cumulative"] = []
  for mut in mutations_in_study:
    gene = mut["gene"].hugoGeneSymbol
    if gene in ner_gene_dict:
      sample_id = mut["sampleId"]
      if sample_id not in ner_gene_dict[gene]:
        ner_gene_dict[gene].append(sample_id)
      if sample_id not in ner_gene_dict["cumulative"]:
        ner_gene_dict["cumulative"].append(sample_id)

  # save ner_gene_dict
  df = pd.DataFrame(ner_gene_dict.values(), index=ner_gene_dict.keys())
  df.to_csv(study_id + "/ner_gene_dict_ercc6_" + study_id + ".txt", sep="\t")

  # calculate mutation percentages for background genes
  cancer_type_summary_ercc6_sample_percentage = {}
  for randomlist, gene_dict in cancer_type_summary_ercc6.items():
    dct = {}
    sample_num = len(gene_dict["samples_with_data"])
    for gene in gene_dict:
      if gene != "samples_with_data":
        dct[gene] = float(format(len(gene_dict[gene]) * 100 / sample_num, ".3f"))
    cancer_type_summary_ercc6_sample_percentage[randomlist] = dct

  # save cancer_type_summary_ercc6_sample_percentage
  df = pd.DataFrame(cancer_type_summary_ercc6_sample_percentage.values(), index = cancer_type_summary_ercc6_sample_percentage.keys())
  df.to_csv(study_id + "/cancer_type_summary_sample_percentage_ercc6_" + study_id + ".txt", sep="\t")

  # calculate mutation percentages for ner
  sample_num = len(ner_gene_dict["samples_with_data"])
  ner_sample_percentages = {}
  for gene in ner_gene_dict:
    if gene != "samples_with_data":
      ner_sample_percentages[gene] = float(format(len(ner_gene_dict[gene]) * 100 / sample_num, ".3f"))

  # save ner_sample_percentages
  df = pd.DataFrame(ner_sample_percentages.values(), index = ner_sample_percentages.keys())
  df.to_csv(study_id + "/ner_sample_percentages_ercc6_" + study_id + ".txt", sep="\t")

if __name__ == "__main__":
  study_id = sys.argv[1]
  knee_index = int(sys.argv[2])  # Convert the second argument to an integer
    
  main(study_id, knee_index) 
