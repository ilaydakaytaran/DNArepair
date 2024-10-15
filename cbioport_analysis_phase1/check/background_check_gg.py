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

  # get df of gg
  df_gg = df_all_repair[(df_all_repair['Pathway'].str.contains('NER')) & (df_all_repair['other_designations']=="GG")]
  df_gg.reset_index(drop=True, inplace=True)

  # get df of other repairs
  df_other_repair = df_all_repair[~(df_all_repair['Pathway'].str.contains("NER"))]
  df_other_repair.reset_index(drop=True, inplace=True)


  # install cbio_py
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
  all_background_genes = {}  # {gene: [backgene1, backgene2, ...] for gene in df_gg["Symbol"]} for study
  for gene in df_gg["GeneID"].values:        # iterate genes in each study
    gene = int(gene)
    background_genes = []
    ner_length = int(df_gg[df_gg["GeneID"] == gene]["protein_length"].iloc[0])

    # set range
    range_min = ner_length * (1-0.05)
    range_max = ner_length * (1+0.05)

    for row in df_other_repair.itertuples(index=False):    # iterate the dataframe for suitable genes
      if int(row.GeneID) not in df_ncbi['GeneID'].values:
        other_length = int(row.protein_length)
        if range_min <= other_length and other_length <= range_max:
          background_genes.append(row.Symbol)

    all_background_genes[gene] = background_genes

  # save all_background_genes
  df = pd.DataFrame(all_background_genes.values(), index=all_background_genes.keys())
  df.to_csv(study_id + "/all_background_genes_gg_" + study_id + ".txt", sep="\t")


  # prepare randomised lists
  import random as r
  r.seed(1234)
  random_backgrounds = {} # {random_list1: [random_gene1...], random_list2: [random_gene1...], ... , random_list30 : [random_gene1...], ner_list: [ner_gene1...]}
  ner_list = []
  for i in range(50):    # prepare 50 random lists
    random_list = []
    for ner_gene, potential_back_genes_list in all_background_genes.items():  # iterate over the potential genes
      if (i == 0):
        ner_list.append(df_gg[df_gg['GeneID']==ner_gene]["Symbol"].values[0])
      random_back_gene = r.choice(potential_back_genes_list)
      while (random_back_gene in random_list):  # don't add the same gene multiple times
        random_back_gene = r.choice(potential_back_genes_list)
      random_list.append(random_back_gene)
    random_backgrounds["random_list"+str(i+1)] = random_list
  random_backgrounds["ner_list"] = ner_list

  # save random_backgrounds
  df = pd.DataFrame(random_backgrounds.values(), index = random_backgrounds.keys())
  df.to_csv(study_id + "/random_backgrounds_gg_50_" + study_id + ".txt", sep="\t") 


  # obtain mutated samples for background genes
  cancer_type_summary_gg = {} # {gene: [] for gene in df_gg["Gene"]} for study_id
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
      cancer_type_summary_gg[randomlist] = gene_dict

  # save cancer_type_summary_gg
  df = pd.DataFrame(cancer_type_summary_gg.values(), index = cancer_type_summary_gg.keys())
  df.to_csv(study_id + "/cancer_type_summary_gg_50_" + study_id + ".txt", sep="\t")


  # obtain mutated samples for ner genes
  ner_gene_dict = {gene: [] for gene in df_gg["Symbol"]}
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
  df.to_csv(study_id + "/ner_gene_dict_gg_" + study_id + ".txt", sep="\t")


  # calculate mutation percentages for background genes
  cancer_type_summary_gg_sample_percentage = {}
  for randomlist, gene_dict in cancer_type_summary_gg.items():
    dct = {}
    sample_num = len(gene_dict["samples_with_data"])
    for gene in gene_dict:
      if gene != "samples_with_data":
        dct[gene] = float(format(len(gene_dict[gene]) * 100 / sample_num, ".3f"))
    cancer_type_summary_gg_sample_percentage[randomlist] = dct

  # save cancer_type_summary_gg_sample_percentage
  df = pd.DataFrame(cancer_type_summary_gg_sample_percentage.values(), index = cancer_type_summary_gg_sample_percentage.keys())
  df.to_csv(study_id + "/cancer_type_summary_sample_percentage_gg_50_" + study_id + ".txt", sep="\t")


  # calculate mutation percentages for ner
  sample_num = len(ner_gene_dict["samples_with_data"])
  ner_sample_percentages = {}
  for gene in ner_gene_dict:
    if gene != "samples_with_data":
      ner_sample_percentages[gene] = float(format(len(ner_gene_dict[gene]) * 100 / sample_num, ".3f"))

  # save ner_sample_percentages
  df = pd.DataFrame(ner_sample_percentages.values(), index = ner_sample_percentages.keys())
  df.to_csv(study_id + "/ner_sample_percentages_gg_" + study_id + ".txt", sep="\t")


  
if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python script.py <study_id> <knee_index>")
    sys.exit(1)
  study_id = sys.argv[1]
  knee_index = int(sys.argv[2])  # Convert the second argument to an integer
    
  main(study_id, knee_index)