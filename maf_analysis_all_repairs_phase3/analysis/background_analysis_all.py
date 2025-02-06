import sys

def main(study_id):
  # import modules
  import pandas as pd
  from kneed import KneeLocator 

  # get df of all protein coding genes
  with open("../maf_analysis_all_repairs/df_all_most_23032024.txt", "r") as f_in:
    df_ncbi = pd.read_csv(f_in, sep='\t')           # initial gene data coming from NCBI with additional swissprot_id and protein_length
  df_ncbi = df_ncbi[df_ncbi['swissprot_id'] != "NotFound"]   # filter our genes without known protein length
  df_ncbi = df_ncbi[['GeneID', 'Symbol', 'protein_length']]

  # get df of all repairs
  with open("../backgrounds/df_all_repair_ner.txt", "r") as f_in:
    df_all_repair = pd.read_csv(f_in, sep="\t")  # df of NER, BER, DSB, MMR and CLR, with swissprotID and protein_length
  df_all_repair.drop("Unnamed: 0", axis=1, inplace=True)

  # get df of non-repair genes
  gene_ids_all_repair = set(df_all_repair["GeneID"])
  df_non_repair = df_ncbi[~df_ncbi["GeneID"].isin(gene_ids_all_repair)]

  # get df of maf file
  with open("../maf_analysis_2/cancer_studies/" + study_id + "/mc3_maf_" + study_id + ".txt") as f:
    maf_df = pd.read_csv(f, sep="\t", low_memory=False)
    maf_df.dropna(how='all', inplace=True)
    maf_df = maf_df.reset_index(drop=True)

  # get all ncbi genes
  mutations_ncbi = {int(gene_id): [] for gene_id in df_ncbi["GeneID"]}   # {gene1: [mut_sample1, mut_sample2...], gene2: [mut_sample1, mut_sample2...]...
                        #  cumulative: [mut_sample1, mut_sample2...]}
  for row in maf_df.itertuples(index=False):
    mut_gene_id = row.Entrez_Gene_Id
    if mut_gene_id in mutations_ncbi:   # skip genes with unknown gene ids
      sample_id = row.sample_id
      if sample_id not in mutations_ncbi[mut_gene_id]:
        mutations_ncbi[mut_gene_id].append(sample_id)

  # get all sample ids in mc3_maf 
  with open("../maf_analysis_2/cancer_studies/" + study_id + "/sample_ids_" + study_id + ".txt") as f:
    line = f.readlines()[0]
    mutations_ncbi["cumulative"] = line.strip().split(",")[:-1]  # has all samples in the study


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
  df_ncbi.to_csv("./cancer_studies_background/" + study_id + "/df_ncbi_" + study_id + ".txt", sep="\t")

  # calculate the knee point
  kneedle = KneeLocator(range(len(df_ncbi["Percentage"])), df_ncbi["Percentage"], S=1, curve="concave", online=True, direction="decreasing", interp_method="interp1d")
  with open("./cancer_studies_background/" + study_id + "/knee_point_" + study_id + ".txt", "w") as f_knee:
    f_knee.write(str(kneedle.knee) + "\t" + str(kneedle.knee_y))

  # get the genes above the knee
  df_ncbi = df_ncbi.head(kneedle.knee+1)


  # filter out the genes below 100 aa length
  df_all_repair = df_all_repair[df_all_repair["protein_length"]>100]


  # get all potential background genes
  all_background_genes = {}  # {gene: [backgene1, backgene2, ...] for gene in df_all_repair["Symbol"]} for study
  for gene in df_all_repair["GeneID"].values:        # iterate genes in each study
    gene = int(gene)
    background_genes = []
    repair_length = int(df_all_repair[df_all_repair["GeneID"] == gene]["protein_length"].iloc[0])

    # set range
    range_min = repair_length * (1-0.05)
    range_max = repair_length * (1+0.05)

    for row in df_non_repair.itertuples(index=False):    # iterate the dataframe for suitable genes
      if int(row.GeneID) not in df_ncbi['GeneID'].values:
        other_length = int(row.protein_length)
        if range_min <= other_length and other_length <= range_max:
          background_genes.append(row.Symbol)

    all_background_genes[gene] = background_genes

  # save all_background_genes
  df = pd.DataFrame(all_background_genes.values(), index=all_background_genes.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/all_background_genes_all_" + study_id + ".txt", sep="\t")

  # prepare randomised lists
  import random as r
  r.seed(1234)
  random_backgrounds = {} # {random_list1: [random_gene1...], random_list2: [random_gene1...], ... , random_list30 : [random_gene1...], repair_list: [repair_gene1...]}
  repair_list = []
  for i in range(1000):    # prepare 1000 random lists
    random_list = []
    for repair_gene, potential_back_genes_list in all_background_genes.items():  # iterate over the potential genes
      if (i == 0):
        repair_list.append(df_all_repair[df_all_repair['GeneID']==repair_gene]["Symbol"].values[0])
      random_back_gene = r.choice(potential_back_genes_list)
      while (random_back_gene in random_list):  # don't add the same gene multiple times
        random_back_gene = r.choice(potential_back_genes_list)
      random_list.append(random_back_gene)
    random_backgrounds["random_list"+str(i+1)] = random_list
  random_backgrounds["repair_list"] = repair_list

  # save random_backgrounds
  df = pd.DataFrame(random_backgrounds.values(), index = random_backgrounds.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/random_backgrounds_all_1000_" + study_id + ".txt", sep="\t") 

  # obtain mutated samples for background genes
  cancer_type_summary_all = {} # {randomlist: genedict} where genedict is dict of mutated sample ids for each gene in randomlist
  for randomlist, genes in random_backgrounds.items():
    gene_dict = {gene: [] for gene in genes}
    gene_dict["samples_with_data"] = maf_df['sample_id'].unique()     # get the unique IDs in study with mutation data

    gene_dict["cumulative"] = []
    for mut in maf_df.itertuples(index=False):
      gene = mut.Hugo_Symbol
      if gene in gene_dict:
        sample_id = mut.sample_id
        if sample_id not in gene_dict[gene]:
          gene_dict[gene].append(sample_id)
        if sample_id not in gene_dict["cumulative"]:
          gene_dict["cumulative"].append(sample_id)
    if randomlist != "repair_list":
      cancer_type_summary_all[randomlist] = gene_dict

  # save cancer_type_summary_all
  df = pd.DataFrame(cancer_type_summary_all.values(), index = cancer_type_summary_all.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/cancer_type_summary_all_1000_" + study_id + ".txt", sep="\t")

  # obtain mutated samples for repair genes
  repair_gene_dict = {gene: [] for gene in df_all_repair["Symbol"]}
  repair_gene_dict["samples_with_data"] = maf_df['sample_id'].unique()     # get the unique IDs in study with mutation data

  repair_gene_dict["cumulative"] = []
  for mut in maf_df.itertuples(index=False):
    gene = mut.Hugo_Symbol
    if gene in repair_gene_dict:
      sample_id = mut.sample_id
      if sample_id not in repair_gene_dict[gene]:
        repair_gene_dict[gene].append(sample_id)
      if sample_id not in repair_gene_dict["cumulative"]:
        repair_gene_dict["cumulative"].append(sample_id)

  # save repair_gene_dict
  df = pd.DataFrame(repair_gene_dict.values(), index=repair_gene_dict.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/repair_gene_dict_all_" + study_id + ".txt", sep="\t")

  # calculate mutation percentages for background genes
  cancer_type_summary_all_sample_percentage = {}
  for randomlist, gene_dict in cancer_type_summary_all.items():
    dct = {}
    sample_num = len(gene_dict["samples_with_data"])
    for gene in gene_dict:
      if gene != "samples_with_data":
        dct[gene] = float(format(len(gene_dict[gene]) * 100 / sample_num, ".3f"))
    cancer_type_summary_all_sample_percentage[randomlist] = dct

  # save cancer_type_summary_all_sample_percentage
  df = pd.DataFrame(cancer_type_summary_all_sample_percentage.values(), index = cancer_type_summary_all_sample_percentage.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/cancer_type_summary_sample_percentage_all_1000_" + study_id + ".txt", sep="\t")

  # calculate mutation percentages for repair genes
  sample_num = len(repair_gene_dict["samples_with_data"])
  repair_sample_percentages = {}
  for gene in repair_gene_dict:
    if gene != "samples_with_data":
      repair_sample_percentages[gene] = float(format(len(repair_gene_dict[gene]) * 100 / sample_num, ".3f"))

  # save repair_sample_percentages
  df = pd.DataFrame(repair_sample_percentages.values(), index = repair_sample_percentages.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/repair_sample_percentages_all_" + study_id + ".txt", sep="\t")

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("Usage: python script.py <study_id>")
    sys.exit(1)
  study_id = sys.argv[1]
    
  main(study_id) 