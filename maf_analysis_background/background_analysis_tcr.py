import sys

def main(study_id):
  # import modules
  import pandas as pd
  from kneed import KneeLocator

  # get df of all protein coding genes
  with open("ncbi_all_gene_result_27122023.txt", "r") as f_in:
    df_ncbi = pd.read_csv(f_in, sep='\t')           # initial gene data coming from NCBI
  df_ncbi = df_ncbi[['GeneID', 'Symbol']]

  # get df of all repairs
  with open("../backgrounds/df_all_repair_ner.txt", "r") as f_in:
    df_all_repair = pd.read_csv(f_in, sep="\t")  # df of NER, BER, DSB, MMR and CLR, with swissprotID and protein_length
  df_all_repair.drop("Unnamed: 0", axis=1, inplace=True)

  # get df of tcr
  df_tcr = df_all_repair[(df_all_repair['Pathway'].str.contains('NER')) & (df_all_repair['other_designations']=="TCR")]
  df_tcr.reset_index(drop=True, inplace=True)

  # get df of other repairs
  df_other_repair = df_all_repair[~(df_all_repair['Pathway'].str.contains("NER"))]
  df_other_repair.reset_index(drop=True, inplace=True)

  # get df of rnap2
  with open("../backgrounds/df_rna_ner.txt", "r") as f_in:
    df_rna = pd.read_csv(f_in, sep="\t")
  df_rna.drop("Unnamed: 0", axis=1, inplace=True)
  df_rna = df_rna[~(df_rna['Pathway'].str.contains("NER"))]

  # get df of maf file
  with open("./cancer_studies/" + study_id + "/mc3_maf_" + study_id + ".txt") as f:
    maf_df = pd.read_csv(f, sep="\t", low_memory=False)
    maf_df.dropna(how='all', inplace=True)
    maf_df = maf_df.reset_index(drop=True)


  # get all ncbi genes
  mutations_ncbi = {int(gene_id): [] for gene_id in df_ncbi["GeneID"]}   # {gene1: [mut_sample1, mut_sample2...], gene2: [mut_sample1, mut_sample2...]...
                        #  cumulative: [mut_sample1, mut_sample2...]}

  for row in maf_df.itertuples(index=False):
    mut_gene_id =  row.Entrez_Gene_Id
    if mut_gene_id in mutations_ncbi:   # skip genes with unknown gene ids
      sample_id = row.sample_id
      if sample_id not in mutations_ncbi[mut_gene_id]:
        mutations_ncbi[mut_gene_id].append(sample_id)

  # get all sample ids in mc3_maf 
  with open("./cancer_studies/" + study_id + "/sample_ids_" + study_id + ".txt") as f:
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
  df_tcr = df_tcr[df_tcr["protein_length"]>100]


  # get all potential background genes
  all_background_genes = {}  # {gene: [backgene1, backgene2, ...] for gene in df_tcr["Symbol"]} for study
  for gene in df_tcr["GeneID"].values:        # iterate genes in each study
    gene = int(gene)
    background_genes = []
    ner_length = int(df_tcr[df_tcr["GeneID"] == gene]["protein_length"].iloc[0])

    # set range
    range_min = ner_length * (1-0.05)
    range_max = ner_length * (1+0.05)

    for row in df_other_repair.itertuples(index=False):    # iterate the dataframe for suitable genes
      if int(row.GeneID) not in df_ncbi['GeneID'].values:
        other_length = int(row.protein_length)
        if range_min <= other_length and other_length <= range_max:
          background_genes.append(row.Symbol)
    if background_genes == []:                 # iterate the rnap2 genes for necessary background
      for row in df_rna.itertuples(index=False):
        if int(row.GeneID) not in df_ncbi["GeneID"].values:
          other_length = int(row.protein_length)
          if range_min <= other_length and other_length <= range_max:
            background_genes.append(row.Symbol)

    all_background_genes[gene] = background_genes

  # save all_background_genes
  df = pd.DataFrame(all_background_genes.values(), index=all_background_genes.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/all_background_genes_tcr_" + study_id + ".txt", sep="\t")

  # prepare randomised lists
  import random as r
  r.seed(1234)
  random_backgrounds = {} # {random_list1: [random_gene1...], random_list2: [random_gene1...], ... , random_list30 : [random_gene1...], ner_list: [ner_gene1...]}
  ner_list = []
  for i in range(1000):    # prepare 1000 random lists
    random_list = []
    for ner_gene, potential_back_genes_list in all_background_genes.items():  # iterate over the potential genes
      if (i == 0):
        ner_list.append(df_tcr[df_tcr['GeneID']==ner_gene]["Symbol"].values[0])
      random_back_gene = r.choice(potential_back_genes_list)
      while (random_back_gene in random_list):  # don't add the same gene multiple times
        random_back_gene = r.choice(potential_back_genes_list)
      random_list.append(random_back_gene)
    random_backgrounds["random_list"+str(i+1)] = random_list
  random_backgrounds["ner_list"] = ner_list

  # save random_backgrounds
  df = pd.DataFrame(random_backgrounds.values(), index = random_backgrounds.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/random_backgrounds_tcr_1000_" + study_id + ".txt", sep="\t") 


  # obtain mutated samples for background genes
  cancer_type_summary_tcr = {} # {randomlist: genedict} where genedict is dict of mutated sample ids for each gene in randomlist
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
    if randomlist != "ner_list":
      cancer_type_summary_tcr[randomlist] = gene_dict

  # save cancer_type_summary_tcr
  df = pd.DataFrame(cancer_type_summary_tcr.values(), index = cancer_type_summary_tcr.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/cancer_type_summary_tcr_1000_" + study_id + ".txt", sep="\t")

  # obtain mutated samples for ner genes
  ner_gene_dict = {gene: [] for gene in df_tcr["Symbol"]}
  ner_gene_dict["samples_with_data"] = maf_df['sample_id'].unique()     # get the unique IDs in study with mutation data

  ner_gene_dict["cumulative"] = []
  for mut in maf_df.itertuples(index=False):
    gene = mut.Hugo_Symbol
    if gene in ner_gene_dict:
      sample_id = mut.sample_id
      if sample_id not in ner_gene_dict[gene]:
        ner_gene_dict[gene].append(sample_id)
      if sample_id not in ner_gene_dict["cumulative"]:
        ner_gene_dict["cumulative"].append(sample_id)

  # save ner_gene_dict
  df = pd.DataFrame(ner_gene_dict.values(), index=ner_gene_dict.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/ner_gene_dict_tcr_" + study_id + ".txt", sep="\t")

  # calculate mutation percentages for background genes
  cancer_type_summary_tcr_sample_percentage = {}
  for randomlist, gene_dict in cancer_type_summary_tcr.items():
    dct = {}
    sample_num = len(gene_dict["samples_with_data"])
    for gene in gene_dict:
      if gene != "samples_with_data":
        dct[gene] = float(format(len(gene_dict[gene]) * 100 / sample_num, ".3f"))
    cancer_type_summary_tcr_sample_percentage[randomlist] = dct

  # save cancer_type_summary_tcr_sample_percentage
  df = pd.DataFrame(cancer_type_summary_tcr_sample_percentage.values(), index = cancer_type_summary_tcr_sample_percentage.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/cancer_type_summary_sample_percentage_tcr_1000_" + study_id + ".txt", sep="\t")

  # calculate mutation percentages for ner
  sample_num = len(ner_gene_dict["samples_with_data"])
  ner_sample_percentages = {}
  for gene in ner_gene_dict:
    if gene != "samples_with_data":
      ner_sample_percentages[gene] = float(format(len(ner_gene_dict[gene]) * 100 / sample_num, ".3f"))

  # save ner_sample_percentages
  df = pd.DataFrame(ner_sample_percentages.values(), index = ner_sample_percentages.keys())
  df.to_csv("./cancer_studies_background/" + study_id + "/ner_sample_percentages_tcr_" + study_id + ".txt", sep="\t")

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("Usage: python script.py <study_id>")
    sys.exit(1)
  study_id = sys.argv[1]
    
  main(study_id) 