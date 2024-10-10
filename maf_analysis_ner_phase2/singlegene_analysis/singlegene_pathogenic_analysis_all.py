import sys
def main(study_id):
  # import modules
  import pandas as pd
  from kneed import KneeLocator
  
  # get df of all protein coding genes
  with open("../ncbi_all_gene_result_27122023.txt", "r") as f_in:
    df_ncbi = pd.read_csv(f_in, sep='\t')           # initial gene data coming from NCBI
  df_ncbi = df_ncbi[['GeneID', 'Symbol']]
  
  # get df of all repairs
  with open("../../backgrounds/df_all_repair_ner.txt", "r") as f_in:
    df_all_repair = pd.read_csv(f_in, sep="\t")  # df of NER, BER, DSB, MMR and CLR, with swissprotID and protein_length
  df_all_repair.drop("Unnamed: 0", axis=1, inplace=True)
  
  # get df of all
  df_all = df_all_repair[(df_all_repair['Pathway'].str.contains('NER')) & ((df_all_repair['other_designations'] == "GG, TCR") | (df_all_repair['other_designations'] == "GG") | (df_all_repair['other_designations'] == "TCR"))]
  df_all = df_all.drop_duplicates(subset='GeneID')  # make it unique based on GeneID --> no duplicate XAB2
  df_all.reset_index(drop=True, inplace=True)
  
  # get df of other repairs
  df_other_repair = df_all_repair[~(df_all_repair['Pathway'].str.contains("NER"))]
  df_other_repair.reset_index(drop=True, inplace=True)

  # get df of rnap2
  with open("../../backgrounds/df_rna_ner.txt", "r") as f_in:
    df_rna = pd.read_csv(f_in, sep="\t")
  df_rna.drop("Unnamed: 0", axis=1, inplace=True)
  df_rna = df_rna[~(df_rna['Pathway'].str.contains("NER"))]

  # get df of maf file
  with open("../cancer_studies/" + study_id + "/mc3_maf_" + study_id + ".txt") as f:
    maf_df = pd.read_csv(f, sep="\t", low_memory=False)
    maf_df.dropna(how='all', inplace=True)
    maf_df = maf_df.reset_index(drop=True)
  
  # filter maf df on clin_sig and polyphen for pathogenicity
  maf_df['CLIN_SIG_list'] = maf_df['CLIN_SIG'].str.split(';')
    # 1 - clin_sig has at least one of this
  condition1 = maf_df['CLIN_SIG_list'].apply(lambda x: any(val in x for val in ['pathogenic', 'likely_pathogenic', 'risk_factor']))
    # 2 - clin_sig has none of this
  condition2 = maf_df['CLIN_SIG_list'].apply(lambda x: not any(val in x for val in ['benign', 'likely_benign', 'uncertain_significance']))
    # 3 - polyphen has at least one of this
  condition3 = maf_df['PolyPhen'].apply(lambda x: any(val in x for val in ['probably_damaging', 'possibly_damaging']))
    # 4 - polyphen has none of this
  condition4 = maf_df['PolyPhen'].apply(lambda x: not any(val in x for val in ['benign', 'unknown']))
    # 5 - a frameshift or nonsense mutation
  condition5 = maf_df["Variant_Classification"].apply(lambda x: any(val in x for val in ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation']))
    # 6 - not a silent mutation
  condition6 = maf_df['Variant_Classification'].apply(lambda x: not any(val in x for val in ['Silent']))
    # combine -> ((1&2)|(3&4)|5)&6
  maf_df_filtered = maf_df[((condition1 & condition2) | (condition3 & condition4) | condition5) & condition6]
  
  # save obtained maf df to a file
  maf_df_filtered.to_csv("../cancer_studies/" + study_id + "/maf_df_filtered_" + study_id + ".txt", sep="\t")

  # get all ncbi genes' mutated sample ids
  mutations_ncbi = {int(gene_id): [] for gene_id in df_ncbi["GeneID"]}   # {gene_id1: [mut_sample1, mut_sample2...], gene_id2: [mut_sample1, mut_sample2...]...
                        #  cumulative: [mut_sample1, mut_sample2...]}
  for row in maf_df_filtered.itertuples(index=False):
    mut_gene_id =  row.Entrez_Gene_Id
    if mut_gene_id in mutations_ncbi:  # skip genes with unknown gene ids
      sample_id = row.sample_id
      if sample_id not in mutations_ncbi[mut_gene_id]:
        mutations_ncbi[mut_gene_id].append(sample_id)

  # get all sample ids in mc3_maf 
  with open("../cancer_studies/" + study_id + "/sample_ids_" + study_id + ".txt") as f:
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
  df_ncbi.to_csv("../cancer_studies/" + study_id + "/df_ncbi_" + study_id + ".txt", sep="\t")
  
  # calculate the knee point
  kneedle = KneeLocator(range(len(df_ncbi["Percentage"])), df_ncbi["Percentage"], S=1, curve="concave", online=True, direction="decreasing", interp_method="interp1d")
  with open("../cancer_studies/" + study_id + "/knee_point_" + study_id + ".txt", "w") as f_knee:
    f_knee.write(str(kneedle.knee) + "\t" + str(kneedle.knee_y))

  # get the genes above the knee
  df_ncbi = df_ncbi.head(kneedle.knee+1)

  # filter out the genes below 100 aa length
  df_all = df_all[df_all["protein_length"]>100]

  for single_gene in df_all["GeneID"].values:
    # get all potential background genes
    all_background_genes = {}  # {gene: [backgene1, backgene2, ...] for single_gene
    single_gene = int(single_gene)
    single_gene_symbol = str(df_all[df_all['GeneID']==single_gene]["Symbol"].iloc[0])
    background_genes = []
    ner_length = int(df_all[df_all["GeneID"] == single_gene]["protein_length"].iloc[0])
    
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
    all_background_genes[single_gene] = background_genes
    
    # if there is not enough genes to create 30 random lists
    if len(all_background_genes[single_gene]) == 0:  
      with open("./all/" + study_id + "/no_random_lists_all_" + study_id + ".txt", "a+") as f_no_random:
        f_no_random.write(single_gene_symbol + "\n") 
    else:
      # save all_background_genes
      with open("./all/" + study_id + "/all_background_genes_all_" + study_id + ".txt", "a+") as f_all_back:
        f_all_back.write(single_gene_symbol)
        for back_gene in all_background_genes[single_gene]:
          f_all_back.write("\t" + back_gene)
        f_all_back.write("\n")
      # prepare randomised lists
      import random as r
      r.seed(1234)
      random_backgrounds = {} # {random_list1: [random_gene1], random_list2: [random_gene2], ... , random_list30 : [random_gene30], ner_list: [single_gene]}
      ner_list = [single_gene_symbol]
      for i in range(30):    # prepare 30 random lists
        random_list = [r.choice(all_background_genes[single_gene])]
        random_backgrounds["random_list"+str(i+1)] = random_list
      random_backgrounds["ner_list"] = ner_list
      
      # save random_backgrounds
      with open("./all/" + study_id + "/random_backgrounds_all_" + study_id + ".txt", "a+") as f_random:
        f_random.write(single_gene_symbol)
        for i in range(30):
          f_random.write("\t" + random_backgrounds["random_list"+str(i+1)][0])
        f_random.write("\n")
    
      # obtain mutated samples for background genes
      cancer_type_summary_all = {} # {randomlist: genedict} where genedict is dict of mutated sample ids for the gene in randomlist
      for randomlist, genes in random_backgrounds.items():
        gene_dict = {genes[0]: []} 
        gene_dict["samples_with_data"] = maf_df['sample_id'].unique()     # get the unique IDs in study with mutation data
        gene_dict["cumulative"] = []
        for mut in maf_df_filtered.itertuples(index=False):
          gene = mut.Hugo_Symbol
          if gene in gene_dict:
            sample_id = mut.sample_id
            if sample_id not in gene_dict[gene]:
              gene_dict[gene].append(sample_id)
            if sample_id not in gene_dict["cumulative"]:
              gene_dict["cumulative"].append(sample_id)
        if randomlist != "ner_list":
          cancer_type_summary_all[randomlist] = gene_dict
      
      # save cancer_type_summary_all
      df = pd.DataFrame(cancer_type_summary_all.values(), index = cancer_type_summary_all.keys())
      df.to_csv("./all/" + study_id + "/cancer_type_summary_all_"+ single_gene_symbol + "_" + study_id + ".txt", sep="\t")
      
      # obtain mutated samples for single gene
      ner_gene_dict = {single_gene_symbol: []}
      ner_gene_dict["samples_with_data"] = maf_df['sample_id'].unique()     # get the unique IDs in study with mutation data
      ner_gene_dict["cumulative"] = []
      for mut in maf_df_filtered.itertuples(index=False):
        gene = mut.Hugo_Symbol
        if gene in ner_gene_dict:
          sample_id = mut.sample_id
          if sample_id not in ner_gene_dict[gene]:
            ner_gene_dict[gene].append(sample_id)
          if sample_id not in ner_gene_dict["cumulative"]:
            ner_gene_dict["cumulative"].append(sample_id)
      
      # save ner_gene_dict
      df = pd.DataFrame(ner_gene_dict.values(), index=ner_gene_dict.keys())
      df.to_csv("./all/" + study_id + "/ner_gene_dict_all_" + single_gene_symbol + "_" + study_id + ".txt", sep="\t")

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
      df.to_csv("./all/" + study_id + "/cancer_type_summary_sample_percentage_all_" + single_gene_symbol + "_" + study_id + ".txt", sep="\t")

      # calculate mutation percentages for ner
      sample_num = len(ner_gene_dict["samples_with_data"])
      ner_sample_percentages = {}
      for gene in ner_gene_dict:
        if gene != "samples_with_data":
          ner_sample_percentages[gene] = float(format(len(ner_gene_dict[gene]) * 100 / sample_num, ".3f"))
      
      # save ner_sample_percentages
      df = pd.DataFrame(ner_sample_percentages.values(), index = ner_sample_percentages.keys())
      df.to_csv("./all/" + study_id + "/ner_sample_percentages_all_" + single_gene_symbol + "_" + study_id + ".txt", sep="\t")

  
if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("Usage: python script.py <study_id>")
    sys.exit(1)
  study_id = sys.argv[1]
    
  main(study_id)
