def get_cosmic_mutations(gene_name:str):

    path = "/home/micnag/bioinformatics/cosmic/Cosmic_GenomeScreensMutant_v98_GRCh38.tsv"
    
    usecols=['GENE_SYMBOL',
         'MUTATION_AA', 'MUTATION_DESCRIPTION',
       'MUTATION_ZYGOSITY', 'LOH', 'CHROMOSOME', 'GENOME_START', 'GENOME_STOP']

    df = dd.read_csv(path, sep="\t", dtype={'CHROMOSOME': 'object',
       'MUTATION_ZYGOSITY': 'object', 'GENOME_START': 'float64',
       'GENOME_STOP': 'float64',
       'LOH': 'object'}, usecols=usecols)

    #we need to switch these tuples and then map the 1letter aa code to 3letter aa for later compatibility.
    lst =  [('Val',"V"), ('Ile',"I"), ('Leu',"L"), ('Glu',"E"), ('Gln',"Q"),
                    ('Asp',"D"), ('Asn',"N"), ('His',"H"), ('Trp',"W"), ('Phe',"F"), ('Tyr',"Y"), 
                    ('Arg',"R"), ('Lys',"K"), ('Ser',"S"), ('Thr',"T"), ('Met',"M"), ('Ala',"A"), 
                    ('Gly',"G"), ('Pro',"P"), ('Cys',"C")]

    lst = [(y, x) for x, y in lst]

    canonical_aas = defaultdict(lambda: "X", lst)

    df_re = df[df["MUTATION_DESCRIPTION"].str.contains("missense")]
    
    df_re = df_re[df_re["GENE_SYMBOL"] == f"{gene_name}"]

    meta = ('Gene name', 'str') 
    df_re['CHROMOSOME'] = df_re['CHROMOSOME'].astype('object')
    df_re['WT_AA'] = df_re['MUTATION_AA'].str[2].apply(lambda x: canonical_aas[x], meta=meta)
    df_re['MUTATION_POSITION'] = df_re['MUTATION_AA'].str[3:-1]
    df_re['MUTATED_AA'] = df_re['MUTATION_AA'].str[-1].apply(lambda x: canonical_aas[x], meta=meta)

    df_re = df_re.drop("MUTATION_AA", axis=1)
    
    cosmic_df = df_re.compute()

    cosmic_df["GENOME_START"] = cosmic_df["GENOME_START"].astype(int)
    cosmic_df["GENOME_STOP"] = cosmic_df["GENOME_STOP"].astype(int)
    
    return cosmic_df

if __name__ == "__main__":
    print("placeholder")
    #your code here guys.
