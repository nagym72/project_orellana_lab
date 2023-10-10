def map_clinvar(Gene_name:str):
    
    path = "/home/micnag/bioinformatics/clinvar_trials/variant_summary.txt"

    use_cols = ["Type", "Name", "GeneSymbol",
           "ClinicalSignificance", "PhenotypeList",
           "Assembly", "ChromosomeAccession", 
           "Chromosome", "Start", "Stop"]

    column_data_types = {
    "Type": str,
    "Name": str,
    "GeneSymbol": str,
    "ClinicalSignificance": str,
    "PhenotypeList": str,
    "Assembly": str,
    "ChromosomeAccession": str,
    "Chromosome": str,
    "Start": int,
    "Stop": int
    }

    df_work = pd.read_csv(path, sep="\t", usecols=use_cols, dtype=column_data_types)
    
    df_work.loc[:, "AA_change"] = df_work["Name"].str.split().str.get(-1)
    df_work.loc[:, "AA_change"] = df_work["AA_change"].str.replace("(", "")
    df_work.loc[:, "AA_change"] = df_work["AA_change"].str.replace(")", "")
    
    df_work.loc[:,"Original_AA"] = df_work["AA_change"].str[2:5]
    df_work.loc[:,"Modified_AA"] = df_work["AA_change"].str[-3:]
    df_work['Position'] = pd.to_numeric(df_work['AA_change'].str[5:-3], errors='coerce')
    
    # Drop rows with NaN values in the 'Position' column
    df_work.dropna(subset=['Position'], inplace=True)
    df_work['Position'] = df_work['Position'].astype(int)
    
    df_work["Genomic_location"] = df_work["Chromosome"] + ":" + df_work["Start"].astype(str)
    df_work["gnomad_aa_change"] = "p." + df_work["Original_AA"] + df_work["Position"].astype(str) + df_work["Modified_AA"]
    
    df_work = df_work.drop("AA_change", axis=1)
    df_work = df_work.drop("Name", axis=1)
    df_work = df_work.drop("Chromosome", axis=1)
    df_work = df_work.drop("Start", axis=1)
    df_work = df_work.drop("Stop", axis=1)
    
    df_work["Allele_count"] = [np.nan] * len(df_work)
    df_work["Allele_number"] = [np.nan] * len(df_work)
    df_work["Allele_frequency"] = [np.nan] * len(df_work)
    
    
    accepted_residues = ["Ala", "Gly", "Ser", "Leu", "Pro",
                    "Ile", "Val", "Phe", "Tyr", "Trp",
                     "His", "Thr", "Asn", "Gln", "Asp", 
                     "Glu","Cys", "Met", "Lys", "Arg"]
    
    
    #filtering based on our Gene name.
    df_filtered = df_work[(df_work["Type"] == "single nucleotide variant") & 
        (df_work["GeneSymbol"] == Gene_name) &
        (df_work["Assembly"] == "GRCh37") & 
        (df_work['Original_AA'].isin(accepted_residues)) &
        (df_work['Modified_AA'].isin(accepted_residues)) ]
    
    return df_filtered

if __name__ == "__main__":
    print("placeholder")
    #your code here guys.
