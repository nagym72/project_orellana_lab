def map_gnomad(Gene_name:str, outpath:str):

    outpath = f"{outpath}/gnomad_datatable.txt"
    #outpath="/home/micnag/bioinformatics/hail_trials/newline_test.tsv"
    
    path = '/home/micnag/bioinformatics/gnomad_raw_data/gnomad_data.mt'
    #path1 = "/home/micnag/bioinformatics/gnomad_raw_data/gnomad_data_2.mt"
    mt = hl.read_matrix_table(path)
    
    substring1 = Gene_name
    substring2 = "missense"
    
    
    mt = mt.annotate_rows(Gene_names=mt.info.vep.map(
        lambda x: x.split("\|")[3]) ,
                      type_of_change = mt.info.vep.map(
        lambda x: x.split("\|")[1]) , 
                      AA_change = mt.info.vep.map(
        lambda x: x.split("\|")[11]) , 
                      ENST_identifier= mt.info.vep.map(
        lambda x: x.split("\|")[6])

    ) 
             
    filtered_mt_2 = mt.filter_rows(
    
    #hl.any(lambda x: hl.str(x).contains(substring3), mt.AA_change)
    hl.any(lambda x: hl.str(x).contains(substring1), mt.info.vep) &
    hl.any(lambda x: hl.str(x).contains(substring2), mt.info.vep)
    
    )
                     
    filtered_mt_3 = filtered_mt_2.annotate_rows(
        Allele_count_int = filtered_mt_2.info.AC,
        Allele_frequency_float = filtered_mt_2.info.AF,
        Allele_number_int = filtered_mt_2.info.AN,
        Gene_name_str = _replace_empty(filtered_mt_2.Gene_names), 
        Mutation_change_str = _replace_empty(filtered_mt_2.AA_change),
        Type_of_change_str = _replace_empty(filtered_mt_2.type_of_change))
    
    
    
    rows_to_keep = ["Gene_name_str", "Mutation_change_str", "Type_of_change_str", "Allele_count_int",
                "Allele_frequency_float", "Allele_number_int"]


    selected_rows = filtered_mt_3.select_rows(
        Allele_count_int=filtered_mt_3.Allele_count_int,
        Allele_frequency_float=filtered_mt_3.Allele_frequency_float,
        Allele_number_int=filtered_mt_3.Allele_number_int,
        Gene_name_str=hl.str(filtered_mt_3.Gene_name_str),
        Mutation_change_str=hl.str(filtered_mt_3.Mutation_change_str),
        Type_of_change_str=hl.str(filtered_mt_3.Type_of_change_str)
            )

    save_buffer = selected_rows.select_rows(*rows_to_keep)
    
    select_rows_out = save_buffer.rows()
    
    select_rows_out.export(outpath)
    
    #this is the location where we save the results.
    return outpath

if __name__ == "__main__":
    print("placeholder")
    #your code here guys.
