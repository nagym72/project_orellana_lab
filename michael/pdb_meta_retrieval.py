#!/bin/env python

import atomium
import os
import sys
from collections import defaultdict
import json

def retrieve_meta_info(pdb_4_digit:str, outpath:str)->dict:

    #print(pdb_files)
    #little helper to deal with the time objects in pdb files.
    def date_encoder(obj):
        if isinstance(obj, date):
            return obj.isoformat()  # Convert date to ISO format


    meta_dict = defaultdict()      
    #pdb_name = os.path.basename(path)[0:4]
    meta_info_storage = os.path.join(outpath, f"{pdb_4_digit}_meta_info_structures.json")
    pdb_codes = pdb_4_digit+".pdb"
    
    sub_dict = defaultdict()
    pdb = atomium.fetch(pdb_4_digit)
    
    sub_dict["title"] = pdb.title
    sub_dict["deposition_date"] = pdb.deposition_date.isoformat()  #isoformat because it is a time object
    sub_dict["technique"] = pdb.technique
    sub_dict["authors"] = pdb.authors
    #sub_dict["file_type"] = pdb.filetype
    #sub_dict["missing_res"] = pdb.missing_residues
    sub_dict["key_words"] = pdb.keywords
    sub_dict["code"] = pdb.code
    sub_dict["classification"] = pdb.classification
    sub_dict["organism"] = pdb.source_organism
    sub_dict["expression_system"] = pdb.expression_system
    sub_dict["resolution"] = pdb.resolution
    sub_dict["r_val"] = pdb.rvalue
    sub_dict["r_free"] = pdb.rfree
    sub_dict["classification"] = pdb.classification
    #print(dir(pdb))
    
    sub_dict['number_of_residues_asymmetric_unit'] = len(pdb.model.residues())
    sub_dict['mass_dalton_asymetric_unit'] = pdb.model.mass 
    
    try:
        pdb1 = atomium.fetch(pdb_4_digit)
    
        assembly = pdb1.generate_assembly(1)
        #print(len(assembly.residues()))
        #print(assembly.mass)
        sub_dict['number_of_residues_biological_unit'] = len(assembly.residues())
        sub_dict['mass_dalton_biological_unit'] = assembly.mass
        
    except Exception as error:
        sub_dict['number_of_residues'] = "Error in computation for length"
        
    with open(meta_info_storage, "w") as file:        
        json.dump(sub_dict, file, default=date_encoder)
            
    #return meta_dict

    return meta_info_storage
    
if __name__ == "__main__":

    print("we start now")
    pdb_4_digit = sys.argv[1]
    outpath = sys.argv[2]
    outlocation = retrieve_meta_info(pdb_4_digit=pdb_4_digit, outpath=outpath)
    print(f"done. Results are saved in: {outlocation}")
