#!/bin/python3

from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB import NeighborSearch
import os
from collections import defaultdict
import sys

# specify path to pdb structures. In the future, this will we done by iterating through os module
# in order to get all files in a given directory.
# path_to_pdb = sys.argv[1]
# name of protein
prot_name = "defaultprot"

# we generate a parser object based on the PDBParser class.
# This object will be used to generate a Structure object from a pdb file.
parser = PDBParser(QUIET=True)


# open pdb file and make a Structure object with our parser class object.

def surr_atoms(inpath, protname, cutoff=8, outpath=os.getcwd()):
    with open(inpath, "r") as pdbfile1:

        # first we need to make extract all atoms from our pdb file.
        structure = parser.get_structure(protname, pdbfile1)

        # Selection.unfold_entities(<structure object>, <level of information that you want>)
        # other levels are "C" for chain, "R" for residue, "A" for atom and so on.
        atom_list = Selection.unfold_entities(structure, "A")

        # lets get the coordinates of all atoms now
        atom_coords = [(atom.get_coord(), atom.get_parent()) for atom in atom_list]
        # parse through atom_list and apply .get_coord() to each retrieved object.
        # we also store for each atom the parent residue

        # we provide as argument here the Selection.unfold.entities object which has all atoms.
        ns = NeighborSearch(atom_list)  # this class object has the .search() method defined in its __init__

        # we will store for each atom all residues that are within 8 A surroundings in this dictionary.
        f"""Keys: atom coordinates
        Values: all residues within {cutoff} A"""

        # we need a counter
        # we will go through all atoms of lets say residue 1: VAL.
        # this has 10 atoms so we need to search 10 times for each atom
        # within 8A cutoff radius, get all surrounding residues, and then
        # we will merge them together into a list and take only unique residues to get
        # rid of redundancy. Come to me Federico if you need explanation in detail given the next section is tricky.
        i = 1

        # we only want aa residue surroundings, excluding solvent and ligands
        aa_lst = ["VAL", "ALA", "GLY", "TRP", "ARG", "LYS", "LEU", "ILE", "ASP", "ASN", "GLN", "GLU", "PRO", "TYR",
                  "PHE",
                  "SER", "THR", "CYS", "MET", "HIS"]

        hits_per_atom_for_surr_residues = defaultdict(list)  # we store all hits in a dictionary
        for atoms in atom_coords:
            # if its a water atom, we are not interested in doing the neighbour search.
            # I will make this even more robust to exclude other ligands by only allowing the 20 aa to be selected.
            if atoms[1].get_resname() not in aa_lst:
                i += 1
                continue
            # if we are no longer in the same residue, we increment by 1 the counter and reset our tmp list.
            if atoms[1].get_id()[1] != i:
                i += 1

            f'''For each atom we will make a search for all surrounding atoms that are within {cutoff} A radius.'''

            proximal_atoms = ns.search(atoms[0], 8, "R")
            # I SET HERE search for atoms[0] because atoms is a tuple containing of coordinates
            # and parent residue name see line 75 + 76 #print(atom_coords[0])

            f"""Synthax: ns.search(<target object>, <Cutoff to be searched for>, 
            <type of information level that should be returned>
            R means we dont want the single atoms that are within {cutoff}A 
            found but instead their corresponding residues. For all atoms we would set <A> instead of <R>"""

            # this function searches through a target (in our case each atom as we loop through all available atoms)
            # and returns a list with all atoms within specified atoms .

            '''Take a look at the following print statement to see whats going on'''
            # print(f"The selected atom has the following coordinates:\nX:{atoms[0]}\nY:{atoms[1]}\nZ:{atoms[2]}\n \
            # These are all Residue ids that are within 8 A vicinity:\n")

            tmp = []  # we store all of them in a temporary list
            for residues in proximal_atoms:  # we go through all residues that were found within cutoff A
                id_x = residues.get_id()[1]
                # get_id gives us a tuple with shape ("", "residue number", "optinal flag").
                # Out of this tuple we want the residue id which is [1]
                # we only want residues that we dont have already in the list.
                # Makes no sense to add stuff that is already in there
                if id_x not in tmp and id_x != i:
                    # we also exclude the residue itself to be added to its neighbours.
                    tmp.append(id_x)
            # if we have all we append the whole list to the dictionary. we take the atoms parent residue name as a key.

            res_name = f"{atoms[1].get_resname()}{atoms[1].get_id()[1]}"
            # this string concatenation is super ugly to look at
            # and very confusing but it does the following:
            f"""{atoms[1].get_resname()} == residue name. In this case its a
                {atoms[1].get_id()[1]} """
            # so the resname would be in this case : VAL1 but you can modify the output as you wish. e.g VAL_1 or 1_VAL
            hits_per_atom_for_surr_residues[res_name].append(tmp)

        # print(hits_per_atom_for_surr_residues)
        # this shows we captures all aa within the protein
        # print(len(hits_per_atom_for_surr_residues))

        # lets make a dictionary containing all residues within 8 A on a residue base instead for all atoms.
        res_dict = defaultdict()
        # we parse through the old values and only add UNIQUE residues to the new dictionary so we dont have duplicates
        """UGLY SOLUTION BUT DOES THE TRICK"""
        for keys, values in hits_per_atom_for_surr_residues.items():
            # print(values)
            tmp_vals = []
            for vals in values:
                # each val is a list corresponding to 1 atom and its surrounded residues
                # which are the list entries. values is a list of lists covering the whole residue.
                for single_res in vals:
                    # now we go through all atoms of the residue.
                    if single_res not in tmp_vals:
                        # we add all atoms that are not already counted from previous atoms
                        # of the same target residue and add their surrounding residues to the tmp_vals list
                        # we don't add the atom itself
                        # print(single_res)
                        tmp_vals.append(single_res)

            # now we got 1 cycle done.
            # This corresponds to going through all atoms of e.g VALINE which has 7 atoms.
            # We group all residues that are neighbouring these 7 atoms and take only the unique ones.
            # This corresponds to the surrounding residues for the whole residue.
            res_dict[keys] = tmp_vals
            # quick look for the results.
            # sorry Federico for a bit mess above.
            # I will maybe refine it but this script works and you can directly implement it
            with open(outpath + protname, "w") as fh_out:
                for keys, values in res_dict.items():
                    full_hit = ""
                    fh_out.write(keys + ",")
                    for single_entries in values:
                        full_hit += str(single_entries) + " "
                    fh_out.write(full_hit + "\n")

                    # for keys, values in hits_per_atom_for_surr_residues.items():
                    #    #take a quick look at the result
                    #    print(f"here comes a list of all residues
                    #    that are in contact with all {len(values)} atoms that {keys} has:")
                    #    for i, vals in enumerate(values):
                    #        print(f"atom {i+1}: {vals}")
                    #    break

if __name__ == "__main__":
    in_path = sys.argv[1]
    prot_name = sys.argv[2]+".csv"
    out_path = sys.argv[3]
    surr_atoms(inpath=in_path, protname=prot_name, outpath=out_path)
    print(f"File containing neighbouring atoms is ready at {out_path + prot_name}")
