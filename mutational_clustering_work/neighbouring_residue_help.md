=============================================================================
Hello, this bash script is intended to be used as a wrapper for my nearest_neighbour_atom_search.py script.
The purpose of this script is to compute neighbours for all atoms in a supplied .pdb file. This works in the following way: We parse through the pdb file and for each single atom we compute the surrounding **residues** that are within 8A vicinity of each atom. The results for each atom will be stored in a list and returned. Additionally, this will also give you as a return a list that is **residue** based, so that the returned list contains each residue in the pdb file and all neighbouring residues within 8A vicinity.

Example how to use:

"$0" [- f <filepath>] or [ -d <dirpath> ] 

in this case, filepath or dirpath should be the location of the pdb file/(s) that will be used to compute for surrounding neighbours.
The cutoff is fixed at 8A (corresponding to ~2 CA-distances) but can be changed in a later version.
The script requires biopython but we will prompt you if you want to install it.
If you already have it installed, still click yes and proceed.
=============================================================================

