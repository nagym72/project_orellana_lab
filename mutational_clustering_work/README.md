#This readme belongs to the nearest_neighbour.sh bashscript and its associated nearest_neighbouring_search.py script.

The intented purpose of the caller2.sh script is to execute the neighbouring_search.py
script either for a single pdb file or a full high throughput pdb files in a directory.


===========================================================================================
**CURRENTLY NOT WORKING. PYTHON SCRIPT NEEDS TO BE ADJUSTED AND THEN IT WORKS AS INTENTED**
===========================================================================================


## Neighbouring_search.py:
We compute for each atom in the supplied pdb file all surrounding atoms
within **8 Angstrom** (currently default but can be modified later by user input).
Additionally we also provide the option to compute the same thing for a residue
base. 

**Exemplary case:**

1)
**all atoms computation:**
VAL has 7 Atoms (not counting hydrogens) and we compute for each of those 7 atoms
all surrounding atoms within 8A vicinity. We retrieve all parent residues in a list where each of those 7 atoms
returns a list of surrounding residues.
In total this will return **7 lists** where each list correspond to the **vicinity of 1 atom** in Val.

2)
**single_residue_based computation**:
VAL is treated as a single residue and we compute all residues surrounding it.
This computation is executed for each atom in VAL but only unique surrounding residues will be kept in a list.
The return will be **a single list** containing all residues that neighbour
 **VAL in any of its 7 atoms**.

**current date: 02-02-2023 [MM-DD-YYYY] author: Michael**
