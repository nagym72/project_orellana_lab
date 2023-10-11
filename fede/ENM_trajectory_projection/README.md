# Elastic Network Model + Normal Modes Trajectory Projection

A bunch of scripts from Domenico which was partly generalized by me (Fede). 
The whole package is launched from the .sh script in this directory.

## Function
Computes ENM normal modes of a structure, projects the input trajectory on the first two normal modes and plots.
Creates a temporary directory in the trajectory folder, copies all the required files in it, calls precompiled fortran codes, creates a shitload of temporary files and deletes the non necessary ones, moves the directory to the final destination.

### Argument inputs: 
==multipdb== trajectory pathname, output name

### Hardcoded inputs: 
5lox structure, ENM_NMA scripts directory, output (parent) directory, chain, normal modes 1 and 2

### Outputs:
NMs, frequencies, projection of trajectory on normal modes

## Manual
The whole package is launched from the .sh script in this directory, like this:
`bash trajectory_projection_on_ENM_NM.sh trajectory.pdb outputname`

## Warnings:
- Fortran scripts in the packages must be previously re-compiled on your machine (not coded in trajectory_projection_on_ENM_NM.sh)

> gfortran ENM_NMA.f90 -o ENM_NMA -lblas -llapack
> gfortran gen_movies_NM.f90 -o gen_movies_NM
> gfortran write_CA.f90 -o write_CA
> gfortran add_chain_label_ali.f90 -o add_chain_label_ali
> gfortran MD_proj_NM.f90 -o MD_proj_NM

> note that lblas and llapack libraries are required
> the execution of ENM is coded in trajectory_projection_on_ENM_NM.sh, but if you want to run it separately, this is an example:

>> ./ENM_NMA 5lox A 8        (name of PDB file + name of chain(s) + cutoff between different chains in Å)
>> ./gen_movies_NM 5lox A 20           (name of PDB file for NMA + chain(s) + scale for motion movie in Å)


- Executables only work with all input files in the same folder, but the file management is done upstream by trajectory_projection_on_ENM_NM.sh
- Fortran executables only accept four-digit names for input structure pdbs, without extensions (e.g. 5lox) 
- This works for single chain and SHOULD support multiple chain (see above) but ask Domenico about it
