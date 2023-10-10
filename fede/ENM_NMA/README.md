# Elastic Network Model + Normal Modes Analaysis

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
`bash trajectory_projection_on_ENM_NM.sh trajectory.pdb nameofchoice`
