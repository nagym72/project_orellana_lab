Compilazione:

gfortran ENM_NMA.f90 -o ENM_NMA -lblas -llapack
gfortran gen_movies_NM.f90 -o gen_movies_NM
gfortran write_CA.f90 -o write_CA
gfortran add_chain_label_ali.f90 -o add_chain_label_ali
gfortran MD_proj_NM.f90 -o MD_proj_NM

Execution

./ENM_NMA 5lox A 8        (name of PDB file + name of chain(s) + cutoff between different chains in Å)
./gen_movies_NM 5lox A 20           (name of PDB file for NMA + chain(s) + scale for motion movie in Å)
bash MD_simulation_projection_NM

Se la compilazione con llapack non funziona, usa comunque i normal modes pre-calcolati in "NMs.txt". Ti basta solo avere il file nello stesso folder in cui hai tutto il resto e in cui lanci l'analisi.

That's it!
Domenico
