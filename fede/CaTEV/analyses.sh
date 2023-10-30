#!/bin/bash

# NB hardcoded for two replicas

workdir=/Volumes/ExFAT/Lab/CaTEV

init="step3_input"
mini_prefix="step4.0_minimization"
equi_prefix="step4.1_equilibration"
prod_prefix="step5_production"



# NON-CALCIUM TRAJECTORIES
dirs=("CKK" "M13")

for dir in "${dirs[@]}"; do

	reps=("r1" "r2")
	mergedir="$workdir/$dir/merged"
	mkdir -p $mergedir

	for rep in "${reps[@]}"; do
	    repdir="$workdir/$dir/$rep"
	    cd "$repdir"

	    # protein-only structure for VMD visualization
		echo "1" | gmx editconf -f "${init}.gro" -o "mdstr_${init}_prot.gro" -ndef

		# PROCESS TRAJECTORY
		echo "1" | gmx trjconv -f "${prod_prefix}.xtc" -s "${equi_prefix}.gro" -o "mdtrj_${rep}_prot_cut.xtc" -pbc whole -dt 1000

		echo "1" | gmx trjconv -f "mdtrj_${rep}_prot_cut.xtc" -s "${equi_prefix}.gro" -o "mdtrj_${rep}_prot_nojump.xtc" -pbc nojump

		echo "1 1" | gmx trjconv -f "mdtrj_${rep}_prot_nojump.xtc" -s "${equi_prefix}.gro" -o "mdtrj_${rep}_prot_cent.xtc" -pbc mol -center

		echo "1 1" | gmx trjconv -f "mdtrj_${rep}_prot_cent.xtc" -s "${equi_prefix}.gro" -o "mdtrj_${rep}_prot_fit.xtc" -fit rot+trans

		rm "mdtrj_${rep}_prot_cut.xtc" "mdtrj_${rep}_prot_nojump.xtc" "mdtrj_${rep}_prot_cent.xtc"
		cp "mdtrj_${rep}_prot_fit.xtc" $mergedir/

		# RMSD RMSF

		echo "4 4" | gmx rms -s "${init}.gro" -f "mdtrj_${rep}_prot_fit.xtc" -o "rmsd_${rep}_init_bb.xvg" -tu ns

		echo "3" | gmx rmsf -s "${init}.gro" -f "mdtrj_${rep}_prot_fit.xtc" -o "rmsf_${rep}_init_ca.xvg" -res

	done

	# MERGED
	rep1="${reps[0]}"
	rep2="${reps[1]}"

	cd $mergedir
	cp ../$rep1/{"${init}.gro","${prod_prefix}.gro","mdstr_${init}_prot.gro"} .

	# Merge
	gmx trjcat -f "mdtrj_${rep1}_prot_fit.xtc" "mdtrj_${rep2}_prot_fit.xtc" -o mdtrj_merged_prot_fit.xtc -cat
	rm "mdtrj_${rep1}_prot_fit.xtc" "mdtrj_${rep2}_prot_fit.xtc"

	# PDB CA-only trajectory for Elastic Network Model Normal Mode Analysis
	
	echo "3" | gmx trjconv -f mdtrj_merged_prot_fit.xtc -s "${init}.gro" -o "mdtrj_merged_ca_fit.pdb"

	# PCA input files
	
	echo "3 3" | gmx covar -f "mdtrj_merged_prot_fit.xtc" -s "${prod_prefix}.gro" -o "eigenval_merged_ca.xvg" -tu ns

	echo "3 3" | gmx anaeig -v eigenvec.trr -f "mdtrj_merged_prot_fit.xtc" -s "${prod_prefix}.gro" -comp "eigcomp_merged.xvg" \
	-rmsf "eigrmsf_merged.xvg" -2d "2d_merged.xvg" -tu ns -first 1 -last 2 -filt "filtered_merged.xtc" -extr "extreme_merged_.pdb" \
	-nframes 10

	

done

mv saltbr* ./salt_bridges
ls ./salt_bridges > ls_saltbridges_${rep}.txt

vmd mdstr_step3_input_prot.gro mdtrj_merged_prot_fit.xtc -e /Users/federicopozzani/Desktop/Lab/Scripts/vmd_saltbridges_hbonds_analysis.tcl -dispdev text