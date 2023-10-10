#!/bin/bash

write_CA=./write_CA
add_chain_label_ali=./add_chain_label_ali
MD_proj_NM=./MD_proj_NM

ref="5lox"											#Name of reference PDB for NMA
chains_ref="A"										#Label of chain(s) for NMA and MD trajectory
name_traj="mdtrj_r1_wt_ca_fit_10frames.pdb"			#Name of MD trajectory file
NM1=1												#Normal mode number to be plotted on X-axis
NM2=2												#Normal mode number to be plotted on Y-axis

#Take only CA atoms from the ref struct
$write_CA $ref $chains_ref

#Read structures in the ensemble and align them to CA-only ref structure 
CA_pdb_suffix="_CA.pdb"
pdb_suffix=".pdb"
app_chains_ref="_$chains_ref"
ref_NM="$ref$app_chains_ref$CA_pdb_suffix"
n_CA=$(wc -l < $ref_NM)
n_MD_frames=0
flag_read_frame=0
while read p; do
	model_label="${p:0:5}"
	if test $model_label == "TER"
	then
		let flag_read_frame=0
	fi
	if test $flag_read_frame == 1
	then
		echo "$p" >> "frame_"$n_MD_frames".pdb"
	fi
	if test $model_label == "MODEL"
	then
		let n_MD_frames=n_MD_frames+1
		let flag_read_frame=1
		touch "frame_"$n_MD_frames".pdb"
	fi
done <$name_traj

for (( f=1; f<=$n_MD_frames; f++ ))
do 
	output_name="frame_"$f"_ali.pdb"
	yes 3 | gmx confrms -f1 $ref_NM -f2 "frame_"$f".pdb" -one -o temp_file.pdb
	grep "CA" temp_file.pdb > $output_name
	rm temp_file.pdb
	$add_chain_label_ali "frame_"$f ${output_name%"$pdb_suffix"}
	rm "frame_"$f".pdb"
done

touch aligned_MD_simulation.pdb
for (( f=1; f<=$n_MD_frames; f++ ))
do 
	echo "MODEL" >> aligned_MD_simulation.pdb
	cat "frame_"$f"_ali_labeled.pdb" >> aligned_MD_simulation.pdb
	echo "TER" >> aligned_MD_simulation.pdb
	echo "ENDMDL" >> aligned_MD_simulation.pdb
done

touch MD_ensemble_proj.txt
for (( f=1; f<=$n_MD_frames; f++ ))
do 
	$MD_proj_NM $n_CA "frame_"$f"_ali_labeled.pdb" $ref_NM
	rm "frame_"$f"_ali_labeled.pdb"
	cat "proj_frame_"$f"_ali_labeled.pdb.txt" >> MD_ensemble_proj.txt
	rm "proj_frame_"$f"_ali_labeled.pdb.txt"
done

bash plot_NM_space.sh $NM1 $NM2