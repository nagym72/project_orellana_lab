#!/bin/bash

#show the user how to use my script
usage() { echo -e "If you require help and information about the purpose of this script type:\n\n$0 -help\n"
	  echo -e "To use this script:\n\n"$0" [ -f <filepath> ] or [ -d <dirpath> ]\n" 1>&2; exit 1; }

while getopts ":f:d:h" fpath; do
	case "${fpath}" in
		f)
			filep="${OPTARG}"
			#check if supplied file is really a file; else exit
			if [ -f "$filep" ];
				then
					echo ""
			else
				echo "File not found"
				usage
			fi
			;;
		d)
			dirp="${OPTARG}"
			#check if supplied file is really a dir; else exit
			if [ -d "$dirp" ];
				then
					echo ""
			else
				echo "Dir not found"
				usage
			fi
			;;
		h)
			mkdownhelp="neighbouring_residue_help.md"
			#read in the description on how to use this file
			while IFS= read -r lines
			do
			echo "$lines"
			done <"$mkdownhelp"
			exit
			;;
		*)
			echo "we are here"
			usage;;
	esac
done

#check if user supplied input

if [ -z "${filep}" ] && [ -z "${dirp}" ];
then
	usage
	echo ""
fi

if [ ! -z "${filep}" ];
then
	echo "File selected: "${filep}""
	echo ""

elif [ ! -z "${dirp}" ];
then
	echo "Dir selected: "${dirp}""
	echo ""
fi


#script requires biopython.
read -rep $'For the Neighbouring script we require Biopython. Shall we install it? [Y/N]\n' answ

#install if the user does not have it already through pip
case "$answ" in
	Y | y)
	printf -- "-%.0s" {1..80}
	pip install biopython
	;;
	N | n)
	echo "If you already have it, please try again and select Y which will work"
	echo "In case you dont want to install it, we cant proceed"
	exit
	;;
	*)
	echo " "$answ" is not a valid input. closing script."
esac


#lets go for the python part
printf -- "-%.0s" {1..80}

if [ ! -d  results_surrounding_residues ];
then
	mkdir results_surrounding_residues; touch "defaultprotname.csv"
elif [ -d results_surrounding_residues ] && [ ! -f "defaultprotname.csv" ];
then
	touch "./results_surrounding_residues/defaultprotname.csv"
fi

python compute_neighbouring_residues.py "$filep" "defaultprotname" "./results_surrounding_residues/"
