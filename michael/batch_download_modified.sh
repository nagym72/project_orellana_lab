#!/bin/bash

# Script to download files from RCSB http file download services.
# Use the -h switch to get help on usage.

if ! command -v curl &> /dev/null
then
    echo "'curl' could not be found. You need to install 'curl' for this script to work."
    exit 1
fi

PROGNAME=$0
BASE_URL="https://files.rcsb.org/download"

usage() {
  cat << EOF >&2
Usage: $PROGNAME -f <file> [-o <dir>] [-c] [-p]

 -f <file>: the input file containing a comma-separated list of PDB ids
 -o  <dir>: the output dir, default: current dir
 -c       : download a cif.gz file for each PDB id
 -p       : download a pdb.gz file for each PDB id (not available for large structures)
 -a       : download a pdb1.gz file (1st bioassembly) for each PDB id (not available for large structures)
 -x       : download a xml.gz file for each PDB id
 -s       : download a sf.cif.gz file for each PDB id (diffraction only)
 -m       : download a mr.gz file for each PDB id (NMR only)
 -r       : download a mr.str.gz for each PDB id (NMR only)
 -k	  : specify chain for output name for later downstream processing (modified michael)
EOF
  exit 1
}

download() {
  retrieve_name="${1:0:4}"".pdb" #i changed it to pdb1! before it was pdb!
  #echo "This is retrieve name:"
  #echo "$retrieve_name"
  url="$BASE_URL/$retrieve_name"       # we need 5ltu.pdb as url target. 
  out=$2/$1      
  #chain="$3" #my mod
  #newout="${1:0:4}" #my mod this is the name e.g 5ltu
  #ending="${1:4:8}" #my mod  this is .pdb
  #new_name=""$newout""_""$chain""$ending"" #my mod
  #out=$2/$new_name #my mod
  #echo "$new_name"
  #echo "Downloading $url to $out"
  
  #if curl does not find a pdb this means it is a cif file! then echo this to stdout and make a new list out of those and rerun script with different paramter e.g -c instead -p   	
  curl -s -f $url -o $out || echo "$1"
  #gzip -d $out || echo "failed to unzip" $out
}

listfile=""
outdir="."
cif=false
pdb=false
pdb1=false
xml=false
sf=false
mr=false
mrstr=false
chain="" #my modification
while getopts f:o:k:cpaxsmr o
do
  case $o in
    (f) listfile=$OPTARG;;
    (o) outdir=$OPTARG;;
    (k) chain=$OPTARG;;  #also my modification
    (c) cif=true;;
    (p) pdb=true;;
    (a) pdb1=true;;
    (x) xml=true;;
    (s) sf=true;;
    (m) mr=true;;
    (r) mrstr=true;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"

#echo "$listfile"
#echo "$outdir"
#echo "$chain"

if [ "$listfile" == "" ]
then
  echo "Parameter -f must be provided"
  exit 1
fi
contents=$(cat $listfile)

# see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
IFS=',' read -ra tokens <<< "$contents"

for token in "${tokens[@]}"
do
  if [ "$cif" == true ]
  then
    download ${token}.cif $outdir
  fi
  if [ "$pdb" == true ]        #my modification
  then
    #echo "we are here"
    #echo "${token}"
    download ${token}.pdb $outdir
  fi
  if [ "$pdb1" == true ]
  then
    download ${token}.pdb1.gz $outdir
  fi
  if [ "$xml" == true ]
  then
    download ${token}.xml.gz $outdir
  fi
  if [ "$sf" == true ]
  then
    download ${token}-sf.cif.gz $outdir
  fi
  if [ "$mr" == true ]
  then
    download ${token}.mr.gz $outdir
  fi
  if [ "$mrstr" == true ]
  then
    download ${token}_mr.str.gz $outdir
  fi

done








