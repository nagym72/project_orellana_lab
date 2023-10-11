#!/bin/bash

nm1=$1
nm2=$2
underscore="_"

gnuplot << EOF

    set key outside
    set key right top 
    set terminal postscript eps solid enhanced color "Helvetica" 16     
    set xlabel "Normal mode $nm1"    
    set ylabel "Normal mode $nm2"

    set output "NM_plot_$nm1$underscore$nm2.jpeg"   
    plot \
    "MD_ensemble_proj.txt" using $nm1:$nm2 notitle with points pt 5 pointsize 1.5 lt -1, \

    replot
    set terminal jpeg
    set output "NM_plot_$nm1$underscore$nm2.jpeg"
    replot