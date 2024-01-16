#This is an exemplary test script for pymol

import csv
import os
from pymol import cmd
from pymol.cgo import *
from math import *
import re

#Automation this time in example script with only EGFR testcase
pdb_sorted_renumbered_path = "/home/michael/bash_cosmic/renumb.pdb"
#load already downloaded and renumbered files from disc
cmd.load(pdb_sorted_renumbered_path)
#display sequence
cmd.set("seq_view",value=1)
#set all b factors to 0 which we will use later for mapping the amount of mutations per residue
cmd.alter("all", "b=0")

#define custom colors for gradient

cmd.set_color("egypt_blue",[16, 52, 166])
cmd.set_color("span_viol",[65, 47, 136])
cmd.set_color("byzantium",[114, 43, 106])
cmd.set_color("amaranth_purple",[162, 38, 75])
cmd.set_color("amaranth_red",[211, 33, 45])
cmd.set_color("deep_carmine_pink",[246, 45, 45])


color_palette = "white gray amaranth_purple amaranth_red " \
                "orange red"

fsavepath = "/home/michael/bash_cosmic/EGFR_b_factor_mapped_mutations.pdb"
fpath = "/home/michael/bash_cosmic/EGFR_mutation_counted_per_position.csv"
with open(fpath, "r") as fh:
    csv_reader = csv.reader(fh)
    for entries in csv_reader:
        #skip header
        if entries[0] == "Position" or entries[1] == "Num_mutations":
            continue
        #print(entries)
        cmd.alter("i. "+entries[0]+" and name CA", "b= entries[1]")
        cmd.alter("i. "+entries[0]+" and name CA", "vdw=4.5*vdw")
        cmd.show("spheres", "i. "+entries[0]+" and name CA")
        #suboptimal solution, size should be normalized based on min max normalization for all mutations
        cmd.set("sphere_scale", 0.1*float(entries[1]), selection='i. '+entries[0])
        #cmd.show("sticks", "i. "+entries[0])

cmd.spectrum("b", palette=color_palette)
cmd.set("sphere_scale", value=5)
#cmd.show("surface", "all")

cmd.save(fsavepath)
