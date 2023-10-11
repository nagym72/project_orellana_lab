#!/usr/bin/env python3
# coding: utf-8

# Checks whether a pdb has gaps or not
# Arg 1 is input pdb structure

import pandas as pd
import numpy as np
import sys

with open (sys.argv[1], 'r') as fin:
#with open ('', 'r') as fin:
#filein= input ('enter PDB file pathname:\n')
#with open (filein, 'r') as fin:
    wholepdb=fin
    
    colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
            (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),
            (78, 80)]

    names = ['ATOM', 'serial', 'name', 'altloc', 'resname', 'chainid', 'resseq',
         'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge']

    pdb = pd.read_fwf(wholepdb, names=names, colspecs=colspecs)
    
    pdb = pdb.drop(['serial', 'resname', 'chainid',
         'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge'], axis=1)
    pdb = pdb.loc[pdb['ATOM'] == 'ATOM']
    pdb = pdb.loc[pdb['name'] == 'CA']
    pdb = pdb.loc[(pdb['altloc'].isna()) | (pdb['altloc']== 'A')]
    pdb = pdb.set_index('resseq')
    pdb.index = pdb.index.astype(int)


res=pdb.index.to_numpy()
lowerBounds = (res+1)[:-1]
upperBounds = (res-1)[1:]
mask = lowerBounds<=upperBounds
upperBounds, lowerBounds = upperBounds[mask], lowerBounds[mask]
ub, lb= tuple(upperBounds), tuple(lowerBounds)
gaps, u, l=[], [], []
c=0
for ele in lb:
    interval=[]
    couple= '{}:{}'.format(lb[c],ub[c])
    gaps.append(couple)
    c += 1

if not a:
    print("No internal gaps")

else:
    print ('NB: Gap bounds are absent residues\n', *gaps, sep='\n')
    
    
    

