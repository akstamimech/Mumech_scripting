#possible solution to centered graph problems.. 

import subprocess, os
import numpy as np 
import matplotlib.pyplot as plt 
import csv
import MDAnalysis as mda
from pathlib import Path
from tqdm import tqdm



PDB_FILE = 'droplet_init.pdb'


current = Path.cwd()
ndx_dir = current / 'cluster_ndx'


for idx, ndx in tqdm(enumerate(sorted(os.listdir('./cluster_ndx/')))):
	out_ndx = ndx_dir / f'maxclust_{idx}.ndx'

	u = mda.Universe(PDB_FILE)
	chainIDs = ['A', 'B']
	with open(out_ndx, 'r') as clust:
		lines = clust.readlines()
		lines = [i.strip() for i in lines[1:]]
		lines = [eval(j) for j in lines[1:]]

	sel = u.select_atoms('all')
	sel2 = sel[lines]
	sel3 = sel2.select_atoms(f'chainID B')

	

