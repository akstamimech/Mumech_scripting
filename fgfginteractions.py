import MDAnalysis as mda
import MDAnalysis.analysis.distances as distances
import numpy as np
import re
from pathlib import Path
import warnings
import subprocess, os
from tqdm import tqdm
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

CUTOFF = 0.7  # nm
NUPS = {
    'Nup116': Path.cwd()
}


def main():
    for nup, path in NUPS.items():
    	FG_occupancy_list = []
    	for idx, pdb_file in tqdm(enumerate(sorted(os.listdir('clusteronly/')))): #mdanalysis universe set, instead of using trajectory, interate over clusteronly pdbs! they have a use now.
            clusterfile = f'clusteronly/clusteronly{idx}_max_clust.pdb'
            u = mda.Universe(clusterfile)
            FG_beads = select_FG_beads(u) #selected specifically FG sequences
            # print(f'# FG-motifs = {len(FG_beads)}')

            FG_occupancy = []
            ndx_file = path / f'cluster_ndx/maxclust_{idx}.ndx'
            in_cluster = u.select_atoms('all') #reading ndx is unnecesary since pdb only contains cluster. We are making this selection just to show the method, or if we want to edit this in the future. 
            FG_in_cluster = in_cluster.intersection(FG_beads)
            # print(f'intersection = {len(FG_in_cluster)}')

            contact_matrix = distances.contact_matrix(FG_in_cluster.positions, cutoff=10*CUTOFF)
            FG_contacts = count_FG_occupancy(contact_matrix)
            FG_occupancy.append(FG_contacts / len(FG_in_cluster))

            if idx == 10: 
                print(FG_occupancy)
            FG_occupancy = FG_occupancy[0]

            FG_occupancy_list.append(FG_occupancy) #what does this value mean?


    	avg = np.mean(np.array(FG_occupancy_list))
    	std = np.std(np.array(FG_occupancy_list))
    	print(FG_occupancy_list)
    	print(f'average occupancy = {avg*100:.4f} % [pm {100*std:.4f} %]\n')


    fig = plt.figure()
    plt.plot([i*10 for i in range(len(FG_occupancy_list))], FG_occupancy_list, label = 'occupancy through time (ns)')
    plt.axvline(x=2500, color = 'green', linestyle = 'dotted')
    plt.xlabel('Time (ns)')
    plt.ylabel('fg-fg interactions (%)')
    plt.legend()
    plt.show()







def count_FG_occupancy(contact_matrix):
    # print('Counting FG beads making FG-FG contact')
    FG_contacts = 0
    for col in range(len(contact_matrix)):
        if np.count_nonzero(contact_matrix[col]) >= 2:  # excluding self-interaction
            FG_contacts += 1
    return FG_contacts

def select_FG_beads(u: mda.Universe) -> mda.Universe:
    # print('Selecting FG beads')
    molecule = u.atoms.segments[0]
    sequence = molecule.residues.sequence(format='string')
    FG_string = ' '.join([str(m.start() + 1) for m in re.finditer('FG', sequence)])
    return u.select_atoms('resid ' + FG_string)



def get_beads_in_cluster(universe, path):
    in_cluster = []
    ndx_file = path / f'cluster_ndx/maxclust_{idx}.ndx'
    with open(ndx_file, 'r') as fid:
        fid.readline()                          # skip first line
        for res in fid.readlines():
            in_cluster.append(int(res) - 1)     # GROMACS ndx file is 1-based, MDAnalysis is 0-based
    return universe.atoms[in_cluster]



if __name__ == '__main__':
    main()