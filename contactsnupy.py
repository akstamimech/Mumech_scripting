import MDAnalysis as mda
import MDAnalysis.analysis.distances as distances
import numpy as np
import re
from pathlib import Path
import warnings
from tqdm import tqdm 

warnings.filterwarnings("ignore")

CUTOFF = 0.7  # nm
NUPS = {
    'Nupy104': Path.cwd(),
}

def main():
    for nup, path in NUPS.items():
        u = load_universe(path)
        
        # print(f'# FG-motifs = {len(FG_beads)}')
        FGKapoccupancy = []
        FG_occupancy = []
        # for ts in u.trajectory[30000:40000:10]:
        counter = 0 
        for ts in u.trajectory:
            # print(int(ts.time), end='\r')
            u.trajectory[counter]
            FG_beads = select_FG_beads(u)
            Kap_beads = select_Kap_beads(u)
            # in_cluster = get_beads_in_cluster(ts.time, u, path)
            FG_in_cluster = FG_beads
            contact_matrix = distances.contact_matrix(FG_in_cluster.positions, cutoff=100*CUTOFF)
            FG_contacts = count_FG_occupancy(contact_matrix)
            FG_occupancy.append(FG_contacts / len(FG_beads + Kap_beads))
            kapnupmatrix = distances.distance_array(Kap_beads.positions, FG_in_cluster.positions)
            # print(kapnupmatrix)
            # kapnupmatrix = distances.contact_matrix(kapnupmatrix, cutoff = 10*CUTOFF)
            # print(kapnupmatrix)
            kapnupmatrix = replace_distances(kapnupmatrix)
            # print(kapnupmatrix)
            # print(kapnupmatrix)
            FGKAPcontacts = count_FG_occupancy(kapnupmatrix)
            print(FGKAPcontacts)
            FGKapoccupancy.append(FGKAPcontacts / len(FG_beads + Kap_beads))

            counter += 1

        avg = np.mean(np.array(FG_occupancy))
        std = np.std(np.array(FG_occupancy))
        print(f'average occupancy = {avg*100:.4f} % [pm {std*100:.4f} %]\n')

        avgkap = np.mean(np.array(FGKapoccupancy))
        stdkap = np.std(np.array(FGKapoccupancy))
        print(f'average Kap occupancy = {avgkap*100:.4f} % [pm {stdkap*100:.4f} %]\n')









def load_universe(path: Path) -> mda.Universe:
    print(f'Loading universe for {path.name}')
    return mda.Universe(path / 'droplet_init.pdb', (path / 'droplet_md.trr').as_posix())


def select_FG_beads(u: mda.Universe) -> mda.Universe:
    # print('Selecting FG beads')
    molecule = u.atoms.segments[0]
    sequence = molecule.residues.sequence(format='string')
    FG_string = ' '.join([str(m.start() + 1) for m in re.finditer('FG', sequence)])
    return u.select_atoms('resid ' + FG_string)

def get_beads_in_cluster(time, universe, path):
    ndx_file = path / 'index.ndx'
    in_cluster = []
    with open(ndx_file, 'r') as fid:
        fid.readline()                          # skip first line
        for res in fid.readlines():
            in_cluster.append(int(res) - 1)     # GROMACS ndx file is 1-based, MDAnalysis is 0-based
    return universe.atoms[in_cluster]


def count_FG_occupancy(contact_matrix):
    # print('Counting FG beads making FG-FG contact')
    FG_contacts = 0
    for col in range(len(contact_matrix)):
        if np.count_nonzero(contact_matrix[col]) >= 2:  # excluding self-interaction
            FG_contacts += 1
    return FG_contacts


def select_Kap_beads(u: mda.Universe) -> mda.Universe: 
     # molecule = u.atoms.segments[0]
     # sequence = molecule.residues.sequence(format='string')
     # FG_string = ' '.join([str(m.start() + 1) for m in re.finditer('FG', sequence)])
     return u.select_atoms('name BX')


def replace_distances(matrix):
    # Create a boolean matrix based on the distance condition
    boolean_matrix = matrix <= (100 * CUTOFF)
    return boolean_matrix

if __name__ == '__main__':
    main()
# , path / 'droplet_eq.trr'