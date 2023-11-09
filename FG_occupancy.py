import MDAnalysis as mda
import MDAnalysis.analysis.distances as distances
import numpy as np
import re
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

CUTOFF = 0.7  # nm
#try adding kap yaml next
NUPS = {
    'Nup116': Path('../../FG_NUP_droplets/NUP116_YEAST-1_725-4.1us')
}


def main():
    for nup, path in NUPS.items():
        u = load_universe(path)
        FG_beads = select_FG_beads(u)
        print(f'# FG-motifs = {len(FG_beads)}')

        FG_occupancy = []
        # for ts in u.trajectory[30000:40000:10]:
        for ts in u.trajectory[20000:40000:100]:
            print(int(ts.time / 1000), end='\r')
            in_cluster = get_beads_in_cluster(ts.time, u, path)
            FG_in_cluster = in_cluster.intersection(FG_beads)
            contact_matrix = distances.contact_matrix(FG_in_cluster.positions, cutoff=10*CUTOFF)
            FG_contacts = count_FG_occupancy(contact_matrix)
            FG_occupancy.append(FG_contacts / len(FG_in_cluster))

        avg = np.mean(np.array(FG_occupancy))
        std = np.std(np.array(FG_occupancy))
        print(f'average occupancy = {avg*100:.4f} % [pm {100*std:.4f} %]\n')









def load_universe(path: Path) -> mda.Universe:
    print(f'Loading universe for {path.name}')
    return mda.Universe(path / 'droplet_init.pdb', (path / 'droplet_eq.trr').as_posix()) #returns string representation 


def select_FG_beads(u: mda.Universe) -> mda.Universe:
    print('Selecting FG beads')
    molecule = u.atoms.segments[0]
    sequence = molecule.residues.sequence(format='string')
    FG_string = ' '.join([str(m.start() + 1) for m in re.finditer('FG', sequence)])
    return u.select_atoms('resid ' + FG_string)

def get_beads_in_cluster(time, universe, path):
    ndx_file = path / f'cluster_ndx/maxclust_{time/10000:.0f}.ndx'
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


if __name__ == '__main__':
    main()
# , path / 'droplet_eq.trr'