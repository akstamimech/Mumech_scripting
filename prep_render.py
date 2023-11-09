import MDAnalysis as mda
import numpy as np
from pathlib import Path
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")

CHAIN = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
KAP_BINDING_SITES = [94, 97, 98, 101, 102, 134, 138, 185, 189, 224, 226, 227, 228,
  230, 262, 264, 265, 268, 269, 272, 324, 327, 331, 334, 373, 377, 385, 453,
  456, 457, 460, 494, 495, 498, 599, 600, 603, 636, 639, 677, 679, 682, 684,
  686, 699, 703, 705, 717, 719, 720, 727, 730, 731, 734, 750]

TPR_FILE = Path('droplet_md.tpr')


def main():
    # select file
    if len(sys.argv) >= 2:
        pdb_in = Path(sys.argv[1])
    else:
        exit('ERROR: Specify PDB file!')
    
    if not pdb_in.is_file():
        exit(f'ERROR: {pdb_in.name} not found!')

    # unwrap coordinates
    pdb_whole = fix_broken_molecules(pdb_in)

    # add chainIDs
    prepare_pdb_file(pdb_whole)
    pdb_whole.unlink()


def fix_broken_molecules(pdb: Path) -> Path:
    # check if tpr file exists
    if TPR_FILE.is_file():
        tpr = TPR_FILE
    else:
        tpr = select_tpr_file()

    # generate GMX call
    pdb_out = Path(pdb.stem + '_whole.pdb')
    gmx_call = [
        'gmx', 'trjconv',
        '-f', pdb,
        '-s', tpr,
        '-pbc', 'whole',
        '-o', pdb_out
    ]

    # run gmx call
    p = subprocess.Popen(gmx_call, stdin=subprocess.PIPE)
    p.communicate(b'0\n')
    p.wait()
    return pdb_out


def select_tpr_file() -> Path:
    print(f'-- Select TPR file --')
    files = list(Path.cwd().rglob("*.tpr"))
    if len(files) == 0:
        exit(f'ERROR: No .tpr files found!')
    elif len(files) == 1:
        fileID = 0
    else:
        [print(f'[{i}] {file.name}') for i, file in enumerate(files)]
        fileID = int(input())

    if fileID >= len(files):
        exit('ERROR: Wrong file ID given!')

    file = files[fileID]
    print(f'Selected the following file: {file.name}\n')
    return file


def prepare_pdb_file(pdb: Path) -> None:
    print('Preparing PDB file...')
    u = mda.Universe('forcefield/topol.top', pdb, topology_format='ITP')
    kap = u.atoms[-1].segid

    # Fix resids (?)
    for seg in u.segments:
        ctr = 1
        for res in seg.residues:
            res.resid = ctr
            ctr += 1

    # update atom properties
    chainIDs = []
    for atom in u.atoms:
        atom.name = 'CA'
        if atom.segid == kap:
            if atom.resid in KAP_BINDING_SITES:
                chainIDs.append('B')
            else:
                chainIDs.append('I')
        else:
            chainIDs.append('F')

    u.add_TopologyAttr('chainID', chainIDs)
    pdb_out = Path(pdb.stem[:-6] + '_CA.pdb')
    u.atoms.write(pdb_out)
    print(f'New PDB file written to {pdb_out.name!r}')


if __name__ == '__main__':
    main()
