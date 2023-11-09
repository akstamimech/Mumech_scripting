""" Add segid labels to PDB file. Useful for visualizing in VMD. """

import os, sys
CHAIN = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

if len(sys.argv) >= 2:
	INPUT = sys.argv[1]
else:
	exit('ERROR: Specify PDB file!')

if not os.path.isfile(INPUT):
	exit(f'ERROR: File {INPUT!r} does not exist!')

if not INPUT.endswith('.pdb'):
	exit('ERROR: Expecting PDB file!')

# Read molecule list --> create chain ID list
freq = []
read_mol = False
with open('forcefield/topol.top') as fid:
	for line in fid.readlines():
		if read_mol:
			freq.append(int(line.split()[1]))
		if '#mols' in line[:-1]:
			read_mol = True
chain_list = [CHAIN[idx] for idx in range(len(freq)) for mol in range(freq[idx]) ]

# Update PDB file
with open(INPUT, 'r') as fid, open('tmp.pdb', 'w+') as out:
	counter = 0
	for line in fid.readlines():
		if line.startswith('ATOM'):
			if int(line[22:26]) == 1:
				counter += 1
			out.write(f'{line[:21]}{chain_list[counter - 1]}{line[22:72]}M{counter:03d}\n')
		else:
			out.write(line)

os.system(f'mv tmp.pdb {INPUT}')
