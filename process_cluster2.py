import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import os, sys
import warnings
warnings.filterwarnings("ignore")

NUP = 'Nup116'
NMOL = 280
Kap = 'KAP95'

SIMULATION_DIR = Path.cwd() / f'{NUP}+{Kap}'

CLUSTER_DIR = 'cluster_ndx'
PDB_FILE = 'droplet_em.pdb'

GREEN = (0, 143/255, 146/255)
PINK = (1, 44/255, 121/255)
BLACK = 'dimgray'

UPDATE = False

def main():
	data_with_Nup116 = process_data_with_Nup116()

	fig = plt.figure() # figsize=(10, 3))
	total_time = plot_data(data_with_Nup116)
	add_labels(total_time)
	# plt.savefig(f'LLPS_{NUP}+Nsp1-{NMOL}.png', dpi=300)
	plt.show()


def process_data_with_Nup116():
	data_file = SIMULATION_DIR / 'cluster_data.npy'
	if UPDATE or not data_file.is_file():
		print('... PROCESSING FILES ...')
		u = mda.Universe(SIMULATION_DIR / PDB_FILE)
		chainIDs = ['A', 'B']		# distinguishes between FG-Nup and Nsp1

		# get chain lengths
		lengths = {}
		for chain in chainIDs:
			sel = u.select_atoms(f'chainID {chain}') #selects atomgroup
			lengths[chain] = len(sel.segments[0].atoms)

		counting_dict = {chain: [] for chain in chainIDs}
		nfiles = len(os.listdir(SIMULATION_DIR / CLUSTER_DIR))
		for ndx in range(nfiles):
			# read maxclust index file
			indices = read_ndx(SIMULATION_DIR / CLUSTER_DIR / f'maxclust_{ndx}.ndx')

			# initialize counter
			for chain in counting_dict.keys():
				counting_dict[chain].append(0)
				
			# process indices
			for i in indices:
				chain = u.atoms[i].chainID
				counting_dict[chain][ndx] += 1

			# get number of molecules
			for chain in counting_dict.keys():
				counting_dict[chain][ndx] = round(counting_dict[chain][ndx] / lengths[chain])

		freq = {'A': NMOL, 'B': NMOL // 2}	# Nup to Nsp1 ratio is 2 : 1
		data = np.array([np.array(counting_dict[chain]) / freq[chain] * 100 for chain in chainIDs])
		np.save(data_file, data)
		return data
	else:
		print('... LOADING DATA FROM FILE ...')
		return np.load(data_file)


def read_ndx(ndx: Path) -> list[int]:
	with open(ndx, 'r') as fid:
		fid.readline()  # skip first line
		return [int(line) - 1 for line in fid.readlines()]


def plot_data(data_with_Nup116):
	Nup = data_with_Nsp1[0]
	Nsp1 = data_with_Nsp1[1]

	plt.plot(0.010 * np.arange(len(Nup)), Nup, color=BLACK, label='Nup49FG')
	plt.plot(0.010 * np.arange(len(Nsp1)), Nsp1, color=PINK, label='Nsp1FG')

	# average with Nsp1
	average = np.mean(data_with_Nsp1[0, -500:])
	total_time = 0.010 * len(Nup)
	plt.axhline(average, ls='dashed', lw=0.8, color=BLACK)
	plt.text(x=total_time*0.01, y=average-1.4, s='avg. after adding Nsp1FG', fontsize=9, color=BLACK, ha='left', va='top')
	
	return total_time


def add_labels(total_time):
	# x-axis
	plt.xlim(0, total_time)
	# plt.xlim(0, 11000)
	plt.xlabel('time (µs)')

	# y-axis
	plt.ylim(0, 100)
	plt.ylabel(r'$N_\mathrm{cp}$ (%)')

	# other
	plt.legend(frameon=False)
	plt.tight_layout()


if __name__ == '__main__':
	main()
