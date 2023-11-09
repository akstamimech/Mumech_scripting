#radial_density
import subprocess, os
import numpy as np 
import matplotlib.pyplot as plt 
import csv
import MDAnalysis as mda
from pathlib import Path
from tqdm import tqdm
import pathlib
import warnings 
warnings.filterwarnings("ignore")

#self-explanatory changeable settings
CUTOFF = 0.7  # nm
START_TIME = 2500  # ns
TIMESTEP = 1  # ns

PDB = 'droplet_init.pdb'
TRR = 'droplet_md.trr'
XTC = 'droplet_md.xtc'
TPR = 'droplet_md.tpr'
FRAMES_PDB = 'droplet_dens-frames.pdb'
FRAMES_DIR = 'frames_centered'

TAG = 'v4'


MASSES = {'ALA' :  71.08, 'ARG' : 156.19, 'ASN' : 114.10, 'ASP' : 115.09, 
 		  'CYS' : 103.14, 'GLN' : 128.13, 'GLU' : 129.12, 'GLY' :  57.05,
 		  'HIS' : 137.14, 'ILE' : 113.16, 'LEU' : 113.16, 'LYS' : 128.17, 
 		  'MET' : 131.19, 'PHE' : 147.18, 'PRO' :  97.12, 'SER' :  87.08,
 		  'THR' : 101.11, 'TRP' : 186.21, 'TYR' : 163.18, 'VAL' :  99.13,
 		  }


fig = plt.figure()
def main():
	# preprocessing
	shrink_trajectory()
	# frames_list = [_, pdb in enumerate(sorted(os.listdir('frames_centered/')))]
	# print(frames_list)
	output_file = initialize_output()
	output_file_kaps = initialize_output_kaps()
	output_file_nups = initialize_output_nups()
	#creating density matrix..?
	density_matrix = []
	density_matrix_kaps = []
	density_matrix_nups = []
	graphcounter = 0 
	# avedenslist = []
	# avedenskaplist = []
	# avedensnuplist = []
	for idx, pdb in tqdm(enumerate(sorted(os.listdir('frames_centered/')))):
		print(pdb) #frame00
		frame_centered = center_frame(pdb)
		os.system(f'python fix_PDB_file.py {frame_centered}')
		u = load_trajectory(frame_centered)
		
		# Determine radial density using spherical shells (from center of geometry of cluster)
		dens = []
		densnup = []
		denskap = []
		box_dim = u.dimensions[0] / 10 	# box size (nm)
		for idx in range(int(box_dim / 2)):
			START = idx * 10
			END   = (idx + 1) * 10

			# volume of the shell
			vol = (4 / 3) * np.pi * (END ** 3 - START ** 3) * 1e-3  # nm ** 3

			# select atoms in shell
			shell = u.select_atoms(f'sphlayer {START} {END} ( all )')
			# if len(shell) == 0: 
			# 	continue
			shellnups = shell.select_atoms(f'chainID A')
			# if len(shellnups) == 0: 
			# 	continue
			shellkaps = shell.select_atoms(f'chainID B')
			# if len(shellkaps) == 0: 
			# 	continue

			# total mass of atoms in shell
			if len(shell) == 0: 
				mass = 0
			else: 
				mass = shell.total_mass() * 1.66053907  # convert to mg
			if len(shellnups) == 0:
				massnups = 0
			else: 
				massnups = shellnups.total_mass() * 1.66053907
			if len(shellkaps) == 0:
				masskaps = 0
			else:
				masskaps = shellkaps.total_mass() * 1.66053907

			# density in shell (mg / mL)
			dens.append(mass / vol) 		# (mg / mL)
			densnup.append(massnups / vol)
			denskap.append(masskaps / vol)
		density_matrix.append(np.array(dens))		# (mg / mL)
		density_matrix_kaps.append(np.array(denskap))
		density_matrix_nups.append(np.array(densnup))
		avedens = np.mean(dens)
		avedenskap = np.mean(denskap)
		avedensnup = np.mean(densnup)
		

		# Write to file
		with open(output_file, 'a+') as out:
			out.write(f"{pdb.split('.')[0]} ")
			for sh in range(len(dens)):
				out.write(f'{dens[sh]:13.8f} ')
			out.write('\n')

		with open(output_file_nups, 'a+') as out:
			out.write(f"{pdb.split('.')[0]} ")
			for sh in range(len(dens)):
				out.write(f'{densnup[sh]:13.8f} ')
			out.write('\n')

		with open(output_file_kaps, 'a+') as out:
			out.write(f"{pdb.split('.')[0]} ")
			for sh in range(len(dens)):
				out.write(f'{denskap[sh]:13.8f} ')
			out.write('\n')

		# if idx == 0:
		# 	plt.plot(range(len(denskap)), denskap, label = 'kap density at 0ns')
		# 	plt.plot(range(len(densnup)), densnup, label = 'nup density at 0ns')

		graphcounter += 1
		# if graphcounter % 1 == 0:
		# 	# plt.plot(range(len(denskap)), denskap, label = f'kap density at {graphcounter}ns')
		# 	# plt.plot(range(len(densnup)), densnup, label = f'nup density at {graphcounter}ns')


	avedenslist = np.mean(density_matrix, axis=0)
	avedenskaplist = np.mean(density_matrix_kaps, axis=0)
	avedensnuplist = np.mean(density_matrix_nups, axis=0)

	plt.plot(range(len(avedenslist)), avedenslist, label = ' average total density')
	plt.plot(range(len(avedenskaplist)), avedenskaplist, label = 'total kap density')
	plt.plot(range(len(avedensnuplist)), avedensnuplist, label = 'average nup density')
	plt.ylabel('density.(mg/ml)')	
	plt.ylim(None, 400)
	plt.xlabel('radial distance(nm)')
	plt.legend()
	plt.title('average densitys for 2500 - 3500 ns')
	plt.savefig("radplotave.png")











def shrink_trajectory(): #shrink and split combined in one, creates all the frames
	# gmx_call = [
	# 	'gmx', 'traj',
	# 	'-f', XTC,
	# 	'-s', PDB,
	# 	'-n', 'index.ndx',
	# 	'-tu', 'ns',
	# 	'-dt', str(TIMESTEP),
	# 	'-b', str(START_TIME),
	# 	'-oxt', FRAMES_PDB
	# ]
	# p0 = subprocess.Popen(gmx_call, stdin=subprocess.PIPE)
	# p0.communicate(b'0\n')
	# p0.wait()

	os.system(f'rm -r {FRAMES_DIR}      &> /dev/null')
	os.system(f'mkdir {FRAMES_DIR}')

	p0 = subprocess.Popen(
		[
			'gmx', 'trjconv',
			'-f', XTC,
			'-s', TPR,
			'-n', 'index.ndx',
			'-b', str(START_TIME),
			'-e', str(3500),
			'-tu', 'ns',
			'-sep',
			'-skip', str(TIMESTEP * 10),  # assuming nstxout = 5000
			'-o', f'./{FRAMES_DIR}/frame.pdb'
		],
		stdin=subprocess.PIPE
	)
	p0.communicate(b'0\n')
	p0.wait()




def initialize_output():
	# Initialize data output file
	matrix_out = pathlib.Path(f'density_matrix-{TAG}.out')
	with open(matrix_out, 'w+') as out:
		out.write('time (ns)   radial density (mg / mL) [using shell of 1 nm]\n')
	return matrix_out



def initialize_output_nups(): 
	# Initialize data output file
	matrix_out_nups = pathlib.Path(f'density_matrix-nups-{TAG}.out')
	with open(matrix_out_nups, 'w+') as out:
		out.write('time (ns)   radial density (mg / mL) [using shell of 1 nm]\n')
	return matrix_out_nups


def initialize_output_kaps(): 
	# Initialize data output file
	matrix_out_kaps = pathlib.Path(f'density_matrix-kaps-{TAG}.out')
	with open(matrix_out_kaps, 'w+') as out:
		out.write('time (ns)   radial density (mg / mL) [using shell of 1 nm]\n')
	return matrix_out_kaps


# def split_trajectory() -> list[str]:
# 	# Determine number of beads
# 	with open(PDB, 'r') as fid:
# 		beads = len(fid.readlines()) - 6
# 		print(f'{beads = }')

# 	# Split frames into separate files and place in frames dir.
# 	os.system(f'rm -r {FRAMES_DIR}      &> /dev/null')
# 	os.system(f'mkdir {FRAMES_DIR}')
# 	# os.system(f'split -l {beads + 6} {FRAMES_PDB} ./{FRAMES_DIR}/frame -d --additional-suffix=.pdb')
# 	# Split the trajectory frames into individual PDBs
# 	p2 = subprocess.Popen(
# 		[
# 			'gmx', 'trjconv',
# 			'-f', trr_file,
# 			'-s', tpr_file,
# 			'-n', index.ndx,
# 			'-sep',
# 			'-skip', str(dt * 10),  # assuming nstxout = 5000
# 			'-nzero', str(nzero),
# 			'-o', frames_dir / 'frame.pdb'
# 		],
# 		stdin=subprocess.PIPE
# 	)
# 	p2.communicate(b'0\n')
# 	p2.wait()
# 	return os.listdir(FRAMES_DIR)



def center_frame(pdb_file): #needs pdb file created inside frames_dir - frames_centered centerframe(frame)
	
	out_ndx = 'maxclust.ndx'
	if os.path.isfile('centered.pdb'):
		os.remove('centered.pdb')

	# get indices of cluster
	gmx_call = ['gmx', 'clustsize', '-f', f'./frames_centered/{pdb_file}',
				'-s', TPR, '-mcn', 'maxclust.ndx', '-mol', '-cut', str(CUTOFF)] #refuses to create maxclust.ndx.. why? 
	p0 = subprocess.Popen(gmx_call)#, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
	p0.wait()
	

	# add index group
	tmp_ndx = 'index_cluster.ndx'
	os.system(f'cat index.ndx maxclust.ndx > {tmp_ndx}')

	# center cluster
	tmp_pdb = 'tmp_clustered.pdb'
	gmx_call = ['gmx', 'trjconv', '-f', f'{FRAMES_DIR}/{pdb_file}',
				'-s', TPR, '-n', tmp_ndx, '-pbc', 'cluster',
				'-center', '-o', tmp_pdb]
	p1 = subprocess.Popen(gmx_call, stdin=subprocess.PIPE) #, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
	p1.communicate(b'7\n7\n0\n')
	p1.wait()

	# put all residues back into the box
	gmx_call = ['gmx', 'trjconv', '-f', tmp_pdb, '-s', TPR, '-n', tmp_ndx, 
				'-pbc', 'atom', '-o', 'centered.pdb']
	p2 = subprocess.Popen(gmx_call, stdin=subprocess.PIPE) #, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
	p2.communicate(b'0\n')
	p2.wait()

	out_pdb = 'centered.pdb'

	# clean up
	temporary_output = [
			"csize.xpm", "csizew.xpm", "nclust.xvg", "maxclust.xvg", out_ndx,
			"avclust.xvg", "histo-clust.xvg", "temp.xvg", tmp_ndx, tmp_pdb
		]
	for filename in temporary_output:
		if os.path.isfile(filename):
			os.remove(filename)
	return out_pdb


def load_trajectory(pdb_file):
	universe = mda.Universe(pdb_file)
	# Add bead masses
	for atom in universe.atoms:
		atom.mass = MASSES.get(atom.resname)
	return universe



if __name__ == '__main__':
	main()
