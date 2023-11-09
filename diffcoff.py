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
import Dfit.Dfit as Dfit
warnings.filterwarnings("ignore") 


PDB = 'droplet_init.pdb'
TRR = 'droplet_md.trr'
XTC = 'droplet_md.xtc'
TPR = 'droplet_em.tpr'
FRAMES_PDB = 'droplet_dens-frames.pdb'
FRAMES_DIR = 'frames_centered'
TOP = './forcefield/topol.top'
DF = './coords-COM.dat'

START_TIME = 0 #in ns
END_TIME = 1000
KAPNUMB = 10

MASSES = {'ALA' :  71.08, 'ARG' : 156.19, 'ASN' : 114.10, 'ASP' : 115.09, 
 		  'CYS' : 103.14, 'GLN' : 128.13, 'GLU' : 129.12, 'GLY' :  57.05,
 		  'HIS' : 137.14, 'ILE' : 113.16, 'LEU' : 113.16, 'LYS' : 128.17, 
 		  'MET' : 131.19, 'PHE' : 147.18, 'PRO' :  97.12, 'SER' :  87.08,
 		  'THR' : 101.11, 'TRP' : 186.21, 'TYR' : 163.18, 'VAL' :  99.13,
 		  }


os.system('rmdir -r kapsonlyndx &> /dev/null')
os.system('rmdir -r diffcoff &> /dev/null')
os.system('mkdir kapsonlyndx &> /dev/null')
os.system('mkdir diffcoff &> /dev/null')

def main(): 
	ndx_file = create_ndx()
	xvg_list = extract_com(ndx_file)
	dat_list = convert_DAT(xvg_list)
#	run_Dfit(dat_list)
	# gmx_msd()


def create_ndx():
	u = mda.Universe(TOP, PDB, topology_format='ITP')
	kaps = u.select_atoms('segid Kap95')
	# print(kaps)
	# print(kaps.segments)
	# print(range(len(kaps.segments)))
	kaplist = list(kaps.segments)
	print(kaplist)

	ndx_file = 'kapsonly.ndx'
	with mda.selections.gromacs.SelectionWriter(ndx_file, mode='w') as ndx: #used segment instead of atom, problem?
		ndx.write(kaps, name='all_kaps')
		for i, kap in enumerate(kaps.segments):
			ndx.write(kap.atoms, name=f'kaps{i}')
	return ndx_file
#created ndx file with all kaps appended

def extract_com(ndx_file):
	# create new tpr file
	tpr_file = 'droplet_kaps.tpr'
	gmx_call = [
		'gmx', 'convert-tpr',
		'-s', TPR,
		'-n', ndx_file,
		'-o', tpr_file
	]
	p0 = subprocess.Popen(gmx_call, stdin=subprocess.PIPE)
	p0.communicate(b'0\n')
	p0.wait()

	# extract com from here on out for each trajectory
	file_list = []
	for i in range(1):
		coords_file = f'coords_kap{i}.xvg'
		file_list.append(coords_file)

		gmx_call = [
			'gmx', 'traj',
			'-f', XTC,
			'-s', tpr_file,
			'-nojump', '-com',
			'-ox', coords_file,
			'-n', ndx_file
		]
		p1 = subprocess.Popen(gmx_call, stdin=subprocess.PIPE)
		p1.communicate(f'{1}\n'.encode())
		p1.wait()

	return file_list

def convert_DAT(file_list):
	out_list = []
	for file in file_list:
		out_file = f'{file[:-4]}-COM.dat'
		out_list.append(out_file)
		with open(file, 'r') as fid, open(out_file, 'w+') as out:
			for line in fid.readlines():
				if line[0] not in '#@':
					tmp = line.split()
					data = [float(tmp[1]), float(tmp[2]), float(tmp[3])]
					out.write(f'{data[0]:8.4f} {data[1]:8.4f} {data[2]:8.4f} \n')
	return out_list

#def run_Dfit(dat_files):
 #   res = Dfit.Dcov(fz = dat_files, m = 20, tmin = 50, tmax = 2500, dt = 0.2, nitmax = 500, imgfmt = 'pdf')
#    res.run_Dfit()
 #   res.analysis(tc = 0.4)
  #  res.finite_size_correction(L = 30.6642, eta = 0.001, tc = 0.4)

	# return ndx_file




""" 
Two problems at the moment: 
- msd is only done for group 0, as in kap0, not the rest. 
- msd conducts diffusion constant in a plane at which height?
"""
# def gmx_msd():
# 	p0 = subprocess.Popen(
# 		[
# 			'gmx', 'msd',
# 			'-f', XTC,
# 			'-s', TPR,
# 			'-n', ndx_file,
# 			'-b', str(START_TIME),
# 			'-e', str(END_TIME),
# 			'-tu', 'ns',
# 			'-beginfit', str(START_TIME), 
# 			'-endfit', str(END_TIME),
# 			'-lateral', 'x', 
# 			# '-skip', str(TIMESTEP * 10),  # assuming nstxout = 5000
# 			'-o', f'diffcoff{k}.xvg'
# 		],
# 		stdin=subprocess.PIPE
# 	)
# 	p0.communicate(b'0\n')
# 	p0.wait()






if __name__ == '__main__':
	main()
