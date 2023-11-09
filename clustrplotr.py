#clustrplotr
import subprocess, os
import numpy as np 
import matplotlib.pyplot as plt 
import csv
import MDAnalysis as mda
from pathlib import Path
from tqdm import tqdm

# from pathlib import path

ndx_dir    = 'cluster_ndx'		# directory with maxclust ndx files


""" 
os compands are for referring to the linux directory and sending commands there
"""
#opening nclust and reading from it

def regular():
	nclust_out = 'cluster.csv'		# output file, made it into a csv earlier
	ndx_dir    = 'cluster_ndx'		# directory with maxclust ndx files


	with open(nclust_out, 'r') as out:
		lines = out.readlines()
		

	linesread = csv.DictReader(lines)

	time = []
	clustersize = []

	for col in linesread: 
		time.append(col['time (ns)'])
		clustersize.append(col['N'])

	time = [eval(i) for i in time]
	clustersize = [eval(i) for i in clustersize]
	print(clustersize)
	print(time)
	# fig = plt.figure()
	# plt.plot(time,clustersize)
	# plt.ylabel('Size (N)')
	# plt.xlabel('time (ns)')
	# plt.show()

regular()
#all this said and done is good, but now we need to differentiate them by members of the cluster, i.e. divide kaps and nups! 
#for this, different files need to be used. The csv and ndx files only include the number of molecules.. 
#So we need to go all the way back to using the frames to create new ndx/out files.. 


#index over each frame, and count the number of nups/kaps in the cluster, then append it to the .out file.
#awk over the .out file and turn it into a csv
#plot the csv



# {}frames.pdb already exist. Time to open them and index over them. Then, find out the number of kaps in the cluster.. 
# first, extract only the cluster as a pdb file. 

os.system('rm -r clusteronly      &> /dev/null')
os.system('mkdir clusteronly')


clusterpath = Path.cwd() / ndx_dir

for idx, pdb_file in enumerate(sorted(os.listdir('frames/'))):
	pdb_file = f'frames/frame{idx}.pdb'
	ndx_file = clusterpath / f'maxclust_{idx}.ndx'
	top_file = 'forcefield/topol.top'
	dt = 10 #ns
	#for each pdb now use the index files maxclust.ndx to extract cluster into its own pdb file
	u = mda.Universe(pdb_file)

	linesread = []

	with open(ndx_file, 'r') as out:
		lines = out.readlines()
		lines = [i.strip() for i in lines[1:]]
		lines = [int(j) for j in lines[1:]]

	if lines[-1] == 121985:
		lines.remove(121985)

	maxgrp = u.select_atoms('all')
	maxgrp = maxgrp[lines]
	# p2 = subprocess.Popen(['gmx', 'extract-cluster', '-f', pdb_file, '-s', top_file, '-clusters', ndx_file,
	# '-o', f'clusteronly{idx}.pdb', '-dt', str(dt)], stdin=subprocess.PIPE)
	# p2.wait()
	print(f'cluster{idx} created!')
	
	#now extract universe as pdb file
	maxgrp.write(f'clusteronly{idx}_max_clust.pdb')

	#mkdir clustersonly 
	#mv

	p1 = subprocess.Popen(['mv', f'clusteronly{idx}_max_clust.pdb', 'clusteronly/' ])
	p1.wait()

	



#hypothetically, lets say the above is true and solved. In that case, we now use the extracted clusters and find out which beads are kaps, which beads are nups. 
#for this, take note that every frame has chains for each nup molecule. we can consider that chains a - whatever indice are nups, unnamed chains are kaps.
#total - chained = number of kap beads. 

clusterdata = [68875, 68875, 69736, 70597, 70597, 72319, 72319, 74041, 74041, 74041, 74041, 75763, 76624, 76624, 77485, 77485, 77485, 77485, 79071, 79071, 79071, 79932, 80793, 81654, 81654, 82515, 82515, 82515, 83376, 83376, 84237, 85098, 84237, 85098, 85098, 85098, 85098, 85098, 85098, 85098, 85098, 85098, 85098, 85098, 85098, 85098, 85959, 85959, 86820, 85959, 86820, 86820, 86820, 87681, 87681, 87681, 87681, 87681, 87681, 87681, 88542, 88542, 89403, 88542, 88542, 88542, 88542, 89403, 87681, 87681, 88542, 88542, 88542, 88542, 89403, 89403, 90264, 90264, 91125, 91986, 92847, 92847, 92847, 93708, 93708, 94569, 94569, 95430, 95430, 94569, 94569, 94569, 95430, 95430, 95430, 95430, 95430, 95430, 95430, 95430, 95430, 96291, 97152, 96291, 96291, 96291, 95430, 95430, 96291, 96291, 97152, 97152, 97152, 98013, 98874, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98013, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 98874, 99735, 99735, 99735, 99735, 99735, 99735, 99735, 99735, 99735, 99735, 99735, 100596, 100596, 100596, 100596, 100596, 100596, 101457, 101457, 101457, 101457, 101457, 101457, 101457, 102318, 102318, 102318, 102318, 103179, 103179, 103179, 102318, 102318, 102318, 102318, 102318, 102318, 102318, 102318, 104040, 104040, 104040, 104901, 104040, 104040, 104040, 104040, 103315, 104040, 104040, 104040, 104901, 104901, 104040, 103315, 103315, 104176, 106487, 106487, 106487, 104176, 104176, 104176, 104176, 104176, 104176, 104176, 106487, 106487, 106487, 106487, 106487, 106487, 106487, 106487, 106487, 107348, 106487, 107348, 106487, 107348, 108209, 107348, 108209, 108209, 108209, 108209, 108209, 108209, 108209, 108209, 108209, 108209, 108209, 108209, 108209, 109070, 109070, 109070, 107484, 107484, 107484, 107484, 107484, 107484, 107484, 107484, 107484, 107484, 107484, 107484, 106623, 107484, 106623, 106623, 107484, 106623, 106623, 106623, 106623, 106623, 106623, 106623, 106623, 106623, 105762, 106623, 106623, 106623, 106623, 106623, 106623, 105762, 105762, 105762, 104901, 104901, 104901, 104901, 104901, 104901, 104901, 105762, 104901, 104901, 104901, 105762, 105762, 105762, 105762, 105762, 105762, 106623, 106623, 106623, 106623, 106623, 106623, 106623, 106623, 107484, 107484, 107484, 107484, 107484, 108345, 108345, 108345, 108345, 109206, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 109206, 109206, 109206, 109206, 109206, 109206, 109206, 109206, 110067, 109206, 109206, 109206, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110067, 110928, 110067, 110928, 110928, 110928, 110928, 109342, 109342, 110203, 109342, 110203, 111064, 110203, 110203, 110203, 110203, 110203, 110203, 110203, 110203, 109342, 110203, 110203, 109342, 107756, 107756, 107756, 107756, 109342, 109342, 109342, 109342, 106170, 106170, 107031, 107031, 107031, 107031, 107031]
time = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000, 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090, 1100, 1110, 1120, 1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1260, 1270, 1280, 1290, 1300, 1310, 1320, 1330, 1340, 1350, 1360, 1370, 1380, 1390, 1400, 1410, 1420, 1430, 1440, 1450, 1460, 1470, 1480, 1490, 1500, 1510, 1520, 1530, 1540, 1550, 1560, 1570, 1580, 1590, 1600, 1610, 1620, 1630, 1640, 1650, 1660, 1670, 1680, 1690, 1700, 1710, 1720, 1730, 1740, 1750, 1760, 1770, 1780, 1790, 1800, 1810, 1820, 1830, 1840, 1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2110, 2120, 2130, 2140, 2150, 2160, 2170, 2180, 2190, 2200, 2210, 2220, 2230, 2240, 2250, 2260, 2270, 2280, 2290, 2300, 2310, 2320, 2330, 2340, 2350, 2360, 2370, 2380, 2390, 2400, 2410, 2420, 2430, 2440, 2450, 2460, 2470, 2480, 2490, 2500, 2510, 2520, 2530, 2540, 2550, 2560, 2570, 2580, 2590, 2600, 2610, 2620, 2630, 2640, 2650, 2660, 2670, 2680, 2690, 2700, 2710, 2720, 2730, 2740, 2750, 2760, 2770, 2780, 2790, 2800, 2810, 2820, 2830, 2840, 2850, 2860, 2870, 2880, 2890, 2900, 2910, 2920, 2930, 2940, 2950, 2960, 2970, 2980, 2990, 3000, 3010, 3020, 3030, 3040, 3050, 3060, 3070, 3080, 3090, 3100, 3110, 3120, 3130, 3140, 3150, 3160, 3170, 3180, 3190, 3200, 3210, 3220, 3230, 3240, 3250, 3260, 3270, 3280, 3290, 3300, 3310, 3320, 3330, 3340, 3350, 3360, 3370, 3380, 3390, 3400, 3410, 3420, 3430, 3440, 3450, 3460, 3470, 3480, 3490, 3500, 3510, 3520, 3530, 3540, 3550, 3560, 3570, 3580, 3590, 3600, 3610, 3620, 3630, 3640, 3650, 3660, 3670, 3680, 3690, 3700, 3710, 3720, 3730, 3740, 3750, 3760, 3770, 3780, 3790, 3800, 3810, 3820, 3830, 3840, 3850, 3860, 3870, 3880, 3890, 3900, 3910, 3920, 3930]
kapbeadsperframe = []


PDB_FILE = 'droplet_init.pdb'


current = Path.cwd()
ndx_dir = current / 'cluster_ndx'

kapbeadsperframe = []
for idx, ndx in tqdm(enumerate(sorted(os.listdir('./cluster_ndx/')))):
	out_ndx = ndx_dir / f'maxclust_{idx}.ndx'

	u = mda.Universe(PDB_FILE)
	chainIDs = ['A', 'B']
	with open(out_ndx, 'r') as clust:
		lines = clust.readlines()
		lines = [i.strip() for i in lines[1:]]
		lines = [eval(j) for j in lines[1:]]

	sel = u.select_atoms('all')
	sel2 = sel[lines[:-1]]
	sel3 = sel2.select_atoms(f'chainID B')
	kapbeadsperframe.append(len(sel3))

# for idx, pdb_file in tqdm(enumerate(sorted(os.listdir('clusteronly/')))):
# 	u = mda.Universe(f'clusteronly/clusteronly{idx}_max_clust.pdb')
# 	chainID = ['X']
# 	kapbeads = 0
# 	sel = u.select_atoms('chainID X') #selects atomgroup
# 	lengths = len(sel.segments[0].atoms)
# 	kapbeads = lengths 	
# 	# print(f'frame{idx} done')
# 	kapbeadsperframe.append(kapbeads)

#must be integers or slices, not str

# then, we consider the kap beads per frame, so that we can plot it in the end. In order to do this run the function first.

kapbeadsperframearr = np.asarray(kapbeadsperframe)
clustersizearr = np.asarray(clusterdata)
nupbeadsperframearr = clustersizearr - kapbeadsperframearr
nupbeadsperframe = nupbeadsperframearr.tolist()


print(nupbeadsperframearr)
print(kapbeadsperframearr)
print(clusterdata)

# finally, we can plot the lists with the rest of the figure. 
fig = plt.figure()
plt.plot(time, kapbeadsperframearr, label = 'Kap capacity')
plt.plot(time, nupbeadsperframearr, label = 'nup capacity')
plt.plot(time, clusterdata, label = 'cluster capacity')
plt.axvline(x=2500, color = 'green', linestyle = 'dotted')
plt.xlabel('Time (ns)')
plt.ylabel('occupation(N)')
plt.legend()
plt.show()


