import subprocess, os
import numpy as np 

nclust_out = 'cluster.out'		# output file
ndx_dir    = 'cluster_ndx'		# directory with maxclust ndx files
cutoff = 0.7 # nm

# Create log files
with open(nclust_out, 'w+') as out:
	out.write('time (ns)    N\n')

# Create directories for files to keep
os.system(f'rm -r {ndx_dir} &> /dev/null')
os.system(f'mkdir {ndx_dir}')

# CREATE chain_md_frames.pdb
##########################################################################################

trr_file = 'droplet_md.xtc'			# input trajectory file
tpr_file = 'droplet_md.tpr'			# input topology file
index = 'index.ndx'
frames_pdb = 'droplet_md_frames.pdb'
dt = 10 # ns

os.system('rm -r frames      &> /dev/null')
os.system('mkdir frames')

# Get the trajectory of the chain for all times in a long .pdb
p2 = subprocess.Popen(['gmx', 'traj', '-f', trr_file, '-s', tpr_file, '-n', index,
	'-tu', 'ns', '-dt', str(dt),
	'-oxt', frames_pdb], stdin=subprocess.PIPE)
p2.communicate(b'0\n')
p2.wait()

# There will now be a long .pdb file with an entry for every frame of the simulation

# SPLIT .PDB FILE INTO INDIVIDUAL FRAMES
##########################################################################################

# Determine number of beads
with open('droplet_em.pdb', 'r') as fid:
	beads = len(fid.readlines()) - 6
	print(f'beads = {beads}')

# Split frames into separate files and place in frames dir.
os.system(f'split --verbose -l {beads + 6} {frames_pdb} ./frames/frame -d --additional-suffix=.pdb')
#for some reason, 90 turns into 9000, 9001 and so on, moreover the beads that are seperated do not have appropriate chainIDs
# RUN CLUSTER ANALYSIS
##########################################################################################

# Loop over frames
for idx, pdb_file in enumerate(sorted(os.listdir('frames/'))): 
	# Create cluster index file
	out_ndx = f'maxclust_{idx}.ndx'
	p1 = subprocess.Popen(['gmx', 'clustsize', '-f', 'frames/'+pdb_file, '-s', tpr_file, '-mcn', out_ndx,
		'-n', 'index.ndx', '-cut', str(cutoff)], stdin=subprocess.PIPE)
	p1.communicate(b'0\n')
	p1.wait()

	# Determine cluster size
	with open(out_ndx, 'r') as fid:
		in_cluster = len(fid.readlines()) - 1

	# Write cluster size
	with open(nclust_out, 'a+') as out:
		out.write(f'{idx*dt:8.1f}    {in_cluster:6}\n')

	# Clean up
	os.system(f'mv {out_ndx} {ndx_dir}')
	os.system('rm avclust.xvg cluster.log clust_size.xvg csizew.xpm csize.xpm histo-clust.xvg '
		'maxclust.xvg nclust.xvg rmsd-clust.xpm rmsd-dist.xvg temp.xvg')
