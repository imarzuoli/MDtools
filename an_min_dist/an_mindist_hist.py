#!/usr/bin/python
# #!/usr/bin/env python3

import MDAnalysis as mda
import argparse
import numpy as np
from MDAnalysis.analysis.distances import distance_array

# Function to compute the distribution of water around an atom type/group of atoms as
# the distribution of the distance between the oxygen of water and the nearest atom in
# the group considered. Histogram averages over the timeframes.
# For an example, see:
# J. Chem. Theory Comput., 2010, 6 (1), pp 325-336 https://pubs.acs.org/doi/10.1021/ct900487a

def histogram_water(pdb, xtc, sel, oxwat = "OW", rmax = None, bins = 1200):
	
	"""
	
	Compute the distribution of the distances between the oxygen of water molecules and
	the nearest atom in the a group, at each time frame (preselected static binning).
	Return the histogram averaged over the timeframes.
	
	For an example, see: J. Chem. Theory Comput., 2010, 6 (1), pp 325-336
	
	
	Parameters
	------------
	pdb   : string
				Path to topology file (any accepted by MDAnalysis)
	trj   : string
				Path to trajectory file (any accepted by MDAnalysis)
	sel   : string
				MDAnalysis style selection for the group to analyse
	oxwat : string
				Atom name for water oxygen. Default "OW"
	rmax  : float
				Distance up to which carry on the analysis.
				Default (None) uses half the length of the shortest box side
	bins  : int
				Nr of bins in the histogram. Default 1200
	
	Returns
	------------
	numpy array of length [bins]
				Average water oxygen-group distances histogram
	
	"""
	
	# Create groups with selections (water oxygen and group of interest)
	ow = u.select_atoms("name " + oxwat)
	grp = u.select_atoms(sel)
	
	nr_waters = ow.positions.shape[0]
	nr_grpatoms = grp.positions.shape[0]
	
	# Set up and fill histogram for every time. If no rmax is given, take half of
	# the smaller box dimension
	minhist = np.zeros(( bins, u.trajectory.n_frames ))
	
	if (rmax == None):
		rmax = np.min(u.coord.dimensions[0:3]/2)
	
	width = np.around([1.0 * rmax/bins],2)
	breaks = np.arange(0, width * (bins + 1), width)
	
	for ts in u.trajectory:
		dist_matrix = distance_array(reference = ow.positions, configuration = grp.positions)
		mindist = dist_matrix.min(axis = 1)
		minhist[ :, ts.frame ] = np.histogram(bins = breaks, a = mindist, density = False)[0]
	
	# Return average histgram (in time)
	return minhist.mean(1)



# Example script to use the function on a xtc trajectory

parser = argparse.ArgumentParser(description='an_mindist_hist')

parser.add_argument('--pdb', dest='pdb', action='store', nargs=1, help='Topology file (any accepted by MDAnalysis)', required=True)
parser.add_argument('--trj', dest='trj', action='store', nargs=1, help='Trajectory file (any accepted by MDAnalysis)', required=True)
parser.add_argument('--selection', dest='sel', action='store', nargs=1, help='MDAnalysis style selection for the group to analyse', required=True)
parser.add_argument('--rmax', dest='rmax', action='store', nargs=1, help='Maximum distance up to which analyse', default=None)
parser.add_argument('--oxwat', dest='oxwat', action='store', nargs=1, default = "OW", help='Atom name for water oxygen. Default \"OW\"')
parser.add_argument('--out', dest='out', action='store', nargs=1, default = "average_water_histo.txt", help='Name for txt output file with histogram. Default "average_water_histo.txt"')
parser.add_argument('--nbins', dest='out', action='store', nargs=1, default = 1200, help='Nr of bins in the histogram. Default 1200')


#  Read parameters
attributes = parser.parse_args()
pdb = attributes.pdb[0]
xtc = attributes.trj[0]
sel = attributes.sel[0]
rmax = attributes.rmax[0]
oxwat = attributes.oxwat
outfile = attributes.out[0]
nbins = attributes.nbins[0]

# Set up universe
u = mda.Universe(pdb, xtc)

average_hist = histogram_water(pdb, xtc, sel, oxwat, rmax, bins = nbins)

np.savetxt(outfile, average_hist, fmt='%s')
