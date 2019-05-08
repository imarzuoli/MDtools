#!/usr/bin/python3

# This script:
# - generate the topology for a complex object (a three branched peptide, designed as three
#   identical sequences attached on a central residues via peptide bonds;
# - starting from the original sequence, it introduces mutations on each of the branches,
#   producing the relative pdb and topology

# SETUP
import numpy
import sys
import argparse
import subprocess
import copy

import MDAnalysis as mda

pwd = subprocess.run(["pwd"], stdout=subprocess.PIPE, encoding = 'utf-8')
mutdir = pwd.stdout.split('\n')[0]+'/MDmutate'
sys.path.insert(0, mutdir)

import auto_top_functions as af
import mutate_functions as mutf
import renumber_function as ren

import os
import Pmw
import itertools

from os.path import splitext


# Create topology for a branched peptide:
# attaching three copies of RRWTWE sequence to
# customised points of RKGB residue.

armtop = "./input/RRWTWE.top"
extop = "./input/RKGB_NH3ter.top"
out_top_tag = "triskelion"

# Load a pdb object for the central residue and create
# a list to host the topology for the 3 branches
topN_orig = af.read_top(extop)
list_topCs = [af.read_top(armtop)]
list_topCs.append(copy.deepcopy(list_topCs[0]))
list_topCs.append(copy.deepcopy(list_topCs[0]))

# N atoms in the central residues to create the peptide bonds
branch_points = [22, 28, 36]
# H atoms to be removed from central residues
branch_remove = [[24, 25], [30, 31], [38, 39]]

# Create a copy to manipulate
topN = topN_orig

for attach in range(len(list_topCs)):
	
	# Pick a branch
	topC = list_topCs[attach]
	atomN = branch_points[attach]
	topN_remove = branch_remove[attach]
	# Assuming ttachment to the last C atom (and a COO terminal)
	atomC = topC.atmnr() - 2
	topC_remove = topC.atmnr()
	
	# Converting to list also in the case there is only one attachment point
	if isinstance(topN_remove, int):
		topN_remove = [topN_remove]
	if isinstance(topC_remove, int):
		topC_remove = [topC_remove]
	
	# Rescale atoms to remove/to which attach a branch if atoms previously
	# listed have been removed (e.g. removing 24 and 25, the next 
	# attachment point becomes 26 and the next removed 28 and 29).
	for next in range(len(branch_points[attach+1:])):
		if branch_points[attach+next+1] > atomN:
			branch_points[attach+next+1] = branch_points[attach+next+1] - len(branch_remove[attach])
			for tip in range(len(branch_remove[attach+next+1])):
				branch_remove[attach+next+1][tip] -= len(branch_remove[attach])
	
	# Include in the topology the peptide bond
	new_top = af.do_peptide_bond(topN, topC, atomN, atomC, topN_remove, topC_remove)
	topN = new_top

# Write an itp file and the corresponding topology (gromos 54a8 standards) which
# read the ipt just generates. A posre fie is generated as well.
# 'merged' is the name given to the molecule in the top, N how many there are.
af.write_itp(new_top, "%s.itp" % out_top_tag)
af.write_54a8top("%s.top" % out_top_tag, "%s.itp" % (out_top_tag), "posreall_%s.itp" % out_top_tag, 'merged', N = 1)

