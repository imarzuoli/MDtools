#!/usr/bin/python

# This script takes a pdb and - calling repeatedly pymol - creates mutations of the original
# sequence according to a scheme of mutations passed by the user.

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


# Load files
orig = "./input/RRWTWE_original.pdb"
nametag = "arm"
seq_length = 6

pwd = os.getcwd()
orig_sequence = mutf.setNames(orig)

# A possible scheme of mutations (mutating positive residues with
# positive ones or maybe polar, keeping other core ones as fixed)
# TAKE CARE, it scales up quickly!
mutations = [[]] * seq_length
mutations[0] = mutf.getAa("plus")
mutations[1] = mutf.getAa("plus")
mutations[3] = mutf.getAa("plus")

# Prepare the list of sequences (all possible combinations)
for res in range(len(mutations)):
	if mutations[res] == []:
		mutations[res] = [mutf.setNames(orig)[res][1]]
sequences = list(itertools.product(*mutations))

# Get all the mutated pdb
for seq in sequences:
	# Mutate, save the pdb and return the new name
	newmutant = mutf.pymutate(orig, seq, nametag)

#	# Optionally prepare the topology with gromacs, if installed
#	bashCommand = ('gmx pdb2gmx -ff gromos54a7 -water spc -ignh -f {} -o {} -p {}'.format(newmutant, newmutant, "%s.top" %newmutant.split(".pdb")[0]))
#	process = subprocess.run(bashCommand.split(), stdin = subprocess.PIPE, stdout = None)
#	subprocess.run(["rm", "posre.itp" ])

# Gather results in a new folder - ATTENTION, if one with the same name exists,
# it will b removed
subprocess.run(["rm", "-r", "%s_mutations/" %nametag ])
subprocess.run(["mkdir", "%s_mutations/" %nametag ]) 
subprocess.Popen(["mv *_%s.pdb %s_mutations/" % (nametag, nametag) ],  shell = True) 

