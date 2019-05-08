#!/usr/bin/env python

import numpy
import sys
import argparse

def isfloat(value):
  try:
    int(value)
    return True
  except:
    return False


parser = argparse.ArgumentParser(description='delete_water')

parser.add_argument('--pdb', dest='pdb_file', action='store', nargs=1, help='Coordinate file - PDB 3.0 name standard', required=True)
parser.add_argument('--out', default='dehydrated.pdb', dest='out_file', action='store', nargs=1, help='Output file - default dehydrated.pdb')
parser.add_argument('--min', dest='MIN', action='store', nargs=1, help='min (A) from which to erase', required=True)
parser.add_argument('--max', dest='MAX', action='store', nargs=1, help='max (A) up which to erase', required=True)
parser.add_argument('--dim', default='z', dest='DIM', action='store', nargs=1, help='perpendicular dimension (x, y or z) - default z')
parser.add_argument('--names', default=['OW','HW1','HW2'], dest='NAMES', action='store', nargs=1, help='string of names of water atoms within a molecule. \
Default \'OW HW1 HW2\' (united atom, GROMOS style)')
parser.add_argument('--resname', default='SOL', dest='RES', action='store', nargs=1, help='residue name of solvent, default \'SOL\'')

attributes = parser.parse_args()
tag =0
line_atoms = []

zmin = float(attributes.MIN[0])
zmax = float(attributes.MAX[0])
dim = attributes.DIM[0]
names = attributes.NAMES[0].split()
res = attributes.RES[0]

tag_limit = len(names)
tag = 0

with open(attributes.pdb_file[0], 'r') as coord:

    out_file = open(attributes.out_file[0],"w")
    
    for line in coord:
        linetag = line.split(None,)[0]
        if (linetag == "ATOM"):
            atomtag = line.split(None,)[3]
            if  (atomtag == res):								# If it is a solvent molecule
                tag = tag + 1
                if (dim == 'x'):
                    zposition = line.split(None,)[5]
                elif (dim == 'y'):
                    zposition = line.split(None,)[6]
                elif (dim == 'z'):
                    zposition = line.split(None,)[7]
                if (float(zposition) < zmin) or (float(zposition) > zmax):
                    line_atoms.append(line)
                if tag == tag_limit:
                    if (len(line_atoms) == tag_limit):
                        for i in range(tag_limit):
                            out_file.write(line_atoms[i])
                    line_atoms = []
                    tag = 0
            else:
                out_file.write(line)
        else:
            out_file.write(line)

    out_file.close()