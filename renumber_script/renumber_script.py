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


parser = argparse.ArgumentParser(description='Reorder pdb according to preexistent topology')

parser.add_argument('--pdb', dest='pdb_file', action='store', nargs=1, help='Coordinate file - PDB 3.0 name standard', required=True)
parser.add_argument('--top', dest='top_file', action='store', nargs=1, help='Topology ot itp file - only the [atoms] section required',required=True)
parser.add_argument('--out', default='renumbered.pdb', dest='out_file', action='store', nargs=1, help='Output file')

attributes = parser.parse_args()

out_name = attributes.out_file[0]
if len(outname) < 4:
	out_name = out_name + '.pdb'
elif out_name[-4:] != '.pdb':
	out_name = out_name + '.pdb'
	
with open(attributes.top_file[0], 'r') as top:
	with open(attributes.pdb_file[0], 'r') as coord:
		
		print( "Reordering " + coord.name + " according to topology file " + top.name )

		# Topology line counter
		top_number = 0
		top_labels = []
		for line in top:
			number = line.split()
			if(number[0:2] == ['[','bonds']):
				break
			if number != []:
				if isfloat(number[0]) == True:
					top_number = top_number+1
					name = number[4]
					top_labels.append(name)

		# Pdb line counter
		pdb_number = 0
		pdb_indeces = []
		pdb_lines = []
		pdb_toswap = []
		for line in coord:
			first = line.split(None,)[0]
			if first == "ATOM" or first == "HETAT" or first == "HETATM":
				pdb_number = pdb_number+1
				pdb_indeces.append(1)
				pdb_toswap.append(1)
				pdb_lines.append(line)
			else:
				pdb_indeces.append(0)

		NMOL = int(pdb_number/top_number)
		if NMOL == pdb_number/top_number:
			print( "I am going to reorder " + str(pdb_number) + " atoms.")
		else:
			print( "Topology cannot match the pdb given (nr. of atoms in pdb is NOT a multiple of number of atoms in top.")
		
		# Print array with lines in new order
		new_order = []
		for nmol in range(0,NMOL):
			for label in top_labels:
				i = 0
				for line in pdb_lines:
					match = line.split()
					if (pdb_toswap[i] == 1) and (match[2] == label):
						new_order.append(line)
						pdb_toswap[i] = 0
						break
					i = i + 1

		print( "Output written on " + str(out_name))
		out_file = open(out_name,"w")
		for line in new_order:
			out_file.write(line)
		out_file.close()


