import numpy
import sys
import argparse

def isfloat(value):
  try:
    int(value)
    return True
  except:
    return False

def isint(value):
  try:
    int(value)
    return True
  except:
    return False


def format_pdb(data):

	''' Adjust spacing of a list of 9 (data) and
	    join them to return a PDB formatted string '''

	# Adjust each field
	# ATOM
	data[0] = data[0].ljust(6)
	# Atom number
	data[1] = data[1].rjust(5)
	# Append extra white space on the right (column 12)
	data[1] = data[1].ljust(6)
	# Atom name
	data[2] = data[2].rjust(4)
	# Append extra white space on the right (column 17)
	data[2] = data[2].ljust(5)
	# Res name
	data[3] = data[3].ljust(4)
	if data[4] == "":
		data[4] = " "
	# Res number
	data[5] = data[5].rjust(4)
	# Append col 27 (iCode) and 3 white spaces
	data[5] = data[5].ljust(8)
	# X
	data[6] = data[6].rjust(8)
	# Y
	data[7] = data[7].rjust(8)
	# Z
	data[8] = data[8].rjust(8)
	data_to_write = "".join(data)
	data_to_write = data_to_write[0:55]
	data_to_write = data_to_write + "\n"
	return data_to_write


def pdb_line(line):

	''' Take a pdb line (line) and split the fields according to the PDB format.
	    Return a list of length 9 '''

	# Adjust each field
	data = []
	# ATOM 6
	data.append(line[0:6])
	# Atom number with extra white space on the right (column 12) 6
	data.append(line[6:12])
	# Atom name with extra white space on the right (column 127) 5
	data.append(line[12:17])
	# Res name with extra white space on the right (column 21) 4
	data.append(line[17:21])
	# Chain name
	data.append(line[21])
	# Res number 4
	data.append(line[22:26])
	# X 8
	data.append(line[30:38])
	# Y 8
	data.append(line[38:46])
	# Z 8
	data.append(line[46:54])

	for i in range(9):
		data[i] = data[i].replace(" ", "")

	return data


def adj_names(pdb_file, out_file):

	''' In a PDB file (pdb_file), if atomname starts with a number
	    e.g. 1HH2, it puts all the figures at the end
	    e.g. HH21. Output to out_file. '''
	
	with open(pdb_file, 'r') as coord:
		
		out_file = open(out_file,"w")
		for line in coord:
			if len(line)>20:
				match = pdb_line(line)
				atm_name = match[2]
				if isint(atm_name[0]) == True:
					tmp = atm_name[1:] + atm_name[0]
					match[2] = tmp
				line_to_write = format_pdb(match)
				out_file.write(line_to_write)
			else:
				out_file.write(line)
		out_file.close()


def reshuffle_top(top_file, pdb_file, out_file, strict = True):

	''' Reshuffle the order of a PDB file (pdb_file) according to an ITP
	    topology file (top_file). If strict is TRUE it does it only if 
	    the number of atoms matches. If FALSE, allows for different
	    numbers taking only the atoms listed in the topology (useful
	    if PDB has extra hydrogens with respect to GROMACS
	    topology): please double check your output! '''

	with open(top_file, 'r') as top:
		with open(pdb_file, 'r') as coord:
		
		
			# Topology line counter
			top_number = 0
			top_labels = []
			top_res = []
			for line in top:
				number = line.split()
				if number != []:
					if isfloat(number[0]) == True:
						top_number = top_number+1
						name = number[4]
						top_labels.append(name)
						top_res.append(number[3])
				if line == "[ bonds ]\n":
					break	

			# Pdb line counter
			pdb_number = 0
			pdb_indeces = []
			pdb_lines = []
			pdb_toswap = []
			for line in coord:
				first = line.split(None,)[0]
				if first == "ATOM" or first == "HETAT":
					pdb_number = pdb_number+1
					pdb_indeces.append(1)
					pdb_toswap.append(1)
					pdb_lines.append(line)
				else:
					pdb_indeces.append(0)
	
			NMOL = int(pdb_number/top_number)

			if strict is True:
				if NMOL == pdb_number/top_number:
					print( "I am going to reorder " + str(pdb_number) + " atoms.")
				else:
					print( "Topology cannot match the pdb given.")
			
			# Print array with lines in new order
			new_order = []
			for nmol in range(0,NMOL):
				for pick in range(len(top_labels)):
					label = top_labels[pick]
					i = 0
					for line in pdb_lines:
						match = pdb_line(line)
						if (pdb_toswap[i] == 1) and (match[2] == label) and (match[3] == top_res[pick]):
							new_order.append(line)
							pdb_toswap[i] = 0
							break
						i = i + 1
	
		out_file = open(out_file,"w")
		for line in new_order:
			out_file.write(line)
		out_file.close()


