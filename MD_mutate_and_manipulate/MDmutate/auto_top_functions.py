import copy

# DICTIONARIES -------------------------------------------------------
angledic = {
  "CH1 C N": "ga_19",
  "N C CH1": "ga_19",
  "C C N": "ga_19",
  "N C C": "ga_19",
  "O C N": "ga_33",
  "N C O": "ga_33",
  "C N H": "ga_32",
  "H N C": "ga_32",
  "C N CH2": "ga_31",
  "CH2 N C": "ga_31"
}

dihdic = {
  "N CH1 C N": ["gd_42", "gd_45"],
  "N C CH1 N": ["gd_42", "gd_45"],
  "N C C N": ["gd_42", "gd_45"],
  "N C C N": ["gd_42", "gd_45"],
  "C N CH2 CH2": ["gd_43", "gd_44"],
  "CH2 CH2 N C": ["gd_43", "gd_44"],
  "CH1 C N CH2": ["gd_14"],
  "CH2 N C CH1": ["gd_14"],
  "C C N CH2": ["gd_14"],
  "CH2 N C C": ["gd_14"]
}


# CLASSES ------------------------------------------------------------

# User classes ---

# Object storing all the topology information
class topology:
	def __init__(info, name, ex, atoms, connect, compound, nmols):
		info.name = name
		info.nrexcl = ex
		# List of atom_top types
		info.atoms = atoms
		# Connect object
		info.connect = connect
		info.compound = compound
		info.nmols = nmols
	
	def atmnr(info):
		return len(info.atoms)
	
	def cgnr(info):
		cgnr = 0
		for i in info.atoms:
			cgnr = max( cgnr, i.cgnr )
		return cgnr


# Internal classes ---

# Object storing a line of the [ atoms ] section of top
class atom_top:
	def __init__(atomstab, nr, type, resnr, resid, atom, cgnr, charge, mass):
		atomstab.nr = nr
		atomstab.type = type
		atomstab.resnr = resnr
		atomstab.resid = resid
		atomstab.atom = atom
		atomstab.cgnr = cgnr
		atomstab.charge = charge
		atomstab.mass = mass

# Object storing the bonds/pairs/angles/dihedrals info
class connect:
	def __init__(bonded, bonds, pairs, angles, dih):
		# List of all bonds, of type bond
		bonded.bonds = bonds
		# List of all pairs, of type pair
		bonded.pairs = pairs
		# ...
		bonded.angles = angles
		bonded.dih = dih

# For each bonded interaction, objects of given shape
class bond:
	def __init__(bb, at1, at2, info):
		bb.at1 = at1
		bb.at2 = at2
		bb.info = info

class pair:
	def __init__(pp, at1, at2):
		pp.at1 = at1
		pp.at2 = at2
		pp.funct = 1
		
class angle:
	def __init__(ang, at1, at2, at3, info):
		ang.at1 = at1
		ang.at2 = at2
		ang.at3 = at3
		ang.funct = 2
		ang.info = info

class dih:
	def __init__(dih, at1, at2, at3, at4, funct, info):
		dih.at1 = at1
		dih.at2 = at2
		dih.at3 = at3
		dih.at4 = at4
		dih.funct = funct
		dih.info = info


# FUNCTIONS ----------------------------------------------------------

#INTERNAL FUNCTIONS ---
def empty_top():
	
	''' Return an empty topology object initializing all the fields to empty lists/NA'''
	
	new = topology("", 0, [], connect([],[],[],[]), 'NA', 0)
	
	return new


def info_reader(stored, mode, line):

	''' Called by read_top(), not to be used by user directly.
	    Reads topology line (line) according to the section it
	    belongs to (mode), and modify the corresponding field
	    of the topology which is going to be written (stored). '''
	
	# No action if mode is exclusions or system

	if mode == "moleculetype":
		stored.name = line.split()[0]
		stored.nrexcl = line.split()[1]
		
	if mode == "atoms":
		specs = line.split()
		atm = atom_top(nr = int(specs[0]), type = specs[1], resnr = int(specs[2]), resid = specs[3], atom = specs[4], cgnr = int(specs[5]), charge = float(specs[6]), mass = float(specs[7]))
		stored.atoms.append(atm)

	if mode == "bonds":
		bnd = bond(at1 = int(line.split()[0]), at2 = int(line.split()[1]), info = '	'.join(line.split()[2:]))
		stored.connect.bonds.append(bnd)

	if mode == "pairs":
		pp = pair(at1 = int(line.split()[0]), at2 = int(line.split()[1]))
		stored.connect.pairs.append(pp)
		
	if mode == "angles":
		ang = angle(at1 = int(line.split()[0]), at2 = int(line.split()[1]), at3 = int(line.split()[2]), info = '	'.join(line.split()[4:]))
		stored.connect.angles.append(ang)
	
	if mode == "dihedrals":
		dihedral = dih(at1 = int(line.split()[0]), at2 = int(line.split()[1]), at3 = int(line.split()[2]), at4 = int(line.split()[3]), funct = int(line.split()[4]), info = '	'.join(line.split()[5:]))
		stored.connect.dih.append(dihedral)
	
	if mode == "molecules":
		stored.nmols = int(line.split()[1])
		stored.compound = line.split()[0]


# IO FUNCTIONS ---
def read_top(top_file):
	
	''' Read topology (top/itp) file and store it in a topology object '''
	
	mode = 'NA'
	stored = empty_top()
	
	with open(top_file, 'r') as top:
		for line in top:
			if len(line.split()) != 0:
				if line.split()[0][0] != ";" and line.split()[0][0] != "#":
					if line.split()[0] == "[":
						mode = line.split()[1]
					else:
						info_reader(stored, mode, line)
	
	return stored


def write_itp(topology, out_itp_file, name = 'merged'):

	''' Write out_itp_file from a topology object. If no name is
	    provided, the molecule name in [ moleculetype ] will be
	    'merged' '''

	itpfile = open(out_itp_file,"w")
	
	itpfile.write("[ moleculetype ]\n")
	itpfile.write("; Name            nrexcl\n")
	itpfile.write("%s             3\n\n" %name)
	
	itpfile.write("[ atoms ]\n")
	for ati in topology.atoms:
		itpfile.write("	" + str(ati.nr) + "	" + ati.type + "	" + str(ati.resnr) + "	" + ati.resid + "	" + ati.atom + "	" + str(ati.cgnr) + "	" + str(ati.charge) + "	" + str(ati.mass) + "\n")
	
	itpfile.write("\n[ bonds ]\n")
	for bnd in topology.connect.bonds:
		itpfile.write("	" + str(bnd.at1) + "	" + str(bnd.at2) + "	" + bnd.info + "\n")
	
	itpfile.write("\n[ pairs ]\n")
	for pp in topology.connect.pairs:
		itpfile.write("	" + str(pp.at1) + "	" + str(pp.at2) + "	" + str(pp.funct) + "\n")
	
	itpfile.write("\n[ angles ]\n")
	for ang in topology.connect.angles:
		itpfile.write("	" + str(ang.at1) + "	" + str(ang.at2) + "	" + str(ang.at3) + "	" + str(ang.funct) + "	" + ang.info + "\n")
	
	itpfile.write("\n[ dihedrals ]\n")
	for dih in topology.connect.dih:
		itpfile.write("	" + str(dih.at1) + "	" + str(dih.at2) + "	" + str(dih.at3) + "	" + str(dih.at4) + "	" + str(dih.funct) + "	" + dih.info + "\n")
	
	itpfile.close() 


def write_54a8top(top_name, itp_name, posre_name = "", name = 'merged', N = 1):
	
	''' Write a 54a8 top file (called top_name), from a topology
	    object itp_file, putting N molecules of type name in the
	    [ molecules ] directory '''
	
	top_out  = open(top_name, "w") 
	top_out.write("#include \"gromos54a8.ff/forcefield.itp\"")
	top_out.write("\n#include \"%s\"" %itp_name)
	
	if posre_name != "":
		top_out.write("\n#ifdef POSRES \n#include \"%s\" \n#endif" % posre_name )
		top_out.write("\n\n#include \"gromos54a8.ff/spc.itp\"")
		top_out.write("\n\n#include \"gromos54a8.ff/ions.itp\"")
		top_out.write("\n\n[ system ]")
		top_out.write("\nmerged")
		top_out.write("\n\n[ molecules ]")
		top_out.write("\n%s       %d\n" % (name, N) )
	
	top_out.close() 


# MANIPULATION FUNCTIONS ---
def top_remove(topology, removendur_atm_nr):

	''' Remove an atom (removendur_atom_nr), identified by atom number,
	    from a topology object, modifying the object directly.
	    Take care of removing bonded interactions not present any more
	    and rescale the atoms number accounting for the removal.
	    Used in do_peptide_bond(). '''
	
	top_final = copy.deepcopy(topology)
	
	if removendur_atm_nr != []:
		removendur_atm_nr = [int(i) for i in removendur_atm_nr]
		for atm in removendur_atm_nr:
			# First remove
			for each in list(top_final.atoms):
				if int(each.nr) == atm:
					top_final.atoms.remove(each)
			for bnd in list(top_final.connect.bonds):
				if bnd.at1 == atm or bnd.at2 == atm:
					top_final.connect.bonds.remove(bnd)
			for pp in list(top_final.connect.pairs):
				if int(pp.at1) == atm or int(pp.at2) == atm:
					top_final.connect.pairs.remove(pp)
			for ang in list(top_final.connect.angles):
				if int(ang.at1) == atm or int(ang.at2) == atm or int(ang.at3) == atm:
					top_final.connect.angles.remove(ang)
			for dd in list(top_final.connect.dih):
				if int(dd.at1) == atm or int(dd.at2) == atm or int(dd.at3) == atm or int(dd.at4) == atm:
					top_final.connect.dih.remove(dd)
		# Then rescale
		adjusting_factor = []
		for number in top_final.atoms:
			adjusting_factor.append(sum([int(number.nr) > atm for atm in removendur_atm_nr ]))
		adjfac_with_gaps = []
		for number in topology.atoms:
			adjfac_with_gaps.append(sum([int(number.nr) > atm for atm in removendur_atm_nr ]))
		for nr in list(removendur_atm_nr):
			adjfac_with_gaps[nr-1] = 0

	if removendur_atm_nr != []:
		for each in range(len(top_final.atoms)):
			top_final.atoms[each].nr -= adjusting_factor[each]
		for each in range(len(top_final.connect.bonds)):
			top_final.connect.bonds[each].at1 -= adjfac_with_gaps[top_final.connect.bonds[each].at1-1]
			top_final.connect.bonds[each].at2 -= adjfac_with_gaps[top_final.connect.bonds[each].at2-1]
		for each in range(len(top_final.connect.pairs)):
			top_final.connect.pairs[each].at1 -= adjfac_with_gaps[top_final.connect.pairs[each].at1-1]
			top_final.connect.pairs[each].at2 -= adjfac_with_gaps[top_final.connect.pairs[each].at2-1]
		for each in range(len(top_final.connect.angles)):
			top_final.connect.angles[each].at1 -= adjfac_with_gaps[top_final.connect.angles[each].at1-1]
			top_final.connect.angles[each].at2 -= adjfac_with_gaps[top_final.connect.angles[each].at2-1]
			top_final.connect.angles[each].at3 -= adjfac_with_gaps[top_final.connect.angles[each].at3-1]
		for each in range(len(top_final.connect.dih)):
			top_final.connect.dih[each].at1 -= adjfac_with_gaps[top_final.connect.dih[each].at1-1]
			top_final.connect.dih[each].at2 -= adjfac_with_gaps[top_final.connect.dih[each].at2-1]
			top_final.connect.dih[each].at3 -= adjfac_with_gaps[top_final.connect.dih[each].at3-1]
			top_final.connect.dih[each].at4 -= adjfac_with_gaps[top_final.connect.dih[each].at4-1]
	
	return top_final, adjusting_factor


def top_rescale(topology, scale_atm = 0, scale_chg = 0):

	''' Rescale object topology renumbering the atoms by adding scale_atm
	    and the partial charged by adding scale_chg.
	    Used in do_peptide_bond(). '''
	
	rescaled_top = copy.deepcopy(topology)
	
	if scale_atm != 0:
	
		for each in range(len(rescaled_top.atoms)):
			rescaled_top.atoms[each].nr += scale_atm
		for each in range(len(rescaled_top.connect.bonds)):
			rescaled_top.connect.bonds[each].at1 += scale_atm
			rescaled_top.connect.bonds[each].at2 += scale_atm
		for each in range(len(rescaled_top.connect.pairs)):
			rescaled_top.connect.pairs[each].at1 += scale_atm
			rescaled_top.connect.pairs[each].at2 += scale_atm
		for each in range(len(rescaled_top.connect.angles)):
			rescaled_top.connect.angles[each].at1 += scale_atm
			rescaled_top.connect.angles[each].at2 += scale_atm
			rescaled_top.connect.angles[each].at3 += scale_atm
		for each in range(len(rescaled_top.connect.dih)):
			rescaled_top.connect.dih[each].at1 += scale_atm
			rescaled_top.connect.dih[each].at2 += scale_atm
			rescaled_top.connect.dih[each].at3 += scale_atm
			rescaled_top.connect.dih[each].at4 += scale_atm
	
	if scale_chg != 0:
		
		for each in range(len(rescaled_top.atoms)):
			rescaled_top.atoms[each].cgnr += scale_chg
	
	return rescaled_top


def merge_scaled_topology(list_top, out_name = "merged"):
	
	''' Merge two topologies in one unique object. Note that
	    the second topology must already be scaled. This function
	    does NOT do that. '''
	
	merg_top = empty_top()
	
	for top in list_top:
		if top != list_top[0]:
			merg_top.name = merg_top.name + " + " + top.name
		else:
			merg_top.name = merg_top.name + top.name
		merg_top.atoms = merg_top.atoms + top.atoms
		merg_top.connect.bonds = merg_top.connect.bonds + top.connect.bonds
		merg_top.connect.pairs = merg_top.connect.pairs + top.connect.pairs
		merg_top.connect.angles = merg_top.connect.angles + top.connect.angles
		merg_top.connect.dih = merg_top.connect.dih + top.connect.dih
		
		merg_top.compound = out_name
		merg_top.nmols = 1
		
	return merg_top


def bonds_peptide_bond(top, aa, bb):

	''' Write to topology top the bonded parameters (bonds, pairs,
	    angles, dihedrals) to create a peptide bond between atom
	    aa and bb, indicated by atom number.
	    Used in do_peptide_bond(). '''

	nn1 = []
	nn2 = []
	for bnd in top.connect.bonds:
		if bnd.at1 == aa:
			nn1.append(bnd.at2)
			if top.atoms[bnd.at2-1].type == "OM":
				top.atoms[bnd.at2-1].type = "O"
		elif bnd.at2 == aa:
			nn1.append(bnd.at1)
			if top.atoms[bnd.at1-1].type == "OM":
				top.atoms[bnd.at1-1].type = "O"
		if bnd.at1 == bb:
			nn2.append(bnd.at2)
			if top.atoms[bnd.at2-1].type == "OM":
				top.atoms[bnd.at2-1].type = "O"
		elif bnd.at2 == bb:
			nn2.append(bnd.at1)
			if top.atoms[bnd.at1-1].type == "OM":
				top.atoms[bnd.at1-1].type = "O"
	
	nnn1 = []
	nnn2 = []
	for i in nn1:
		pre = []
		for bnd in top.connect.bonds:
			if bnd.at1 == i and bnd.at2 != aa and bnd.at2 != bb:
				pre.append(bnd.at2)
			elif bnd.at2 == i and bnd.at1 != aa and bnd.at1 != bb:
				pre.append(bnd.at1)
		nnn1.append(pre)
	for i in nn2:
		pre = []
		for bnd in top.connect.bonds:
			if bnd.at1 == i and bnd.at2 != aa and bnd.at2 != bb:
				pre.append(bnd.at2)
			elif bnd.at2 == i and bnd.at1 != aa and bnd.at1 != bb:
				pre.append(bnd.at1)
		nnn2.append(pre)
	
	# Add bond
	top.connect.bonds.append( bond(aa, bb, "2	gb_10") )
	# Add angles
	for i in nn1:
		info = [top.atoms[i-1].type, top.atoms[aa-1].type, top.atoms[bb-1].type]
		top.connect.angles.append( angle(i, aa, bb, angledic[" ".join(info)]) )
	for i in nn2:
		info = [top.atoms[i-1].type, top.atoms[bb-1].type, top.atoms[aa-1].type]
		top.connect.angles.append( angle(i, bb, aa, angledic[" ".join(info)]) )
	# Add pairs and dihedrals
	# N-C-...-...
	for i in range(len(nn2)):
		if nnn2[i] != []:
			for j in range(len(nnn2[i])):
				top.connect.pairs.append( pair(aa, nnn2[i][j]) )
				if top.atoms[nnn2[i][j]-1].type != "H" and top.atoms[nnn2[i][j]-1].atom != "CB" and top.atoms[nnn2[i][j]-1].atom != "O":
					info = [ top.atoms[aa-1].type, top.atoms[bb-1].type, top.atoms[nn2[i]-1].type, top.atoms[nnn2[i][j]-1].type ]
					toadd = len(dihdic[" ".join(info)])
					for dihnr in range(toadd):
						top.connect.dih.append( dih(aa, bb, nn2[i], nnn2[i][j], 1, dihdic[" ".join(info)][dihnr]) )
	# ...-N-C-...
	for i in nn2:
		for j in nn1:
			top.connect.pairs.append( pair(j, i) )
			if top.atoms[i-1].type != "H" and top.atoms[j-1].type != "H" and top.atoms[i-1].atom != "CB" and top.atoms[j-1].atom != "CB" and top.atoms[i-1].type != "O" and top.atoms[j-1].type != "O":
					info = [top.atoms[j-1].type, top.atoms[aa-1].type, top.atoms[bb-1].type, top.atoms[i-1].type]
					toadd = len(dihdic[" ".join(info)])
					for dihnr in range(toadd):
						top.connect.dih.append( dih(j, aa, bb, i, 1, dihdic[" ".join(info)][dihnr]) )
	# ...-...-N-C
	for i in range(len(nn1)):
		if nnn1[i] != []:
			for j in range(len(nnn1[i])):
				top.connect.pairs.append( pair(nnn1[i][j], bb) )
				if top.atoms[nnn1[i][j]-1].type != "H" and top.atoms[nnn1[i][j]-1].atom != "CB" and top.atoms[nnn1[i][j]-1].type != "O":
					info = [top.atoms[nnn1[i][j]-1].type, top.atoms[nn1[i]-1].type, top.atoms[aa-1].type, top.atoms[bb-1].type]
					toadd = len(dihdic[" ".join(info)])
					for dihnr in range(toadd):
						top.connect.dih.append( dih(nnn1[i][j], nn1[i], aa, bb, 1, dihdic[" ".join(info)][dihnr]) )
	# Add improper dihedrals
	if top.atoms[aa-1].type == "N":
		CH2_N = 0
		H_N = 0
		O_C = 0
		CA_C = 0
		for i in nn1:
			if top.atoms[i-1].type == "CH2":
				CH2_N = i
			if top.atoms[i-1].type == "H":
				H_N = i
		for i in nn2:
			if top.atoms[i-1].type == "O":
				O_C = i
			if top.atoms[i-1].atom == "CA":
				CA_C = i
		top.connect.dih.append( dih(aa, bb, CH2_N, H_N, 2, "gi_1" ))
		top.connect.dih.append( dih(bb, CA_C, aa, O_C, 2, "gi_1" ))
	elif top.atoms[bb-1].type == "N":
		CH2_N = 0
		H_N = 0
		O_C = 0
		CA_C = 0
		for i in nn2:
			if top.atoms[i-1].type == "CH2":
				CH2_N = i
			if top.atoms[i-1].type == "H":
				H_N = i
		for i in nn1:
			if top.atoms[i-1].type == "O":
				O_C = i
			if top.atoms[i-1].atom == "CA":
				CA_C = i
		top.connect.dih.append( dih(bb, aa, CH2_N, H_N, 2, "gi_1" ))
		top.connect.dih.append( dih(aa, CA_C, bb, O_C, 2, "gi_1" ))


def do_peptide_bond(topN, topC, atomN, atomC, topN_remove, topC_remove):

	''' Merge two topologies via a peptide bond: topN is the topology object
	    hosting the atom to be used as N-ter, topC the one as C-ter, atomN and
	    atomC are their atom numbers in the respective topologies.
	    topN_remove the extra hydrogen/atoms to be removed when creating the
	    bonds, indicated by atom number in topN. Similar for topC_remove in
	    topC.
	    It removes extra atoms, adjust terminal charges, rescale topology accordingly,
	    merge topologies, add bonded paramters (bonds, pairs, angles, dihedrals) '''
	
	# To be integrated in one function
	topN_rem, rmN = top_remove(topN, topN_remove)
	mapped_atomN = atomN - rmN[atomN-1]
	topC_rem, rmC = top_remove(topC, topC_remove)
	mapped_atomC = atomC - rmC[atomC-1]
	# Fix final charge of C-ter, assuming oxygen is after C
	topC_rem.atoms[topC_rem.atoms[mapped_atomC].nr-2].charge = 0.45
	topC_rem.atoms[topC_rem.atoms[mapped_atomC].nr-1].charge = -0.45
	
	# Fix final charge of N-ter, assuming H is after and fishing for nearby C
	N_neighbours = []
	for bnd in topN.connect.bonds:
		if bnd.at1 == mapped_atomN and topN_rem.atoms[bnd.at2-1].type[0] == "C":
			N_neighbours.append(bnd.at2)
		elif bnd.at2 == mapped_atomN:
			N_neighbours.append(bnd.at1)
	#topN.atoms[topN.atoms[N_neighbours[0]].nr-1].charge = 0
	topN_rem.atoms[topN_rem.atoms[mapped_atomN].nr-1].charge = -0.31
	topN_rem.atoms[topN_rem.atoms[mapped_atomN].nr-2].charge = 0.31
	
	# Prepare: first topN, then topC
	topC_rescaled = top_rescale(topC_rem, topN_rem.atmnr(), topN_rem.cgnr())
	
	new_top = merge_scaled_topology([topN_rem, topC_rescaled])
	
	bonds_peptide_bond(new_top, mapped_atomN, topN_rem.atmnr() + mapped_atomC)

	return new_top

# def identify_bonds: to map a pair of atom to a gb value (?)
# def missing_bonds: to flag a pair of atom without gb value.
#		     Put bond the same but throw warning and
#		     leave a blank/standard value?
