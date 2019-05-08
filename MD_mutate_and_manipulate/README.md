Collection of scripts to manipulate peptide files, in order to easily create topology files
when a peptide bond is inserted, to manipulate (reorder) a pdb file if the topology
describes the molecule with a different ordering of the atoms, and to introduce mutations
in an automated way.

-----------------------------------------------------------------------------------
Scripts:

example_create_peptide_bonds.py
The script create the topology for a multi branched peptide inserting peptide bonds
at the required positions from the topology of the single pieces.
It requires the index of the atoms to be joined, and the ones to be erased (e.g. extra
hydrogens).
It takes care of:
- rescaling the index in one topology and append it to the first one;
- erase the extra atoms and rescale the remaining ones to have a continuous
  numbering in all the topology sections (i.e. erase their atom information, but also
  all the bonds, angles, ... where they appear;
- add the bond, angles, ... information for the peptide bond created. For regular atom
  types, it retrieves automatically the appropriate bond, angle, ... details in a gromos
  54a8 like fashion (gb*, ga*, ...)
  
[The default run calls RKGB_NH3ter.top and RRWTWE.top from the input folder, 
and produces a arm_mutations folder (see example_results)]

example_mutate_peptide.py
[Backend pymol]
Given a list with the mutations proposed for each amino acids in the sequence,
creates all the possible combinations of sequence from these mutations and
generated the pd. files.
!!! At present, no search for the best rotamer is implemented !!!

[The default run calls RRWTWE_original.pdb from the input folder, and produces
triskelion.itp and triskelion.top (see example_results)]

-----------------------------------------------------------------------------------
MDmutate/renumber_functions

Requirements: numpy

# Functions

format_pdb(data)
	''' Adjust spacing of a list of 9 (data) and
	    join them to return a PDB formatted string '''


pdb_line(line) [not intended as a user function]
	''' Take a pdb line (line) and split the fields according to the PDB format.
	    Return a list of length 9 '''


adj_names(pdb_file, out_file)
	''' In a PDB file (pdb_file), if atomname starts with a number
	    e.g. 1HH2, it puts all the figures at the end
	    e.g. HH21. Output to out_file. '''


reshuffle_top(top_file, pdb_file, out_file, strict = True)
	''' Reshuffle the order of a PDB file (pdb_file) according to an ITP
	    topology file (top_file). If strict is TRUE it does it only if 
	    the number of atoms matches. If FALSE, allows for different
	    numbers taking only the atoms listed in the topology (useful
	    if PDB has extra hydrogens with respect to GROMACS
	    topology): please double check your output! '''


-----------------------------------------------------------------------------------
MDmutate/auto_top_functions

Requirements: copy

# Objects

topology
	Attributes:
	info.name	[string]
	info.nrexcl	[int]
	info.atoms	[atom_top object]
	info.connect	[connect object]
	info.compound	[string]
	info.nmols	[int]
	
	Methods:
	atmnr(info): returns number of atoms in topology
	cgnr(info): returnsnumber of charge groups in topology

atom_top [not intended as a user class]
	Attributes:
	atomstab.nr	[int]
	atomstab.type	[string]
	atomstab.resnr	[int]
	atomstab.resid	[string]
	atomstab.atom	[int]
	atomstab.cgnr	[int]
	atomstab.charge	[float]
	atomstab.mass	[float]

connect [not intended as a user class]
	Attributes:
	bonded.bonds	[bond object]
	bonded.pairs	[pair object]
	bonded.angles	[ngle object]
	bonded.dih	[dih object]

bond [not intended as a user class]
	Attributes:
	bb.at1		[int]
	bb.at2		[int]
	bb.info		[string]

pair [not intended as a user class]
	Attributes:
	pp.at1		[int]
	pp.at2		[int]
	pp.funct	[string]

angle [not intended as a user class]
	Attributes:
	ang.at1		[int]
	ang.at2		[int]
	ang.at3		[int]
	ang.funct	[int = 2]
	ang.info	[string]

dih [not intended as a user class]
	Attributes:
	dih.at1		[int]
	dih.at2		[int]
	dih.at3		[int]
	dih.at4		[int]
	dih.funct	[int]
	dih.info	[tring]


# Functions 

# IO functions

empty_top() [not intended as a user function]
	''' Return an empty topology object initializing all the fields to empty lists/NA'''


info_reader(stored, mode, line) [not intended as a user function]
	''' Called by read_top(), not to be used by user directly.
	    Reads topology line (line) according to the section it
	    belongs to (mode), and modify the corresponding field
	    of the topology which is going to be written (stored). '''


read_top(top_file)	
	''' Read topology (top/itp) file and store it in a topology object '''


write_itp(topology, out_itp_file, name = 'merged')
	''' Write out_itp_file from a topology object. If no name is
	    provided, the molecule name in [ moleculetype ] will be
	    'merged' '''


write_54a8top(top_name, itp_name, posre_name = "", name = 'merged', N = 1)
	''' Write a 54a8 top file (called top_name), from a topology
	    object itp_file, putting N molecules of type name in the
	    [ molecules ] directory '''


# MANIPULATION functions

top_remove(topology, removendur_atm_nr)
	''' Remove an atom (removendur_atom_nr), identified by atom number,
	    from a topology object, modifying the object directly.
	    Take care of removing bonded interactions not present any more
	    and rescale the atoms number accounting for the removal.
	    Used in do_peptide_bond(). '''


top_rescale(topology, scale_atm = 0, scale_chg = 0):
	''' Rescale object topology renumbering the atoms by adding scale_atm
	    and the partial charged by adding scale_chg.
	    Used in do_peptide_bond(). '''


merge_scaled_topology(list_top, out_name = "merged"):
	''' Merge two topologies in one unique object. Note that
	    the second topology must already be scaled. This function
	    does NOT do that. '''


bonds_peptide_bond(top, aa, bb) [not intended as a user function]
	''' Write to topology top the bonded parameters (bonds, pairs,
	    angles, dihedrals) to create a peptide bond between atom
	    aa and bb, indicated by atom number.
	    Used in do_peptide_bond(). '''


def do_peptide_bond(topN, topC, atomN, atomC, topN_remove, topC_remove)
	''' Merge two topologies via a peptide bond: topN is the topology object
	    hosting the atom to be used as N-ter, topC the one as C-ter, atomN and
	    atomC are their atom numbers in the respective topologies.
	    topN_remove the extra hydrogen/atoms to be removed when creating the
	    bonds, indicated by atom number in topN. Similar for topC_remove in
	    topC.
	    It removes extra atoms, adjust terminal charges, rescale topology accordingly,
	    merge topologies, add bonded paramters (bonds, pairs, angles, dihedrals) '''


-----------------------------------------------------------------------------------
MDmutate/mutate_functions

Requirements: numpy, copy, pymol

# Functions

setNames(obj)
	''' Retrieve the list of amino acids from a pdb pymol object '''


getAa(cathegory)
	''' Return list of all amino acids ("all"), positively charged
	    ones ("plus"), negatively charged ones ("minus"),
	    polar ("polar") and hydrophobic ("hypho") '''


getOne(three)
	''' Returns the one letter amino acid code from the three letters one '''


minimize(selection='all', forcefield='MMFF94s', method='Conjugate Gradients', nsteps0= 500, conv=0.0001, cutoff=False, cut_vdw=6.0, cut_elec=8.0)
	''' Energy minimize a structure in pymol '''


pymutate(pdb, new_seq, nametag, ter = [-1])
	''' Mutate a pymol pdb object according to a new sequence.
	    Terminals are charged by default (-1 value) '''



