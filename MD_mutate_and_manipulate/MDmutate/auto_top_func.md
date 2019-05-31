# MDmutate/auto_top_functions

Requirements: copy

---

## Objects

```
topology
	Attributes:
	info.name		[string]
	info.nrexcl		[int]
	info.atoms		[atom_top object]
	info.connect	[connect object]
	info.compound	[string]
	info.nmols		[int]
	
	Methods:
	atmnr(info): returns number of atoms in topology
	cgnr(info): returnsnumber of charge groups in topology

atom_top [not intended as a user class]
	Attributes:
	atomstab.nr		[int]
	atomstab.type		[string]
	atomstab.resnr		[int]
	atomstab.resid		[string]
	atomstab.atom		[int]
	atomstab.cgnr		[int]
	atomstab.charge	[float]
	atomstab.mass		[float]

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
	pp.funct		[string]

angle [not intended as a user class]
	Attributes:
	ang.at1		[int]
	ang.at2		[int]
	ang.at3		[int]
	ang.funct		[int = 2]
	ang.info		[string]

dih [not intended as a user class]
	Attributes:
	dih.at1		[int]
	dih.at2		[int]
	dih.at3		[int]
	dih.at4		[int]
	dih.funct		[int]
	dih.info		[string]
```

---

## Functions 

## IO functions

```
empty_top() [not intended as a user function]

	''' Return an empty topology object initializing all the fields to empty lists/NA'''
```

```
info_reader(stored, mode, line) [not intended as a user function]

	''' Called by read_top(), not to be used by user directly.
	    Reads topology line (line) according to the section it
	    belongs to (mode), and modify the corresponding field
	    of the topology which is going to be written (stored). '''
```

```
read_top(top_file)	

	''' Read topology (top/itp) file and store it in a topology object '''
```

```
write_itp(topology, out_itp_file, name = 'merged')

	''' Write out_itp_file from a topology object. If no name is
	    provided, the molecule name in [ moleculetype ] will be
	    'merged' '''
```

```
write_54a8top(top_name, itp_name, posre_name = "", name = 'merged', N = 1)

	''' Write a 54a8 top file (called top_name), from a topology
	    object itp_file, putting N molecules of type name in the
	    [ molecules ] directory '''
```


## MANIPULATION functions

```
top_remove(topology, removendur_atm_nr)

	''' Remove an atom (removendur_atom_nr), identified by atom number,
	    from a topology object, modifying the object directly.
	    Take care of removing bonded interactions not present any more
	    and rescale the atoms number accounting for the removal.
	    Used in do_peptide_bond(). '''
```

```
top_rescale(topology, scale_atm = 0, scale_chg = 0):

	''' Rescale object topology renumbering the atoms by adding scale_atm
	    and the partial charged by adding scale_chg.
	    Used in do_peptide_bond(). '''
```

```
merge_scaled_topology(list_top, out_name = "merged"):

	''' Merge two topologies in one unique object. Note that
	    the second topology must already be scaled. This function
	    does NOT do that. '''
```

```
bonds_peptide_bond(top, aa, bb) [not intended as a user function]

	''' Write to topology top the bonded parameters (bonds, pairs,
	    angles, dihedrals) to create a peptide bond between atom
	    aa and bb, indicated by atom number.
	    Used in do_peptide_bond(). '''
```

```
do_peptide_bond(topN, topC, atomN, atomC, topN_remove, topC_remove)

	''' Merge two topologies via a peptide bond: topN is the topology object
	    hosting the atom to be used as N-ter, topC the one as C-ter, atomN and
	    atomC are their atom numbers in the respective topologies.
	    topN_remove the extra hydrogen/atoms to be removed when creating the
	    bonds, indicated by atom number in topN. Similar for topC_remove in
	    topC.
	    It removes extra atoms, adjust terminal charges, rescale topology accordingly,
	    merge topologies, add bonded paramters (bonds, pairs, angles, dihedrals) '''
```
