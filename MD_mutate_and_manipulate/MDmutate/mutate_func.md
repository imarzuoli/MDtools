# MDmutate/mutate_functions

Requirements: numpy, copy, pymol

---

## Functions

```
setNames(obj)

	''' Retrieve the list of amino acids from a pdb pymol object '''
```

```
getAa(cathegory)

	''' Return list of all amino acids ("all"), positively charged
	    ones ("plus"), negatively charged ones ("minus"),
	    polar ("polar") and hydrophobic ("hypho") '''
```

```
getOne(three)
	''' Returns the one letter amino acid code from the three letters one '''
```

```
minimize(selection='all', forcefield='MMFF94s', method='Conjugate Gradients', nsteps0= 500, conv=0.0001, cutoff=False, cut_vdw=6.0, cut_elec=8.0)

	''' Energy minimize a structure in pymol '''
```

```
pymutate(pdb, new_seq, nametag, ter = [-1])

	''' Mutate a pymol pdb object according to a new sequence.
	    Terminals are charged by default (-1 value) '''
```
