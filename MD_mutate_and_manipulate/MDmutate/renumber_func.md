# MDmutate/renumber_functions

Requirements: numpy

---

## Functions

```
format_pdb(data)

	''' Adjust spacing of a list of 9 (data) and
	    join them to return a PDB formatted string '''
```

```
pdb_line(line) [not intended as a user function]

	''' Take a pdb line (line) and split the fields according to the PDB format.
	    Return a list of length 9 '''
```

```
adj_names(pdb_file, out_file)

	''' In a PDB file (pdb_file), if atomname starts with a number
	    e.g. 1HH2, it puts all the figures at the end
	    e.g. HH21. Output to out_file. '''
```

```
reshuffle_top(top_file, pdb_file, out_file, strict = True)

	''' Reshuffle the order of a PDB file (pdb_file) according to an ITP
	    topology file (top_file). If strict is TRUE it does it only if 
	    the number of atoms matches. If FALSE, allows for different
	    numbers taking only the atoms listed in the topology (useful
	    if PDB has extra hydrogens with respect to GROMACS
	    topology): please double check your output! '''
```
