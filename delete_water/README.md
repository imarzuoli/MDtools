# delete_water

Script
Delete water molecules in a given slice of space (perpendicular to either x, y or z). A molecule is erased if at least one of its atom is inside the selected region. Default is for atomistic description, options allows to change it for CG descriptions. 

```
usage: delete_water.py [-h] --pdb PDB_FILE [--out OUT_FILE] --min MIN --max MAX [--dim DIM] [--names NAMES] [--resname RES]

  -h, --help      show this help message and exit
  --pdb  PDB_FILE 	Coordinate file - PDB 3.0 name standard
  --out   OUT_FILE 	Output file - default dehydrated.pdb
  --min  MIN            	min (A) from which to erase
  --max MAX          	max (A) up which to erase
  --dim  DIM        	perpendicular dimension (x, y or z) - default z
  --names NAMES     string of names of water atoms within a molecule. Default 'OW HW1 HW2' (united atom, GROMOS style)
  --resname RES       string with residue name of solvent, default 'SOL'
```
  
Requirements: numpy
