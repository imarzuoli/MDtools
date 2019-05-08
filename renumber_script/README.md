# renumber_script

Script
Reshuffle the atoms of a pdb. file according to a topology file (reading the [atoms] section of a top or tip file). at present it does NOT renumber them afterwards - this can easily be obtained via vmx genconf for example. If the pdb file contains multiple copies of the same molecule it loops through the topology information until all the pd. lines are processed. It generates an error and exit if the nr of atoms in the pd. are not a multiple of the ones in the topology.

usage: renumber_script.py [-h] --pdb PDB_FILE --top TOP_FILE [--out OUT_FILE]

  -h, --help       		show this help message and exit
  --pdb PDB_FILE  	Coordinate file - PDB 3.0 name standard
  --top TOP_FILE  	Topology ot itp file - only the [atoms] section required
  --out OUT_FILE  	Output file
  
  Requirements: numpy  
 