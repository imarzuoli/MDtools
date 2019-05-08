MDtools - Collections of scripts for MD simulations set up and analysis


# an_min_hist
Relying on a MDAnalysis backend, the script computes the histogram of the distances between the oxygen of water molecules (default name OW, customisable) and a group defined in any MDAnalysis compatible way (given by the selection option). The histogram is computed at each frame of the trajectory and the results averaged. Options allow to control the number of bins of the histograms and the maximum distance to analyse (to avoid spending time in the bulk water regions).
[Requirements: numpy, MDAnalysis]


# contacts_analysis.py
Collection of functions to effectively analyse the number of contacts between different subset of atoms and store the structural information for each contact for further analysis.
[Requirements: numpy, MDAnalysis, panda]


# delete_water
Delete water molecules in a given slice of space (perpendicular to either x, y or z). A molecule is erased if at least one of its atom is inside the selected region. Default is for atomistic description, options allows to change it for CG descriptions. 
[Requirements: numpy]


# MD_mutate_and_manipulate
Collection of scripts to manipulate peptide files, in order to easily create topology files
when a peptide bond is inserted, to manipulate (reorder) a pdb file if the topology
describes the molecule with a different ordering of the atoms, and to introduce mutations
in an automated way.
[Requirements: numpy, copy, pymol (for the mutation part only)]


# renumber_script
Reshuffle the atoms of a pdb. file according to a topology file (reading the [atoms] section of a top or tip file). at present it does NOT renumber them afterwards - this can easily be obtained via vmx genconf for example. If the pdb file contains multiple copies of the same molecule it loops through the topology information until all the pd. lines are processed. It generates an error and exit if the nr of atoms in the pd. are not a multiple of the ones in the topology.
[Requirements: numpy]

