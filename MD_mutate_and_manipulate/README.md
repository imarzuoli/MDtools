# MD_mutate_and_manipulate

Collection of scripts and functions to manipulate pdb and topology files.

In particular:

* [MDmutate/auto_top_functions] easily creates new topology files when an extra
   peptide bond is inserted between two proteins;
* [MDmutate/renumber_functions] to manipulate (reorder) a pdb file if the topology
   describes the molecule with a different ordering of the atoms, and convert atom names
   with mixed letters and numbers e.g. from 2HH1 to HH12;
* [MDmutate/mutate_functions] introduces mutations in a structures in an automated way
   generating one pdb file for each ne sequence.

-----------------------------------------------------------------------------------
### Scripts:

```example_create_peptide_bonds.py```

The script create the topology for a multi branched peptide inserting peptide bonds
at the required positions from the topology of the single pieces.
It requires the index of the atoms to be joined, and the ones to be erased (e.g. extra
hydrogens).

It takes care of:
- rescaling the index in one topology and append it to the first one;
- erase the extra atoms and rescale the remaining ones to have a continuous
  numbering in all the topology sections (i.e. erase their atom information, but also
  all the bonds/angles... where they appear;
- add the bonds/angles... information for the peptide bond created. For regular atom
  types, it retrieves automatically the appropriate bonds/angles/... label in a gromos
  54a8 like fashion (gb*/ga*/...)
  
The default run calls input/RKGB_NH3ter.top and input/RRWTWE.top and produces
a arm_mutations folder (as in the example_results folder).

-----------------------------------------------------------------------------------

```example_mutate_peptide.py```

[Backend pymol]

Given a list with the mutations proposed for each amino acids in the sequence,
creates all the possible combinations of sequence from these mutations and
generated the pd. files.

!!! At present, no search for the best rotamer is implemented !!!

The default run calls input/RRWTWE_original.pdb and produces
triskelion.itp and triskelion.top (as in the example_results folder).

-----------------------------------------------------------------------------------

## For details on the functions of each section, see the md files in the MDmutate folder.
