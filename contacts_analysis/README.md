# contacts_analysis

Collection of functions to effectively analyse the number of contacts between different subset of atoms and store the structural information for each contact for further analysis.

Requirements: numpy, MDAnalysis, panda

### List of user intended functions

computeContacts(universe, sel1 = "name CA", sel2 = None, cutoff = 6, lifetime = 0.5, delta_frames = 1)	
```	
	Compute contacts within cutoff occurring between atoms within a set selected by sel1
	or between atoms of sel1 and toms of sel2 (if given). Contacts are filtered for the
	ones present more than a given fraction (lifetime).
	Trajectory is analysed every delta_frames.
	Returns a list contaiing the number of contacts in time and a list with information
	for every contact found.
	
	
	Parameters
	------------
	universe    : MDAnalysis.universe
				Path to topology file (any accepted by MDAnalysis)
	sel1        : string
				String to select subset 1 of atoms (default "name CA")
	sel2        : string
				String to select subset 2 of atoms (default None)
	cutoff      : float
				Cutoff distance for contacts, in Angstrom (default 6)
	lifetime    : float
				Threshold to filter contacts appearing less than a fraction of time frames
				analysed. Values between 0 and 1, default 0.5
	delta_frames: int
				Stride on the trajectory frames. Default 1 (all frames are analysed)
	
	Returns
	------------
	list [np.array, shape (n_frames_analysed,2) with time and number of contacts,
	list shape (nr_all_contacts - surviving more than the lifetime threshold, 12)
	listing for each contact
			[0] Pair_nr, 		[1] A1_pdb_id, 		[2] A1_resname, [3] A1_resid, \n \
			[4] A1_resnumber, 	[5] A1_mol, 		[6] A2_pdb_id, 	[7] A2_resname, \n \
			[8] A2_resid, 		[9] A2_resnumber, 	[10] A2_mol, 	[11] occupational_time")
```
	

getContactsByRes(info_contact_table)
```
	Given a list of contacts (as per computeContacts() output), computes the number of
	contacts between different amino acids, breaking them by amino acid type.
	Returns a panda data frame of size NxN, with N the different types of amino acid
	present in the input list of contacts.
	
	Parameters
	------------
	info_contact_table: list of contacts as per computeContacts() output
	
	Returns
	------------
	list [np.array, shape (n_frames_analysed,2) with time and number of contacts,
	list shape (nr_all_contacts - surviving more than the lifetime threshold, 12)
	listing for each contact
			[0] Pair_nr, 		[1] A1_pdb_id, 		[2] A1_resname, [3] A1_resid, \n \
			[4] A1_resnumber, 	[5] A1_mol, 		[6] A2_pdb_id, 	[7] A2_resname, \n \
			[8] A2_resid, 		[9] A2_resnumber, 	[10] A2_mol, 	[11] occupational_time")
```
	