#!/usr/bin/env python3

import MDAnalysis as mda
import argparse
import sys
import numpy as np
import pandas as pd
from MDAnalysis.analysis.distances import contact_matrix


# chainsToMols                if repeated molecules of fixed aa length but no info in pdb, assign
#                             chain nr based on breaks by such length
# computeContacts             according to cutoff, between atoms in same trajectory;
#                             produce list with time series and list with all contact info
# computeMixedContacts        as above but between atoms of two different trajectories
# filter_by_type              filter contact list for type of res1 res2
# plot_contacts_bars          plot bars of contact per restype, broken by restype
# readable_output             transforms list of contact info in nice data frame to print
# sanitizeContacts            remove contacts in same molecule occurring along the chain, 1
#                             and 2 contacts apart
# table_by_residueType        from list of all contact and info, get a plottable table of contacts
#                             by residue


def computeContacts(universe, sel1 = "name CA", sel2 = None, cutoff = 6, lifetime = 0.5, delta_frames = 1, resmol = None, collapse = True):
	
	"""
	
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
	resmol      : int
				Nr of residues in a molecule, to identify the molecule ID. Useful when
				residues are listed sequentially, not restarting from 1 at each molecule.
				If none, the number of residues with different ID in selection 1 is taken.
				Molecule ID is meaningful only if the system is made of N copies of the
				same molecule. Default None
	collapse    : bool
				If True, removes contacts linking the same pair of residues, keeping only
				the first one detected
	
	Returns
	------------
	list [np.array, shape (n_frames_analysed,2) with time and number of contacts,
	list shape (nr_all_contacts - surviving more than the lifetime threshold, 12)
	listing for each contact
			[0] Pair_nr, 		[1] A1_pdb_id, 		[2] A1_resname, [3] A1_resid, \n \
			[4] A1_resnumber, 	[5] A1_mol, 		[6] A2_pdb_id, 	[7] A2_resname, \n \
			[8] A2_resid, 		[9] A2_resnumber, 	[10] A2_mol, 	[11] occupational_time")
	
	"""
	
	if type(sel1) is str:
		grp1 = universe.select_atoms(str(sel1))
	else:
		grp1 = sel1
	
	if (sel2 is None):
		[contacts_series, info_contact_table] = _computeSameGroupContacts(universe, grp1, cutoff, lifetime, delta_frames, resmol, collapse)
		
	else:
		if type(sel2) is str:
		
			grp2 = universe.select_atoms(str(sel2))
		else:
			grp2 = sel2
			
		if grp1 == grp2:
			[contacts_series, info_contact_table] = _computeSameGroupContacts(universe, grp1, cutoff, lifetime, delta_frames, resmol, collapse)
		else:
			[contacts_series, info_contact_table] = _computeTwoGroupsContacts(universe, grp1, grp2, cutoff, lifetime, delta_frames, resmol, collapse)
		
	return [contacts_series, info_contact_table]



def getContactsByRes(info_contact_table):
	
	"""
	
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
	
	"""
	restype_atm1 = [el[2] for el in info_contact_table] 
	restype_atm2 = [el[7] for el in info_contact_table] 
	
	aa_list = set(restype_atm1 + restype_atm2)
	
	contacts_by_res = pd.DataFrame(data = 0, index = aa_list, columns = aa_list)
	
	for item in info_contact_table:
		contacts_by_res.loc[item[2]][item[7]] += 1
		
	return contacts_by_res



def _computeSameGroupContacts(universe, grp1, cutoff, lifetime, delta_frames, resmol, collapse):
	
	print("\nComputing contacts in trajectory:\n%s\nwithin group of atoms %s (%d instances), using a %f A cutoff" \
						% (universe.trajectory.filename, list(set(grp1.names)), len(grp1.names), cutoff) )
	print("")
	
	sliced_traj = universe.trajectory[::delta_frames]
	steps = int(np.ceil(universe.trajectory.n_frames/(delta_frames*1.0)))
	
	# Compute distances
	contacts = []
	contacts_series = []
	
	distmat = np.zeros(( steps, len(grp1), len(grp1) ))
	
	for ts in sliced_traj:
		sys.stdout.write("Computing distances in frame %d/%d   \r" % (ts.frame, universe.trajectory.n_frames) )
		sys.stdout.flush()
		distmat[ ts.frame/delta_frames ] = contact_matrix(coord = grp1.positions, cutoff = cutoff)
	
	sys.stdout.write("Computing distances in frame %d/%d\n" % (universe.trajectory.n_frames, universe.trajectory.n_frames) )
	
	# Filtering for occupancy < lifetime
	print("Filtering out pairs with lifetime < " + str(lifetime*100) + "% frames analysed..." )
	
	occupancymat = distmat.sum(0) / ( 1.0 * steps )
	
	# Computing contacts timeseries
	contacts_series = np.zeros(( steps ))
	for ts in universe.trajectory[::delta_frames]:
		sys.stdout.write("Computing total contacts in frame %d/%d     \r" % (ts.frame, universe.trajectory.n_frames) )
		sys.stdout.flush()
		contacts_series[ts.frame/delta_frames] = np.sum(distmat[ts.frame/delta_frames])		# sistema
	sys.stdout.write("Computing total contacts in frame %d/%d     \n" % (universe.trajectory.n_frames, universe.trajectory.n_frames) )
	
	times = [ts.time for ts in universe.trajectory[::delta_frames]]
	contacts_series = np.vstack((times, contacts_series)).transpose()
	
	# Assign pair information
	
	bool_occupancy = occupancymat > lifetime
	nr_kept = np.sum(bool_occupancy)
	
	info_contact_table = []
	if (resmol is None):
		resinmol = max(grp1.atoms.resids)
	else:
		resinmol = resmol
	
	pair_nr = 0
	
	for atom1 in range(len(grp1)):
		for atom2 in range(atom1+1,len(grp1)):
			if bool_occupancy[atom1,atom2]:
				
				pair_nr = pair_nr + 1
				sys.stdout.write("Assign pairs information: %d/%d     \r" % (pair_nr, nr_kept) )
				sys.stdout.flush()
				
				tmp_pair = [ pair_nr , grp1.atoms.ids[int(atom1)], grp1.atoms.resnames[int(atom1)], grp1.atoms.resids[int(atom1)], \
							grp1.atoms.resindices[int(atom1)] + 1, (grp1.atoms.resindices[int(atom1)] + 1)/resinmol + 1, \
							grp1.atoms.ids[int(atom2)], grp1.atoms.resnames[int(atom2)], grp1.atoms.resids[int(atom2)], \
							grp1.atoms.resindices[int(atom2)] + 1, (grp1.atoms.resindices[int(atom2)] + 1)/resinmol + 1, \
							occupancymat[atom1,atom2] ]
				
				# Excluding subsequent or two residues apart in a chain
				if ( tmp_pair[5] != tmp_pair[10] or tmp_pair[3] + 2 < tmp_pair[8] ):
					
					# If collapse is True, checking if the contact was already inserted
					if collapse:
						# Check if pair already exists
						new = 1
						for item in info_contact_table:
							if (tmp_pair[4] == item[4] and tmp_pair[9] == item[9]):
								new = 0
						if new == 1:
							info_contact_table.append(tmp_pair)
					else:
						info_contact_table.append(tmp_pair)				
	
	sys.stdout.write("Assigning pairs information: %d/%d     \n" % (nr_kept, nr_kept) )
	
	print("\n \
OUTPUT [contacts_series, info_contacts_table] \n \
	contacts_series:		np.array, shape (n_frames,2) \n \
	info_contacts_table:	list, shape (nr_all_contacts, 12) \n \
		for each contact: \n \
			[0] Pair_nr, 		[1] A1_pdb_id, 		[2] A1_resname, [3] A1_resid, \n \
			[4] A1_resnumber, 	[5] A1_mol, 		[6] A2_pdb_id, 	[7] A2_resname, \n \
			[8] A2_resid, 		[9] A2_resnumber, 	[10] A2_mol, 	[11] occupational_time")
	
	return [contacts_series, info_contact_table]



def _computeTwoGroupsContacts(universe, grp1, grp2, cutoff, lifetime, delta_frames, resmol, collapse):
	
	print("\nComputing contacts in trajectory\n%s,\nbetween group of atoms %s (%d instances) and atoms %s (%d instances), using a %f A cutoff" \
						% (universe.trajectory.filename, list(set(grp1.names)), len(grp1.names), list(set(grp2.names)), len(grp2.names), cutoff) )
	print("")
	
	sliced_traj = universe.trajectory[::delta_frames]
	steps = int(np.ceil(universe.trajectory.n_frames/(delta_frames*1.0)))
	
	# Compute distances
	contacts = []
	contacts_series = []
	
	total = grp1 + grp2
	
	distmat = np.zeros(( steps, len(grp1), len(grp2) ))
	
	for ts in sliced_traj:
		sys.stdout.write("Computing distances in frame %d/%d   \r" % (ts.frame, universe.trajectory.n_frames) )
		sys.stdout.flush()
		distmat[ ts.frame/delta_frames ] = contact_matrix(coord = total.positions, cutoff = cutoff)[0:len(grp1), len(grp1):]
	
	sys.stdout.write("Computing distances in frame %d/%d\n" % (universe.trajectory.n_frames, universe.trajectory.n_frames) )
	
	# Filtering for occupancy < lifetime
	print("Filtering out pairs with lifetime < " + str(lifetime*100) + "% frames analysed..." )
	
	occupancymat = distmat.sum(0) / ( 1.0 * steps )
	
	# Computing contacts timeseries
	contacts_series = np.zeros(( steps ))
	for ts in universe.trajectory[::delta_frames]:
		sys.stdout.write("Computing total contacts in frame %d/%d     \r" % (ts.frame, universe.trajectory.n_frames) )
		sys.stdout.flush()
		contacts_series[ts.frame/delta_frames] = np.sum(distmat[ ts.frame/delta_frames, :, : ])
	sys.stdout.write("Computing total contacts in frame %d/%d     \n" % (universe.trajectory.n_frames, universe.trajectory.n_frames) )
	
	times = [ts.time for ts in universe.trajectory[::delta_frames]]
	contacts_series = np.vstack((times, contacts_series)).transpose()
	
	# Assign pair information
	
	bool_occupancy = occupancymat > lifetime
	nr_kept = np.sum(bool_occupancy)
	
	info_contact_table = []
	if (resmol is None):
		resinmol = max(total.atoms.resids)
	else:
		resinmol = resmol
	
	pair_nr = 0
	
	for atom1 in range(len(grp1)):
		for atom2 in range(len(grp2)):
			if bool_occupancy[atom1,atom2]:
				
				pair_nr = pair_nr + 1
				sys.stdout.write("Assign pairs information: %d/%d     \r" % (pair_nr, nr_kept) )
				sys.stdout.flush()
				
				tmp_pair = [ pair_nr , total.atoms.ids[int(atom1)], total.atoms.resnames[int(atom1)], total.atoms.resids[int(atom1)], \
							total.atoms.resindices[int(atom1)] + 1, (total.atoms.resindices[int(atom1)] + 1)/resinmol + 1, \
							total.atoms.ids[int(atom2)+len(grp1)], total.atoms.resnames[int(atom2)+len(grp1)], total.atoms.resids[int(atom2)+len(grp1)], \
							total.atoms.resindices[int(atom2)+len(grp1)] + 1, (total.atoms.resindices[int(atom2)+len(grp1)] + 1)/resinmol + 1, \
							occupancymat[atom1,atom2] ]
				
				# Excluding subsequent residues in a chain
				if ( tmp_pair[5] != tmp_pair[10] or tmp_pair[3] + 2 < tmp_pair[8] ):
					
					# If collapse is True, checking if the contact was already inserted
					if collapse:
						# Check if pair already exists
						new = 1
						for item in info_contact_table:
							if (tmp_pair[4] == item[4] and tmp_pair[9] == item[9]):
								new = 0
						if new == 1:
							info_contact_table.append(tmp_pair)
					else:
						info_contact_table.append(tmp_pair)
	
	sys.stdout.write("Assigning pairs information: %d/%d     \n" % (nr_kept, nr_kept) )
	
	print("\n \
OUTPUT [contacts_series, info_contacts_table] \n \
	contacts_series:		np.array, shape (n_frames,2) \n \
	info_contacts_table:	list, shape (nr_all_contacts, 12) \n \
		for each contact: \n \
			[0] Pair_nr, 		[1] A1_pdb_id, 		[2] A1_resname, [3] A1_resid, \n \
			[4] A1_resnumber, 	[5] A1_mol, 		[6] A2_pdb_id, 	[7] A2_resname, \n \
			[8] A2_resid, 		[9] A2_resnumber, 	[10] A2_mol, 	[11] occupational_time")
	
	return [contacts_series, info_contact_table]
