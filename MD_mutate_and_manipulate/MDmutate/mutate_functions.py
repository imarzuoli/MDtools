from pymol import cmd
from pymol import stored

# Functions
# Get original names
def setNames(obj):
    orig_sequence = {}
    cmd.load(obj) 
    cmd.select("prot", "name ca")
    cmd.do("stored.names = []")
    cmd.do("iterate (prot), stored.names.append((resi, resn))")
    for i in stored.names:
        orig_sequence[i[0]] = i[1] 
    return stored.names
    cmd.do('delete all') 

# Get one letter aa code
def getOne(three):
    trans = { 
       'ALA':'A',
       'ARG':'R',
       'ASN':'N',
       'ASP':'D',
       'CYS':'C',
       'GLU':'E',
       'GLN':'Q',
       'GLY':'G',
       'HIS':'H',
       'ILE':'I',
       'LEU':'L',
       'LYS':'K',
       'MET':'M',
       'PHE':'F',
       'PRO':'P',
       'SER':'S',
       'THR':'T',
       'TRP':'W',
       'TYR':'Y',
       'VAL':'V'
       } 
    return trans[three]

# Dictionary of AA
def getAa(cathegory):
	''' Return list of all amino acids ("all"), positively charged
	    ones ("plus"), negatively charged ones ("minus"),
	    polar ("polar") and hydrophobic ("hypho") '''
	
	aa_dic = {
	  "all": ['ALA', 'ARG', 'ASN', 'ASP',
	                        'CYS', 'GLU', 'GLN', 'GLY',
	                        'HIS', 'ILE', 'LEU', 'LYS',
	                        'MET', 'PHE', 'PRO', 'SER',
	                        'THR', 'TRP', 'TYR', 'VAL' ],
	  "plus": ['ARG', 'HIS', 'LYS' ],
	  "minus": ['ASP', 'GLU' ],
	  "polar": ['ASN', 'CYS', 'GLN', 'SER',
	                        'THR', 'TYR' ],
	  "hypho": ['ALA', 'GLY', 'ILE', 'LEU',
	                        'MET', 'PHE', 'PRO', 'TRP',
	                        'VAL' ]
	}
	return aa_dic[cathegory]

                        
# Minimization
def minimize(selection='all', forcefield='MMFF94s', method='Conjugate Gradients', nsteps0= 500, conv=0.0001, cutoff=False, cut_vdw=6.0, cut_elec=8.0):
    ''' Energy minimize a structure in pymol '''
    
    name = cmd.get_legal_name(selection)
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('pdb', 'pdb')
    mol = ob.OBMol()
    obconversion.ReadString(mol, pdb_string)
    mol.AddHydrogens()
    ff = ob.OBForceField.FindForceField(forcefield) ## GAFF, MMFF94s, MMFF94, UFF, Ghemical
    ff.Setup(mol)
    if cutoff == True:
        ff.EnableCutOff(True)
        ff.SetVDWCutOff(cut_vdw)
        ff.SetElectrostaticCutOff(cut_elec)
    if method == 'Conjugate Gradients':
        ff.ConjugateGradients(nsteps0, conv)
    else:
        ff.SteepestDescent(nsteps0, conv)
    ff.GetCoordinates(mol)
    nrg = ff.Energy()
    pdb_string = obconversion.WriteString(mol)
    cmd.delete(name)
    if name == 'all':
        name = 'all_'
    cmd.read_pdbstr(pdb_string, name)
    print('#########################################')
    print('The Energy of %s is %8.2f %s       '  % (name, nrg, ff.GetUnit()))
    print('#########################################')
    return nrg


# Perform mutation in pymol
def pymutate(pdb, new_seq, nametag, ter = [-1]):
    ''' Mutate a pymol pdb object according to a new sequence.
        Terminals are charged by default (-1 value) '''
    
    cmd.load(pdb)
    orig_sequence = setNames(pdb)
   
    if len(new_seq) != len(orig_sequence):
        raise ValueError('Length of mutated sequence not matching original one')
        
    cmd.do('wizard mutagenesis')
    cmd.do('refresh_wizard')

    for new in range(len(orig_sequence)):
    
        if new_seq[new] != orig_sequence[new][1]:
            if new+1 in ter:
                cmd.get_wizard().set_n_cap("posi")
            cmd.get_wizard().set_mode(new_seq[new])
            cmd.get_wizard().do_select(str(orig_sequence[new][0])+'/')
            cmd.frame(1)
            cmd.get_wizard().apply()
    
    # Protonation of the N.
    #cmd.edit("n%d" % int(site), None, None, None, pkresi=0, pkbond=0)
    
     #           # Protonation of the C.
      #          cmd.do("select c%d, name c and %d/" % (int(site), int(site)))
       #         cmd.edit("c%d" % int(site), None, None, None, pkresi=0, pkbond=0)
        #        cmd.do("h_fill") 

    seq_onelett = ''
    for lett in range(len(orig_sequence)):
        seq_onelett = seq_onelett + getOne(new_seq[lett])

    name = "./%s_%s.pdb" % (seq_onelett, nametag)
    cmd.save(name)
    cmd.do('delete all') 
    cmd.set_wizard('done')
    
    return name
