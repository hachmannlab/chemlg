from __future__ import print_function
import sys
import pybel
from openbabel import OBAtomAtomIter
#import subprocess
import time
import scipy
from collections import defaultdict
import operator
#import cStringIO
from mpi4py import MPI
import os
import threading
from itertools import islice,chain


###################################################################################################

class OutputGrabber(object):
    """
    Class used to grab standard output or another stream.
    """
    def __init__(self, stream=None, threaded=False):
        self.origstream = stream
        self.origstream3 = stream
        self.origstream2 = sys.stdout
        self.origstreamfd = self.origstream.fileno()
        self.capturedtext = ""
        # Create a pipe so the stream can be captured:
        self.pipe_out, self.pipe_in = os.pipe()
        pass

    def start(self):
        """
        Start capturing the stream data.
        """
        self.capturedtext = ""
        # Save a copy of the stream:
        self.streamfd = os.dup(self.origstreamfd)
        # Replace the Original stream with our write pipe
        os.dup2(self.pipe_in, self.origstreamfd)
        pass

    def stop(self):
        """
        Stop capturing the stream data and save the text in `capturedtext`.
        """
        # Flush the stream to make sure all our data goes in before
        # the escape character.
        #self.origstream.flush()
        os.dup2(self.streamfd, self.origstreamfd)
        os.close(self.pipe_out)
        pass

## Openbabel sometime prints out an error related to sterochemistry. 
## This is usully the case when the molecule do not have any sterechemistry information.
## Because this error is from OpenBabel which is written in C, it is quite challenging 
## to capture error statements. 
## OutputGrabber is a class which capture system errors from any program (including C0. 
## We will have to be careful with this class because it might also capture system errors
## not related to stereochemistry and we would never know what the system error was. 

out = OutputGrabber(sys.stderr)

def molecule(smiles, code):
    """Set docstring here.

    Parameters
    ----------
    smiles: string,
                smiles string of molecule
    code: string, 
                code for molecule

    Returns
    -------
    mol: dict, 
                dict with keys - smiles, code, can_smiles, reverse_smiles, mol_ob 
    """
    mol = {'smiles': smiles, 'code': code}
    obm = pybel.readstring("smi", smiles)
    mol['can_smiles'] = obm.write("can")
    obm_can = pybel.readstring("smi", mol['can_smiles'])
    mol['reverse_smiles'] = reverse_mol(obm_can, list(obm_can.atoms))
    return mol

def get_atom_pair_list(smiles, type1or2):
    mol=pybel.readstring('smi',smiles)
    atoms=list(mol)
    atom_pair_list=[]
    for atom in atoms:
        if atom.OBAtom.GetAtomicNum() ==88:            
            index=atom.OBAtom.GetIdx()
            # OBAtomAtomIter for iterating over neighboring atoms
            for atom1 in OBAtomAtomIter(mol.OBMol.GetAtom(index)):
                if atom1.GetAtomicNum() !=6: 
                    continue
                atom1_idx=atom1.GetIdx()
                atom_atoms=[]
                if type1or2==1:
                    # iterate over C's neighboring atoms (other C's in ring)
                    for atom2 in OBAtomAtomIter(atom1):
                        if atom2.GetAtomicNum() !=6 or atom2.IsInRing()==False: 
                            continue
                        hcount = atom2.ExplicitHydrogenCount() + atom2.ImplicitHydrogenCount()
                        if hcount==0:
                            continue
                        atom_atoms.append(atom2.GetIdx())
                    if len(atom_atoms)==2:
                        for aa in atom_atoms:
                            atom_pair=[atom1_idx,aa,index]      # [a1-index, a2-index, Ra-index]
                            atom_pair_list.append(atom_pair)
                    if len(atom_atoms)==1:
                        atom_pair_list.append([atom1_idx,atom_atoms[0],index])
                if type1or2==2:
                    # iterate over C's neighboring atoms (other C's in ring)
                    for atom2 in OBAtomAtomIter(atom1):
                        if atom2.GetAtomicNum() ==88 or atom2.IsInRing()==False: 
                            continue
                        hcount = atom2.ExplicitHydrogenCount() + atom2.ImplicitHydrogenCount()
                        if hcount==0:
                            continue
                        atom_atoms.append(atom2.GetIdx())
                    for idx in atom_atoms:
                        for atom2 in OBAtomAtomIter(atom1):
                            In_ring=False
                            atom_pair_list_tmp=[]
                            atom2_idx=atom2.GetIdx()
                            if atom2.GetAtomicNum()!=6 or atom2_idx==idx: 
                                continue
                            for atom3 in OBAtomAtomIter(atom2):
                                atom3_idx=atom3.GetIdx()
                                if atom3_idx==atom1_idx:
                                    continue
                                if atom3.IsInRing()==False:
                                    In_ring=True
                                    continue
                                hcount = atom3.ExplicitHydrogenCount() + atom3.ImplicitHydrogenCount()
                                if hcount!=0:
                                    atom_pair_list_tmp.append([atom1_idx,atom2_idx,idx,atom3_idx,index])
                            if In_ring==False:
                                atom_pair_list += atom_pair_list_tmp
                        
                            
    return atom_pair_list, len(atoms)

def get_fused(mol1, mol2, rules):
    """Fusion

    Parameters
    ----------
    mol1: molecule class object
    mol2: molecule class object
    rules: 

    Returns
    -------

    """
    list1, size1 = get_atom_pair_list(mol1['can_smiles'],1)
    list2, size2 = get_atom_pair_list(mol2['can_smiles'],2)
    smiles_combi = mol1['can_smiles'] + '.' + mol2['can_smiles']
    lib_can, lib_can_nRa = [], []
    code = mol1.code + ':' + mol2.code
    for item1 in list1:
        for item2 in list2:
            mol_combi= pybel.readstring("smi",smiles_combi)
            a1_to_set=mol_combi.OBMol.GetAtom(item1[0])
            a2_to_set=mol_combi.OBMol.GetAtom(item1[1])
            a3_to_set=mol_combi.OBMol.GetAtom(size1+item2[2])
            a4_to_set=mol_combi.OBMol.GetAtom(size1+item2[3])

            chng_arom=False
            if a3_to_set.IsAromatic()==True and a4_to_set.IsAromatic()==True:
                chng_arom=True
            bond_to_del=mol_combi.OBMol.GetBond(size1+item2[4],size1+item2[0])
            mol_combi.OBMol.DeleteBond(bond_to_del)
            bond_to_del=mol_combi.OBMol.GetBond(item1[0],item1[2])
            mol_combi.OBMol.DeleteBond(bond_to_del)
            bond_to_del=mol_combi.OBMol.GetBond(item2[0]+size1,item2[1]+size1)
            mol_combi.OBMol.DeleteBond(bond_to_del)
            bond_to_del=mol_combi.OBMol.GetBond(item2[0]+size1,item2[2]+size1)
            mol_combi.OBMol.DeleteBond(bond_to_del)
            bond_to_del=mol_combi.OBMol.GetBond(item2[1]+size1,item2[3]+size1)
            mol_combi.OBMol.DeleteBond(bond_to_del)
            atoms=list(mol_combi.atoms)
            
            mol_combi.OBMol.AddBond(item1[1],item2[2]+size1,1,0,-1)
            mol_combi.OBMol.AddBond(item1[0],item2[3]+size1,1,0,-1)
            
            
            mol_combi_new= pybel.readstring("smi",str(mol_combi))
            atoms_new=list(mol_combi_new)
            
            for atoms in atoms_new:
                index=atoms.OBAtom.GetIdx()
                neigh_atm=False
                for atom in OBAtomAtomIter(mol_combi_new.OBMol.GetAtom(index)):
                    neigh_atm=True
                if neigh_atm==False:
                    mol_combi_new.OBMol.DeleteAtom(mol_combi_new.OBMol.GetAtom(index))
            
            if a1_to_set.GetHeteroValence()==1 or a2_to_set.GetHeteroValence()==1:
                continue
                
            can_mol_combi = mol_combi_new.write("can")
            
            if can_mol_combi not in lib_can_nRa:
                mol_wt=str(int(mol_combi_new.OBMol.GetMolWt()))
                atoms=list(mol_combi_new.atoms)
                shd_add=True
                for atom in atoms:                
                    a =atom.OBAtom.CountBondsOfOrder(3)
                    b =atom.OBAtom.CountBondsOfOrder(2)
                    c =atom.OBAtom.CountBondsOfOrder(1)
                    tot_bnds= a*3+b*2+c*1
                    
                    if tot_bnds>4:
                        shd_add=False
                        
                if if_add(mol_combi_new,rules,code)==True and shd_add==True:
                    lib_can.append([str(can_mol_combi),mol_wt])

                for atom in atoms:
                    if atom.OBAtom.GetAtomicNum()==88:
                        atom.OBAtom.SetAtomicNum(1)
                
                lib_can_nRa.append(str(can_mol_combi))
    
                
            if chng_arom==True:
                
                if not a1_to_set.IsAromatic():
                    a1_to_set.SetAromatic()
                    a2_to_set.SetAromatic()
                mol_combi_new= pybel.readstring("smi",str(mol_combi))
                atoms_new=list(mol_combi_new)
                
                for atoms in atoms_new:
                    index=atoms.OBAtom.GetIdx()
                    neigh_atm=False
                    for atom in OBAtomAtomIter(mol_combi_new.OBMol.GetAtom(index)):
                        neigh_atm=True
                    if neigh_atm==False:
                        mol_combi_new.OBMol.DeleteAtom(mol_combi_new.OBMol.GetAtom(index))

                can_mol_combi = mol_combi_new.write("can")

                if can_mol_combi not in lib_can_nRa:
                    mol_wt=str(int(mol_combi_new.OBMol.GetMolWt()))
                    
                    atoms=list(mol_combi_new.atoms)
                    #print mol_combi_new,'before'
                    shd_add=True
                    for atom in atoms:
                        a =atom.OBAtom.CountBondsOfOrder(3)
                        b =atom.OBAtom.CountBondsOfOrder(2)
                        c =atom.OBAtom.CountBondsOfOrder(1)
                        tot_bnds= a*3+b*2+c*1
                        
                        if tot_bnds>4:
                            shd_add=False
                    #print if_add(mol_combi_new,mol_wt,rules,code,'f')
                    if if_add(mol_combi_new,rules,code)==True and shd_add==True:
                        lib_can.append([str(can_mol_combi),mol_wt])

                        #print mol_combi_new,'after'
                    for atom in atoms:
                        if atom.OBAtom.GetAtomicNum()==88:
                            atom.OBAtom.SetAtomicNum(1)
                    
                    lib_can_nRa.append(str(can_mol_combi))
    lib_can_c=[]
    for item in lib_can:
        lib_can_c.append([item[0][:-2], item[1], mol1['code'] + ':' + mol2[['code']]])
            
    return lib_can

def lipinski(mol):
    """Return the values of the Lipinski descriptors.

    Parameters
    ----------
    mol: 

    Returns
    -------

    """
   desc = {'molwt': mol.molwt,

      'HBD': len(HBD.findall(mol)),

      'HBA': len(HBA.findall(mol)),

      'logP': mol.calcdesc(['logP']) ['logP']}

   return desc

def if_add(molc, rules, code):
    """Determines if a molecule should be included in the library based on the generation rules given.

    Parameters
    ----------
    molc: string, 
            reversed canonical smiles of molecule whose validity is to be determined
    rules: dict, 
            dict of rules for library generation 
    code: string,
            molecule code

    Returns
    -------
    boolean: True/False
    """
    # check validity of molecule
    # pybel returns true even for wrong molecules by converting them to the correct strings.
    if molc == '':  return False
    if not pybel.readstring("smi", molc):
        print("invalid", pybel.readstring("smi", molc))
        return False
    # create a new copy of the dict to avoid changing the original dict    
    mol = molecule(molc, code)
    mol_ob = pybel.readstring("smi", mol['reverse_smiles'])
    add = True

    for atom in list(mol_ob.atoms):
        ## Removing Francium and Radium atoms
        if atom.OBAtom.GetAtomicNum()==87 or atom.OBAtom.GetAtomicNum()==88:
            atom.OBAtom.SetAtomicNum(1)

    if 'include_bb' in rules:
        for item in rules['include_bb']:
            if item not in code:
                add = False
                del mol
                return add
    
    if '2' in rules:
        bonds = mol_ob.OBMol.NumBonds()
        if not rules['2'][0] < bonds < rules['2'][1]:
            add = False
            del mol
            return add
    
    if '3' in rules:
        if not rules['3'][0] < len(list(mol_ob.atoms)) < rules['3'][1]:
            add = False
            del mol
            return add

    if '4' in rules:
        if not rules['4'][0] < int(mol_ob.OBMol.GetMolWt()) < rules['4'][1]:
            add = False
            del mol
            return add

    # Calculating no.of rings
    if '5' in rules:
        rings = len(mol_ob.OBMol.GetSSSR())
        if not rules['5'][0] < rings < rules['5'][1]:
            add = False
            del mol
            return add

    # Calculating no.of aromatic and non-aromatic rings
    if '6' in rules or '7' in rules:
        no_ar, no_non_ar= 0, 0
        for r in mol_ob.OBMol.GetSSSR():
            if r.IsAromatic():
                no_ar += 1
            else:
                no_non_ar += 1
        if '6' in rules:
            if not rules['6'][0] < no_ar < rules['6'][1]:
                add = False
                del mol
                return add
        if '7' in rules:
            if not rules['7'][0] < no_non_ar < rules['7'][1]:
                add = False
                del mol
                return add

    if '8' in rules:
        no_s_bonds=unique_structs(mol_ob,"*-*")
        if not rules['8'][0] < no_s_bonds < rules['8'][1]:
            add = False
            del mol
            return add

    if '9' in rules:
        no_d_bonds=unique_structs(mol_ob,"*=*")
        if not rules['9'][0] < no_d_bonds < rules['9'][1]:
            add = False
            del mol
            return add

    if '10' in rules:
        no_t_bonds=unique_structs(mol_ob,"*#*")
        if not rules['10'][0] < no_t_bonds < rules['10'][1]:
            add = False
            del mol
            return add
    
    if 'heteroatoms' in rules:
        for item in rules['heteroatoms']:
            no_at=0
            if item[0]=='C' or item[0]=='S' or item[0]=='N' or item[0]=='O' or item[0]=='c' or item[0]=='s' or item[0]=='n' or item[0]=='o':
                no_at=unique_structs(mol_ob, item[0].lower()) + unique_structs(mol_ob, item[0].upper())
            else:
                no_at=unique_structs(mol_ob,item[0])
            if no_at > item[1]:
                add = False
                del mol
                return add
    
    if 'lipinski' in rules:
        descriptors = lipinski(mol_ob)
        if not passes_all_rules(descriptors):
            add = False
            del mol
            return add

    if 'fingerprint' in rules:
        for mol_to_comp in rules['fingerprint']:
            mol2 = pybel.readstring('smi',mol_to_comp[0])
            tanimoto = mol_ob.calcfp()|mol2.calcfp()            
            if tanimoto < float(mol_to_comp[1]):
                add = False
                del mol
                return add

    if '14' in rules:
        for item in rules['14']:
            no_occ=unique_structs(mol_ob,item)
            if no_occ == 0 :
                add = False
                del mol
                return add 
    
    if '15' in rules:
        for item in rules['15']:
            no_occ=unique_structs(mol_ob,item)
            if no_occ>0 :
                add = False
                del mol
                return add    

    if 'bb_final_lib' in rules:                             # if this is present, then the value is false. 
        if mol['code'].count('-') == 0 and mol['code'].count(':') == 0:
            add = False
            del mol
            return add

    del mol
    return add

def create_link(mol1, mol2, rules):
    """This function creates all possible links between two given molecules.

    Parameters
    ----------
    mol1: molecule dict
    mol2: molecule dict
    rules: rules dict

    Returns
    -------
    library_two: list, list of molecule class objects
    """
    library_two = []
    mol1_ob = pybel.readstring("smi", mol1['reverse_smiles'])
    mol2_ob = pybel.readstring("smi", mol2['reverse_smiles'])
    mol1_index_list = get_index_list(mol1, list(mol1_ob.atoms))
    mol2_index_list = get_index_list(mol2, list(mol2_ob.atoms))
    
    smiles_combi = mol1['reverse_smiles'] + '.' + mol2['reverse_smiles']
    
    for index1 in mol1_index_list:
        for index2 in mol2_index_list:
            mol_combi= pybel.readstring("smi",smiles_combi)
            mol_combi.OBMol.AddBond(index1, index2 + len(mol1_ob.atoms), 1, 0, -1)
            can_mol_combi = mol_combi.write("can")
            code = mol1['code'] + '-' + mol2['code']
            if if_add(can_mol_combi, rules, code):
                temp = molecule(can_mol_combi, code)
                library_two.append(temp)
    return library_two

def get_index_list(mol, atoms):
    """Returns the list of index numbers of atoms that can be reacted in a molecule while making sure no duplicate sites are returned.

    Parameters
    ----------
    mol: molecule dict
    atoms: 

    Returns
    -------
    atoms_index: list
    """
    
    can_smiles_list=[]
    atoms_index=[]

    for atom in atoms:
        
        # Counting the number of hydrogens attached to the atom
        hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
        if hcount==0:
            continue # Do not do anything if there are no hydrogens attached
        
        newmol = pybel.readstring("smi", mol['reverse_smiles'])

        # Attach Francium atom. Makes it easy to remove duplicates
        newmol.OBMol.InsertAtom(Fratom) 
        # Get index of atom
        index=atom.OBAtom.GetIdx() 
        # Create a bond between Fr atom and the curent atom
        newmol.OBMol.AddBond(index, len(atoms)+1, 1, 0, -1) 
        
        out.start() # Capturing stderr from openbabel (C program) 
        # Making use of canonical smiles to remove duplicates
        can_smiles = newmol.write("can")
        out.stop() # Closing capture of stderr

        if can_smiles not in can_smiles_list:
            can_smiles_list.append(can_smiles)
            index=atom.OBAtom.GetIdx()
            atoms_index.append(index)
    
    return atoms_index


passes_all_rules = lambda desc: (desc ['molwt'] <= 500 and

         desc ['HBD'] <= 5 and desc ['HBA'] <= 10 and

         desc ['logP'] <= 5)

HBD = pybel.Smarts("[#7,#8;!H0]")

HBA = pybel.Smarts("[#7,#8]")


def get_rules(config_file):
    """ Library generator config.
        
    Parameters
    ----------
    config_file: 

    Returns
    -------
    rules_dict: dict
    lib_args: list
    """
    rules_dict, lib_args = {}, []
    print_l('Provided rules')

    for i,line in enumerate(config_file):
        if i == 0:
            continue
        print_l(line[:-1])

        if '==' in line:
            words = line.split('==')
            value = words[1].strip()

            if value == 'None':
                continue
            
            elif i == 1:
                in_frags = []
                if not isinstance(eval(value), tuple):
                    tmp_str = "ERROR: Wrong generation rule provided for "+line
                    print_le(tmp_str,"Aborting due to wrong generation rule.")
                for item in eval(value):
                    if 'F' in item:
                        in_frags.append(item.strip())
                    else:
                        tmp_str = "ERROR: Wrong generation rule provided for "+line
                        print_le(tmp_str,"Aborting due to wrong generation rule.")
                rules_dict['include_bb'] = in_frags
                continue
                
            elif i == 11:
                if isinstance(eval(value), tuple):
                    rules_dict['heteroatoms'] = eval(value)
                else:
                    tmp_str = "ERROR: Wrong generation rule provided for "+line
                    print_le(tmp_str,"Aborting due to wrong generation rule.")
                continue
        
            elif i == 12:
                if value == "True":
                    rules_dict['lipinski'] = True
                elif value == "False":
                    continue
                else:
                    tmp_str = "ERROR: Wrong generation rule provided for "+line
                    print_le(tmp_str,"Aborting due to wrong generation rule.")
                continue

            elif i == 13: # This rule is for fingerprint matching
                target_mols = eval(value)
                smiles_to_comp = []
                if not isinstance(target_mols, tuple) or scipy.shape(target_mols)[1] != 2:
                    tmp_str = "ERROR: Wrong generation rule provided for "+line
                    tmp_str = tmp_str+"Provide the molecule as SMILES/InChI followed by the Tanimoto index. For example, '[(c1ccccc1, 20), ]'. \n"
                    print_le(tmp_str,"Aborting due to wrong fingerprint rule.")
                for target_smiles, tanimoto_index in target_mols:
                    smiles = check_building_blocks(target_smiles,i+1,rulesFile)
                    smiles_to_comp.append([smiles,tanimoto_index])
                rules_dict['fingerprint'] = smiles_to_comp            
                continue

            elif i == 14 or i == 15:  # This rule for substructure inclusion and exclusion
                smiles_l = []
                if not isinstance(eval(value), tuple):
                    tmp_str = "ERROR: Wrong generation rule provided for "+line
                    tmp_str = tmp_str+"Provide molecule(s) as a list. \n"
                    print_le(tmp_str,"Aborting due to wrong substructure rule.")
                for item in eval(value):
                    smiles = check_building_blocks(item,i+1,rulesFile)
                    smiles_l.append(smiles)
                rules_dict[str(i)] = tuple(smiles_l)
                continue

            elif i == 16: # This rule is for inclusion of building blocks in the final library
                if value != 'True' and value != 'False':
                    tmp_str = "ERROR: Wrong generation rule provided for "+line
                    tmp_str = tmp_str+"Provide either True or False. \n"
                    print_le(tmp_str,"Aborting due to wrong generation rule.")
                
                if value == 'False':
                    rules_dict['bb_final_lib'] = False
                continue
            
            elif value != 'None':
                if not isinstance(eval(value), tuple) or len(eval(value)) != 2:
                    tmp_str = "ERROR: Wrong generation rule provided for "+line
                    tmp_str = tmp_str+"Provide the range in tuple format (min, max). \n"
                    print_le(tmp_str,"Aborting due to wrong generation rule.")
                
                else:
                    rules_dict[str(i)] = eval(value)
            
        elif '::' in line:
            words = line.split('::')
            value = words[1].strip()
            lib_args.append(value)

    return rules_dict, lib_args

def check_building_blocks(smiles, line, file_name):
    """Validate the building blocks input (smiles or inchi) and return the smiles of the molecule

    Parameters
    ----------
    smiles: string
        one molecule in 2d format (smiles/inchi)
    line: integer 
        line number in input file
    file_name: string 
        file name

    Returns
    -------
    smiles: string
        molecule in smiles format
    """
    inchi_bb, smiles_bb = True, True
    try:
        mol = pybel.readstring("inchi",smiles)
        smiles = str(mol)
        smiles = smiles.strip()
    except:
        inchi_bb = False
    
    try:
        mol = pybel.readstring("smi",smiles)
    except:
        smiles_bb = False

    if inchi_bb == False and smiles_bb == False:
        tmp_str = 'Error: The SMILES/InChI string(\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct SMILES/InChI.'.format(smiles,line,file_name)
        print_le(tmp_str,"Aborting due to wrong molecule description.")
    else:
        return smiles

    return smiles

def unique_structs(mol, smarts):
    smarts = pybel.Smarts(smarts)
    smarts.obsmarts.Match(mol.OBMol)
    num_unique_matches = len(smarts.findall(mol))
    return num_unique_matches
    
def generator(combi_type, init_mol_list):
    """
    This function creates a new generation with the initial building blocks provided and the current generation of molecules

    Parameters
    ----------
    init_mol_list: list, 
                    input molecules with duplicates removed
    combi_type: string, 
                    link/fusion
    
    Returns
    -------
    library[-1]: list,
                    Final molecular library of dict objects after all the generations
    """
    library = []
    library.append(init_mol_list)
    for gen in range(gen_len):
        print("library length: ", len(library), "gen: ", gen)
        prev_gen_mol_list = library[gen]
        lib_temp, new_gen_mol_list = [], []
        
        ## Now individual processors will generate molecules based on the list of molecules in that processor list.
        chunks_list = scipy.array_split(range(len(prev_gen_mol_list)), mpisize)
        # chunks = comm.scatter(chunks_list,root=0)
        if chunks_list[rank].any():
            for i in chunks_list[rank]:                       
                for mol2 in init_mol_list:
                    if combi_type.lower() == 'link':
                        new_gen_mol_list += create_link(prev_gen_mol_list[i],mol2,rules_dict)
                    if combi_type.lower() == 'fusion':
                        new_gen_mol_list += get_fused(prev_gen_mol_list[i],mol2,rules_dict)
        new_gen_mol_list = comm.gather(new_gen_mol_list, root=0)
        if rank == 0:
            new_gen_mol_list = list(chain.from_iterable(new_gen_mol_list))      # flatten out the list
        new_gen_mol_list = comm.bcast(new_gen_mol_list, root=0)

        # runs only for the last generation
        if gen == gen_len-1:
            list_to_scatter = list(chain.from_iterable(library))            # flatten out the list
            new_gen_mol_list = list_to_scatter + new_gen_mol_list           # add the last generation to library
            chunks_list = scipy.array_split(range(len(new_gen_mol_list)), mpisize)
            # chunks = comm.scatter(chunks_list,root=0)
            for i in chunks_list[rank]:
                mol_ob = pybel.readstring("smi", new_gen_mol_list[i]['reverse_smiles'])
                for atom in list(mol_ob.atoms):
                    ## Removing Francium and Radium atoms. It is easy to convert Francium atom to hydrogen than deleting the atom
                    if atom.OBAtom.GetAtomicNum() == 87 or atom.OBAtom.GetAtomicNum() == 88:
                        atom.OBAtom.SetAtomicNum(1)
                new_gen_mol_list[i]['reverse_smiles'] = reverse_mol(mol_ob, list(mol_ob.atoms))
        # Creating a dictionary of molecules. This is a faster way to prevent duplicates
        ## SMILES dictionary might have duplicates. Since the duplicates will only be in one Mol Wt
        ## list of the dictionary, we can divide Mol Wts into available processors.
        
        smiles_dict = defaultdict(list) 
        for l2 in new_gen_mol_list:
            mol_ob = pybel.readstring("smi", l2['reverse_smiles'])
            smiles_dict[int(mol_ob.OBMol.GetMolWt())].append(l2)                 # appending dicts of molecule in dictionary
        items = list(smiles_dict.items())
        chunks_list = scipy.array_split(range(len(items)), mpisize)
        # chunks = comm.scatter(chunks_list,root=0)
        for items_ind in chunks_list[rank]:
            mol_wt, mols = items[items_ind][0], items[items_ind][1]                     # mols --> list of molecules in that mol_wt category
            tmp_list = []
            for i in mols:
                mol_ob = pybel.readstring("smi", i['reverse_smiles'])
                if mol_ob.write("can") not in tmp_list:
                    tmp_list.append(mol_ob.write("can"))
                    lib_temp.append(i)
        lib_temp = comm.gather(lib_temp, root=0)
        if rank == 0:
            lib_temp = list(chain.from_iterable(lib_temp))
            library.append(lib_temp)
            print_l('Total molecules generated in generation number '+str(gen+1)+' is '+str(len(library[gen+1])))
        library = comm.bcast(library, root=0)
    
    
    wt2 = MPI.Wtime()
    print_l('Total time taken in generation number '+str(gen+1)+' is '+str('%.3g'%(wt2-wt1))+'\n')
    return library[-1]
    
def reverse_mol(mol, atoms):
    """
    If a molecule has no Ra atom, all available H atom sites are considered for reaction.
    If a molecule has Ra atoms, only those sites will be considered for reaction. For the reaction, the H atoms 
    are replaced with Fr atoms so that they are not considered as potential sites, and the Ra atoms are converted to H atoms. 

    Parameters
    ----------
    mol: OpenBabel object of a molecule
    atoms: list of atoms in mol 

    Returns
    -------
    smiles: string, canonical smiles
    """
    
    atom_num = [] # to append atomic number of atoms present in molecule

    ## Make a list of atomic numbers for all atoms that are in the molecule
    for atom in atoms:
        atom_num.append(atom.OBAtom.GetAtomicNum())
    
    
    if 88 in atom_num:
        ## generate a new molecule for each loop so that the old one is not changed
        #newmol=pybel.readstring("smi",smiles)
        
        for atom in atoms:
            ## check the number of hydrogens attached to the atom
            hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount() # zero for H atom
            index = atom.OBAtom.GetIdx()
            
            ## This is replacing hydrogen atoms with Francium atom
            while hcount != 0:
                size = len(list(mol.atoms))
                mol.OBMol.InsertAtom(Fratom)
                mol.OBMol.AddBond(index,size+1,1,0,-1)
                hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
                
        ## As the atoms are now changed in molecule, we will have to define atom list again
        atoms = list(mol.atoms)
        
        ## Replace all Radium atoms with hydrogen, which makes them site points.
        
        for atom in atoms:
            if atom.OBAtom.GetAtomicNum() == 88:
                atom.OBAtom.SetAtomicNum(1)
                   
    smiles = mol.write("can")[:-2]
    return smiles

def print_l(sentence):
    if rank == 0:
        print(sentence)
        logfile.write(str(sentence)+"\n")

def print_le(sentence, msg="Aborting the run"):
    if rank == 0:
        print(sentence)
        logfile.write(sentence+"\n")
        error_file.write(sentence+"\n")
        sys.exit(msg)
    else:
        sys.exit()

def print_e(sentence):
    if rank == 0:
        print(sentence)
        error_file.write(sentence+"\n")
        
# def library_generator():              # for import calls from other files
if __name__ == "__main__":

    ## initializing MPI to time, to check the MPI efficiency
    wt1 = MPI.Wtime()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpisize = comm.Get_size()
    
    if rank == 0:
        logfile = open('logfile.txt','a')
        error_file = open('error_file.txt','a')
        
    ## Defining Ra atom and Francium atom
    myRa = pybel.readstring('smi',"[Ra]")
    Raatom = myRa.OBMol.GetAtom(1)
    myFr = pybel.readstring('smi',"[Fr]")
    Fratom = myFr.OBMol.GetAtom(1)
    
    try :
        rulesFile = open("config.dat")
    except:
        tmp_str = "Generation rules file "+rulesFile+" does not exist. "
        tmp_str = tmp_str+"Please provide correct generation rule file.\n"
        print_le(tmp_str,"Aborting due to wrong file.")
    
    print_l("Reading generation rules \n")
    rules_dict, args = get_rules(rulesFile)
    print_l(rules_dict)
    BB_file, combi_type, gen_len, outfile_type, max_fpf, lib_name = args
    gen_len, max_fpf = int(gen_len), int(max_fpf)

    ## Reading the building blocks from the input file
    initial_mols = []
    print_l("Reading building blocks from the file \'"+BB_file+'\'\n')

    try :
        infile = open(BB_file)
    except:
        tmp_str = "Building blocks file "+BB_file+" does not exist. "
        tmp_str = tmp_str+"Please provide correct building blocks file.\n"
        print_le(tmp_str,"Aborting due to wrong file.")

    i_smi_list = []
    for i,line in enumerate(infile):
        smiles = line.strip()
        if smiles.isspace() or len(smiles)==0 or smiles[0]=='#':
            continue
        
        if '[X]' in smiles:
            smiles = smiles.replace('[X]','[Ra]')

        smiles = check_building_blocks(smiles,i+1,BB_file)
        # removing duplicates in the input list based on canonical smiles
        temp = molecule(smiles, 'F'+str(len(initial_mols)+1))
        is_duplicate = False
        for z in initial_mols:
            if temp['can_smiles'] not in z['can_smiles']:
                continue
            is_duplicate = True
        if not is_duplicate:
            initial_mols.append(temp)
            i_smi_list.append(temp['can_smiles'])

    print_l('Number of buidling blocks provided = '+str(len(initial_mols))+'\n')
    
    print_l('unique SMILES: \n')
    
    print_l(i_smi_list)
    print_l('=============================================================================\n')
    
    final_list = generator(combi_type, initial_mols)
    print_l('Total number of molecules generated = '+str(len(final_list))+'\n')

    output_dest = os.getcwd() + '/'
    outfile_name = output_dest + "Final_smiles_output.dat"
    outfile = open(outfile_name, "w")
    print_l('Writing SMILES to file \''+outfile_name+'\' along with corresponding molecule code.\n')

    if rank == 0:
        outfile.write('Molecule_Smiles, Combination_Code\n')
        for i in final_list:
            mol_ob = pybel.readstring("smi", i['reverse_smiles'])
            outfile.write(mol_ob.write("can") + ',' + i['code'] + '\n')
        
    ## Creating a seperate output file for each Molecule.
    ## The files are written to folder with specified no. of files per folder.

    else:
        
        print_l('Writing molecules with molecule type '+str(outfile_type)+'\n')
    
        smiles_to_scatter = []
        if rank == 0:
            if not os.path.exists(output_dest + lib_name + outfile_type):
                os.makedirs(output_dest + lib_name + outfile_type)
            smiles_to_scatter=[]
            for i in range(mpisize):
                start = int(i*(len(final_list))/mpisize)
                end = int((i+1)*(len(final_list))/mpisize)-1
                list_to_add = final_list[start:end+1]
                list_to_add = list_to_add+[len(final_list),start,end]
                smiles_to_scatter.append(list_to_add)
        else:
            smiles_to_scatter = []
    
        ## Dividing the list into processors
        smiles_list = comm.scatter(smiles_to_scatter,root=0)
        
        final_list_len = smiles_list[-3]
        start = smiles_list[-2]
        end = smiles_list[-1]
        smiles_list = smiles_list[0:-3]
        
        ratio_s = int(start/max_fpf)
        ratio_e = int(end/max_fpf)
        
        if end+1 == final_list_len:
            ratio_e = ratio_e+1
        
        for i in range(ratio_s,ratio_e):
            
            if not os.path.exists(output_dest + lib_name + outfile_type+"/"+str(i+1)+"_"+str(max_fpf)):
                os.makedirs(output_dest + lib_name + outfile_type+"/"+str(i+1)+"_"+str(max_fpf))

        folder_no = ratio_s+1
        for i, val in enumerate(range(start,end+1)):
            mol_ob = pybel.readstring("smi", smiles_list[i]['reverse_smiles'])
            mymol = pybel.readstring("smi", mol_ob.write("can"))
            mymol.make3D(forcefield='mmff94', steps=50)
            mymol.write(outfile_type, output_dest + lib_name +outfile_type+"/"+str(folder_no)+"_"+str(max_fpf)+"/"+str(val+1)+"."+outfile_type,overwrite=True)

            if (val+1)%max_fpf == 0:
                folder_no = folder_no+1
        
    print_l('File writing terminated successfully'+'\n')
    wt2 = MPI.Wtime()
    print_l('Total time_taken '+str('%.3g'%(wt2-wt1))+'\n')
    sys.stderr.close()
    sys.exit()
