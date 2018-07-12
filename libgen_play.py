#!/usr/bin/env python

_SCRIPT_NAME = "Library_Generator"
_SCRIPT_VERSION = "v0.1.20"
_REVISION_DATE = "01/19/2018"
_AUTHOR = "Mohammad Atif Faiz Afzal (m27@buffalo.edu) and Johannes Hachmann (hachmann@buffalo.edu) "
_DESCRIPTION = "This is a package for generating molecular libraries."

# Version history timeline:
# v0.1.1 (7/12/2015): Made sure the Jchem (reactor) suite works  
# v0.1.2 (7/24/2015): Use Openbabel to create the library and react
# v0.1.3 (8/01/2015): link only carbons which have atleast one hydrogen bond 
# v0.1.4 (8/03/2015): Getting rid of jchem dependency 
# v0.1.5 (8/05/2015): Parse arguments 
# v0.1.6 (8/08/2015): Get only those reactions specified by the user
# v0.1.7 (8/11/2015): Implementing parallel algorithm
# v0.1.8 (8/15/2015): Including generation rules
# v0.1.9 (8/21/2015): Further reducing the computation time (efficient parallel code)
# v0.1.10 (8/24/2015): Include fusion 
# v0.1.11 (8/28/2015): Debugging fusion. Removing duplicates in an effective manner
# v0.1.12 (9/11/2015): Added detailed comments
# v0.1.13 (11/13/2015): Saving molecules in different formats 
# v0.1.14 (12/14/2015): Specifying maximum files per folder 
# v0.1.15 (12/23/2015): Tracking the order of combination in liking 
# v0.1.15 (12/29/2015): Included generation rules
# v0.1.16 (01/03/2016): Included log file and error file
# v0.1.17 (02/25/2016): Added more options, changed to Ra and Fr and included lower limit
# v0.1.18 (01/02/2017): Further reducing the computation time (efficient parallel code), changes to the molecule data structure, included Fingerprint matching as well as substructure inclusion and exclusion
# v0.1.19 (12/26/2017): Changed the first rule in generation rules file. This rule is to include specified building blocks in the final library.
# v0.1.20 (01/19/2018): Added a new rule. Fixed a bug in rule number of specific atoms rule.

###################################################################################################
# TASKS OF THIS SCRIPT:
# -Molecular library generation
###################################################################################################

###################################################################################################

# TODO:
# -Get the reaction without reactor-Done
# -Try to run reactor from python instead subprocess -Don't need this anymore
# -Try to parallelize the code - Parallelized the core part of the code
# -Import SMILES from data file -Done
# -Capture the stderr from C program ie. openbabel. -Done (this took a very long time to debug)
# -Increase scalability of parallel - Done
# -Include Fr removals in an efficient manner- Done
# -Track the link or fusion order while building - Done
# -Include binding rules - Done
# -Inlude max no. of specific atoms - Done
# -Include log file - Done
# -Include Lipinski's rule- Done
# -Improve output style - Partially done
# -Include Fingerprint match criteria - Done 
# -Include substructure inclusion and removal criteria - Done
# -The input can be mixture of smiles and inchi - Done
# -Include range instead of Max value for generation rules - Done
# -Fix the SMILES recognition error due to split by '=' sign - Done

# -Specify order of fusion and link - This requires a large chunk of rewriting of the code, not done yet
# -Detailed comments- Partially done
# -Implement a smart duplicate removal system - Using Molwt for now, will have to look for better options
# -Use error handling properly - Mostly done
# -Include stopping criteria - Mostly done
# -Merge genetic algorithm -  Currently it is a separate code
# -Parallelize the reading of input SMILES incase the input is large
# -Include a polymer builder

###################################################################################################

import sys
sys.path.insert(0, "/user/m27/pkg/openbabel/2.3.2/lib")
import pybel
from openbabel import OBAtomAtomIter
#import subprocess
import time
import argparse
#import scipy
from collections import defaultdict
import operator


###################################################################################################

#import cStringIO
from mpi4py import MPI
import os
#import threading
from lib_modules import libgen_classes
from itertools import islice,chain


'''
This is to remove the stereo chemistry information from smiles provided.
'''
def remove_stereochemistry(smiles):
    return smiles.replace("@", "").replace("/", "-").replace("\\", "-")



def get_num_struc(mol,smarts):
    smarts = pybel.Smarts(smarts)
 
    smarts.obsmarts.Match(mol.OBMol)
    num_matches = sum(1 for indicies in smarts.obsmarts.GetMapList())
    num_unique_matches = len(smarts.findall(mol))
    return num_unique_matches


def if_del(mol,rules,code):
    
    #print mol,'after'

    if rules[1][0] != 0:
        bonds = mol.OBMol.NumBonds()
        #print rules[1][1]
        
        if bonds<rules[1][0]:
            #print 'bonds'
            return False
    
    if rules[2][0] != 0:
        no_atoms = len(mol.atoms)
        
        if no_atoms<rules[2][0]:
            #print 'no_atoms'
            return False

    if rules[3][0] != 0:
        
        mol_wt = mol.OBMol.GetMolWt()
        
        if int(mol_wt)<rules[3][0]:
            #print 'mol_wt'
            return False

    if isinstance(rules[4][0],(int)):
        if rules[4][0] != 0:
            rings = len(mol.OBMol.GetSSSR())
            if rings<int(rules[4][0]):
                #print 'rings'
                return False

    if isinstance(rules[5][0],(int)) or isinstance(rules[6][0],(int)):
        
        if rules[5][0] != 0 or rules[6][0] != 0:
            no_ar = 0
            no_non_ar = 0
            for r in mol.OBMol.GetSSSR():
                if r.IsAromatic():
                    no_ar = no_ar+1
                else:
                    no_non_ar = no_non_ar+1
            if isinstance(rules[5][0],(int)):
                if no_ar<int(rules[5][0]) :
                    return False
            if isinstance(rules[6][0],(int)):
                if no_ar<int(rules[6][0]) :
                    return False

    if isinstance(rules[7][0],(int)):
        if rules[7][0] != 0:
            no_s_bonds = get_num_struc(mol,"*-*")
            if no_s_bonds<int(rules[7][0]) :
                return False

    if isinstance(rules[8][0],(int)):
        if rules[8][0] != 0:
            no_d_bonds = get_num_struc(mol,"*=*")
            if no_d_bonds<int(rules[8][0]) :
                return False

    if isinstance(rules[9][0],(int)):
        if rules[9][0] != 0:
            no_t_bonds = get_num_struc(mol,"*#*")
            if no_t_bonds<int(rules[9][0]) :
                return False    
                
    if rules[12] != 'None':
        for mol_to_comp in rules[12]:
            mol2 = pybel.readstring('smi',mol_to_comp[0])
            tanimoto = mol.calcfp()|mol2.calcfp()            
            if tanimoto<float(mol_to_comp[1]):
                return False
           
    if rules[13] != 'None':
        for item in rules[13]:
            
            no_occ = get_num_struc(mol,item)
            if no_occ == 0 :
                return False            

    if rules[15] == 'False':
        if code.count('-') == 0 and code.count(':') == 0:
            return False

    return True
    
    
'''
This function creates a new generation with the building blocks provided
and the current generation of molecules
'''
#@profile
def create_gen_lev(smiles_list_gen,ini_list,combi_type,gen):

    library_can, library_full, library = [], [], []
    
    # Creating a dictionary of molecules. This is a faster way to prevent duplicates
    
    smiles_dict = defaultdict(list) 

    ## Adding all previous generation molecules to dictionary 
    if rank == 0:
        smiles_to_scatter = []
        for i in xrange(mpisize):
            list_to_scatter = smiles_list_gen
            start = int(i*(len(list_to_scatter))/mpisize)
            end = int((i+1)*(len(list_to_scatter))/mpisize)-1
            smiles_to_scatter.append(list_to_scatter[start:end+1])
    else:
        smiles_to_scatter = []
    
    ## Dividing the list into processors
    scattered_gen = comm.scatter(smiles_to_scatter,root=0)
                        

    ## Now individual processors will generate molecules based on the list of molecules 
    ## in that processor list.

    for smiles1 in scattered_gen:
        for smiles2 in ini_list:
            ## Check for the combination type and send the smiles to respective functions
            ## create_link function (included in libgen_classes) creates linked molecules.
            ## Whereas, get_fused (included in libgen_classes) creates fused molecule.
            ## The return type for these functions is a list with first element the SMILES
            ## of the generated molecule and second element its molecular weight.

            if combi_type == 'link':
                library_full = library_full+libgen_classes.create_link_c([smiles1[0],smiles1[2]],[smiles2[0],smiles2[2]],rules_l)
            if combi_type == 'fusion':
                library_full = library_full+libgen_classes.get_fused_c([smiles1[0],smiles1[2]],[smiles2[0],smiles2[2]],rules_l)

    ## Now we have to delete the Fr atoms for linked atoms and Ra atoms for Fused atoms
    ## This can be easily done parallely as the jobs are independent of each other
    ## Making the list ready to scatter between processors

    if gen == gen_len-1:
        if rank == 0:
            smiles_to_scatter = []
            for i in xrange(mpisize):
                list_to_scatter = list(chain.from_iterable(Global_list))
                start = int(i*(len(list_to_scatter))/mpisize)
                end = int((i+1)*(len(list_to_scatter))/mpisize)-1
                smiles_to_scatter.append(list_to_scatter[start:end+1])
        else:
            smiles_to_scatter = []
    
        ## Dividing the list into processors
        scattered_list = comm.scatter(smiles_to_scatter,root=0)
        
        library_full = scattered_list+library_full

        atm_del_l = []
        for pos,smiles in enumerate(library_full):
            
            mol_combi = pybel.readstring("smi",smiles[0])
            
            atoms = list(mol_combi.atoms)
            del_idx = []
            ## iterating over all atoms of the molecule
            for atom in atoms:
                ## Removing Francium atoms
                if atom.OBAtom.GetAtomicNum() == 87:
                    ## it is easy to convert Francium atom to hydrogen than deleting the atom
                    atom.OBAtom.SetAtomicNum(1)
                
                ## Removing Radium atoms
                if atom.OBAtom.GetAtomicNum() == 88:
                    atom.OBAtom.SetAtomicNum(1)
                    #print mol_combi
            ## mark the index of the molecules that do not the lower limit in the rules list
            if if_del(mol_combi,rules_l,smiles[2]) == False:
                
                #print 'mol_combi',mol_combi
                atm_del_l.append(pos)
            
            ## After removing Ra atoms, Fusion molecules list might contain duplicates.
            ## So convert all the SMILES to canonical so that it can be deleted later
            
            if combi_type == 'fusion':
                can_mol_combi = mol_combi.write("can")
                library_full[pos][0] = str(can_mol_combi)[:-2]
                library_full[pos][1] = int(mol_combi.OBMol.GetMolWt())
            ## If linked type, then there are no duplicates so no need for canonical
            if combi_type=='link':       
                library_full[pos][0] = str(mol_combi)[:-2]
                library_full[pos][1] = int(mol_combi.OBMol.GetMolWt())

        ## deleted the atoms that were marked to be deleted
        for pos_del in atm_del_l[::-1]:
            del library_full[pos_del]
    ## once the new generation is ready in processor, collect all the lists into the master
    ## processor

    library_gather = comm.gather(library_full,root=0)

    # ## As there can be duplicates in fused molecules we have to remove the duplicates
    # if combi_type=='fusion':
    #     lib_smiles_list=list(set(lib_smiles_list))
    
    #print gen
    #print library_gather

    ## concatinating the gathered list into SMILES dictionary based on Mol Wt
    if rank == 0:
        for l1 in library_gather:
            for l2 in l1:
                mol_wt = l2[1]
                smiles_dict[mol_wt].append([l2[0],l2[2]])
            
    ## SMILES dictionary might have duplicates. Since the duplicates will only be in one Mol Wt
    ## list of the disctionary, we can divide Mol Wts into available processors.
    
    ## Preparing the dictionary to scatter. This is similar the way described above for dividing 
    ## SMILES list
    if rank == 0:
        smiles_dict_scatter = []
        items = smiles_dict.items()
        for i in xrange(mpisize):
            start = int(i*(len(smiles_dict))/mpisize)
            end = int((i+1)*(len(smiles_dict))/mpisize)-1
            smiles_dict_scatter.append(dict(item for item in items[start:end+1]))
    else:
        smiles_dict_scatter = []
    
    ## Scatterring the list
    smiles_dict_indi = comm.scatter(smiles_dict_scatter,root=0)
    
    ## Now the duplicates in each Mol Wt can be removed by using 'set'.

    library_indi = []
    for mol_wt in smiles_dict_indi:
        a = smiles_dict_indi[mol_wt]
        # b= sorted(a, key=operator.itemgetter(0))
        # #print a
        # #print b
        # b_zip=zip(*b)
        # #print b_zip
        # b_dup=[[],[]]
        # for i,item in enumerate(b_zip[0]):
        #     if item not in b_dup[0]:
        #         b_dup[0].append(item)
        #         b_dup[1].append(b_zip[1][i])
        #         library_indi.append([item,str(mol_wt),b_zip[1][i]])

        tmp_list = []
        for i, item in enumerate(a):
            if item[0] not in tmp_list:
                tmp_list.append(item[0])
                library_indi.append([item[0],str(mol_wt),item[1]])

    library_g = comm.gather(library_indi,root=0)

    if rank == 0:
        for item in library_g:
            library = library+item

    ## printing out smiles after every generation. This mus not be included in the final version
    ## only for test purpose
    
    # if rank==0:
    #     outdata="gen_"+str(gen+1)+"_smiles.dat"
    #     outfile = open(outdata, "w")
        
    #     print_l('Writing molecules SMILES to file \''+outdata+'\' along with corresponding code.\n')
        
    #     outfile.write('Sl.No,Molecule_Smiles,Combination_Code\n')
        
    #     for i, smiles in enumerate(library):
    #         outfile.write(str(i+1)+','+smiles[0]+','+smiles[1]+'\n')
    
    ## Return the list with all molecules until this generation number
    #print library
    if rank == 0:
        return library
    else:
        return []       
'''
For a given molecule, this fucntion converts all the atoms that are not connected to 
Ra atom to fully satisfied (no hydrogens). This is done by attaching Fr atoms to these atoms.
''' 
def reverse_mol(smiles):
    
    mol = pybel.readstring("smi",smiles)            
    atoms = list(mol.atoms)
    atom_num = [] # to append atomic number of atoms present in molecule

    ## Make a list of atomic numbers for all atoms that are in the molecule
    for atom in atoms:
        atom_num.append(atom.OBAtom.GetAtomicNum())
    
    ## Changes will done only if there is a Radium atom in the molecule. If no Radium atom is
    ## present then every atom is considered as a site point for combination. 
    
    if 88 in atom_num:
        ## generate a new molecule for each loop so that the old one is not changed
        #newmol=pybel.readstring("smi",smiles)
        
        for atom in atoms:
            ## check the number of hydrogens attached to the atom
            hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
            
            #if hcount==0:
            #    continue
            
            index = atom.OBAtom.GetIdx()
            
            ### Below section (commented) can be used to the include site points containing Ra to have multiple 
            ### combinations at the site point
            
            #this_atom=mol.OBMol.GetAtom(index)
            #atom_num_n=[]
            
            ## Check all the atoms that attached to this atom
            ## and append the index to atom_num_n list
            #for atom_n in OBAtomAtomIter(this_atom):
            #    check_88=atom_n.GetAtomicNum()
            #    atom_num_n.append(check_88)
            ## if Ra present in the attached atoms, then do not do anything.
            #if 88 in atom_num_n:
            #    continue
            ##
            
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

    
def generation_test(combi_type):
    lib_smiles_list = []
    lib_can = []
    
    ## It is important to remove any duplicates that are in the provided SMILES list. 
    ## If not removed then the algorithm makes many more duplicates.
    
    for smiles in smiles_list_c:
        ## using canonical smiles, duplicates can be removed
        mol_combi = pybel.readstring("smi",smiles[0]) 
        mol_wt = str(int(mol_combi.OBMol.GetMolWt()))
        can_mol_combi = mol_combi.write("can")
        if can_mol_combi in lib_can:
            continue
        lib_can.append([can_mol_combi,mol_wt,smiles[1]])

    ## If the combination type is link, modify the molecule with Francium atoms
    if combi_type == 'link':
        for i in xrange(len(lib_can)):
            lib_smiles_list.append([reverse_mol(lib_can[i][0]),lib_can[i][1],lib_can[i][2]])
    ## If the combination type is fusion, do not do anything. 
    if combi_type == 'fusion':
        for i in xrange(len(lib_can)):
            lib_smiles_list.append([lib_can[i][0][:-2],lib_can[i][1],lib_can[i][2]]) # Note: Blank spaces are removed
    
    Global_list.append(lib_smiles_list)
    for gen in xrange(gen_len):
        
        ## lib_smiles_list includes all the molecules list up until the current generation 
        Global_list.append(create_gen_lev(Global_list[gen],Global_list[0],combi_type,gen))

        if gen<gen_len-1:
            print_l('Total molecules generated in generation number '+str(gen+1)+' is '+str(len(Global_list[gen+1])))
        else:
            a = sum([len(a) for a in Global_list[:-1]])

            if rules_l[15] == 'False':
                a=a-len(Global_list[0])

            length = len(Global_list[-1])-a
            if length<0:
                length = 0
            print_l('Total molecules generated in generation number '+str(gen+1)+' is '+str(length))
            
        ## printing out time after each generation
        wt2 = MPI.Wtime()

        print_l('Total time taken in generation number '+str(gen+1)+' is '+str('%.3g'%(wt2-wt1))+'\n')
    
    print Global_list
    # ## Now we have to delete the Fr atoms for linked atoms and Ra atoms for Fused atoms
    # ## This can be easily done parallely as the jobs are independet of each other
    # ## Making the list ready to scatter between processors

    # if rank ==0:
    #     smiles_to_scatter=[]
    #     for i in xrange(mpisize):
    #         start=int(i*(len(lib_smiles_list))/mpisize)
    #         end=int((i+1)*(len(lib_smiles_list))/mpisize)-1
    #         smiles_to_scatter.append(lib_smiles_list[start:end+1])
    # else:
    #     smiles_to_scatter=[]
    
    # lib_smiles_list=comm.scatter(smiles_to_scatter,root=0)
    
    # for pos,smiles in enumerate(lib_smiles_list):
        
    #     mol_combi= pybel.readstring("smi",smiles)
    #     atoms=list(mol_combi.atoms)
    #     del_idx=[]
    #     ## iterating over all atoms of the molecule
    #     for atom in atoms:
    #         ## Removing Francium atoms
    #         if atom.OBAtom.GetAtomicNum()==87:
    #             ## it is easy to convert Francium atom to hydrogen than deleting the atom
    #             atom.OBAtom.SetAtomicNum(1)
                
    #         ## Removing Radium atoms
    #         if atom.OBAtom.GetAtomicNum()==88:
    #             atom.OBAtom.SetAtomicNum(1)
    #     ## After removing Ra atoms, Fusion molecules list might contain duplicates.
    #     ## So convert all the SMILES to canonical tso that it can be deleted later
    #     if combi_type=='fusion':
    #         can_mol_combi = mol_combi.write("can")
    #         lib_smiles_list[pos]=str(can_mol_combi)[:-2]
    #     ## If linked type, then there are no duplicates so no need for canonical
    #     if combi_type=='link':       
    #         lib_smiles_list[pos]=str(mol_combi)[:-2]
    
    # ## Gather the list from all processors. Note that this is list of lists
    # library_gather=comm.gather(lib_smiles_list,root=0)

    # ## Re initiate the list
    # lib_smiles_list=[]
    
    # ## Convert gathered list of one smiles list
    # if rank==0:
    #     for list_s in library_gather:
    #         lib_smiles_list=lib_smiles_list+list_s
    #     print len(lib_smiles_list)
        
    #     ## As there can be duplicates in fused molecules we have to remove the duplicates
    #     if combi_type=='fusion':
    #         lib_smiles_list=list(set(lib_smiles_list))
    
    ## priniting out the time after getting the final list of molecules
    # wt2 = MPI.Wtime()
    
    #if rank==0:    
    #    print 'Total time_taken',wt2-wt1    

    ## test to see if there are any duplicates
    ## Uncomment the below part to test if any duplicates in final list
    
    #library_can=[]
    
    # for smiles in lib_smiles_list:
        
    #     mol_combi= pybel.readstring("smi",smiles)
    #     can_mol_combi = mol_combi.write("can")
    #     pos=pos+1
    #     if can_mol_combi in library_can:
    #         print can_mol_combi
    #         continue
    #     library_can.append(can_mol_combi)
    # print len(library_can)
    
    # if rank==0:
    #     pass
    #     print len(library_can),'after checking for duplicates'
    ###
    
    return Global_list[-1]


def check_if_mol(smiles,line,file_name):
    if check_if_inchi(smiles) == True:            
        this_mol = pybel.readstring("inchi",smiles)
        smiles = str(this_mol)
        smiles = smiles.strip()
    elif check_if_smiles(smiles) == False:    
        ## check if smiles
        tmp_str = 'Error: The SMILES/InChI string(\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct SMILES/InChI.'.format(smiles,line,file_name)
        print_le(tmp_str,"Aborting due to wrong molecule description.")
    return smiles

## This function is to check if the provided SMILES string are valid    
def check_if_smiles(smiles):
    data = True
    if rank == 0:
        try:
            mol = pybel.readstring("smi",smiles)
        except:
            data = False
    data = comm.bcast(data, root=0)
    
    return data

## This function is to check if the provided InChI string are valid    
def check_if_inchi(inchi):
    try:
        mol = pybel.readstring("inchi",inchi)
    except:
        return False
    return True

def print_l(sentence):
    if rank == 0:
        print sentence
        logfile.write(str(sentence)+"\n")

def print_le(sentence,msg="Aborting the run"):
    if rank == 0:
        print sentence
        logfile.write(sentence+"\n")
        error_file.write(sentence+"\n")
        sys.exit(msg)
    else:
        sys.exit()
def print_e(sentence):
    if rank == 0:
        print sentence
        error_file.write(sentence+"\n")
        
        

def get_rules(rulesFile):
    rules_l = []
    print_l('Provided rules')

    for i,lines in enumerate(rulesFile):
        if i == 0:
            continue
        print_l(lines[:-1])

        if '==' not in lines:
            tmp_str = "ERROR: Wrong generation rule provided for "+lines
            print_le(tmp_str,"Aborting due to wrong generation rule.")

        words = lines.split('==')
        value = words[1].strip()
        if i == 1:
            in_frags = []
            for item in words[1][:-1].split(','):
                in_frags.append(item.strip())
            if 'F' not in words[1]:
                in_frags = []
            rules_l.append(in_frags)
            #print words[1],rules_l,'rules_l'
            continue
            
            ## below chuck won't work because the combinations are graph
            # ex_combis=words[1][:-1].split(',')
            # for i in xrange(len(ex_combis)):
            #     item=ex_combis[i]
            #     if '-' in item:
            #         shuffle=item.split('-')
            #         ex_combis.append(shuffle[1]+'-'+shuffle[0])
            #     if ':' in item:
            #         shuffle=item.split(':')
            #         ex_combis.append(shuffle[1]+':'+shuffle[0])
            # rules_l.append(ex_combis)
            # continue

        if i == 11:
            atomsg = value
            #print atomsg,'atomsg'
            atoms_l = []
            if atomsg == 'None':
                pass
            else:
                atoms = atomsg.split(',')
                for atom in atoms:
                    atomi = atom.split('-')
                    atoms_l.append([atomi[0].strip(),int(atomi[1])])
            
            rules_l.append(atoms_l)
            continue
       
        if i == 12:
            #print words,i
            rules_l.append(value)
            continue

        if i == 13: # This rule is for FP matching
            #print words,i
            smiles_to_comp = value
            if smiles_to_comp == 'None':
                rules_l.append('None')            
                continue
            smiles_to_comp = smiles_to_comp.split(',')
            smiles_to_comp_l = []

            for item in smiles_to_comp:
                if '-' not in item:
                    tmp_str = "ERROR: Wrong generation rule provided for "+lines
                    tmp_str = tmp_str+"Privide the molecule as SMILES/InChI followed by the Tanimoto index. For example, 'c1ccccc1-20'. \n More molecules can be provided separated by a comma. \n"
                    print_le(tmp_str,"Aborting due to wrong fingerprint rule.")
                item_split = item.split('-')
                smiles = check_if_mol(item_split[0][2:-1],i+1,rule_file)
                smiles_to_comp_l.append([smiles,item_split[1][:-1]])
            rules_l.append(smiles_to_comp_l)            
            continue

        if i == 14 or i == 15:  # This rule for substructure inclusion and exclusion
            if value == 'None':
                rules_l.append('None')            
                continue
            smiles_l = []
            for item in value.split(','):
                smiles = check_if_mol(item[1:-1],i+1,rule_file)
                smiles_l.append(smiles)
            rules_l.append(smiles_l)
            print value.split(',')
            continue

        if i == 16: # This rule is for inclusion of building blocks in the final library
            if value != 'True' and value != 'False':
                tmp_str = "ERROR: Wrong generation rule provided for "+lines
                tmp_str = tmp_str+"Provide either True or False. \n"
                print_le(tmp_str,"Aborting due to wrong generation rule.")
            
            if value == 'False':
                rules_l.append('False')
            continue
        
        if value != 'None':
            if '-' not in words[1]:
                tmp_str = "ERROR: Wrong generation rule provided for "+lines
                tmp_str = tmp_str+"Privide the range of number required. For example, 10-20. \n"
                print_le(tmp_str,"Aborting due to wrong generation rule.")
            
            values = words[1].split('-')
            
            try:
                val_l = int(values[0])
                val_u = int(values[1])

            except ValueError:
                tmp_str = "ERROR: Wrong generation rule provided for "+lines
                tmp_str = tmp_str+"Please provide either a number or None as input.\n"
                tmp_str = tmp_str+"Use - to separate the range of numbers. For example, 10-20.\n"
                print_le(tmp_str,"Aborting due to wrong generation rule.")
            val = [val_l,val_u]
            rules_l.append(val)
            
        else:
            rules_l.append([0,'None'])
            
        # if i==2 or i==3 or i==4:
        #     #print words,i
        #     if value=='None':
        #         words[1]=0

        #     rules_l.append(int(words[1]))
        #     continue


        
        #rules_l.append(value)
        
    print_l("\nRules list \n"+str(rules_l) )
    return rules_l


def lower_limit(smiles_l,rules_l):
    
    return smiles_l


if __name__ == "__main__":
    
    ## initializing MPI to time, to check the MPI efficiency
    wt1 = MPI.Wtime()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpisize = comm.Get_size()


    if rank == 0:
        logfile = open('logfile.txt','a',0)
        error_file = open('error_file.txt','a',0)
        
        libgen_classes.banner(logfile, _SCRIPT_NAME, _SCRIPT_VERSION, _REVISION_DATE, _AUTHOR, _DESCRIPTION)

    ##Initializing number of generations value
    gen_len = 1
    
    ## Defining Ra atom and Francium atom
    myRa = pybel.readstring('smi',"[Ra]")
    Raatom = myRa.OBMol.GetAtom(1)
    myFr = pybel.readstring('smi',"[Fr]")
    Fratom = myFr.OBMol.GetAtom(1)
    
    ## Argument parser desription
    parser = argparse.ArgumentParser(description='This is a pacakge to generate a combinatorial library of molecules based on the building blocks provided. Please provide the building blocks in the a file in either SMILES form or InChi.')
    
    parser.add_argument('-i',"--input_file", action='store', dest='file_name', default='building_blocks.dat', help="provide the building block data file. Default is building_blocks.dat file.")
    
    parser.add_argument('-t',"--molecule_type", action='store', dest='mol_type', default='SMILES', help="Mention the molecule type in this section. Default is SMILES")
    
    parser.add_argument('-c',"--combination_type", action='store', dest='combi_type', default='Link', help="Mention the combination type in this section. Link for linking and Fusion for fusing the molecules. Default is link.")

    parser.add_argument('-g',"--generation_levels", action='store', dest='gen_len', default='1', help="Give the number of maximum combinations required in each molecule. Default is 1.")

    parser.add_argument('-output',"--output_type", action='store', dest='oft', default='smi', help="Give the type of molecule format for the generated library. Default is SMILES format.")

    parser.add_argument('-max_fpf',"--max_files_per_folder", action='store', dest='max_fpf', default=1000, help="Maximum number of files that are in a single folder. Having a large number of files in a single folder may hinder performace. Default is 10000 files per folder.")

    parser.add_argument('-rf',"--rule_file", action='store', dest='rule_file', default='generation_rules.dat', help="Specified file should contain the generation rules. Order of the rules is fixed. If the order is changed then the program runs into error . Default is generation_rules.dat.")

    parser.add_argument('-n', "--lib_name", action='store', dest='lib_name', default='new_library_', help="The name of the new library to be created. Default is new_library")

    parser.add_argument('-htps',"--ChemHTPS", action='store_true', dest='chemhtps', default=False, help="Tells the library generator whether it is being run from inside the ChemHTPS program or not for purposes of output file destination. Default is False.")

    ## defining arguments
    args = parser.parse_args()
    mol_type = args.mol_type.lower()
    combi_type = args.combi_type.lower()
    max_fpf = int(args.max_fpf)
    gen_len = int(args.gen_len)
    rule_file = args.rule_file
    BB_file = args.file_name
    oft = args.oft.lower()

    ## setting output directory
    if args.chemhtps:
        output_dest = os.getcwd() + '/screeninglib/'
    else:
        output_dest = os.getcwd() + '/'
    
    print_l("Total generation levels provided = "+str(gen_len)+'\n')
    print_l("Combination type is "+str(combi_type)+'\n')
    print_l("Building blocks type is "+str(mol_type)+'\n')
    print_l("Output molecule type is "+str(args.oft).upper()+'\n')

    
    smiles_list = []
    
    
    ## Reading the building rules from generation_rules.dat
    print_l("Reading generation rules from the file \'"+rule_file+'\'\n')

    try :
        rulesFile = open(rule_file)
    except:
        tmp_str = "Generation rules file "+rule_file+" does not exist. "
        tmp_str = tmp_str+"Please provide correct generation rule file.\n"
        print_le(tmp_str,"Aborting due to wrong file.")
        
    # rules list order
    '''
    1. Exclude combinations
    2. No.of bonds
    3. Max no.of atomS
    4. Max mol.weight
    5. Max no.of rings
    6. Max no.of aromatic rings
    7. Max no.of non aromatic rings
    8. Max no.of single bonds
    9. Max no.of double bonds
    10. Max no.of triple bonds
    11. Max no.of specific atoms
    12. Lipinski's rule
    13. Combination order
    14. symmetry

    '''

    rules_l = get_rules(rulesFile)


    ## Reading the building blocks from the input file
    print_l("Reading building blocks from the file \'"+BB_file+'\'\n')

    try :
        infile = open(BB_file)
    except:
        tmp_str = "Building blocks file "+BB_file+" does not exist. "
        tmp_str = tmp_str+"Please provide correct building blocks file.\n"
        print_le(tmp_str,"Aborting due to wrong file.")

    
    ## Read molecules provided in the input 
    for i,line in enumerate(infile):
        smiles = line.strip()
        if smiles.isspace() or len(smiles)==0 or smiles[0]=='#':
            continue
        ## if the input is InChI, then convert all into SMILES
        # if mol_type == 'inchi':
        #     continue
        if 'X' in smiles:
            smiles = smiles.replace('[x]','[Ra]')
            smiles = smiles.replace('[X]','[Ra]')

        smiles = check_if_mol(smiles,i+1,args.file_name)
            #kill_me=comm.gather(tmp_str,root=0)

        #if check_if_smiles(smiles):
        smiles_list.append(smiles)

    print_l('Number of buidling blocks provided = '+str(len(smiles_list))+'\n')
    
    print_l('SMILES of the building blocks\n')
    
    print_l(smiles_list)
    print_l('=============================================================================\n')
    
    ##assigning the code for each building block
    
    smiles_list_c = []
    for i, item in enumerate(smiles_list):
        
        smiles_list_c.append([smiles_list[i],'F'+str(i+1)])

    Global_list = []
    ## generation_test funtion generates combinatorial molecules
    final_list = generation_test(combi_type)

    ## Removing the molecules according to the lower limit in generation rules
    
    final_list = lower_limit(final_list,rules_l)
    

    final_list_len=len(final_list)

    print_l('Total number of molecules generated = '+str(final_list_len)+'\n')

    #if rank==0:
    F_smi = output_dest + "Final_smiles_output.dat"
    outfile = open(F_smi, "w")

    print_l('Writing molecules SMILES to file \''+F_smi+'\' along with corresponding code.\n')


    if rank == 0:
        outfile.write('Sl.No,Molecule_Smiles,Combination_Code\n')
    
    for i, smiles in enumerate(final_list):
        outfile.write(str(i+1)+','+smiles[0]+','+smiles[2]+'\n')
              

    if oft == 'smi':
        if rank == 0:
            if not os.path.exists(output_dest + args.lib_name + oft):
                os.makedirs(output_dest + args.lib_name + oft)
        outdata = output_dest + args.lib_name + oft + "/Final_smiles_output.smi"
        outfile = open(outdata, "w")
        
        print_l('Writing molecules SMILES to file \''+outdata+'\'\n')
        
        if rank == 0:
            for i, smiles in enumerate(final_list):
                outfile.write(smiles[0]+'\n')
	    os.system('cp '+BB_file+' '+rule_file+' '+F_smi+' '+output_dest+args.lib_name+oft+'/.')
        
    ## Creating a seperate output file for each Molecule.
    ## The files are written to folder with specified no. of files per folder.

    if oft != 'smi':
        
        print_l('Writing molecules with molecule type '+str(oft)+'\n')
    
        smiles_to_scatter = []
        if rank == 0:
            if not os.path.exists(output_dest + args.lib_name + oft):
                os.makedirs(output_dest + args.lib_name + oft)
            smiles_to_scatter=[]
            for i in xrange(mpisize):
                start = int(i*(final_list_len)/mpisize)
                end = int((i+1)*(final_list_len)/mpisize)-1
                list_to_add = final_list[start:end+1]
                list_to_add = list_to_add+[final_list_len,start,end]
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
        
        for i in xrange(ratio_s,ratio_e):
            
            if not os.path.exists(output_dest + args.lib_name + oft+"/"+str(i+1)+"_"+str(max_fpf)):
                os.makedirs(output_dest + args.lib_name + oft+"/"+str(i+1)+"_"+str(max_fpf))

        folder_no = ratio_s+1
        for i, val in enumerate(xrange(start,end+1)):
            mymol = pybel.readstring("smi",smiles_list[i][0])
            mymol.make3D(forcefield='mmff94', steps=50)
            mymol.write(oft, output_dest + args.lib_name +oft+"/"+str(folder_no)+"_"+str(max_fpf)+"/"+str(val+1)+"."+oft,overwrite=True)

            if (val+1)%max_fpf == 0:
                folder_no = folder_no+1
	if rank == 0:
	    os.system('cp '+BB_file+' '+rule_file+' '+F_smi+' '+output_dest+args.lib_name+oft+'/.')
            
            

    print_l('File writing terminated successfully'+'\n')
                
    wt2 = MPI.Wtime()
    
    print_l('Total time_taken '+str('%.3g'%(wt2-wt1))+'\n')
    
    
    sys.stderr.close()
    
