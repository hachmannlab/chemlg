#!/usr/bin/env python

_SCRIPT_NAME = "Library_Generator"
_SCRIPT_VERSION = "v0.1.5"
_REVISION_DATE = "8/05/2015"
_AUTHOR = "M. Atif Afzal (m27@buffalo.edu)"
_DESCRIPTION = "This is the a library generating molecular libraries."

# Version history timeline:
# v0.1.1 (7/12/2015): Made sure the Jchem (reactor) suite works  
# v0.1.2 (7/24/2015): Use Openbabel to create the library and react
# v0.1.3 (8/01/2015): link only carbons which have atleast one hydrogen bond 
# v0.1.4 (8/03/2015): Getting rid of jchem dependency 
# v0.1.5 (8/05/2015): Parse arguments and implement a smart duplicate removal system

###################################################################################################
# TASKS OF THIS SCRIPT:
# -Molecular library generation
###################################################################################################

###################################################################################################
# TODO:
# -Get the reaction without reactor-Done
# -Try to run reactor from python instead subprocess (create jchem python library) -Don't need this anymore
# -Try to parallelize the code
# -Include stopping creteria
# -Import SMILES from data file -Done
# -Use error handling properly -Partially done

###################################################################################################

import sys
#sys.path.insert(0, "/user/m27/pkg/openbabel/2.3.2/lib")
import pybel
import openbabel
import subprocess
import time
import argparse

###################################################################################################

import cStringIO
import sys



def remove_stereochemistry(smiles):
    return smiles.replace("@", "").replace("/", "-").replace("\\", "-")


def get_index_list(smiles):
    mol=pybel.readstring("smi",smiles)            
    atoms=list(mol.atoms)
    size=len(atoms)
    can_smiles_list=[]
    atoms_index=[]
    atom_num=[]
    for atom in atoms:
        atom_num.append(atom.OBAtom.GetAtomicNum())
    if 12 in atom_num:
        print 'yes'
        newmol=pybel.readstring("smi",smiles)
        for atom in atoms:
            
            if atom.OBAtom.GetAtomicNum()==12:
                index=atom.OBAtom.GetIdx()
                mol.OBMol.DeleteAtom(mol.OBMol.GetAtom(index))
                atoms_index.append(index-1)
                size=len(list(mol.atoms))
        smiles=str(mol)
    else:
        for atom in atoms:
            hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
            print atom.OBAtom.GetAtomicNum()
            if hcount==0:
                continue
            ## remove duplicates
            newmol=pybel.readstring("smi",smiles)
            newmol.OBMol.InsertAtom(Mgatom)
            index=atom.OBAtom.GetIdx()
            newmol.OBMol.AddBond(index,size+1,1,0,-1)
            
            can_smiles = newmol.write("can")
            if can_smiles in can_smiles_list:
                continue
            can_smiles_list.append(can_smiles)
            index=atom.OBAtom.GetIdx()
            atoms_index.append(index)
    return atoms_index, size, smiles
           
def create_gen_lev(smiles_list_gen,ini_list):

    #full_size=len(smiles_list_gen)
    #print full_size
    library_can=[]
    library=[]
    for smiles in smiles_list_gen:
        mol_combi= pybel.readstring("smi",smiles)
        can_mol_combi = mol_combi.write("can")
        if can_mol_combi in library_can:
            continue
        library_can.append(can_mol_combi)
        
    # for i in xrange(full_size):
    #     for j in xrange(i,full_size):
                    
    for smiles1 in smiles_list_gen:
        for smiles2 in ini_list:
            
            #smiles1=smiles_list_gen[i]
            #smiles2=smiles_list_gen[j]
            
            mol1_index_list, size1=get_index_list(smiles1)
            mol2_index_list, size2=get_index_list(smiles2)
            #print mol1_index_list
            
            smiles_combi=smiles1+'.'+smiles2
            
            for index1 in mol1_index_list:
                for index2 in mol2_index_list:
                    
                    mol_combi= pybel.readstring("smi",smiles_combi)
                    
                    mol_combi.OBMol.AddBond(index1,index2+size1,1,0,-1)
                    can_mol_combi = mol_combi.write("can")
                    if can_mol_combi in library_can:
                        continue
                    #print 'smile1',index1,'smile2',index2
                    library_can.append(can_mol_combi)
                    library.append(str(mol_combi)[:-2])
                    #print 'newmol',mol_combi
    #print library_can
    #print library
    full_library_list=smiles_list_gen+library
    return library

def generation_test():
    lib_smiles_list=smiles_list
    

    for gen in xrange(gen_len):
        lib_smiles_list=lib_smiles_list+create_gen_lev(lib_smiles_list,smiles_list)
        #print len(lib_smiles_list)
    
    
    print len(lib_smiles_list)
    ## test to see if there are any duplicates
    library_can=[]
    for smiles in lib_smiles_list:
        mol_combi= pybel.readstring("smi",smiles)
        can_mol_combi = mol_combi.write("can")
        if can_mol_combi in library_can:
            continue
        library_can.append(can_mol_combi)
    
    print len(library_can)
    return lib_smiles_list
    
    
    
def check_if_smiles(smiles):
    try:
        mol= pybel.readstring("smi",smiles)
    except:
        return False
    return True

def check_if_inchi(inchi):
    try:
        mol= pybel.readstring("inchi",inchi)
    except:
        return False
    return True

if __name__ == "__main__":
    
    smiles1="c1(C)ccccc1" 
    smiles2="c1occc1"
    gen_len=1
    myMg=pybel.readstring('smi',"[Mg]")
    Mgatom=myMg.OBMol.GetAtom(1)

    #smiles_list=[smiles1,smiles2]

    parser = argparse.ArgumentParser(description='This is a pacakge to generate a combinatorial library of molecules based on the building blocks provided. Please provide the building blocks in the a file in either SMILES form or InChi.')
    
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')

    parser.add_argument('-i',"--input_file", action='store', dest='file_name', default='building_blocks.dat', help="provide the building block data file. Default is building_blocks.dat file.")
    
    parser.add_argument('-t',"--molecule_type", action='store', dest='mol_type', default='SMILES', help="Mention the molecule type in this section.")
    
    args = parser.parse_args()
    mol_type=args.mol_type.lower()
    #print mol_type
    smiles_list=[]
    
    #infile=open(args.file_name)
    
    #num_lines = sum(1 for line in infile)
    
    infile=open(args.file_name)

    # capture_output = cStringIO.StringIO()
    # sys.stdout = capture_output
    
    for i,line in enumerate(infile):
        smiles= line.strip()
        #smiles=smiles[0].strip()
        if smiles.isspace() or len(smiles)==0 or smiles[0]=='#':
            continue
        if mol_type == 'inchi':
            if check_if_inchi(smiles)==False:
                message='Error: The InChI string(\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct InChI.'.format(smiles,i+1,args.file_name)
            this_mol=pybel.readstring("inchi",smiles)
            smiles=str(this_mol)
            smiles=smiles.strip()

        ## check if smiles
        if check_if_smiles(smiles)==False:
            message='Error: The SMILES string(\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct SMILES.'.format(smiles,i+1,args.file_name)
            sys.exit(message)
            
        if check_if_smiles(smiles):
            smiles_list.append(smiles)
    
    print 'Number of smiles provided =',len(smiles_list)
    
    print smiles_list

    final_list=generation_test()
    
    print final_list
    #sys.stdout = sys.__stdout_
