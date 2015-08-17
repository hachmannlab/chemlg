#!/usr/bin/env python

_LIB_NAME = "Library_Generator"
_LIB_VERSION = "v0.1.3"
_REVISION_DATE = "8/01/2015"
_AUTHOR = "M. Atif Afzal (m27@buffalo.edu)"
_DESCRIPTION = "This is the a library generating molecular libraries."

# Version history timeline:
# v0.1.1 (7/12/2015): Made sure the Jchem (reactor) suite works  
# v0.1.2 (7/24/2015): Use Openbabel to create the library and react
# v0.1.3 (8/01/2015): link only carbons which have atleast one hydrogen bond 

###################################################################################################
# TASKS OF THIS LIBRARY:
# -Molecular library generation
###################################################################################################

###################################################################################################
# TODO:
# -Make sure to remove duplicates for building block replicates- Done
# -Enter the library elements to a list
# -Try to run reactor from python instead subprocess (create jchem python library)
# -Try to parallelize the code
# -Include stopping creteria

###################################################################################################

#import sys
#sys.path.insert(0, "/user/m27/pkg/openbabel/2.3.2/lib")
import pybel
import openbabel
import subprocess
import time

###################################################################################################




def remove_stereochemistry(smiles):
    return smiles.replace("@", "").replace("/", "-").replace("\\", "-")


smiles1="c1(C)ccccc1" 
smiles2="c1ccccc1"

myMg=pybel.readstring('smi',"[Mg]")
Mgatom=myMg.OBMol.GetAtom(1)
print Mgatom

smiles_list=[smiles1,smiles2]

reaction= "[Mg:3][c,C:1].[Mg:4][c,C:2]>>[c,C:1][c,C:2].[Mg:3][Mg:4]"

bashCommand = "react -r "+reaction+" "+smiles1+" "+smiles2

#subprocess.check_call(bashCommand.split())

def attach_Mg(atom):
    # A general purpose SMILES writer would check to see if the atom
    # can be written without the []s and if so, only show the symbol.
    # This version is easier, though incomplete (no isotope support)
    # The goal is to get the point across and not make the best SMILES
    smiles=""
    if atom.OBAtom.IsAromatic():
        smiles += atom.OBAtom.GetType().lower()
    else:
        smiles += atom.OBAtom.GetType()

    smiles += "([Mg])"
    
    return smiles


def get_smiles_permutation(smiles):
    smiles_extended=[] 
    can_smiles_list=[]

    mol= pybel.readstring("smi",smiles)
    print mol
    atoms=list(mol.atoms)

    
    size=len(list(mol))
    

    for atom in atoms:
        hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
        if hcount==0:
            continue
        newmol=pybel.readstring("smi",smiles)
        newmol.OBMol.InsertAtom(Mgatom)
        index=atom.OBAtom.GetIdx()
        newmol.OBMol.AddBond(index,size+1,1,0,-1)
        ## remove duplicates
        can_smiles = newmol.write("can")
        print 'can_smiles',can_smiles
        if can_smiles in can_smiles_list:
            continue
        can_smiles_list.append(can_smiles)
        #print can_smiles
        
        print 'newmol',newmol
        smiles_extended.append(str(newmol)[:-2])
    #print smiles_extended
    return smiles_extended

    # for atom in atoms:
    #     hcount = atom.GetNumExplicitHs() + atom.GetNumImplicitHs()
    #     if hcount==0:
    #         continue
    #     else:
    #         new_smiles=
    #         smiles_extended.append
        

def generation_test():
    smiles_list_extended=[]
    for smiles in smiles_list:
        smiles = remove_stereochemistry(smiles)
        
        mol= pybel.readstring("smi",smiles)

        smiles_list_extended = smiles_list_extended+get_smiles_permutation(smiles)

    print smiles_list_extended
    full_size=len(smiles_list_extended)
    library_can=[]
    library=[]
    for smile1 in xrange(full_size):
        for smile2 in xrange(smile1,full_size):
            bashCommand = "react -r \""+reaction+"\" \""+smiles_list_extended[smile1]+"\" \""+smiles_list_extended[smile2]+"\" -x 1"
            print bashCommand

            tstart=time.time()
            print tstart

            p= subprocess.Popen(bashCommand, shell=True, stdout=subprocess.PIPE)
            p = p.stdout.read()
            
            #tend=time.clock()
            elapsed =(time.time()-tstart)
            print elapsed
            
            #req_smiles=p.split('\n')
            #req_smiles=req_smiles[0]
            req_smiles=p[:-1]

            newmol=pybel.readstring("smi",req_smiles)
            print newmol
            can_req_smiles = newmol.write("can")
            
            print can_req_smiles
            if can_req_smiles in library_can:
                continue
            print 'smile1',smile1,'smile2',smile2
            library_can.append(can_req_smiles)
            library.append(req_smiles)
            
            #subprocess.check_call(bashCommand.split(), stdout=subprocess.PIPE)
            #p = p.stdout.read()
            #print p
    print library_can
    print library
        
if __name__ == "__main__":
    generation_test()
