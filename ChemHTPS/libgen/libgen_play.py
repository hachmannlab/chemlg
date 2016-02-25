#!/usr/bin/env python


_SCRIPT_NAME = "Library_Generator"
_SCRIPT_VERSION = "v0.1.16"
_REVISION_DATE = "01/03/2015"
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

###################################################################################################
# TASKS OF THIS SCRIPT:
# -Molecular library generation
###################################################################################################

###################################################################################################
# TODO:
# -Get the reaction without reactor-Done
# -Try to run reactor from python instead subprocess -Don't need this anymore
# -Try to parallelize the code - Parallelized the core part of the code
# -Include stopping creteria
# -Import SMILES from data file -Done
# -Implement a smart duplicate removal system - Using Molwt for now, will have to look for better options
# -Capture the stderr from C program ie. openbabel. -Done (this took a very long time to debug)
# -Increase scalability of parallel - Done
# -Detailed comments- Partially done
# -Include He removals in an efficient manner- Done
# -Track the link or fusion order while building - Done for linking and fusion
# -Include binding rules - Done
# -Inlude max no. of specific atoms- Done
# -Include log file - Done
# -Include Lipinski's rule- Done
# -Use error handling properly -Partially done
# -Improve output style - Partially done
# -Specify order of fusion and link - This requires a large chunk of rewriting of the code, not done yet
# -Include range instead of Max value for generation rules - Upper range taken care, lower range not done

###################################################################################################

import sys
#sys.path.insert(0, "/user/m27/pkg/openbabel/2.3.2/lib")
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
from itertools import islice


'''
This is to remove the stereo chemistry information from smiles provided.
'''
def remove_stereochemistry(smiles):
    return smiles.replace("@", "").replace("/", "-").replace("\\", "-")

'''
This function creates a new generation with the building blocks provided
and the current of molecules
'''
def create_gen_lev(smiles_list_gen,ini_list,combi_type,gen):

    library_can=[]
    library_full=[]
    library=[]
    
    # Creating a dictionary of molecules. This is a faster way to prevent duplicates
    
    smiles_dict = defaultdict(list) 

    ## Adding all previous generation molecules to dictionary 
    ini_lib_full=[]
    for smiles in smiles_list_gen+ini_list:
        
        mol_combi= pybel.readstring("smi",smiles[0])
        
        ## Here the dictionary is based on molecular weight of compund. 
        ## All same molecular weight compunds will be appended to the same list. 
        ## This data structures helps in indentifying duplicates fast.
        ## This way of structure also makes it possible to delete duplicates in 
        ## a parallel fashion. Otherwise, it is extremely challenging to parallelize
        ## duplicates removal algorithm.

        mol_wt=str(int(mol_combi.OBMol.GetMolWt()))
        
        #smiles_dict[mol_wt].append(str(mol_combi))
        
        can_mol_combi = mol_combi.write("can")
        ini_lib_full.append([str(can_mol_combi),mol_wt,smiles[1]])
        
                        
    ## It is better to initialize a list to be divided among processors in the master 
    ## processor. Also, put the list as empty in other processors.
    
    if rank!=0:
        smiles_gen_list=None
    
    ## Scatter, function of comm, is used to divide the list among the processors. But, 
    ## scatter only works when the size of the list is equal to the number of processors
    ## available. Therefore, we will have to restructure the list into a list of lists with 
    ## size equal to number of processors.
    
    if rank ==0:
        smiles_to_scatter=[]
        for i in xrange(mpisize):
            start=int(i*(len(smiles_list_gen))/mpisize)
            end=int((i+1)*(len(smiles_list_gen))/mpisize)-1
            smiles_to_scatter.append(smiles_list_gen[start:end+1])
    else:
        smiles_to_scatter=[]
    
    ## Dividing the list into processors
    smiles_list_gen=comm.scatter(smiles_to_scatter,root=0)

    ## Now individual processors will generate molecules based on the list of molecules 
    ## in that processor list.

    for smiles1 in smiles_list_gen:
        for smiles2 in ini_list:
            ## Check for the combination type and send the smiles to respective functions
            ## create_link function (included in libgen_classes) creates linked molecules.
            ## Whereas, get_fused (included in libgen_classes) creates fused molecule.
            ## The return type for these functions is a list with first element the SMILES
            ## of the generated molecule and second element its molecular weight.
            
            if combi_type=='link':
                library_full=library_full+libgen_classes.create_link_c(smiles1,smiles2,rules_l)
            if combi_type=='fusion':
                library_full=library_full+libgen_classes.get_fused_c(smiles1,smiles2,rules_l)

    ## Now we have to delete the He atoms for linked atoms and Mg atoms for Fused atoms
    ## This can be easily done parallely as the jobs are independet of each other
    ## Making the list ready to scatter between processors

    if gen==gen_len-1:

        library_full=ini_lib_full+library_full
        for pos,smiles in enumerate(library_full):
            mol_combi= pybel.readstring("smi",smiles[0])
            atoms=list(mol_combi.atoms)
            del_idx=[]
            ## iterating over all atoms of the molecule
            for atom in atoms:
                ## Removing Helium atoms
                if atom.OBAtom.GetAtomicNum()==2:
                    ## it is easy to convert Helium atom to hydrogen than deleting the atom
                    atom.OBAtom.SetAtomicNum(1)
                
                ## Removing Magnesium atoms
                if atom.OBAtom.GetAtomicNum()==12:
                    atom.OBAtom.SetAtomicNum(1)
            ## After removing Mg atoms, Fusion molecules list might contain duplicates.
            ## So convert all the SMILES to canonical so that it can be deleted later
            if combi_type=='fusion':
                can_mol_combi = mol_combi.write("can")
                library_full[pos][0]=str(can_mol_combi)
                library_full[pos][1]=int(mol_combi.OBMol.GetMolWt())
            ## If linked type, then there are no duplicates so no need for canonical
            if combi_type=='link':       
                library_full[pos][0]=str(mol_combi)
                library_full[pos][1]=int(mol_combi.OBMol.GetMolWt())
    ## once the new generation is ready in processor, collect all the lists into the master
    ## processor
    library_gather=comm.gather(library_full,root=0)

    # ## As there can be duplicates in fused molecules we have to remove the duplicates
    # if combi_type=='fusion':
    #     lib_smiles_list=list(set(lib_smiles_list))
    
    #print gen
    #print library_gather

    ## concatinating the gathered list into SMILES dictionary based on Mol Wt
    if rank ==0:
        for l1 in library_gather:
            for l2 in l1:
                mol_wt=l2[1]
                smiles_dict[mol_wt].append([l2[0],l2[2]])
        if gen!=gen_len-1:
         
            for l2 in ini_lib_full:
                mol_wt=l2[1]
                smiles_dict[mol_wt].append([l2[0],l2[2]])
            
    ## SMILES dictionary might have duplicates. Since the duplicates will only be in one Mol Wt
    ## list of the disctionary, we can divide Mol Wts into available processors.
    
    ## Preparing the dictionary to scatter. This is similar the way described above for dividing 
    ## SMILES list
    if rank ==0:
        smiles_dict_scatter=[]
        items=smiles_dict.items()
        for i in xrange(mpisize):
            start=int(i*(len(smiles_dict))/mpisize)
            end=int((i+1)*(len(smiles_dict))/mpisize)-1
            smiles_dict_scatter.append(dict(item for item in items[start:end+1]))
    else:
        smiles_dict_scatter=[]
    
    ## Scatterring the list
    smiles_dict_indi=comm.scatter(smiles_dict_scatter,root=0)
    
    ## Now the duplicates in each Mol Wt can be removed by using 'set'.

    for mol_wt in smiles_dict_indi:
        # smiles_dict_indi[mol_wt]=list(set(smiles_dict_indi[mol_wt]))
        a=smiles_dict_indi[mol_wt]
        b= sorted(a, key=operator.itemgetter(0))
        
        b_zip=zip(*b)
        
        b_dup=[[],[]]
        for i,item in enumerate(b_zip[0]):
            if item not in b_dup[0]:
                b_dup[0].append(item)
                b_dup[1].append(b_zip[1][i])
                
        smiles_dict_indi[mol_wt]=zip(b_dup[0],b_dup[1])
    
    ## This is to see if computation time can be further decreased
    # for mol_wt in smiles_dict_indi:
        
    #     for smiles_indi in smiles_dict_indi[mol_wt]:
    #         library.append(smiles_indi[:-2])
           
    
    ## Gather the list of SMILES from each processor to the master processor
    smiles_dict_g=comm.gather(smiles_dict_indi,root=0)
   
    ## This is to see if computation time can be further decreased
    #library_g=comm.gather(library,root=0)
   
    # if rank ==0:
    #     for list_g_indi in library_g:
    #         library=library+list_g_indi
    
    ## Now that the duplicates are removed we can make a list. This is because it 
    ## is easy to deal with list and easy to parallelize.
    if rank==0:
        for dict_indi in smiles_dict_g:
            for mol_wt in dict_indi:
                for smiles_indi in dict_indi[mol_wt]:
                    library.append([smiles_indi[0][:-2],smiles_indi[1]])

    ## printing out smiles after every generation. This mus not be included in the final version
    ## only for test purpose
    
    if rank==0:
        outdata="gen_"+str(gen+1)+"_smiles.dat"
        outfile = open(outdata, "w")
        
        print_l('Writing molecules SMILES to file \''+outdata+'\' along with corresponding code.\n')
        
        outfile.write('Sl.No,Molecule_Smiles,Combination_Code\n')
        
        for i, smiles in enumerate(library):
            outfile.write(str(i+1)+','+smiles[0]+','+smiles[1]+'\n')
    
    ## Return the list with all molecules until this generation number
    #print library
    if rank==0:
        return library
    else:
        return []       
'''
For a given molecule, this fucntion converts all the atoms that are not connected to 
Mg atom to fully satisfied (no hydrogens). This is done by attaching He atoms to these atoms.
''' 
def reverse_mol(smiles):
    
    mol=pybel.readstring("smi",smiles)            
    atoms=list(mol.atoms)
    atom_num=[] # to append atomic number of atoms present in molecule

    ## Make a list of atomic numbers for all atoms that are in the molecule
    for atom in atoms:
        atom_num.append(atom.OBAtom.GetAtomicNum())
    
    ## Changes will done only if there is a Magnesium atom in the molecule. If no Magnesium atom is
    ## present then every atom is considered as a site point for combination. 
    
    if 12 in atom_num:
        ## generate a new molecule for each loop so that the old one is not changed
        #newmol=pybel.readstring("smi",smiles)
        
        for atom in atoms:
            ## check the number of hydrogens attached to the atom
            hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
            
            #if hcount==0:
            #    continue
            
            index=atom.OBAtom.GetIdx()
            
            ### Below section (commented) can be used to the include site points containing Mg to have multiple 
            ### combinations at the site point
            
            #this_atom=mol.OBMol.GetAtom(index)
            #atom_num_n=[]
            
            ## Check all the atoms that attached to this atom
            ## and append the index to atom_num_n list
            #for atom_n in OBAtomAtomIter(this_atom):
            #    check_12=atom_n.GetAtomicNum()
            #    atom_num_n.append(check_12)
            ## if Mg present in the attached atoms, then do not do anything.
            #if 12 in atom_num_n:
            #    continue
            ##
            
            ## This is replacing hydrogen atoms with Helium atom
            while hcount!=0:
                size=len(list(mol.atoms))
                mol.OBMol.InsertAtom(Heatom)
                mol.OBMol.AddBond(index,size+1,1,0,-1)
                hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
                
        ## As the atoms are now changed in molecule, we will have to define atom list again
        atoms=list(mol.atoms)
        
        ## Replace all Magnesium atoms with hydrogen, which makes them site points.
        
        for atom in atoms:
            if atom.OBAtom.GetAtomicNum()==12:
                atom.OBAtom.SetAtomicNum(1)
                   
    smiles=str(mol)[:-2]
    return smiles

    
def generation_test(combi_type):
    lib_smiles_list=[]
    lib_can=[]
    
    ## It is important to remove any duplicates that are in the provided SMILES list. 
    ## If not removed then the algorithm makes many more duplicates.
    
    for smiles in smiles_list_c:
        ## using canonical smiles, duplicates can be removed
        mol_combi= pybel.readstring("smi",smiles[0]) 
        #mol_wt=str(int(mol_combi.OBMol.GetMolWt()))
        can_mol_combi = mol_combi.write("can")
        if can_mol_combi in lib_can:
            continue
        lib_can.append([can_mol_combi,smiles[1]])

    ## If the combination type is link, modify the molecule with Helium atoms
    if combi_type=='link':
        for i in xrange(len(lib_can)):
            lib_smiles_list.append([reverse_mol(lib_can[i][0]),lib_can[i][1]])
    ## If the combination type is fusion, do not do anything. 
    if combi_type=='fusion':
        for i in xrange(len(lib_can)):
            lib_smiles_list.append([lib_can[i][0][:-2],lib_can[i][1]]) # Note: Blank spaces are removed
    
    ini_list=lib_smiles_list
    
    
    
    ## Generating molecules at each generation until the required length

    for gen in xrange(gen_len):
        
        ## lib_smiles_list includes all the molecules list up until the current generation 
        lib_smiles_list=create_gen_lev(lib_smiles_list,ini_list,combi_type,gen)
        
        print_l('Total molecules generated in generation number '+str(gen+1)+' is '+str(len(lib_smiles_list)))
            
        ## printing out time after each generation
        wt2 = MPI.Wtime()

        print_l('Total time taken in generation number '+str(gen+1)+' is '+str('%.3g'%(wt2-wt1))+'\n')
    
    # ## Now we have to delete the He atoms for linked atoms and Mg atoms for Fused atoms
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
    #         ## Removing Helium atoms
    #         if atom.OBAtom.GetAtomicNum()==2:
    #             ## it is easy to convert Helium atom to hydrogen than deleting the atom
    #             atom.OBAtom.SetAtomicNum(1)
                
    #         ## Removing Magnesium atoms
    #         if atom.OBAtom.GetAtomicNum()==12:
    #             atom.OBAtom.SetAtomicNum(1)
    #     ## After removing Mg atoms, Fusion molecules list might contain duplicates.
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
    
    return lib_smiles_list


## This function is to check if the provided SMILES string are valid    
def check_if_smiles(smiles):
    data=True
    if rank==0:
        try:
            mol= pybel.readstring("smi",smiles)
        except:
            data= False
    data = comm.bcast(data, root=0)
    
    return data

## This function is to check if the provided InChI string are valid    
def check_if_inchi(inchi):
    try:
        mol= pybel.readstring("inchi",inchi)
    except:
        return False
    return True

def print_l(sentence):
    if rank ==0:
        print sentence
        logfile.write(str(sentence)+"\n")

def print_le(sentence,msg="Aborting the run"):
    if rank ==0:
        print sentence
        logfile.write(sentence+"\n")
        error_file.write(sentence+"\n")
        sys.exit(msg)
    else:
        sys.exit()
def print_e(sentence):
    if rank ==0:
        print sentence
        error_file.write(sentence+"\n")
        
        

def get_rules(rulesFile):
    rules_l=[]
    print_l('Provided rules')

    for i,lines in enumerate(rulesFile):
        if i==0:
            continue
        print_l(lines[:-1])

        if '=' not in lines:
            tmp_str = "ERROR: Wrong generation rule provided for "+lines
            print_le(tmp_str,"Aborting due to wrong generation rule.")

        words=lines.split('=')

        if i==1:
            ex_combis=words[1][:-1].split(',')
            for i in xrange(len(ex_combis)):
                item=ex_combis[i]
                if '-' in item:
                    shuffle=item.split('-')
                    ex_combis.append(shuffle[1]+'-'+shuffle[0])
                if ':' in item:
                    shuffle=item.split(':')
                    ex_combis.append(shuffle[1]+':'+shuffle[0])
            rules_l.append(ex_combis)
            continue

        if i==11:
            atomsg=words[1].strip()
            #print atomsg,'atomsg'
            atoms_l=[]
            if atomsg=='None':
                pass
            else:
                atoms=atomsg.split(',')
                for atom in atoms:
                    atomi=atom.split('-')
                    atoms_l.append([atomi[0].strip(),int(atomi[1])])
            
            rules_l.append(atoms_l)
            continue
       
        if i==12:
            #print words,i
            rules_l.append(words[1].strip())
            continue

        
        if words[1].strip()!='None':
            if '-' not in words[1]:
                tmp_str = "ERROR: Wrong generation rule provided for "+lines
                tmp_str=tmp_str+"Privide the range of number required. For example, 10-20. \n"
                print_le(tmp_str,"Aborting due to wrong generation rule.")
            
            values=words[1].split('-')
            
            try:
                val_l = int(values[0])
                val_u = int(values[1])

            except ValueError:
                tmp_str = "ERROR: Wrong generation rule provided for "+lines
                tmp_str=tmp_str+"Please provide either a number or None as input.\n"
                tmp_str=tmp_str+"Use - to separate the range of numbers. For example, 10-20.\n"
                print_le(tmp_str,"Aborting due to wrong generation rule.")
            val=[val_l,val_u]
            rules_l.append(val[1])
            
        else:
            rules_l.append('None')
            
        # if i==2 or i==3 or i==4:
        #     #print words,i
        #     if words[1].strip()=='None':
        #         words[1]=0

        #     rules_l.append(int(words[1]))
        #     continue


        
        #rules_l.append(words[1].strip())
        
    print_l("\nRules list \n"+str(rules_l) )
    return rules_l


if __name__ == "__main__":
    
    ## initializing MPI to time, to check the MPI efficiency
    wt1 = MPI.Wtime()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpisize = comm.Get_size()


    if rank==0:
        logfile = open('logfile.txt','a',0)
        error_file = open('error_file.txt','a',0)
        
        libgen_classes.banner(logfile, _SCRIPT_NAME, _SCRIPT_VERSION, _REVISION_DATE, _AUTHOR, _DESCRIPTION)

    ##Initializing number of generations value
    gen_len=1
    
    ## Defining Mg atom and Helium atom
    myMg=pybel.readstring('smi',"[Mg]")
    Mgatom=myMg.OBMol.GetAtom(1)
    myHe=pybel.readstring('smi',"[He]")
    Heatom=myHe.OBMol.GetAtom(1)
    
    ## Argument parser desription
    parser = argparse.ArgumentParser(description='This is a pacakge to generate a combinatorial library of molecules based on the building blocks provided. Please provide the building blocks in the a file in either SMILES form or InChi.')
    
    parser.add_argument('-i',"--input_file", action='store', dest='file_name', default='building_blocks.dat', help="provide the building block data file. Default is building_blocks.dat file.")
    
    parser.add_argument('-t',"--molecule_type", action='store', dest='mol_type', default='SMILES', help="Mention the molecule type in this section. Default is SMILES")
    
    parser.add_argument('-c',"--combination_type", action='store', dest='combi_type', default='Link', help="Mention the combination type in this section. Link for linking and Fusion for fusing the molecules. Default is link.")

    parser.add_argument('-g',"--generation_levels", action='store', dest='gen_len', default='1', help="Give the number of maximum combinations required in each molecule. Default is 1.")

    parser.add_argument('-output',"--output_type", action='store', dest='oft', default='smi', help="Give the type of molecule format for the generated library. Default is SMILES format.")

    parser.add_argument('-max_fpf',"--max_files_per_folder", action='store', dest='max_fpf', default=1000, help="Maximu number of files that are in a single folder. Having a large number of files in a single folder may hinder performace. Default is 10000 files per folder.")

    parser.add_argument('-rf',"--rule_file", action='store', dest='rule_file', default='generation_rules.dat', help="Specified file should contain the generation rules. Order of the rules is fixed. If the order is changed then the program runs into error . Default is generation_rules.dat.")

    parser.add_argument('-n', "--lib_name", action='store', dest='lib_name', default='new_library_', help="The name of the new library to be created. Default is new_library")

    parser.add_argument('-htps',"--ChemHTPS", action='store_true', dest='chemhtps', default=False, help="Tells the library generator whether it is being run from inside the ChemHTPS program or not for purposes of output file destination. Default is False.")

    ## defining arguments
    args = parser.parse_args()
    mol_type=args.mol_type.lower()
    combi_type=args.combi_type.lower()
    max_fpf=int(args.max_fpf)
    gen_len=int(args.gen_len)
    rule_file=args.rule_file
    BB_file=args.file_name
    oft=args.oft.lower()

    ## setting output directory
    if args.chemhtps:
        output_dest = os.getcwd() + '/screeninglib/'
    else:
        output_dest = os.getcwd() + '/'
    
    print_l("Total generation levels provided = "+str(gen_len)+'\n')
    print_l("Combination type is "+str(combi_type)+'\n')
    print_l("Building blocks type is "+str(mol_type)+'\n')
    print_l("Output molecule type is "+str(args.oft).upper()+'\n')

    
    smiles_list=[]
    
    
    ## Reading the building rules from generation_rules.dat
    print_l("Reading generation rules from the file \'"+rule_file+'\'\n')

    try :
        rulesFile=open(rule_file)
    except:
        tmp_str = "Generation rules file "+rule_file+" does not exist. "
        tmp_str=tmp_str+"Please provide correct generation rule file.\n"
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
    rules_l=get_rules(rulesFile)


    ## Reading the building blocks from the input file
    print_l("Reading building blocks from the file \'"+BB_file+'\'\n')

    try :
        infile=open(BB_file)
    except:
        tmp_str = "Building blocks file "+BB_file+" does not exist. "
        tmp_str=tmp_str+"Please provide correct building blocks file.\n"
        print_le(tmp_str,"Aborting due to wrong file.")

    
    ## Read molecules provided in the input 
    for i,line in enumerate(infile):
        smiles= line.strip()
        if smiles.isspace() or len(smiles)==0 or smiles[0]=='#':
            continue
        ## if the inpur is InChI, then convert all into SMILES
        if mol_type == 'inchi':
            if check_if_inchi(smiles)==False:
                tmp_str='Error: The InChI string (\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct InChI.'.format(smiles,i+1,args.file_name)
                print_le(tmp_str,"Aborting due to wrong file.")
                
            this_mol=pybel.readstring("inchi",smiles)
            smiles=str(this_mol)
            smiles=smiles.strip()

        ## check if smiles
        if check_if_smiles(smiles)==False:
            tmp_str='Error: The SMILES string(\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct SMILES.'.format(smiles,i+1,args.file_name)
            print_le(tmp_str,"Aborting due to wrong file.")
            kill_me=comm.gather(tmp_str,root=0)

        if check_if_smiles(smiles):
            smiles_list.append(smiles)
    print_l('Number of buidling blocks provided = '+str(len(smiles_list))+'\n')
    
    print_l('SMILES of the building blocks\n')
    
    print_l(smiles_list)
    print_l('=============================================================================\n')
    
    ##assigning the code for each building block
    
    smiles_list_c=[]
    for i,item in enumerate(smiles_list):
        
        smiles_list_c.append([smiles_list[i],'F'+str(i+1)])

    ## generation_test funtion generates combinatorial molecules
    final_list=generation_test(combi_type)

    final_list_len=len(final_list)

    print_l('Total number of molecules generated = '+str(final_list_len)+'\n')
    
    outdata=output_dest + "Final_smiles_output.dat"
    outfile = open(outdata, "w")

    print_l('Writing molecules SMILES to file \''+outdata+'\' along with corresponding code.\n')

    outfile.write('Sl.No,Molecule_Smiles,Combination_Code\n')
    
    for i, smiles in enumerate(final_list):
        outfile.write(str(i+1)+','+smiles[0]+','+smiles[1]+'\n')
              

    if oft=='smi':
        outdata=output_dest + "Final_smiles_output.smi"
        outfile = open(outdata, "w")
        
        print_l('Writing molecules SMILES to file \''+outdata+'\'\n')
        
        for i, smiles in enumerate(final_list):
            outfile.write(smiles[0]+'\n')
        
    ## Creating a seperate output file for each Molecule.
    ## The files are written to folder with specified no. of files per folder.

    if oft!='smi':
        
        print_l('Writing molecules with molecule type '+str(oft)+'\n')
    
        if not os.path.exists(output_dest + args.lib_name + oft):
            os.makedirs(output_dest + args.lib_name + oft)
        smiles_to_scatter=[]
        if rank ==0:
            smiles_to_scatter=[]
            for i in xrange(mpisize):
                start=int(i*(final_list_len)/mpisize)
                end=int((i+1)*(final_list_len)/mpisize)-1
                list_to_add=final_list[start:end+1]
                list_to_add=list_to_add+[final_list_len,start,end]
                smiles_to_scatter.append(list_to_add)
        else:
            smiles_to_scatter=[]
    
        ## Dividing the list into processors
        smiles_list=comm.scatter(smiles_to_scatter,root=0)
        
        final_list_len=smiles_list[-3]
        start=smiles_list[-2]
        end=smiles_list[-1]
        smiles_list=smiles_list[0:-3]
        
        ratio_s=int(start/max_fpf)
        ratio_e=int(end/max_fpf)
        
        if end+1==final_list_len:
            ratio_e=ratio_e+1
        
        for i in xrange(ratio_s,ratio_e):
            
            if not os.path.exists(output_dest + args.lib_name + oft+"/"+str(i+1)+"_"+str(max_fpf)):
                os.makedirs(output_dest + args.lib_name + oft+"/"+str(i+1)+"_"+str(max_fpf))

        folder_no=ratio_s+1
        for i,val in enumerate(xrange(start,end+1)):
            mymol= pybel.readstring("smi",smiles_list[i][0])
            mymol.make3D(forcefield='mmff94', steps=50)
            mymol.write(oft, output_dest + args.lib_name +oft+"/"+str(folder_no)+"_"+str(max_fpf)+"/"+str(val+1)+"."+oft,overwrite=True)

            if (val+1)%max_fpf==0:
                folder_no=folder_no+1
            
            

    print_l('File writing terminated successfully'+'\n')
                
    wt2 = MPI.Wtime()
    
    print_l('Total time_taken '+str('%.3g'%(wt2-wt1))+'\n')
    
    
    sys.stderr.close()
    
