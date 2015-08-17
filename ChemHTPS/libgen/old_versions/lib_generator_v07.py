#!/usr/bin/env python

_SCRIPT_NAME = "Library_Generator"
_SCRIPT_VERSION = "v0.1.7"
_REVISION_DATE = "8/11/2015"
_AUTHOR = "M. Atif Afzal (m27@buffalo.edu)"
_DESCRIPTION = "This is the a library generating molecular libraries."

# Version history timeline:
# v0.1.1 (7/12/2015): Made sure the Jchem (reactor) suite works  
# v0.1.2 (7/24/2015): Use Openbabel to create the library and react
# v0.1.3 (8/01/2015): link only carbons which have atleast one hydrogen bond 
# v0.1.4 (8/03/2015): Getting rid of jchem dependency 
# v0.1.5 (8/05/2015): Parse arguments 
# v0.1.6 (8/08/2015): Get only those reactions specified by the user
# v0.1.7 (8/11/2015): Implementing parallel algorithm

###################################################################################################
# TASKS OF THIS SCRIPT:
# -Molecular library generation
###################################################################################################

###################################################################################################
# TODO:
# -Get the reaction without reactor-Done
# -Try to run reactor from python instead subprocess -Don't need this anymore
# -Try to parallelize the code
# -Include stopping creteria
# -Import SMILES from data file -Done
# -Use error handling properly -Partially done
# -Implement a smart duplicate removal system
# -Capture the stderr from C program ie. openbabel. -Done (this took a very long time to debug)
# 

###################################################################################################

import sys
#sys.path.insert(0, "/user/m27/pkg/openbabel/2.3.2/lib")
import pybel
import openbabel
import subprocess
import time
import argparse
#import scipy
from collections import defaultdict
#from my_classes import OutputGrabber

###################################################################################################

#import cStringIO
from mpi4py import MPI
import os
import threading

#out = OutputGrabber(sys.stderr)


# class OutputGrabber(object):
#     """
#     Class used to grab standard output or another stream.
#     """
#     escape_char = "\b"

#     def __init__(self, stream=None, threaded=False):
#         self.origstream = stream
#         self.threaded = threaded
#         if self.origstream is None:
#             self.origstream = sys.stdout
#         self.origstreamfd = self.origstream.fileno()
#         self.capturedtext = ""
#         #print 'origstreamfd',self.origstreamfd
#         # Create a pipe so the stream can be captured:
#         self.pipe_out, self.pipe_in = os.pipe()
#         #print self.pipe_in
#         pass

#     def start(self):
#         """
#         Start capturing the stream data.
#         """
#         #self.capturedtext = ""
#         # Save a copy of the stream:
#         #print os.read(self.origstreamfd,100)
#         self.streamfd = os.dup(self.origstreamfd)
#         #print 'origstreamfd',self.origstreamfd
#         #print 'streamfd',self.streamfd
#         # Replace the Original stream with our write pipe
#         os.dup2(self.pipe_in, self.origstreamfd)
#         #print 'pipe', self.pipe_in
#         if self.threaded:
#             # Start thread that will read the stream:
#             self.workerThread = threading.Thread(target=self.readOutput)
#             self.workerThread.start()
#             # Make sure that the thread is running and os.read is executed:
#             time.sleep(0.01)
#         pass

#     def stop(self):
#         """
#         Stop capturing the stream data and save the text in `capturedtext`.
#         """
#         # Flush the stream to make sure all our data goes in before
#         # the escape character.
        
#         self.origstream.flush()
#         #os.fsync(self.origstreamfd)
#         #print self.origstreamfd
#         #os.fsync(self.origstreamfd)
#         #self.pipe_out_file=os.fdopen(self.pipe_out,'r')
#         #print type(self.pipe_in)
#         #os.fsync(self.pipe_out)
        
#         #print os.read(self.origstream,1)
#         # Print the escape character to make the readOutput method stop:
# #        self.origstream.write(self.escape_char)
#         # if self.threaded:
#         #     # wait until the thread finishes so we are sure that
#         #     # we have until the last character:
#         #     self.workerThread.join()
#         # else:
#         #     self.readOutput()
#         # Close the pipe:
#         # Restore the original stream:
        
#         os.dup2(self.streamfd, self.origstreamfd)
        
#         os.close(self.pipe_out)
        
#         #os.close(self.origstreamfd)
#         #os.close(self.streamfd)
#         pass

#     # def readOutput(self):
#     #     """
#     #     Read the stream data (one byte at a time)
#     #     and save the text in `capturedtext`.
#     #     """
#     #     while True:
#     #         data = os.read(self.pipe_out, 1)  # Read One Byte Only
#     #         if self.escape_char in data:
#     #             break
#     #         if not data:
#     #             break
#     #         self.capturedtext += data
#     #     pass


class OutputGrabber(object):
    """
    Class used to grab standard output or another stream.
    """
    escape_char = "z"

    def __init__(self, stream=None, threaded=False):
        self.origstream = stream
        self.origstream3 = stream
        self.threaded = threaded
        if self.origstream is None:
            self.origstream = sys.stdout
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
        if self.threaded:
            # Start thread that will read the stream:
            self.workerThread = threading.Thread(target=self.readOutput)
            self.workerThread.start()
            # Make sure that the thread is running and os.read is executed:
            time.sleep(0.01)
        pass

    def stop(self):
        """
        Stop capturing the stream data and save the text in `capturedtext`.
        """
        # Flush the stream to make sure all our data goes in before
        # the escape character.
        #self.origstream.flush()
 #       self.origstream2.flush()
  #      self.origstream3.flush()
    
        # Print the escape character to make the readOutput method stop:
        #self.origstream.write('z')
        if self.threaded:
            # wait until the thread finishes so we are sure that
            # we have until the last character:
            self.workerThread.join()
        #else:
            #self.readOutput()
        # Close the pipe:
        # Restore the original stream:
        os.dup2(self.streamfd, self.origstreamfd)
        os.close(self.pipe_out)
        pass


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
        hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
        if hcount==0:
            continue
        ## remove duplicates
        newmol=pybel.readstring("smi",smiles)
        newmol.OBMol.InsertAtom(Heatom)
        index=atom.OBAtom.GetIdx()
        newmol.OBMol.AddBond(index,size+1,1,0,-1)
        
        out.start()
        can_smiles = newmol.write("can")
        out.stop()
#        print 'me'
        #print str(sys.stderr)
        #print can_smiles
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
    library_full=[]
    library=[]
    #library_gather=[]
    smiles_dict = defaultdict(list) 
    for smiles in smiles_list_gen+ini_list:
        mol_combi= pybel.readstring("smi",smiles)
        mol_wt=str(int(mol_combi.OBMol.GetMolWt()))
        can_mol_combi = mol_combi.write("can")
                
        if can_mol_combi in smiles_dict[mol_wt]:
            continue
        smiles_dict[mol_wt].append(can_mol_combi)
        
        # if can_mol_combi in library_can:
        #     continue
        # library_can.append(can_mol_combi)
                        

    if rank!=0:
        smiles_gen_list=None
    
    if rank ==0:
        smiles_to_scatter=[]
        for i in xrange(mpisize):
            start=int(i*(len(smiles_list_gen))/mpisize)
            end=int((i+1)*(len(smiles_list_gen))/mpisize)-1
            smiles_to_scatter.append(smiles_list_gen[start:end+1])
    else:
        smiles_to_scatter=[]
    
    smiles_list_gen=comm.scatter(smiles_to_scatter,root=0)

    for smiles1 in smiles_list_gen:
        for smiles2 in ini_list:
            
            mol1_index_list, size1, smiles1=get_index_list(smiles1)
            
            mol2_index_list, size2, smiles2=get_index_list(smiles2)
            
            smiles_combi=smiles1+'.'+smiles2
            
            for index1 in mol1_index_list:
                for index2 in mol2_index_list:
                    
                    
                    mol_combi= pybel.readstring("smi",smiles_combi)
                    mol_combi.OBMol.AddBond(index1,index2+size1,1,0,-1)
                    library_full.append(str(mol_combi)[:-2])
            
    library_gather=comm.gather(library_full,root=0)
    
    if rank ==0:
        for l1 in library_gather:
            for l2 in l1:
                mol_combi= pybel.readstring("smi",l2)
                mol_wt=str(int(mol_combi.OBMol.GetMolWt()))
                can_mol_combi = mol_combi.write("can")
                    
                if can_mol_combi in smiles_dict[mol_wt]:
                    continue
                smiles_dict[mol_wt].append(can_mol_combi)

                # if can_mol_combi in library_can:
                #     continue
                # library_can.append(can_mol_combi)
                
                library.append(str(mol_combi)[:-2])

        return library
    else:
        return []       
        #print 'library_gatherd_length',len(library_gather)
    
    #library_gather=comm.bcast(library_gather,root=0)
    
    # if rank==0:
    #     for combi in library_gather:
    #         mol_combi= pybel.readstring("smi",combi)
    #         mol_wt=str(int(mol.OBMol.GetMolWt()))
    #         can_mol_combi = mol_combi.write("can")
    #         if can_mol_combi in library_can:
    #             continue
    #         library_can.append(can_mol_combi)
    #         library.append(str(mol_combi)[:-2])
            
    #     full_library_list=smiles_list_gen+library
    #     return library
    # else:
    #     return []

def reverse_mol(smiles):
    
    mol=pybel.readstring("smi",smiles)            
    atoms=list(mol.atoms)
    atom_num=[]
    for atom in atoms:
        atom_num.append(atom.OBAtom.GetAtomicNum())
    if 12 in atom_num:
        newmol=pybel.readstring("smi",smiles)
        del_index=[]
        for atom in atoms:
            hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()

            if hcount==0:
                continue
            index=atom.OBAtom.GetIdx()
            this_atom=mol.OBMol.GetAtom(index)
            atom_num_n=[]
            
            for atom_n in openbabel.OBAtomAtomIter(this_atom):
                check_12=atom_n.GetAtomicNum()
                atom_num_n.append(check_12)
            if 12 in atom_num_n:
                continue
                        
            while hcount!=0:
                size=len(list(mol.atoms))
                mol.OBMol.InsertAtom(Heatom)
                mol.OBMol.AddBond(index,size+1,1,0,-1)
                hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
        
        atoms=list(mol.atoms)
        del_idx=[]
        for atom in atoms:
            if atom.OBAtom.GetAtomicNum()==12:
                del_idx.append(atom.OBAtom.GetIdx())
        for idx in reversed(del_idx):
            mol.OBMol.DeleteAtom(mol.OBMol.GetAtom(idx))
                   
    smiles=str(mol)[:-2]
    return smiles
    
def generation_test():
    lib_smiles_list=[]
    lib_can=[]
    
    for smiles in smiles_list:
        mol_combi= pybel.readstring("smi",smiles)
        mol_wt=str(int(mol_combi.OBMol.GetMolWt()))
        can_mol_combi = mol_combi.write("can")
        if can_mol_combi in lib_can:
            continue
        lib_can.append(can_mol_combi)
        lib_smiles_list.append(smiles)
   
    for i in xrange(len(lib_smiles_list)):
        lib_smiles_list[i]=reverse_mol(lib_smiles_list[i])

    #print lib_smiles_list
    for gen in xrange(gen_len):
        lib_smiles_list=lib_smiles_list+create_gen_lev(lib_smiles_list,smiles_list)
        #print len(lib_smiles_list)
    
    if rank ==0:
        print len(lib_smiles_list)
    ## test to see if there are any duplicates
    # library_can=[]
    # pos=0
    # for smiles in lib_smiles_list:
    #     mol_combi= pybel.readstring("smi",smiles)
    #     atoms=list(mol_combi.atoms)
    #     del_idx=[]
    #     for atom in atoms:
    #         if atom.OBAtom.GetAtomicNum()==2:
    #             del_idx.append(atom.OBAtom.GetIdx())
    #     for idx in reversed(del_idx):
    #         mol_combi.OBMol.DeleteAtom(mol_combi.OBMol.GetAtom(idx))
        
    #     lib_smiles_list[pos]=str(mol_combi)[:-2]
    #     can_mol_combi = mol_combi.write("can")
    #     pos=pos+1
    #     if can_mol_combi in library_can:
    #         print can_mol_combi
    #         continue
    #     library_can.append(can_mol_combi)
    
    # if rank==0:
    #     print len(library_can)
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
    
    wt1 = MPI.Wtime()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpisize = comm.Get_size()
    #print rank,mpisize

    smiles1="c1(C)ccccc1" 
    smiles2="c1occc1"
    gen_len=1
    myMg=pybel.readstring('smi',"[Mg]")
    Mgatom=myMg.OBMol.GetAtom(1)
    myHe=pybel.readstring('smi',"[He]")
    Heatom=myHe.OBMol.GetAtom(1)

    out = OutputGrabber(sys.stderr)

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
            if check_if_inchi(smiles)==False and rank==0:
                message='Error: The InChI string(\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct InChI.'.format(smiles,i+1,args.file_name)
                if rank==0:
                    print message
                sys.exit()
            this_mol=pybel.readstring("inchi",smiles)
            smiles=str(this_mol)
            smiles=smiles.strip()

        ## check if smiles
        if check_if_smiles(smiles)==False:
            message='Error: The SMILES string(\'{}\') provided in line {} of data file \'{}\' is not valid. Please provide correct SMILES.'.format(smiles,i+1,args.file_name)
            if rank==0:
                print message
            sys.exit()
            
        if check_if_smiles(smiles):
            smiles_list.append(smiles)
    
    if rank==0:
        print 'Number of smiles provided =',len(smiles_list)
    
        print smiles_list


    final_list=generation_test()

    wt2 = MPI.Wtime()
    if rank==0:    
        print 'time_taken',wt2-wt1
    
    #print final_list
    
    # file = open("Final_smiles_output.dat", "w")
    # k=1
    # file.write('Sl.No,Molecule_Smiles\n')

    # for smiles in final_list:
    #     file.write(str(k)+','+smiles+'\n')
    #     k=k+1
                
    
    sys.stderr.close()
