import os
import sys
import threading
import pybel
from openbabel import OBAtomAtomIter




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


## 
def get_atom_pair_list(smiles,type1or2):
    mol=pybel.readstring('smi',smiles)
    atoms=list(mol)
    atom_pair_list=[]
    for atom in atoms:
        if atom.OBAtom.GetAtomicNum() ==12:            
            index=atom.OBAtom.GetIdx()
            for atom1 in OBAtomAtomIter(mol.OBMol.GetAtom(index)):
                if atom1.GetAtomicNum() !=6: 
                    continue
                atom1_idx=atom1.GetIdx()
                atom_atoms=[]
                if type1or2==1:

                    for atom2 in OBAtomAtomIter(atom1):
                        if atom2.GetAtomicNum() !=6 or atom2.IsInRing()==False: 
                            continue
                        hcount = atom2.ExplicitHydrogenCount() + atom2.ImplicitHydrogenCount()
                        if hcount==0:
                            continue
                        atom_atoms.append(atom2.GetIdx())
                    if len(atom_atoms)==2:
                        for aa in atom_atoms:
                            atom_pair=[atom1_idx,aa,index]
                            atom_pair_list.append(atom_pair)
                    if len(atom_atoms)==1:
                        atom_pair_list.append([atom1_idx,atom_atoms[0],index])
                if type1or2==2:
                    for atom2 in OBAtomAtomIter(atom1):
                    
                        if atom2.GetAtomicNum() ==12 or atom2.IsInRing()==False: 
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
                                atom_pair_list= atom_pair_list+ atom_pair_list_tmp
                        
                            
    return atom_pair_list, len(atoms)


def get_fused_mol(smiles1,smiles2):
    list1,size1= get_atom_pair_list(smiles1,1)
    #print list1,'list1'
    list2,size2=get_atom_pair_list(smiles2,2)
    #print list2,'list2'
    smiles_combi=smiles1+'.'+smiles2
    lib_can=[]
    lib_can_nMg=[]
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
            
            if can_mol_combi not in lib_can_nMg:
                mol_wt=str(int(mol_combi_new.OBMol.GetMolWt()))
                lib_can.append([str(can_mol_combi),mol_wt])
                atoms=list(mol_combi_new.atoms)
                for atom in atoms:
                    if atom.OBAtom.GetAtomicNum()==12:
                        atom.OBAtom.SetAtomicNum(1)
                lib_can_nMg.append(str(can_mol_combi))
    
                
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
                
                if can_mol_combi not in lib_can_nMg:
                    mol_wt=str(int(mol_combi_new.OBMol.GetMolWt()))
                    lib_can.append([str(can_mol_combi),mol_wt])
                    atoms=list(mol_combi_new.atoms)
                    for atom in atoms:
                        if atom.OBAtom.GetAtomicNum()==12:
                            atom.OBAtom.SetAtomicNum(1)
                    lib_can_nMg.append(str(can_mol_combi))
    
    return lib_can
        
def get_fused(smiles1,smiles2):
    #smiles1=smiles1[:-2]
    #smiles2=smiles2[:-2]
    #print smiles1,smiles2
    lib_can1=get_fused_mol(smiles1,smiles2)
    #lib_can2=get_fused_mol(smiles2,smiles1)
    # for list1 in lib_can1:
    #    if list1[0] in
    #lib_can=lib_can1+lib_can2
    lib_can=lib_can1
    #print lib_can
    return lib_can


## This function creates all possible links when two smiles are provided
def create_link(smiles1,smiles2):
    library_full=[]
    
    mol1_index_list, size1, smiles1=get_index_list(smiles1)
            
    mol2_index_list, size2, smiles2=get_index_list(smiles2)
    
    smiles_combi=smiles1+'.'+smiles2
    
    for index1 in mol1_index_list:
        for index2 in mol2_index_list:
            mol_combi= pybel.readstring("smi",smiles_combi)
            mol_combi.OBMol.AddBond(index1,index2+size1,1,0,-1)
            can_mol_combi = mol_combi.write("can")
            mol_wt=str(int(mol_combi.OBMol.GetMolWt()))
            library_full.append([str(can_mol_combi),mol_wt])
    return library_full
    

## Defining Mg atom and Helium atom
myMg=pybel.readstring('smi',"[Mg]")
Mgatom=myMg.OBMol.GetAtom(1)
myHe=pybel.readstring('smi',"[He]")
Heatom=myHe.OBMol.GetAtom(1)

## Openbabel sometime prints out an error related to sterochemistry. 
## This is usully the case when the molecule do not have any sterechemistry information.
## Because this error is from OpenBabel which is written in C, it is quite challenging 
## to capture error statements. 
## OutputGrabber is a class which capture system errors from any program (including C0. 
## We will have to be careful with this class because it might also capture system errors
## not related to stereochemistry and we would never know what the system error was. 

out = OutputGrabber(sys.stderr)

'''
 This function returns the list with index numbers of atoms that can be reacted. It also returns the corresponsing smiles string.
'''
def get_index_list(smiles):
    
    mol=pybel.readstring("smi",smiles)            
    atoms=list(mol.atoms)
    size=len(atoms)

    can_smiles_list=[]
    atoms_index=[]
    atom_num=[]

    for atom in atoms:
        
        # Counting the number of hydrogens attached to the atom
        hcount = atom.OBAtom.ExplicitHydrogenCount() + atom.OBAtom.ImplicitHydrogenCount()
        if hcount==0:
            continue # Do not do anything if there are no hydrogens attached
        
        newmol=pybel.readstring("smi",smiles)

        # Attach Helium atom. Makes it easy to remove duplicates
        newmol.OBMol.InsertAtom(Heatom) 
        # Get index of atom
        index=atom.OBAtom.GetIdx() 
        # Create a bond between He atom and the curent atom
        newmol.OBMol.AddBond(index,size+1,1,0,-1) 
        
        out.start() # Capturing stderr from openbabel (C program) 
        # Making use of canonical smiles to remove duplicates
        can_smiles = newmol.write("can")
        out.stop() # Closing capture of stderr

        if can_smiles in can_smiles_list:
            continue
        can_smiles_list.append(can_smiles)
        index=atom.OBAtom.GetIdx()
        atoms_index.append(index)
        
    return atoms_index, size, smiles


'''
This is to remove the stereo chemistry information from smiles provided.
'''
def remove_stereochemistry(smiles):
    return smiles.replace("@", "").replace("/", "-").replace("\\", "-")

