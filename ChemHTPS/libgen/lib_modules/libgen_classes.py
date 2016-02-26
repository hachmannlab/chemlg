import os
import sys
import threading
import pybel
from openbabel import OBAtomAtomIter
import time



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
        if atom.OBAtom.GetAtomicNum() ==88:            
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
                                atom_pair_list= atom_pair_list+ atom_pair_list_tmp
                        
                            
    return atom_pair_list, len(atoms)


def get_fused_mol(smiles1,smiles2):

    list1,size1= get_atom_pair_list(smiles1,1)
    # print list1,'list1'
    list2,size2=get_atom_pair_list(smiles2,2)
    # print list2,'list2'
    smiles_combi=smiles1+'.'+smiles2
    lib_can=[]
    lib_can_nRa=[]
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
                lib_can.append([str(can_mol_combi),mol_wt])
                atoms=list(mol_combi_new.atoms)
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
                    lib_can.append([str(can_mol_combi),mol_wt])
                    atoms=list(mol_combi_new.atoms)
                    for atom in atoms:
                        if atom.OBAtom.GetAtomicNum()==88:
                            atom.OBAtom.SetAtomicNum(1)
                    lib_can_nRa.append(str(can_mol_combi))
    
    return lib_can

def get_fused_mol_c(smiles1,smiles2,rules):
    list1,size1= get_atom_pair_list(smiles1[0],1)
    #print list1,'list1'
    list2,size2=get_atom_pair_list(smiles2[0],2)
    #print list2,'list2'
    smiles_combi=smiles1[0]+'.'+smiles2[0]
    lib_can=[]
    lib_can_nRa=[]
    code=smiles1[1]+':'+smiles2[1]
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
                        
                if if_add(mol_combi_new,mol_wt,rules,code,'f')==True and shd_add==True:
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
                    if if_add(mol_combi_new,mol_wt,rules,code,'f')==True and shd_add==True:
                        lib_can.append([str(can_mol_combi),mol_wt])

                        #print mol_combi_new,'after'
                    for atom in atoms:
                        if atom.OBAtom.GetAtomicNum()==88:
                            atom.OBAtom.SetAtomicNum(1)
                    
                    lib_can_nRa.append(str(can_mol_combi))
    lib_can_c=[]
    for item in lib_can:
        lib_can_c.append([item[0],item[1],smiles1[1]+':'+smiles2[1]])
            
    
    return lib_can_c
        
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

def get_fused_c(smiles1,smiles2,rules):
    #smiles1=smiles1[:-2]
    #smiles2=smiles2[:-2]
    #print smiles1,smiles2
    lib_can1=get_fused_mol_c(smiles1,smiles2,rules)
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

def get_num_struc(mol,smarts):
    smarts = pybel.Smarts(smarts)
 
    smarts.obsmarts.Match(mol.OBMol)
    num_matches = sum(1 for indicies in smarts.obsmarts.GetMapList())
    num_unique_matches = len(smarts.findall(mol))
    return num_unique_matches

HBD = pybel.Smarts("[#7,#8;!H0]")

HBA = pybel.Smarts("[#7,#8]")

def lipinski(mol):

   """Return the values of the Lipinski descriptors."""

   desc = {'molwt': mol.molwt,

      'HBD': len(HBD.findall(mol)),

      'HBA': len(HBA.findall(mol)),

      'logP': mol.calcdesc(['logP']) ['logP']}

   return desc

passes_all_rules = lambda desc: (desc ['molwt'] <= 500 and

         desc ['HBD'] <= 5 and desc ['HBA'] <= 10 and

         desc ['logP'] <= 5)


def if_add(mol,mol_wt,rules,code,c_type='l'):
    
    #print rules
    if c_type=='f':
        mol=pybel.readstring('smi',str(mol))

    #print mol,'before'
    atoms=list(mol.atoms)
    del_idx=[]

    ## iterating over all atoms of the molecule
    for atom in atoms:
        ## Removing Francium atoms
        if atom.OBAtom.GetAtomicNum()==87:
            ## it is easy to convert Francium atom to hydrogen than deleting the atom
            atom.OBAtom.SetAtomicNum(1)
            
            ## Removing Radium atoms
        if atom.OBAtom.GetAtomicNum()==88:
            atom.OBAtom.SetAtomicNum(1)

    mol_wt=mol.OBMol.GetMolWt()
    #print mol,'after'

    if code in rules[0]:
        
        #print 'code',rules[3]
        return False

    # Calculating no.of rings
    bonds = mol.OBMol.NumBonds()
    #print rules[1][1]
    
    #print rules[1][1]
    if bonds>rules[1][1] and rules[1][1]!=0:
        #print 'bonds'
        return False
    no_atoms=len(mol.atoms)

    if no_atoms>rules[2][1] and rules[2][1]!=0:
        #print 'no_atoms'
        return False

    if int(mol_wt)>rules[3][1] and rules[3][1]!=0:
        #print 'mol_wt'
        return False

    if isinstance(rules[4][1],(int)):
        rings = len(mol.OBMol.GetSSSR())
        if rings>int(rules[4][1]):
            #print 'rings'
            return False

    if isinstance(rules[5][1],(int)) or isinstance(rules[6][1],(int)):
        no_ar=0
        no_non_ar=0
        for r in mol.OBMol.GetSSSR():
            if r.IsAromatic():
                no_ar=no_ar+1
            else:
                no_non_ar=no_non_ar+1
        if isinstance(rules[5][1],(int)):
            if no_ar>int(rules[5][1]) :
                return False
        if isinstance(rules[6][1],(int)):
            if no_ar>int(rules[6][1]) :
                return False

    if isinstance(rules[7][1],(int)):
        no_s_bonds=get_num_struc(mol,"*-*")
        if no_s_bonds>int(rules[7][1]) :
            return False

    if isinstance(rules[8][1],(int)):
        no_d_bonds=get_num_struc(mol,"*=*")
        if no_d_bonds>int(rules[8][1]) :
            return False

    if isinstance(rules[9][1],(int)):
        no_t_bonds=get_num_struc(mol,"*#*")
        if no_t_bonds>int(rules[9][1]) :
            return False
    
    
    if isinstance(rules[10],list):
        for item in rules[10]:
            no_at=0
            if item[0]=='C'or'S'or'N'or'O'or'c'or's'or'n'or'o':
                no_at=get_num_struc(mol,item[0].lower())
                no_at=no_at+get_num_struc(mol,item[0].upper())
            else:
                no_at=get_num_struc(mol,item[0])
            #print no_at,'no_atom'
            if no_at>item[1]:
                #print no_atom,"hello"
                return False
    
    if rules[11].lower()=='true':

        descriptors = lipinski(mol)
        
        if not passes_all_rules(descriptors):
            #print "hello"
            return False

    return True
    
    

def create_link_c(smiles1,smiles2,rules):
    library_full=[]
    
    mol1_index_list, size1, smiles1=get_index_list_c(smiles1)
            
    mol2_index_list, size2, smiles2=get_index_list_c(smiles2)
    
    smiles_combi=smiles1[0]+'.'+smiles2[0]
    
    for index1 in mol1_index_list:
        for index2 in mol2_index_list:
            mol_combi= pybel.readstring("smi",smiles_combi)
            mol_combi.OBMol.AddBond(index1,index2+size1,1,0,-1)
            can_mol_combi = mol_combi.write("can")
            mol_wt=str(int(mol_combi.OBMol.GetMolWt()))
            code=smiles1[1]+'-'+smiles2[1]
            if if_add(mol_combi,mol_wt,rules,code)==True:
                library_full.append([str(can_mol_combi),mol_wt,code])
    return library_full
    

## Defining Ra atom and Francium atom
myRa=pybel.readstring('smi',"[Ra]")
Raatom=myRa.OBMol.GetAtom(1)
myFr=pybel.readstring('smi',"[Fr]")
Fratom=myFr.OBMol.GetAtom(1)

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

        # Attach Francium atom. Makes it easy to remove duplicates
        newmol.OBMol.InsertAtom(Fratom) 
        # Get index of atom
        index=atom.OBAtom.GetIdx() 
        # Create a bond between Fr atom and the curent atom
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

def get_index_list_c(smiles):
    
    mol=pybel.readstring("smi",smiles[0])
            
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
        
        newmol=pybel.readstring("smi",smiles[0])

        # Attach Francium atom. Makes it easy to remove duplicates
        newmol.OBMol.InsertAtom(Fratom) 
        # Get index of atom
        index=atom.OBAtom.GetIdx() 
        # Create a bond between Fr atom and the curent atom
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


## from lib_jcode.py
def banner(logfile, SCRIPT_NAME, SCRIPT_VERSION, REVISION_DATE, AUTHOR, DESCRIPTION,):
    """(banner):
        Banner for this little script.
    """
    str = []
    str.append("============================================================================== ")
    str.append(SCRIPT_NAME + " " + SCRIPT_VERSION + " (" + REVISION_DATE + ")")
    str.append(AUTHOR)
    str.append("============================================================================== ")
    str.append(time.ctime())
    str.append("")    
    str.append(DESCRIPTION)
    str.append("")

    print 
    for line in str:
        print line
        logfile.write(line + '\n')
