from __future__ import print_function
from ipywidgets import interact, Layout
import ipywidgets as widgets
#import ipywidgets.Layout as ly
from IPython.display import display, clear_output
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import sys
from IPython.display import Javascript
sys.path.insert(0, "/user/m27/pkg/openbabel/2.3.2/lib")
import pybel


style = {'description_width': 'initial'}
count =0
BB_list=[]

"""def GUI_1():
    global BB_list
    style = {'description_width': 'initial'}

    BB = widgets.Text(
        placeholder='Type the SMILES of building block',
        description='Building block:',style=style,
         )
    display(BB)
    button = widgets.Button(description="Add the building block")
    display(button)
    BB_dis = widgets.HTML(
        description='Building blocks added:',style=style,
        )
    display(BB_dis)
    
    button_BB_dis = widgets.Button(description="Visualize building block")
    display(button_BB_dis)

    def on_button_clicked(b):
        global count, BB_list
        count=count+1
        display(Javascript('IPython.notebook.execute_cells([3])'))
        if count>0:
            BB_list.append(BB.value)
            building_blocks2=open("building_blocks.txt","a+")
            building_blocks2.write(BB.value + '\n')
            building_blocks2.close()
            BB_dis.value = str(BB_list)

        else:
            building_blocks1 = open("building_blocks.txt", "w+")
            building_blocks1.write('Building blocks are:' + '\n' + BB.value + '\n')
            building_blocks1.close()
            count=2

    def on_dis_bb_clicked(b):

        global BB_list
        mol_list = [Chem.MolFromSmiles(smiles) for smiles in BB_list]
        [mol.SetProp('_Name', 'B'+str(i)) for i,mol in enumerate(mol_list)]
        ibu=Chem.Draw.MolsToGridImage(mol_list,legends=[mol.GetProp('_Name') for mol in mol_list])
        display(ibu)
        #BB_viz()

                        
            
    button.on_click(on_button_clicked)
    button_BB_dis.on_click(on_dis_bb_clicked)

def BB_viz():
    
    global BB_list
    button_BB_dis = widgets.Button(description="Visualize building block")
    display(button_BB_dis)

    def on_dis_bb_clicked(b):

        global BB_list
        clear_output()
        button_BB_dis = widgets.Button(description="Visualize building block")
        display(button_BB_dis)
        mol_list = [Chem.MolFromSmiles(smiles) for smiles in BB_list]
        [mol.SetProp('_Name', 'B'+str(i)) for i,mol in enumerate(mol_list)]
        ibu=Chem.Draw.MolsToGridImage(mol_list,legends=[mol.GetProp('_Name') for mol in mol_list])
        display(ibu)
            
    button_BB_dis.on_click(on_dis_bb_clicked)"""

def GUI_2():
    style = {'description_width': 'initial','font_weight':'bold'}
    space_box = widgets.Box(layout=widgets.Layout(height ='55px', width='90%')) 
    second=widgets.Button(description='Proceed to the next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    
    third=widgets.Button(description='Proceed to the next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    
    display (second)      
    def second_section(q):
        display(widgets.HTML(value="<font color=crimson><font size=5><b><u>SECTION 2: GENERATION RULES</font>"))

        style = {'description_width': 'initial','font_weight':'bold'}
        building_blocks = widgets.Text(value='None',
                                   placeholder='Type the building blocks to be used',
                                   description='Building Blocks:',
                                   style=style,
                                   )    
        bb_intro=widgets.HTML("""Specify the building blocks which must be present in all the molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        bb=widgets.VBox(children=[bb_intro,space_box,building_blocks])
                                
        style = {'description_width': 'initial'}
        bondmin = widgets.Text(description='Minimum no of bonds', value='None', style=style)
        bondmax = widgets.Text(description='Maximum no of bonds', value='None', style=style)
        bonds_intro=widgets.HTML("""This the minimum and maximum bonds section wherein the user can enter the minimum and maximum number of bonds requirements""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        x = widgets.HBox([bondmin, bondmax])
        no_bonds=widgets.VBox(children=[bonds_intro,space_box,x
                                ])

        atommin = widgets.Text(description='Minimum No of Atoms', value='None', style=style)
        atommax = widgets.Text(description='Maximum No of Atoms', value='None', style=style)
        y = widgets.HBox([atommin, atommax])
        atom_intro=widgets.HTML("""Specify the minimum and maximum number of bonds that must be present in each molecule""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        atom=widgets.VBox(children=[atom_intro,space_box,y])

        molmin = widgets.Text(description='Minimum Molecular Weight', value='None', style=style)
        molmax = widgets.Text(description='Maximum Molecular Weight', value='None', style=style)
        z = widgets.HBox([molmin, molmax])
        mol_intro=widgets.HTML("""Specify the range of the molecular weight of the molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        mol=widgets.VBox(children=[mol_intro,space_box,z])
        
        ringmin = widgets.Text(description='Minimum no of Rings', value='None', style=style)
        ringmax = widgets.Text(description='Maximum no of Rings', value='None', style=style)
        a = widgets.HBox([ringmin, ringmax])
        ring_intro=widgets.HTML("""Specify the range of the number of rings present in the  molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        ring=widgets.VBox(children=[ring_intro,space_box,a])

        armin = widgets.Text(description='Minimum no of Aromatic rings', value='None', style=style)
        armax = widgets.Text(description='Maximum no of Aromatic Rings', value='None', style=style)
        b = widgets.HBox([armin, armax])
        aring_intro=widgets.HTML("""Specify the range of the number of aromatic rings present in the  molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        aring=widgets.VBox(children=[aring_intro,space_box,b])

        narmin = widgets.Text(description='Minimum no of Non Aromatic rings', value='None', style=style)
        narmax = widgets.Text(description='Maximum no of Non Aromatic Rings', value='None', style=style)
        c = widgets.HBox([narmin, narmax])
        naring_intro=widgets.HTML("""Specify the range of the number of non aromatic rings present in the  molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        naring=widgets.VBox(children=[naring_intro,space_box,c])
        
        smin = widgets.Text(description='Minimum no of Single Bonds', value='None', style=style)
        smax = widgets.Text(description='Maximum no of Single Bonds', value='None', style=style)
        d = widgets.HBox([smin, smax])
        smin_intro=widgets.HTML("""Specify the range of the number of single bonds present in the  molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        sbond=widgets.VBox(children=[smin_intro,space_box,d])

        dmin = widgets.Text(description='Minimum no of Double Bonds', value='None', style=style)
        dmax = widgets.Text(description='Maximum no of Double Bonds', value='None', style=style)
        e = widgets.HBox([dmin, dmax])
        dmin_intro=widgets.HTML("""Specify the range of the number of double bonds present in the  molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        dbond=widgets.VBox(children=[dmin_intro,space_box,e])

        tmin = widgets.Text(description='Minimum no of Tiple Bonds', value='None', style=style)
        tmax = widgets.Text(description='Maximum no of Triple Bonds', value='None', style=style)
        f = widgets.HBox([tmin, tmax])
        tmin_intro=widgets.HTML("""Specify the range of the number of triple bonds present in the  molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        tbond=widgets.VBox(children=[tmin_intro,space_box,f])

        specific_atoms = widgets.Text(
            value='None',
            placeholder='Type the number of specific atoms',
            description='Maximum No. of Specific Atoms:',
            style=style)
        satom_intro=widgets.HTML("""Specify the specific atoms that must be present in the  molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        satom=widgets.VBox(children=[satom_intro,space_box,specific_atoms])

        lipinski_rule = widgets.RadioButtons(
            options=['True', 'False'],
            value='False',
            description='Lipinski rule:', style=style,
            disabled=False)
        lr_intro=widgets.HTML("""Lipinski rule to be considered or not for the generation of molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        lr=widgets.VBox(children=[lr_intro,space_box,lipinski_rule])
        
        
        fingerprint_matching = widgets.Text(
            value='None',
            placeholder='Finger print match',
            description='Finger print matching:',
            style=style)
        fm_intro=widgets.HTML("""Fingerprint matching to be considered or not for the generation of molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        fm=widgets.VBox(children=[fm_intro,space_box,fingerprint_matching])

        substructure = widgets.Text(
            value='None',
            placeholder='Type the substructure SMILES',
            description='Substructure:',
            style=style)
        sub_intro=widgets.HTML("""Define the substructure which must be included in the molecule during generation""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        sub=widgets.VBox(children=[sub_intro,space_box,substructure])

        substructure_exclusion = widgets.Text(
            value='None',
            placeholder='Type the substructure exclusion SMILES',
            description='Substructure exclusion:',
            style=style)
        sube_intro=widgets.HTML("""Define the substructure which must be excluded in the molecule during generation""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        sube=widgets.VBox(children=[sube_intro,space_box,substructure_exclusion])

        Include_BB = widgets.RadioButtons(
            options=['Yes', 'No'],
            value='No',
            description='Include BB:', style=style,
            disabled=False
        )
        ibb_intro=widgets.HTML("""Define the building blocks to be included in all the molecules""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        ibb=widgets.VBox(children=[ibb_intro,space_box,Include_BB])


        """Symmetry = widgets.RadioButtons(
            options=['Yes', 'No'],
            value='No',
            description='Symmetry:', style=style,
            disabled=False
        )"""
        accordion = widgets.Accordion(
            children=[bb, no_bonds, atom, mol, ring, aring, naring, sbond, dbond, tbond, satom, lr, fm,
                      sub, sube, ibb])
        accordion.set_title(0, 'Building Blocks')
        accordion.set_title(1, 'No of Bonds')
        accordion.set_title(2, 'No of Atoms')
        accordion.set_title(3, 'Molecular Weight')
        accordion.set_title(4, 'No of Rings')
        accordion.set_title(5, 'No of Aromatic Rings')
        accordion.set_title(6, 'No of Non Aromatic Rings')
        accordion.set_title(7, 'No of Single Bonds')
        accordion.set_title(8, 'No of Double Bonds')
        accordion.set_title(9, 'No of Triple Bonds')
        accordion.set_title(10, 'Specific Atoms')
        accordion.set_title(11, 'Lipinski Rule')
        accordion.set_title(12, 'Fingerprint Matching')
        accordion.set_title(13, 'Substructure')
        accordion.set_title(14, 'Substructure Exclusion')
        accordion.set_title(15, 'Include BB')
        #accordion.set_title(16, 'Symmetry')
        display(accordion)
        def generation_file():
            if bondmin.value=="":
                bondmin.value='None'
            elif building_blocks.value=="":
                building_blocks.value='None'
            elif bondmax.value=="":
                bodmax.value='None'
            elif atommin.value=="":
                atommin.value='None'
            elif atommax.value=="":
                atommax.vlaue='None'
            elif molmin.value=="":
                molmin.value='None'
            elif molmax.value=="":
                molmax.value='None'
            elif ringmin.value=="":
                ringmin.value='None'
            elif armin.value=="":
                armin.value='None'
            elif armax.value=="":
                armax.value='None'
            elif narmin.value=="":
                narmin.value='None'
            elif narmax.value=="":
                narmax.value='None'
            elif smin.value=="":
                smin.value='None'
            elif smax.value=="":
                smax.value='None'
            elif dmin.value=="":
                dmin.value='None'
            elif dmax.value=="":
                dmax.value='None'
            elif tmin.value=="":
                tmin.value='None'
            elif tmin.value=="":
                tmin.value='None'
            elif tmax.value=="":
                tmax.value='None'
            elif specific_atoms.value=="":
                specific_atoms.value='None'
            elif fingerprint_matching.value=="":
                fingerprint_matching.value='None'
            elif susbtructure.value=="":
                substructure.value='None'
            elif substructure_exclusion.value=="":
                substructure_exclusion.value=='None'
                
                
                
            generation = open("generation_rules.dat", "w+")
            generation.write(
                "Please input generation rules below. Do not change the order of the options" + '\n' + "1. Include building blocks ==" + building_blocks.value + '\n' +
                "2. Min and max no. of bonds == " + bondmin.value + "," + bondmax.value + '\n' + "3. Min and max no. of atoms == " + atommin.value + "," + atommax.value + '\n' +
                "4. Min and max mol. weight == " + molmin.value + "," + molmax.value + '\n' + "5. Min and max no. of rings == " + ringmin.value + "," + ringmax.value + '\n' +
                "6. Min and max no. of aromatic rings == " + armin.value + "," + armax.value + '\n' + "7. Min and max no. of non aromatic rings ==" + narmin.value + "," + narmax.value + '\n' +
                "8. Min and max no. of single bonds == " + smin.value + "," + smax.value + '\n' + "9. Min and max no. of double bonds == " + dmin.value + "," + dmax.value + '\n' +
                "10. Min and max no. of triple bonds == " + tmin.value + "," + tmax.value + '\n' + "11. Max no. of specific atoms == " + specific_atoms.value + '\n' +
                "12. Lipinski's rule == " + lipinski_rule.value + '\n' + "13. Fingerprint matching ('c1ccccc1'-0.1), ('C1CCCC1'-0.1) == " + fingerprint_matching.value + '\n'
                                                                                                                                                                         "14. Substructure == " + substructure.value + '\n' + "15. Substructure exclusion == " + substructure_exclusion.value + '\n' +
                "15. Include_BB == " + Include_BB.value + '\n')
            generation.close()

        button1 = widgets.Button(description="Generate rules file",layout=widgets.Layout(border='solid 1px black'),style=style)
        display(button1)

        def on_button_clicked(b):
            generation_file()
            
            
        display(third)    
        button1.on_click(on_button_clicked)
    second.on_click(second_section)
    
    def third_section(v):
        style = {'description_width': 'initial','font_weight':'bold'}
        display(widgets.HTML(value="<font color=crimson><font size=5><b><u>SECTION 3: GENERATION OF COMMAND LINE</font>"))

        input_file = widgets.Text(value='building_blocks.dat',
                              placeholder='Type the name of input file',
                              description='Input file name:',
                              style=style,
                              )
        input_intro=widgets.HTML("""Type the name of the input file""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        i_f=widgets.VBox(children=[input_intro,space_box,input_file])

        molecule_type = widgets.Dropdown(
            options=['SMILES', 'SMARTS', 'INCHI'],
            value='SMILES',
            description='Molecule Type:',
            style=style,
            disabled=False)
        mt_intro=widgets.HTML("""Enter the type of representation of the molecule""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        mt=widgets.VBox(children=[mt_intro,space_box,molecule_type])
        
        combination_type = widgets.RadioButtons(
            options=['Fuse', 'Link'],
            value='Link',
            description='Combination type:',
            style=style,
            disabled=False
        )
        comb_intro=widgets.HTML("""Enter the combination type required for generation of molecule""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        comb=widgets.VBox(children=[comb_intro,space_box,combination_type])
        
        generation_level = widgets.BoundedIntText(
            value=1,
            min=1,
            max=100,
            step=1,
            description='Generation level:',
            style=style,
            disabled=False)
        gl_intro=widgets.HTML("""Enter the genration level required""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        gl=widgets.VBox(children=[gl_intro,space_box,generation_level])
        

        output_type = widgets.Dropdown(
            options=['SMILES', 'SMARTS', 'INCHI'],
            value='SMILES',
            description='Output Type:',
            style=style,
            disabled=False, )
        op_intro=widgets.HTML("""Type the name of the input file""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        op=widgets.VBox(children=[op_intro,space_box,output_type])
        

        max_files = widgets.IntText(
            value=10000,
            description='Maximum files per folder:',
            style=style,
            disabled=False)
        mf_intro=widgets.HTML("""Enter the maximum number of files required in one folder""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        mf=widgets.VBox(children=[mf_intro,space_box,max_files])
        

        library_name = widgets.Text(
            value='new_library_',
            placeholder='Type the nameof the library',
            description='Library Name:',
            style=style)
        ln_intro=widgets.HTML("""Type the name of the library""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        ln=widgets.VBox(children=[ln_intro,space_box,library_name])
        

        style = {'description_width': 'initial', }
        chemHTPS = widgets.RadioButtons(
            options=['Yes', 'No'],
            value='No',
            description='Run with ChemHTPS:', style=style,
            disabled=False
        )
        ch_intro=widgets.HTML("""Enter whether ChemHTPS shall be included during the generation""",
        layout = widgets.Layout(height = '45px', width = '90%',
                    size = '20'))
        cH=widgets.VBox(children=[ch_intro,space_box,chemHTPS])
        
        children = [i_f, mt, comb, gl, op, mf, ln,
                    cH]
        tab = widgets.Tab()
        tab.children = children
        tab.set_title(0, 'Input File')
        tab.set_title(1, 'Molecule Type')
        tab.set_title(2, 'Combination Type')
        tab.set_title(3, 'Generation Level')
        tab.set_title(4, 'Output File')
        tab.set_title(5, 'Maximum No of Files')
        tab.set_title(6, 'Library Name')
        tab.set_title(7, 'ChemHTPS')

        display(tab)
        button = widgets.Button(description="Generate Command line", layout= Layout(width= 'auto',border='solid 1px black'),style=style)
        display(button)

        def on_button_clicked(b):
            opt = "--input_file %s --molecule_type %s --combination_type %s --generation_levels %s --output_type %s--max_files_per_folder%s --rule_file generation.dat --lib_name %s --ChemHTPS %s" % (
            input_file.value, molecule_type.value, combination_type.value, generation_level.value, output_type.value,
            max_files.value, library_name.value, chemHTPS.value)
            print(opt)

        button.on_click(on_button_clicked)
    third.on_click(third_section)
