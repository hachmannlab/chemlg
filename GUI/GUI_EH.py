
from __future__ import print_function
from ipywidgets import interact, Layout
import ipywidgets as widgets
from IPython.display import display
import sys
sys.path.insert(0, "/user/m27/pkg/openbabel/2.3.2/lib")
import pybel

style = {'description_width': 'initial','font_weight':'bold'}
count =0

def GUI():
    style = {'description_width': 'initial','font_weight':'bold'}
    space_box = widgets.Box(layout=widgets.Layout(height ='55px', width='90%')) 
    display(widgets.HTML(value="<font color=crimson><font size=5><b><u>SECTION 1: BUILDING BLOCKS</font>"))
    second=widgets.Button(description='Proceed to the next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    
    third=widgets.Button(description='Proceed to the next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    
         

    name=widgets.Text(description="File name",placeholder="Type the file to be imported",style=style)
     
    enter=widgets.Button(description='Upload file',layout= Layout(width='30%',border='solid 1px black'),style=style)
    BB=widgets.Text(description='Type building blocks',style=style) 
    another=widgets. Button(description='Add building block', layout= Layout(width= '45%',border='solid 1px black'),style=style)
    existing_file_intro = widgets.HTML("""<font size=3>This the Building blocks section wherein the user can enter the building blocks by uploading an existing building blocks file or create a building blocks file by entering the SMILES of each building block""",
        layout = widgets.Layout(height = '60px', width = '90%',
                    size = '20'))
    existing_filebox = widgets.VBox(children=[existing_file_intro,space_box,name,space_box,enter]
                                )
    type_smiles=widgets.Button(description='Type smiles',layout= Layout(width= '45%',border='solid 1px black'),style=style)
    type_smiles_intro = widgets.HTML("""<font size=3>The user can type in single SMILES of a molecules or can enter multiple SMILES separated by a <b> ',' </b> in one attempt and click on the add building blocks button and after all the building blocks have been added the user can create the building blocks file by clicking and the create building block button and proceed to the next section""",
        layout = widgets.Layout(height = '85px', width = '90%',
                    size = '20'))
    def on_proceed(e):
        o=name.value
        try:
            existing= open(o,"r")
            c=existing.readlines()
            for line in c:
                try:
                    mol=pybel.readstring("smi",line)
                except:
                    incorrect= widgets.HTML(line +" in file is a incorrect SMILE")
                    display(incorrect)
        except:
            print("The Building Blocks file does not exist")
             
                    
        
        
    enter.on_click(on_proceed)    
  
    BB_list=[]
    
   
        
        
    def on_another_clicked(i):
        z=BB.value
        z=str(z)
        f=z.split(',')
        for element in f:
            try:
                mol=pybel.readstring("smi",element)
                BB_list.append(element)
            except:
                
                print("You typed a wrong SMILES "+element)
            
        BB.value=""
        
   
    final_BB=widgets.Button( description='Create Building Blocks file',layout= Layout(width= '45%',border='solid 1px black'),style=style)
    
    def on_button_clicked(d):
        building_blocks= open("building_blocks.txt", "w+")
    
        building_blocks.write('Building blocks are:'+'\n')
        for h in range(len(BB_list)):
            building_blocks.write(BB_list[h] +'\n')
        building_blocks.close()
    
    final_BB.on_click(on_button_clicked)
    another.on_click(on_another_clicked)
    type_smilesbox = widgets.VBox(children=[type_smiles_intro,space_box,BB,space_box,another,final_BB
                                ])
    tab1 = widgets.Tab(
            children=[existing_filebox, type_smilesbox],style=style)
    tab1.set_title(0, 'Upload file')
    tab1.set_title(1, 'Type smiles')
    display (tab1)
    display (second)    
    def second_section(q):
        display(widgets.HTML(value="<font color=crimson><font size=5><b><u>SECTION 2: GENERATION RULES</font>"))

        style = {'description_width': 'initial','font_weight':'bold'}
        building_blocks = widgets.Text(value='F1',
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
            try:
                mol=pybel.readstring("smi",building_blocks.value)
            except:
                print("You entered a wrong SMILE for the building block")
            
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



