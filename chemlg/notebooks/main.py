import rdkit
#from __future__ import print_function
from ipywidgets import interact, Layout
import ipywidgets as widgets
from IPython.display import display
import os
import pybel
import sys
from IPython.display import Javascript
from ipywidgets import interact, Layout
import ipywidgets as widgets
from IPython.display import display, clear_output, Javascript
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import os

style = {'description_width': 'initial','font_weight':'bold'}
count =0
BB_list=[]
smiles=[]


# run this function to start the GUI builder.
def config_builder():
        display(widgets.HTML(value="<font color=crimson><font size=5><b><u>CHEMLG</font>"))
        """This function allows the user to provide the GA function required for library generation if the user wants to provide constraints in building the library"""
        def constraints(x):
            GA_func=widgets.Textarea(description="GA function",placeholder="Type the GA function to be incorporated for library selection", layout=Layout(width='80%', height='300px'))
            if x== 'with constraints':
                display(GA_func)
                zero=widgets.Button(description="Next")
                display(zero)
                def on_button_click(a):
                    if GA_func.value=="":
                        print("GA function empty")
                    else:
                        gen_func = open("GA_script.txt", "w+")
                        gen_func.write(GA_func.value)
                        gen_func.close()
                        building_blocks()
                zero.on_click(on_button_click)
            else:
                building_blocks()


        interact(constraints, x=widgets.RadioButtons(
                        options=['with constraints', 'without constraints'],
                        value='without constraints',
                        description='Run ChemLG ',
                        disabled=False))

def building_blocks():
    """This function allows the user to provide the building blocks required for library generation"""
    style = {"description_width": "initial","font_weight":"bold"}
    space_box = widgets.Box(layout=widgets.Layout(height ='20px', width='90%')) 
    display(widgets.HTML(value="<font color=crimson><font size=5><b><u>BUILDING BLOCKS</font>"))
    second=widgets.Button(description='Next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    
    third=widgets.Button(description='Next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    


    name=widgets.Text(description="File name",placeholder="Type the file to be imported",style=style)

    enter=widgets.Button(description='Upload file',layout= Layout(width='30%',border='solid 1px black'),style=style)
    BB=widgets.Text(description='Building Blocks',style=style) 
    another=widgets.Button(description='Add building block', layout= Layout(width= '45%',border='solid 1px black'),style=style)
    existing_file_intro = widgets.HTML("""<font size=3>Upload a file containing all the building blocks""",
        layout = widgets.Layout(height = '60px', width = '90%',
                        size = '20'))
    existing_filebox = widgets.VBox(children=[existing_file_intro,space_box,name,space_box,enter] )
    type_smiles=widgets.Button(description='Individual SMILES input',layout= Layout(width= '120%',border='solid 1px black'),style=style)
    type_smiles_intro = widgets.HTML("""<font size=3>Create a building blocks file by entering individual/comma separated SMILES of each building block.""",
            layout = widgets.Layout(height = '55px', width = '90%',
                        size = '20'))

    def on_proceed(e):
        file_name=name.value


        try:
            existing= open(file_name,"r")
            if os.stat(file_name).st_size == 0:
                print("Building blocks file is empty")
            else:
                lines=existing.readlines()[1:]

                for line in lines:

                    smiles.append(line)

                    try:
                        mol=pybel.readstring("smi",line)
                        correct=widgets.HTML("""<font size=3>The file has been uploaded successfully""",
                layout = widgets.Layout(height = '85px', width = '90%',
                            size = '20'))

                        #|display(Javascript('IPython.notebook.execute_cells([2])'))

                    except:
                        incorrect= widgets.HTML(line +" in file is a incorrect SMILES")
                        display(incorrect)
                if len(smiles)>0:
                            mol_lists = [Chem.MolFromSmiles(smile) for smile in smiles]
                            [mol.SetProp('_Name','B'+str(i)) for i,mol in enumerate(mol_lists)]
                            ibu1=Chem.Draw.MolsToGridImage(mol_lists)
                            display(ibu1)
                            generation_rules()
        except:
            print("The building blocks file does not exist")
        #clear_output()



    enter.on_click(on_proceed)    

    BB_list=[]




    def on_another_clicked(i):
        global f
        clear_output()

        for element in f:

            try:
                mol=pybel.readstring("smi",element)
                BB_list.append(element)

            except:

                print("You typed a wrong SMILES "+element)

        BB.value=""

        display(widgets.HTML(value="<font color=crimson><font size=5><b><u>BUILDING BLOCKS</font>"))
        display(tab1)
    final_BB=widgets.Button( description='Create building blocks file',layout= Layout(width= '45%',border='solid 1px black'),style=style)
    
    def on_button_clicked(d):


        #display(Javascript('IPython.notebook.execute_cells([1])'))
        building_blocks= open("building_blocks.dat", "w")

        
        for h in range(len(BB_list)):

            building_blocks.write(BB_list[h] +'\n')

        building_blocks.close()
        generation_rules()

    def on_visualization_clicked(t):
        #This function visualizes the smiles provided by the user
        global f
        z=BB.value
        z=str(z)
        f=z.split(',')
        try:
            mol_list = [Chem.MolFromSmiles(element) for element in f]
            [mol.SetProp('_Name','B'+str(i)) for i,mol in enumerate(mol_list)]
            ibu=Chem.Draw.MolsToGridImage(mol_list)
            #,legends=[mol.GetProp('_Name') for mol in mol_list]
            display(ibu)
        except:
            print("You typed a wrong SMILES. Unable to visualize")

        #display(Javascript('IPython.notebook.execute_cells([1])'))

    visualize=widgets.Button(description='Visualize',layout= Layout(width= '45%',border='solid 1px black'),style=style)


    visualize.on_click(on_visualization_clicked)    
    final_BB.on_click(on_button_clicked)
    another.on_click(on_another_clicked)
    type_smilesbox = widgets.VBox(children=[type_smiles_intro,BB,space_box,visualize,another,final_BB
                                    ])
    tab1 = widgets.Tab(
            children=[type_smilesbox,existing_filebox],style=style)
    tab1.set_title(0, ' Individual SMILES')
    tab1.set_title(1, 'Upload file')
    display (tab1)
    

def generation_rules():
    """This function is used to take input for the generation rules from user."""
    style = {'description_width': 'initial','font_weight':'bold'}
    space_box = widgets.Box(layout=widgets.Layout(height ='55px', width='90%')) 
    second=widgets.Button(description='Next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    
    third=widgets.Button(description='Next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    
    display (second)      
    def second_section(q):
        
        if os.stat("building_blocks.dat").st_size <= 0:
            print("building blocks file is empty")
            
            
        else:
            display(widgets.HTML(value="<font color=crimson><font size=5><b><u>GENERATION RULES</font>"))

            style = {'description_width': 'initial','font_weight':'bold'}
            building_blocks = widgets.Text(value='None',
                                       placeholder='Enter building blocks',
                                       description=' Specific Building Blocks:',
                                       style=style,
                                       )    
            bb_intro=widgets.HTML("""Specify the building blocks which must be present in all the molecules in the final library""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            bb=widgets.VBox(children=[bb_intro,space_box,building_blocks])

            style = {'description_width': 'initial'}
            bondmin = widgets.Text(description='Minimum', value='None', style=style)
            bondmax = widgets.Text(description='Maximum', value='None', style=style)
            bonds_intro=widgets.HTML("""Enter the minimum and maximum values of the total number of bonds for all the molecules in the final library (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            x = widgets.HBox([bondmin, bondmax])
            no_bonds=widgets.VBox(children=[bonds_intro,space_box,x
                                    ])
           

            atommin = widgets.Text(description='Minimum', value='None', style=style)
            atommax = widgets.Text(description='Maximum', value='None', style=style)
            y = widgets.HBox([atommin, atommax])
            atom_intro=widgets.HTML("""Specify the minimum and maximum number of atoms that must be present in each molecule in the generation library  (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            atom=widgets.VBox(children=[atom_intro,space_box,y])
           

            molmin = widgets.Text(description='Minimum', value='None', style=style)
            molmax = widgets.Text(description='Maximum', value='None', style=style)
            z = widgets.HBox([molmin, molmax])
            mol_intro=widgets.HTML("""Specify the range of the molecular weight of the molecules in the generation library  (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            mol=widgets.VBox(children=[mol_intro,space_box,z])
           

            ringmin = widgets.Text(description='Minimum', value='None', style=style)
            ringmax = widgets.Text(description='Maximum', value='None', style=style)
            a = widgets.HBox([ringmin, ringmax])
            ring_intro=widgets.HTML("""Specify the range of the number of rings present in the  molecules in the generation library  (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            ring=widgets.VBox(children=[ring_intro,space_box,a])
           
            armin = widgets.Text(description='Minimum', value='None', style=style)
            armax = widgets.Text(description='Maximum', value='None', style=style)
            b = widgets.HBox([armin, armax])
            aring_intro=widgets.HTML("""Specify the range of the number of aromatic rings present in the  molecules in the generation library  (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            aring=widgets.VBox(children=[aring_intro,space_box,b])
            

            narmin = widgets.Text(description='Minimum', value='None', style=style)
            narmax = widgets.Text(description='Maximum', value='None', style=style)
            c = widgets.HBox([narmin, narmax])
            naring_intro=widgets.HTML("""Specify the range of the number of non aromatic rings present in the  molecules in the generation library  (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            naring=widgets.VBox(children=[naring_intro,space_box,c])
            

            smin = widgets.Text(description='Minimum', value='None', style=style)
            smax = widgets.Text(description='Maximum', value='None', style=style)
            d = widgets.HBox([smin, smax])
            smin_intro=widgets.HTML("""Specify the range of the number of single bonds present in the  molecules in the generation library  (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            sbond=widgets.VBox(children=[smin_intro,space_box,d])
            

            dmin = widgets.Text(description='Minimum', value='None', style=style)
            dmax = widgets.Text(description='Maximum', value='None', style=style)
            e = widgets.HBox([dmin, dmax])
            dmin_intro=widgets.HTML("""Specify the range of the number of double bonds present in the  molecules in the generation library  (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            dbond=widgets.VBox(children=[dmin_intro,space_box,e])
           

            tmin = widgets.Text(description='Minimum', value='None', style=style)
            tmax = widgets.Text(description='Maximum', value='None', style=style)
            f = widgets.HBox([tmin, tmax])
            tmin_intro=widgets.HTML("""Specify the range of the number of triple bonds present in the  molecules in the generation library  (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            tbond=widgets.VBox(children=[tmin_intro,space_box,f])
            

            specific_atoms = widgets.Text(
                value='None',
                description='Element:',
                style=style)
            
            

            satom_intro=widgets.HTML("""Specify the maximum number of heteroatoms that must be present in the  molecules in the final library.Expected type: tuple of tuple(s) (('Cl', 10), )""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            spec=widgets.HBox(children=[specific_atoms])
            satom=widgets.VBox(children=[satom_intro,space_box,spec])

            lipinski_rule = widgets.RadioButtons(
                options=['True', 'False'],
                value='False',
                description='Lipinski rule:', style=style,
                disabled=False)
            lr_intro=widgets.HTML("""Lipinski's rule of five is a rule of thumb to evaluate druglikeness or determine if a chemical compound with a certain pharmacological or biological activity has chemical properties and physical properties that would make it a likely orally active drug in humans. Select the choice to incorporate Linpinski rule for generating molecules in the final library""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            lr=widgets.VBox(children=[lr_intro,space_box,lipinski_rule])


            fingerprint_matching = widgets.Text(
                value='None',
                placeholder='Finger print match',
                description='Finger print matching:',
                style=style)
            fm_intro=widgets.HTML("""Molecular fingerprints encode molecular structure in a series of binary digits that represent the presence or absence of particular substructures in the molecule. Comparing fingerprints will allow you to determine the similarity between two molecules. Type the target molecule and the Tanimoto index.example: ('c1ccccc1'-0.1), ('C1CCCC1'-0.1)""",
            layout = widgets.Layout(height = '60px', width = '90%',
                        size = '20'))
            fm=widgets.VBox(children=[fm_intro,space_box,fingerprint_matching])

            substructure = widgets.Text(
                value='None',
                placeholder='Type the substructure SMILES',
                description='Substructure:',
                style=style)
            sub_intro=widgets.HTML("""Enter substructures in SMARTS format which must be included in all molecules in the final library
            """,
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            sub=widgets.VBox(children=[sub_intro,space_box,substructure])

            substructure_exclusion = widgets.Text(
                value='None',
                placeholder='Type the substructure exclusion SMILES',
                description='Substructure exclusion:',
                style=style)
            sube_intro=widgets.HTML("""Enter substructures in SMARTS format which must be excluded in all molecules in the final library
            """,
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            sube=widgets.VBox(children=[sube_intro,space_box,substructure_exclusion])

            Include_BB = widgets.RadioButtons(
                options=['True', 'False'],
                value='True',
                description='Include BB:', style=style,
                disabled=False
            )
            ibb_intro=widgets.HTML("""Should the initial building blocks be included in the final library""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            ibb=widgets.VBox(children=[ibb_intro,space_box,Include_BB])


            """Symmetry = widgets.RadioButtons(
                options=['Yes', 'No'],
                value='No',
                description='Symmetry:', style=style,
                disabled=False
            )"""
            input_file = widgets.Text(value='building_blocks.dat',
                                  placeholder='Type the name of input file',
                                  description='Enter building blocks file name:',
                                  style=style,layout = widgets.Layout(height = '45px', width = '40%',
                        size = '20')
                                  )
            input_intro=widgets.HTML("""Type the name of the input file""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            i_f=widgets.VBox(children=[input_file])

            """molecule_type = widgets.Dropdown(
                options=['SMILES', 'SMARTS', 'INCHI'],
                value='SMILES',
                description='Molecule Type:',
                style=style,
                disabled=False)
            mt_intro=widgets.HTML('Enter the type of representation of the molecule',
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            mt=widgets.VBox(children=[mt_intro,space_box,molecule_type])"""

            combination_type = widgets.RadioButtons(
                options=['Fusion', 'Link'],
                value='Link',
                description='Combination type:',
                style=style,
                disabled=False
            )
            comb_intro=widgets.HTML("""Select the combination type for generating the molecules""",
            layout = widgets.Layout(height = '50px', width = '90%',
                        size = '20'))
            comb=widgets.VBox(children=[comb_intro,combination_type])

            generation_level = widgets.BoundedIntText(
                value=1,
                min=1,
                max=100,
                step=1,
                description='Generation level:',
                style=style,
                disabled=False)
            gl_intro=widgets.HTML("""Enter the number of generations to run the library generator""",
            layout = widgets.Layout(height = '50px', width = '90%',
                        size = '20'))
            gl=widgets.VBox(children=[gl_intro,generation_level])


            output_type = widgets.Dropdown(
                options=['smi', 'xyz',],
                value='smi',
                description='Output Type:',
                style=style,
                disabled=False, )
            op_intro=widgets.HTML("""Specify the output file format""",
            layout = widgets.Layout(height = '50px', width = '90%',
                        size = '20'))
            op=widgets.VBox(children=[op_intro,output_type])


            max_files = widgets.IntText(
                value=10000,
                description='Maximum files per folder:',
                style=style,
                disabled=False)
            """mf_intro=widgets.HTML('Enter the maximum number of files required in one folder',
            layout = widgets.Layout(height = '45px', width = 'auto',
                        size = '20'))"""
            mf=widgets.VBox(children=[max_files])


            library_name = widgets.Text(
                value='new_library_',
                placeholder='Type the nameof the library',
                description='Library Name:',
                style=style)
            """ln_intro=widgets.HTML('Type the name of the library',
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))"""
            ln=widgets.VBox(children=[library_name])

            accordion = widgets.Accordion(
                children=[bb, no_bonds, atom, mol, ring, aring, naring, sbond, dbond, tbond, satom, lr, fm,
                          sub, sube, ibb, i_f, comb, gl, op, mf, ln])
            accordion.set_title(0, 'Include Building Blocks')
            accordion.set_title(1, 'Number of Bonds')
            accordion.set_title(2, 'Number of Atoms')
            accordion.set_title(3, 'Molecular Weight Range')
            accordion.set_title(4, 'Number of Rings')
            accordion.set_title(5, 'Number of Aromatic Rings')
            accordion.set_title(6, 'Number of Non-Aromatic Rings')
            accordion.set_title(7, 'Number of Single Bonds')
            accordion.set_title(8, 'Number of Double Bonds')
            accordion.set_title(9, 'Number of Triple Bonds')
            accordion.set_title(10, 'Heteroatoms')
            accordion.set_title(11, 'Lipinski Rule')
            accordion.set_title(12, 'Fingerprint Matching')
            accordion.set_title(13, 'Substructure Inclusion')
            accordion.set_title(14, 'Substructure Exclusion')
            accordion.set_title(15, 'Include initial Building Blocks')
            accordion.set_title(16, 'Building Blocks')
            accordion.set_title(17, 'Combination Type')
            accordion.set_title(18, 'No of generations')
            accordion.set_title(19, 'Output File Format')
            accordion.set_title(20, 'Maximum No of Files')
            accordion.set_title(21, 'Library Name')


            #accordion.set_title(16, 'Symmetry')
            display(accordion)
            def generation_file():
                #This function generates the config.dat file based on the user input for the generation rules

                display(widgets.HTML(value="""<font size=3>Configuration file created""",layout = widgets.Layout(height = '60px', width = '90%',size = '20')))
                

                generation = open("config.dat", "w+")
                generation.write("Please input generation rules below. Do not change the order of the options" + '\n' + "1. Include building blocks == " + building_blocks.value + '\n')
                if bondmin.value==bondmax.value=="None":
                    generation.write("2. Min and max no. of bonds == None"+'\n')
                else:
                    tbonds=(int(bondmin.value),int(bondmax.value))
                    generation.write("2. Min and max no. of bonds =="+str(tbonds)+'\n')
                if atommin.value==atommax.value=="None":
                    generation.write("3. Min and max no. of atoms == None"+ '\n')
                else:
                    tatoms=(int(atommin.value),int(atommax.value))
                    generation.write("3. Min and max no. of atoms == "+ str(tatoms)+ '\n')
                if molmin.value==molmax.value=="None":
                    generation.write("4. Min and max mol. weight == None" + '\n')
                else:
                    tmol=(int(molmin.value),int(molmax.value))
                    generation.write("4. Min and max mol. weight == "+str(tmol) + '\n')
                if ringmin.value==ringmax.value=="None":
                    generation.write("5. Min and max no. of rings == None"+ '\n')
                else:
                    tring=(int(ringmin.value),int(ringmax.value))
                    generation.write("5. Min and max no. of rings == "+str(tring)+ '\n')
                if armin.value==armax.value=="None":
                    generation.write("6. Min and max no. of aromatic rings == None" + '\n')
                else:
                    taro=(int(armin.value),int(armax.value))
                    generation.write("6. Min and max no. of aromatic rings == " +str(taro)+ '\n')
                if narmin.value==narmax.value=="None":
                    generation.write("7. Min and max no. of non aromatic rings == None" + '\n')
                else:
                    tnar=(int(narmin.value),int(narmax.value))
                    generation.write("7. Min and max no. of non aromatic rings == "+str(tnar) + '\n')
                if smin.value==smax.value=="None":
                    generation.write("8. Min and max no. of single bonds == None" + '\n')
                else:
                    tsbond=(int(smin.value),int(smax.value))
                    generation.write("8. Min and max no. of single bonds == "+str(tsbond) + '\n')
                if dmin.value==dmax.value=="None":
                    generation.write("9. Min and max no. of double bonds == None" + '\n')
                else:
                    tdbond=(int(dmin.value),int(dmax.value))
                    generation.write("9. Min and max no. of double bonds == "+str(tdbond) + '\n')
                if tmin.value==tmax.value=="None":
                    generation.write("10. Min and max no. of triple bonds == None"+'\n')
                else:
                    ttbond=(int(tmin.value),int(tmax.value))
                    generation.write("10. Min and max no. of triple bonds == "+str(ttbond)+'\n') 
                generation.write("11. Max no. of specific atoms == " + specific_atoms.value + '\n' +
                    "12. Lipinski's rule == " + lipinski_rule.value + '\n' + "13. Fingerprint matching == " + fingerprint_matching.value + '\n'+"14. Substructure == " + substructure.value + '\n' + "15. Substructure exclusion == " + substructure_exclusion.value + '\n' +
                    "15. Include_BB == " + Include_BB.value + '\n'+'\n'+'\n'+'\n'+'\n'+'\n'+
                    #"Building blocks file :: "+input_file.value+'\n'+
                    "Combination type for molecules :: "+combination_type.value+'\n'+
                    "Number of generations :: "+str(generation_level.value)+'\n'+
                    "Molecule format in output file ::"+output_type.value+'\n'+
                    "Maximum files per folder :: "+str(max_files.value)+'\n'+"Library name :: "+library_name.value+'\n' )

                generation.close()

            button1 = widgets.Button(description="Generate configuration file",layout=widgets.Layout(width='20%',border='solid 1px black'),style=style)

            display(button1)

            def on_button_clicked(b):
                #This function takes the environment name to run chemlg as the input from the user
                clear_output()
                generation_file()
                global envi
                envi_intro=widgets.HTML("""Specify the virual environment in which you want to run the library generator""",
            layout = widgets.Layout(height = '50px', width = '90%',
                        size = '20'))
                envi=widgets.Text(
                value='None',
                placeholder='Enter Environment name',
                description='Environment name',
                style=style)
                environment=widgets.VBox(children=[envi_intro,envi])
                display(environment)
                run=widgets.Button(description='Run Chemlg',layout= Layout(width= 'auto',border='solid 1px black'),style=style)
                display(run)
                run.on_click(on_run_clicked)
            def on_run_clicked(c):
                if envi.value=="":
                    envi.value=='None'
                    os.system("chemlgshell -i config.dat -b building_blocks.dat -o ./")
                else:
                    try:
                        os.system("activate "+envi.value)
                    except:
                        os.system("source activate"+envi.value)
                    os.system("chemlgshell -i config.dat -b building_blocks.dat -o ./") 




            button1.on_click(on_button_clicked)
    second.on_click(second_section)
config_builder()
