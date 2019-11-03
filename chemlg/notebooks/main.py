from ipywidgets import interact, Layout, Label
import ipywidgets as widgets
import os
import pybel
import sys
from IPython.display import display, clear_output, Javascript
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
import os.path
import csv
import pandas as pd
from collections import Counter

style = {'description_width': 'initial','font_weight':'bold'}
count =0
BB_list=[]
smiles=[]


# run this function to start the GUI builder.
def config_builder():
    
        display(widgets.HTML(value="<font color=blue><font size=5><b><u>CHEMLG</font>"))
        """This function allows the user to provide the GA function required for library generation if the user wants to provide constraints in building the library"""
        global input_direc
        input_direc=widgets.Text(description="Input file directory",placeholder="",style=style)
        input_direc_intro=widgets.HTML(description="Provide the directory to save all the input files",style=style)
        input_direc_button=widgets.Button(description='Continue',layout= Layout(width= '20%',border='solid 1px black'),style=style)
        inp=widgets.VBox(children=[input_direc_intro,input_direc,input_direc_button])
        display(inp)

        def on_input_clicked(o):
            display(widgets.HTML(value="<font color=blue><font size=5><b><u>BUILDING BLOCKS</font>"))
            building_blocks()
        input_direc_button.on_click(on_input_clicked)

        

def building_blocks():
    """This function allows the user to provide the building blocks required for library generation"""
    style = {"description_width": "initial","font_weight":"bold"}
    space_box = widgets.Box(layout=widgets.Layout(height ='20px', width='90%')) 
    
        
#Section for providing SMILES input

    name=widgets.Text(description="File name",placeholder="Type the file to be imported",style=style)
    BB=widgets.Text(description='Building Blocks',style=style)  #Textbox for SMILES input
    type_smiles=widgets.Button(description='Individual SMILES input',layout= Layout(width= '120%',border='solid 1px black'),style=style)
    type_smiles_intro = widgets.HTML("""<font size=3>Create a building blocks file by entering individual/comma separated SMILES of each building block.""",
            layout = widgets.Layout(height = '55px', width = '90%',
                        size = '20'))
    
    another=widgets.Button(description='Add building block', layout= Layout(width= '45%',border='solid 1px black'),style=style) #Button for adding another building block
    
#Section for uploading building blocks file

    existing_file_intro = widgets.HTML("""<font size=3>Upload a file containing all the building blocks""",
        layout = widgets.Layout(height = '60px', width = '90%',
                        size = '20'))
    enter=widgets.Button(description='Upload file',layout= Layout(width='30%',border='solid 1px black'),style=style)
      
    existing_filebox = widgets.VBox(children=[existing_file_intro,space_box,name,space_box,enter] )
    

    def upload(e):
        
        #This function checks if the building blocks file uploaded is empty and also visualizes the building blocks present in the file
        
        file_name=name.value
        try:
            existing= open(file_name,"r")
            if os.stat(file_name).st_size == 0:
                print("Building blocks file is empty")
            else:
                lines=existing.readlines()
                for line in lines:
                    if line.startswith(('#','\n','\t')):
                        pass
                    else:
                        try:
                            mol=pybel.readstring("smi",line)
                            correct=widgets.HTML("""<font size=3>The file has been uploaded successfully""",
                    layout = widgets.Layout(height = '85px', width = '90%',
                                size = '20'))
                            smiles.append(line)

                            #|display(Javascript('IPython.notebook.execute_cells([2])'))
                        except:
                            incorrect= widgets.HTML(line +" in file is a incorrect SMILES") #checking the SMILES in the building blocks file
                            display(incorrect)
                # visualizing the SMILES from the building blocks file
                if len(smiles)>0:
                    
                            mol_lists = [Chem.MolFromSmiles(smile) for smile in smiles]
                            [mol.SetProp('_Name','B'+str(i)) for i,mol in enumerate(mol_lists)]
                            ibu1=Chem.Draw.MolsToGridImage(mol_lists)
                            display(ibu1)
                            generation_rules()
        except:
            # Checking if the building blocks file exists
            print("The building blocks file does not exist") # Executing upload function on click of enter button
           
    enter.on_click(upload) 
    BB_list=[]




    def on_another_clicked(i):
        # This function is for adding another building block
        global f
        clear_output()

        for element in f:

            try:
                mol=pybel.readstring("smi",element)
                BB_list.append(element)

            except:

                print("You typed a wrong SMILES "+element)

        BB.value=""

        display(widgets.HTML(value="<font color=blue><font size=5><b><u>BUILDING BLOCKS</font>"))
        display(tab1)
    final_BB=widgets.Button( description='Create building blocks file',layout= Layout(width= '45%',border='solid 1px black'),style=style)
    
    def on_button_clicked(d):
        # This function creates the building_blocks.dat file
        if input_direc.value=="":
            building_blocks= open("building_blocks.dat", "w")
        else:
            building_blocks= open(os.path.join(input_direc.value+"/building_blocks.dat"), "w")
            
        
        for h in range(len(BB_list)):
            building_blocks.write(BB_list[h] +'\n')

        building_blocks.close()
        generation_rules()

    def on_visualization_clicked(t):
        
        #This function visualizes the SMILES
        
        global f
        z=BB.value
        z=str(z)
        f=z.split(',')
        try:
            mol_list = [Chem.MolFromSmiles(element) for element in f]
            [mol.SetProp('_Name','B'+str(i)) for i,mol in enumerate(mol_list)]
            ibu=Chem.Draw.MolsToGridImage(mol_list)
            display(ibu)
        except:
            print("You typed a wrong SMILES. Unable to visualize")

    visualize=widgets.Button(description='Visualize',layout= Layout(width= '45%',border='solid 1px black'),style=style)
    visualize.on_click(on_visualization_clicked)    
    final_BB.on_click(on_button_clicked)
    another.on_click(on_another_clicked)
    type_smilesbox = widgets.VBox(children=[type_smiles_intro,BB,space_box,visualize,another,final_BB])
    tab1 = widgets.Tab(children=[type_smilesbox,existing_filebox],style=style)
    tab1.set_title(0, ' Individual SMILES')
    tab1.set_title(1, 'Upload file')
    display (tab1)
    

def generation_rules():
    """This function is used to take input for the generation rules from user."""
    
    style = {'description_width': 'initial','font_weight':'bold'}
    space_box = widgets.Box(layout=widgets.Layout(height ='55px', width='90%')) 
    second=widgets.Button(description='Next section',layout= Layout(width= 'auto',border='solid 1px black'),style=style)    
    display (second) 
    
    def second_section(q):
        
        if os.stat("building_blocks.dat").st_size <= 0:
            print("building blocks file is empty")
            
            
        else:
            display(widgets.HTML(value="<font color=blue><font size=5><b><u>USER DEFINED CONSTRAINTS</font>"))

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
            
           # Constraints
        
            style = {'description_width': 'initial'}
            bondmin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            bondmax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            bonds_intro=widgets.HTML("""Enter the minimum and maximum values of the following (integers)""",
            layout = widgets.Layout(height = '45px', width = '90%',
                        size = '20'))
            bonds = widgets.HBox([bondmin, bondmax])

            atommin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            atommax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            atoms = widgets.HBox([atommin, atommax])
            
            molmin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            molmax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            molecular_weight = widgets.HBox([molmin, molmax])

            ringmin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            ringmax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            rings = widgets.HBox([ringmin, ringmax])


            armin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            armax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            aromatic_rings = widgets.HBox([armin, armax])

            narmin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            narmax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            nonaromatic_rings = widgets.HBox([narmin, narmax])

            smin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            smax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            single_bonds=widgets.HBox([smin,smax])


            dmin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            dmax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            double_bonds = widgets.HBox([dmin, dmax])


            tmin = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            tmax = widgets.Text(description='', value='None',layout= Layout(width='30%'), style=style)
            triple_bonds = widgets.HBox([tmin, tmax])

            lmin = widgets.HTML(description='', value='Minumum',layout= Layout(width= '30%'), style=style)
            lmax = widgets.HTML(description='', value='Maximum',layout= Layout(width= '30%'), style=style)
            hl=widgets.HBox([lmin,lmax])
            left_box=widgets.VBox([Label(''),Label('Bonds'),Label('Atoms'),Label('Molecular Weight'),Label('Rings'),Label('Aromatic Rings'),Label('Non Aromatic Rings'),Label('Single Bonds'),Label('Double Bonds'),Label('Triple Bonds)')])
            Right_box=widgets.VBox([hl,bonds,atoms,molecular_weight,rings,aromatic_rings,nonaromatic_rings,single_bonds,double_bonds,triple_bonds])
            Hbox=widgets.HBox([left_box,Right_box],layout=Layout(border='solid 1px',width='60%', align_items='stretch'))
            min_max_intro=widgets.HTML("""Specify the minimum and maximum values of the following:""",
            layout = widgets.Layout(height = '45px', width = '90%',size = '20'))
            min_max=widgets.VBox(children=[min_max_intro,space_box,Hbox])

            specific_atoms = widgets.Text(
                value='None',
                description='Element:',
                style=style)

            satom_intro=widgets.HTML("""Specify the maximum number of heteroatoms that must be present in the  molecules in the final library.Expected type: tuple of tuple(s) (('Cl', 10), )""",
            layout = widgets.Layout(height = '45px', width = '90%',size = '20'))
            spec=widgets.HBox(children=[specific_atoms])
            satom=widgets.VBox(children=[satom_intro,space_box,spec])

            lipinski_rule = widgets.RadioButtons(options=['True', 'False'],value='False',description='Lipinski rule:', style=style, disabled=False)
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
                description='Include BB:', style={'description_width': 'initial'},
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
            
            accordion = widgets.Tab(
                children=[bb, min_max,satom, lr, fm,
                          sub, sube, ibb],style=style)
            accordion.set_title(0, 'Include Building Blocks')
            accordion.set_title(1, 'Min max')
            accordion.set_title(2, 'Heteroatoms')
            accordion.set_title(3, 'Lipinski Rule')
            accordion.set_title(4, 'Fingerprint Matching')
            accordion.set_title(5, 'Substructure Inclusion')
            accordion.set_title(6, 'Substructure Exclusion')
            accordion.set_title(7, 'Include initial Building Blocks')
            display(accordion)
            submit=widgets.Button(description="Submit",layout=widgets.Layout(width='20%',border='solid 1px black'),style=style)

            display(submit)
            button1 = widgets.Button(description="Generate configuration file",layout=widgets.Layout(width='20%',border='solid 1px black'),style=style)

        
        def command_line_arguments(b):
      #This function is for taking the command line arguments to run ChemLG from user  
    
            display(widgets.HTML(value="<font color=blue><font size=5><b><u>COMMAND LINE ARGUMENTS</font>"))


            """input_file = widgets.Text(value='building_blocks.dat',
                                  placeholder='Type the name of input file',
                                  description='Enter building blocks file name:',
                                  style=style,layout = widgets.Layout(height = '45px', width = '40%',
                        size = '20')
                                  )
            

            molecule_type = widgets.Dropdown(
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
            
            global generation_level
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
            
            global directory
            directory_intro=widgets.HTML("""Specify the directory in which you want to save the output files. Default value- current directory""",
            layout = widgets.Layout(height = '50px', width = '90%',
                        size = '20'))
            directory=widgets.Text(value="",placeholder='Enter path',description='Directory path',style=style)
            direc=widgets.VBox(children=[directory_intro,space_box,directory])
            

            arguments=widgets.Tab(children=[comb, gl, op, mf, ln, direc])
            
            arguments.set_title(0, 'Combination Type')
            arguments.set_title(1, 'No of generations')
            arguments.set_title(2, 'Output File Format')
            arguments.set_title(3, 'Maximum No of Files')
            arguments.set_title(4, 'Library Name')
            arguments.set_title(5, 'Output Directory')
            display(arguments)

            display(button1)
            

            
            def generation_file():
                #This function generates the config.dat file based on the user input for the generation rules

                display(widgets.HTML(value="""<font size=3>Configuration file created""",layout = widgets.Layout(height = '60px', width = '90%',size = '20')))
                if input_direc.value=="":
                    generation = open("config.dat", "w+")
                else:
                    generation=open(os.path.join(input_direc.value+"/config.dat"),"w")
                 
                generation.write("Please input generation rules below. Do not change the order of the options" + '\n' + "1. Include building blocks == " + building_blocks.value + '\n')
                if bondmin.value==bondmax.value=="None":
                    generation.write("2. Min and max no. of bonds == None"+'\n')
                else:
                    print("Changed values of No of bonds are Min "+bondmin.value+" Max "+bondmax.value)
                    tbonds=(int(bondmin.value),int(bondmax.value))
                    generation.write("2. Min and max no. of bonds =="+str(tbonds)+'\n')
                if atommin.value==atommax.value=="None":
                    generation.write("3. Min and max no. of atoms == None"+ '\n')
                else:
                    print("Changed values of No of atoms are Min "+atommin.value+" Max No "+atommax.value)
                    tatoms=(int(atommin.value),int(atommax.value))
                    generation.write("3. Min and max no. of atoms == "+ str(tatoms)+ '\n')
                if molmin.value==molmax.value=="None":
                    generation.write("4. Min and max mol. weight == None" + '\n')
                else:
                    print("Changed values of molecular weight are min "+molmin.value+" Max No of bonds "+molmax.value)
                    tmol=(int(molmin.value),int(molmax.value))
                    generation.write("4. Min and max mol. weight == "+str(tmol) + '\n')
                if ringmin.value==ringmax.value=="None":
                    generation.write("5. Min and max no. of rings == None"+ '\n')
                else:
                    print("Changed values of No of rings are Min "+ringmin.value+" Max No "+ringmax.value)
                    tring=(int(ringmin.value),int(ringmax.value))
                    generation.write("5. Min and max no. of rings == "+str(tring)+ '\n')
                if armin.value==armax.value=="None":
                    generation.write("6. Min and max no. of aromatic rings == None" + '\n')
                else:
                    print("Changed values of No of aromativ rings are Min "+armin.value+" Max "+armax.value)
                    taro=(int(armin.value),int(armax.value))
                    generation.write("6. Min and max no. of aromatic rings == " +str(taro)+ '\n')
                if narmin.value==narmax.value=="None":
                    generation.write("7. Min and max no. of non aromatic rings == None" + '\n')
                else:
                    print("Changed values of No of non-aromatic rings are Min "+narmin.value+" Max "+narmax.value)
                    tnar=(int(narmin.value),int(narmax.value))
                    generation.write("7. Min and max no. of non aromatic rings == "+str(tnar) + '\n')
                if smin.value==smax.value=="None":
                    generation.write("8. Min and max no. of single bonds == None" + '\n')
                else:
                    print("Changed values of No of single bonds are Min "+smin.value+" Max "+smax.value)
                    tsbond=(int(smin.value),int(smax.value))
                    generation.write("8. Min and max no. of single bonds == "+str(tsbond) + '\n')
                if dmin.value==dmax.value=="None":
                    generation.write("9. Min and max no. of double bonds == None" + '\n')
                else:
                    print("Changed values of No of double bonds are Min "+dmin.value+" Max "+dmax.value)
                    tdbond=(int(dmin.value),int(dmax.value))
                    generation.write("9. Min and max no. of double bonds == "+str(tdbond) + '\n')
                if tmin.value==tmax.value=="None":
                    generation.write("10. Min and max no. of triple bonds == None"+'\n')
                else:
                    print("Changed values of No of triple bonds are Min "+tmin.value+" Max "+tmax.value)
                    ttbond=(int(tmin.value),int(tmax.value))
                    generation.write("10. Min and max no. of triple bonds == "+str(ttbond)+'\n') 
                generation.write("11. Max no. of specific atoms == " + specific_atoms.value + '\n' +
                    "12. Lipinski's rule == " + lipinski_rule.value + '\n' + "13. Fingerprint matching == " + fingerprint_matching.value + '\n'+"14. Substructure == " + substructure.value + '\n' + "15. Substructure exclusion == " + substructure_exclusion.value + '\n' + '\n'+'\n'+'\n'+'\n'+'\n'+'\n'+
                    
                    "Combination type for molecules :: "+combination_type.value+'\n'+
                    "Number of generations :: "+str(generation_level.value)+'\n'+
                    "Molecule format in output file ::"+output_type.value+'\n'+
                    "Maximum files per folder :: "+str(max_files.value)+'\n'+"Library name :: "+library_name.value+'\n' )

                generation.close()

               

            def on_button_clicked(b):
                #This function takes the environment name to run chemlg as the input from the user
                generation_file()
                global envi
                global environment
                global run
                global fitness
                
                envi_intro=widgets.HTML("""Specify the virual environment in which you want to run the library generator""",
            layout = widgets.Layout(height = '50px', width = '90%',
                        size = '20'))
                envi=widgets.Text(
                value='None',
                placeholder='Enter Environment name',
                description='Environment name',
                style=style)
                environment=widgets.VBox(children=[envi_intro,envi],border='solid 2px black')
                #display(environment)
                
                run=widgets.Button(description='Run Chemlg',layout= Layout(width= 'auto',border='solid 1px black'),style=style)
                #display(run)
                #run.on_click(on_run_clicked)
                
                def constraints(x):
                    GA_func=widgets.Textarea(description="GA function",placeholder="Type the GA function to be incorporated for library selection", layout=Layout(width='80%', height='300px'))
                    if x== 'with constraints':
                        display(GA_func)
                        zero=widgets.Button(description="Next")
                        display(zero)
                        def on_ga_button_click(a):
                            
                                
                            if GA_func.value=="":
                                print("GA function empty")
                            else:
                                try:
                                    code=compile(GA_func.value,"", "exec")
                                    exec(code)
                                    f = open('objective_func.py', 'w+')
                                    f.write(GA_func.value)
                                    f.close()
                                    fitness = widgets.Text(
                                        value='None',
                                        placeholder='Type fitness value as a tuple',
                                        description='Fitness Value',
                                        style=style)
                                    fitness_intro=widgets.HTML("""Enter the fitness values. Expected type: tuple of tuple(s): (("max",0.7),)
                                    """,layout = widgets.Layout(height = '45px', width = '90%',size = '20'))
                                    fit=widgets.VBox(children=[fitness_intro,fitness])
                                    
                                    display(environment)
                                    display(run)
                                    
                                    run.on_click(use_ga)
                                except:
                                    print("Error in objective function")
                                

                        zero.on_click(on_ga_button_click)
                    else:
                        
                        display(environment)
                        display(run)
                        run.on_click(on_run_clicked)


                interact(constraints, x=widgets.RadioButtons(options=['with constraints', 'without constraints'],value='without constraints',description='Run ChemLG ',disabled=False))
                display(interact)
            def use_ga(j):
                
                from chemlg.constrained import GeneticAlgorithm                     # import genetic algorithm
                from objective_func import objective
                ga_library = GeneticAlgorithm(evaluate=objective,
                                              fitness     =  fitness.value,
                                              bb_file     =  'building_blocks.dat',
                                              config_file =  'config.dat',
                                              output_dir  =  directory.value)

                best_ind_df, best_individual = ga_library.search(n_generations=generation_level.value)
                
                
            def on_run_clicked(c):
            
                if envi.value=="":
                    envi.value=='None'
                    if directory.value=="":
                        os.system("chemlgshell -i config.dat -b building_blocks.dat -o ./output/")
                    else:
                         os.system("chemlgshell -i config.dat -b building_blocks.dat -o "+directory.value)
                else:
                    try:
                        os.system("activate "+envi.value)
                    except:
                        os.system("source activate"+envi.value)
                    if directory.value=="":
                        os.system("chemlgshell -i config.dat -b building_blocks.dat -o ./")
                    else:
                        os.system("chemlgshell -i config.dat -b building_blocks.dat -o "+directory.value)
                
                statistics_heading=widgets.HTML(value="""<font color=blue><font size=5><b><u>Statistics</font>""")
                display(statistics_heading)
                statistics_intro=widgets.HTML(""" Generate Statistics of the generated library
                                    """,layout = widgets.Layout(height = '45px', width = '90%',size = '20'))

                stats=widgets.Button(description='Run Statistics',layout= Layout(width= '20%',border='solid 1px black'),style=style)
                statbox=widgets.VBox(children=[statistics_intro,stats],border='solid 2px black')
                display(statbox)
                stats.on_click(on_stats_clicked)

            def on_stats_clicked(d):
                """This function gives a count of each building block in a molecule"""
                from chemlg.feasibility import generate_statistics
                generate_statistics()
                feasibility_heading=widgets.HTML(value="""<font color=blue><font size=5><b><u>Synthetic Feasibility</font>""")
                display(feasibility_heading)
                feasibility_intro=widgets.HTML(""" Test the synthetic feasibility of the generated library molecules
                                    """,layout = widgets.Layout(height = '45px', width = '90%',size = '20'))
                feasibility=widgets.Button(description='Run Feasibility Analysis',layout= Layout(width='20%',border='solid 1px black'),style=style)
                feasi=widgets.VBox(children=[feasibility_intro,feasibility],border='solid 2px black')
                display(feasi)
                feasibility.on_click(on_feasibility_clicked)
            def on_feasibility_clicked(c):
                """This function is used to calculate the feasibility of generated library molecules"""
                from chemlg.feasibility import feasibility
                feasibility()




            button1.on_click(on_button_clicked)
        submit.on_click(command_line_arguments)
    
    second.on_click(second_section)

config_builder()
