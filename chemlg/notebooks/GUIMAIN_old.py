from __future__ import print_function
from ipywidgets import interact, Layout
import ipywidgets as widgets
from IPython.display import display
style = {'description_width': 'initial'}
count =0
def GUI():
    style = {'description_width': 'initial'}
    
   
    style = {'description_width': 'initial'}
    BB = widgets.Text(
        placeholder='Type the SMILES of building block',
        description='Building block:',style=style,
         )
    display(BB)
    button = widgets.Button(description="Add the building block")
    display(button)
 
    def on_button_clicked(b):
        global count
        count=count+1
        if count>1:
            building_blocks2=open("building_blocks.txt","a+")
            building_blocks2.write(BB.value + '\n')
            building_blocks2.close()
            
            
        else:
            building_blocks1 = open("building_blocks.txt", "w+")
            building_blocks1.write('Building blocks are:' + '\n' + BB.value + '\n')
            building_blocks1.close()
            count=2
            
            
    button.on_click(on_button_clicked)
    style = {'description_width': 'initial'}
    building_blocks = widgets.Text(value='F1',
                                   placeholder='Type the building blocks to be used',
                                   description='Building Blocks:',
                                   style=style,
                                   )    
    
    style = {'description_width': 'initial'}
    bondmin = widgets.Text(description='Minimum no of bonds', value='None', style=style)
    bondmax = widgets.Text(description='Maximum no of bonds', value='None', style=style)
    x = widgets.HBox([bondmin, bondmax])

    atommin = widgets.Text(description='Minimum No of Atoms', value='None', style=style)
    atommax = widgets.Text(description='Maximum No of Atoms', value='None', style=style)
    y = widgets.HBox([atommin, atommax])

    molmin = widgets.Text(description='Minimum Molecular Weight', value='None', style=style)
    molmax = widgets.Text(description='Maximum Molecular Weight', value='None', style=style)
    z = widgets.HBox([molmin, molmax])

    ringmin = widgets.Text(description='Minimum no of Rings', value='None', style=style)
    ringmax = widgets.Text(description='Maximum no of Rings', value='None', style=style)
    a = widgets.HBox([ringmin, ringmax])

    armin = widgets.Text(description='Minimum no of Aromatic rings', value='None', style=style)
    armax = widgets.Text(description='Maximum no of Aromatic Rings', value='None', style=style)
    b = widgets.HBox([armin, armax])

    narmin = widgets.Text(description='Minimum no of Non Aromatic rings', value='None', style=style)
    narmax = widgets.Text(description='Maximum no of Non Aromatic Rings', value='None', style=style)
    c = widgets.HBox([narmin, narmax])

    smin = widgets.Text(description='Minimum no of Single Bonds', value='None', style=style)
    smax = widgets.Text(description='Maximum no of Single Bonds', value='None', style=style)
    d = widgets.HBox([smin, smax])

    dmin = widgets.Text(description='Minimum no of Double Bonds', value='None', style=style)
    dmax = widgets.Text(description='Maximum no of Double Bonds', value='None', style=style)
    e = widgets.HBox([dmin, dmax])

    tmin = widgets.Text(description='Minimum no of Tiple Bonds', value='None', style=style)
    tmax = widgets.Text(description='Maximum no of Triple Bonds', value='None', style=style)
    f = widgets.HBox([tmin, tmax])

    specific_atoms = widgets.Text(
        value='None',
        placeholder='Type the number of specific atoms',
        description='Maximum No. of Specific Atoms:',
        style=style)

    lipinski_rule = widgets.RadioButtons(
        options=['True', 'False'],
        value='False',
        description='Lipinski rule:', style=style,
        disabled=False
    )

    fingerprint_matching = widgets.Text(
        value='None',
        placeholder='Finger print match',
        description='Finger print matching:',
        style=style)

    substructure = widgets.Text(
        value='None',
        placeholder='Type the substructure SMILES',
        description='Substructure:',
        style=style)

    substructure_exclusion = widgets.Text(
        value='None',
        placeholder='Type the substructure exclusion SMILES',
        description='Substructure exclusion:',
        style=style)

    Include_BB = widgets.RadioButtons(
        options=['Yes', 'No'],
        value='No',
        description='Include BB:', style=style,
        disabled=False
    )

    """Symmetry = widgets.RadioButtons(
        options=['Yes', 'No'],
        value='No',
        description='Symmetry:', style=style,
        disabled=False
    )"""
    accordion = widgets.Accordion(
        children=[building_blocks, x, y, z, a, b, c, d, e, f, specific_atoms, lipinski_rule, fingerprint_matching,
                  substructure, substructure_exclusion, Include_BB])
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

    button1 = widgets.Button(description="Generate rules file")
    display(button1)

    def on_button_clicked(b):
        generation_file()

    button1.on_click(on_button_clicked)
    input_file = widgets.Text(value='building_blocks.dat',
                              placeholder='Type the name of input file',
                              description='Input file name:',
                              style=style,
                              )
    molecule_type = widgets.Dropdown(
        options=['SMILES', 'SMARTS', 'INCHI'],
        value='SMILES',
        description='Molecule Type:',
        style=style,
        disabled=False,
    )
    combination_type = widgets.RadioButtons(
        options=['Fuse', 'Link'],
        value='Link',
        description='Combination type:',
        style=style,
        disabled=False
    )
    generation_level = widgets.BoundedIntText(
        value=1,
        min=1,
        max=100,
        step=1,
        description='Generation level:',
        style=style,
        disabled=False)

    output_type = widgets.Dropdown(
        options=['SMILES', 'SMARTS', 'INCHI'],
        value='SMILES',
        description='Output Type:',
        style=style,
        disabled=False, )

    max_files = widgets.IntText(
        value=10000,
        description='Maximum files per folder:',
        style=style,
        disabled=False)

    library_name = widgets.Text(
        value='new_library_',
        placeholder='Type the nameof the library',
        description='Library Name:',
        style=style)

    style = {'description_width': 'initial', }
    chemHTPS = widgets.RadioButtons(
        options=['Yes', 'No'],
        value='No',
        description='Run with ChemHTPS:', style=style,
        disabled=False
    )
    children = [input_file, molecule_type, combination_type, generation_level, output_type, max_files, library_name,
                chemHTPS]
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
    button = widgets.Button(description="Generate Command line", layout= Layout(width= 'auto'))
    display(button)

    def on_button_clicked(b):
        opt = "--input_file %s --molecule_type %s --combination_type %s --generation_levels %s --output_type %s--max_files_per_folder%s --rule_file generation.dat --lib_name %s --ChemHTPS %s" % (
        input_file.value, molecule_type.value, combination_type.value, generation_level.value, output_type.value,
        max_files.value, library_name.value, chemHTPS.value)
        print(opt)

    button.on_click(on_button_clicked)







