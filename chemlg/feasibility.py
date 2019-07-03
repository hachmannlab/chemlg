import csv
import pandas as pd
from collections import Counter
import csv
import pybel
from ipywidgets import interact, Layout, Label
import ipywidgets as widgets
from IPython.display import display

style = {'description_width': 'initial','font_weight':'bold'}
def feasibility():
    with open('database.txt') as csv_file:         #database.txt is the database file from MolPort
        csv_reader = csv.reader(csv_file, delimiter='\t')
        line_count = 0
        commercial=[]
        for row in csv_reader:
            if line_count == 0:
               
                line_count += 1
            else:
                commercial.append(row[0])
               
                line_count += 1
    with open('final_smiles.csv') as csv_file:    # final_smiles.csv is the file of generated molecules
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        library=[]
        for row in csv_reader:
            library.append(row[0])
    cutoff=widgets.Text(description="Cut-off Value",placeholder="Enter the cut off value of Tanimoto index",style=style)
    cutoff_button=widgets.Button(description='Confirm cut-off',layout= Layout(width= '35%',border='solid 1px black'))
    display(cutoff)
    display(cutoff_button)
    
    #cutoff=float(input("Enter the cut off value of tanimoto index"))

    for mol1 in library:
        final_library=pybel.readstring('smi',mol1) #reading generated library molecule
    final=[]
    tanimoto1=[]
    fll=[]
    
    def cutoff_clicked(m):
        for i in library:
            fl=pybel.readstring('smi',i)
            for mol in commercial:
                smiles=pybel.readstring('smi',mol) #reading library molecule
                tanimoto=smiles.calcfp()|fl.calcfp()

                if tanimoto>=cutoff.value:
                    tanimoto1.append(tanimoto)
                    final.append(mol)
                    fll.append(i)
                else:
                    pass


           # print(final)
            #print(tanimoto1)
        data={'Molecule':final, 'Tanimoto':tanimoto1,'Final library':fll}
        df=pd.DataFrame(data)
        pd.set_option('display.max_colwidth', -1)
        #print(df)
        grouped=df.groupby(['Final library']).count()[['Molecule']]
        sorted_df = grouped.sort_values(by='Molecule', ascending=False)
        sorted_df.to_excel("dataframe.xlsx")

        df.to_excel("output.xlsx")
    cutoff_button.on_click(cutoff_clicked)
        
    
def generate_statistics():
    """This function gives a count of each building block in a molecule"""

    with open('final_library.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        library=[]
        smiles=[]
        for row in csv_reader:
            library.append(row[0])
            smiles.append(row[1])
        #print(library)

        #print(smiles)

    building_blocks=[]
    for bb in library:
        building_blocks.append(bb.split('-'))
    #print(building_blocks)
    number_of_building_blocks=[]
    for block in building_blocks:
        number_of_building_blocks.append(dict(Counter(block)))

    df=pd.DataFrame(number_of_building_blocks,index=smiles)
    df.fillna(0,inplace= True)
    df['Total building blocks'] = df.sum(axis=1, numeric_only=True)
    pd.set_option('display.max_colwidth', -1)
    #print(df)
    df.to_excel("statistics.xlsx")
#generate_statistics()
#feasibility()


#user needs to isntall openpyxl