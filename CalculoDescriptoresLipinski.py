from rdkit import Chem
from rdkit.Chem import Lipinski
import pandas as pd

# Exportando base de datos de smiles a una lista

smiles = pd.read_csv('./SMILES.csv', delimiter=',')  
smiles = smiles.transpose()
smiles = smiles.values.tolist()

# Recorriendo la lista para calcular descriptores

DescriptoresLipinski = {}
for x in smiles:
    for molecula in x:  
        Lipinski_Diccionario = {}
        m = Chem.MolFromSmiles(molecula)
        DescriptoresLipinski_Diccionario = {}

        DescriptoresLipinski_Diccionario['FractionCSP3'] = Lipinski.FractionCSP3(m)
        DescriptoresLipinski_Diccionario['HeavyAtomCount'] = Lipinski.HeavyAtomCount(m)
        DescriptoresLipinski_Diccionario['NHOHCount'] = Lipinski.NHOHCount(m)
        DescriptoresLipinski_Diccionario['NOCount'] = Lipinski.NOCount(m)
        DescriptoresLipinski_Diccionario['NumAliphaticCarbocycles'] = Lipinski.NumAliphaticCarbocycles(m)
        DescriptoresLipinski_Diccionario['NumAliphaticHeterocycles'] = Lipinski.NumAliphaticHeterocycles(m)
        DescriptoresLipinski_Diccionario['NumAromaticCarbocycles'] = Lipinski.NumAromaticCarbocycles(m)
        DescriptoresLipinski_Diccionario['NumHAcceptors'] = Lipinski.NumHAcceptors(m)
        DescriptoresLipinski_Diccionario['NumHDonors'] = Lipinski.NumHDonors(m)
        DescriptoresLipinski_Diccionario['NumHeteroatoms'] = Lipinski.NumHeteroatoms(m)
        DescriptoresLipinski_Diccionario['NumRotatableBonds'] = Lipinski.NumRotatableBonds(m)
        DescriptoresLipinski_Diccionario['NumSaturatedCarbocycles'] = Lipinski.NumSaturatedCarbocycles(m)
        DescriptoresLipinski_Diccionario['NumSaturatedHeterocycles'] = Lipinski.NumSaturatedHeterocycles(m)
        DescriptoresLipinski_Diccionario['NumSaturatedRings'] = Lipinski.NumSaturatedRings(m)
        DescriptoresLipinski_Diccionario['RingCount'] = Lipinski.RingCount(m)
        
        DescriptoresLipinski[molecula] = DescriptoresLipinski_Diccionario

# Los descriptores almacenados en un diccionario se convierten a dataframes, se guardan como csv
DescriptoresLipinski_Dataframe = pd.DataFrame.from_dict(DescriptoresLipinski) 
DescriptoresLipinski_Dataframe = DescriptoresLipinski_Dataframe.transpose()

DescriptoresLipinski_Dataframe.to_csv("./DescriptorsCSV/DescriptoresLipinski_Dataframe.csv")
