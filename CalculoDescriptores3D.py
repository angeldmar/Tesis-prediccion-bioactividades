from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors3D
from rdkit.Chem import rdMolDescriptors
import pandas as pd

# Exportando base de datos de smiles a una lista
smiles = pd.read_csv('./SMILES.csv', delimiter=',')
smiles = smiles.transpose()
smiles = smiles.values.tolist()

# Recorriendo la lista para calcular descriptores

i = 0

Descriptores3D = {}
for x in smiles:
    for molecula in x:
        i += 1
        print("Generando conformación 3D de {}, molecula número {}".format(molecula, i))
        Descriptores3D_diccionario = {}
        m = Chem.MolFromSmiles(molecula)
        m2 = Chem.AddHs(m)
        AllChem.EmbedMolecule(m2)
        AllChem.MMFFOptimizeMolecule(m2)

        print("Iniciando calculos de descriptores3D de {}".format(molecula))
        Descriptores3D_diccionario['CalcPBF'] = rdMolDescriptors.CalcPBF(m2)
        Descriptores3D_diccionario['CalcAUTOCORR3D'] = rdMolDescriptors.CalcAUTOCORR3D(m2)
        Descriptores3D_diccionario['CalcRDF'] = rdMolDescriptors.CalcRDF(m2)
        Descriptores3D_diccionario['CalcMORSE'] = rdMolDescriptors.CalcMORSE(m2)
        Descriptores3D_diccionario['CalcWHIM'] = rdMolDescriptors.CalcWHIM(m2)
        Descriptores3D_diccionario['CalcGETAWAY'] = rdMolDescriptors.CalcGETAWAY(m2)
        Descriptores3D_diccionario['Asphericity'] = Descriptors3D.Asphericity(m2)
        Descriptores3D_diccionario['Eccentricity'] = Descriptors3D.Eccentricity(m2)
        Descriptores3D_diccionario['InertialShapeFactor'] = Descriptors3D.InertialShapeFactor(m2)
        Descriptores3D_diccionario['NPR1'] = Descriptors3D.NPR1(m2)
        Descriptores3D_diccionario['NPR2'] = Descriptors3D.NPR2(m2)
        Descriptores3D_diccionario['PMI1'] = Descriptors3D.PMI1(m2)
        Descriptores3D_diccionario['PMI2'] = Descriptors3D.PMI2(m2)
        Descriptores3D_diccionario['PMI3'] = Descriptors3D.PMI3(m2)
        Descriptores3D_diccionario['RadiusOfGyration'] = Descriptors3D.RadiusOfGyration(m2)
        Descriptores3D_diccionario['SpherocityIndex'] = Descriptors3D.SpherocityIndex(m2)

        Descriptores3D[molecula] = Descriptores3D_diccionario

        print("Finalizando calculos de descriptores3D de {}".format(molecula))

# Los descriptores almacenados en un diccionario se convierten a dataframes
# Se realizan correcciones en el formato del dataframe
# Se guardan como csv

Descriptores3D_Dataframe = pd.DataFrame.from_dict(Descriptores3D)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.transpose()

Descriptores3D_Dataframe_subset1 = Descriptores3D_Dataframe.astype('string')
Descriptores3D_Dataframe_subset1['CalcAUTOCORR3D'] = Descriptores3D_Dataframe_subset1['CalcAUTOCORR3D'].str.replace("[", "").str.replace("]", "")
Descriptores3D_Dataframe_subset1 = Descriptores3D_Dataframe_subset1['CalcAUTOCORR3D'].str.split(pat=",", expand=True).add_prefix('CalcAUTOCORR3D_value_').fillna(0)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.join(Descriptores3D_Dataframe_subset1)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.drop('CalcAUTOCORR3D', axis=1)

Descriptores3D_Dataframe_subset2 = Descriptores3D_Dataframe.astype('string')
Descriptores3D_Dataframe_subset2['CalcRDF'] = Descriptores3D_Dataframe_subset2['CalcRDF'].str.replace("[", "").str.replace("]", "")
Descriptores3D_Dataframe_subset2 = Descriptores3D_Dataframe_subset2['CalcRDF'].str.split(pat=",", expand=True).add_prefix('CalcRDF_value_').fillna(0)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.join(Descriptores3D_Dataframe_subset2)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.drop('CalcRDF', axis=1)

Descriptores3D_Dataframe_subset3 = Descriptores3D_Dataframe.astype('string')
Descriptores3D_Dataframe_subset3['CalcMORSE'] = Descriptores3D_Dataframe_subset3['CalcMORSE'].str.replace("[", "").str.replace("]", "")
Descriptores3D_Dataframe_subset3 = Descriptores3D_Dataframe_subset3['CalcMORSE'].str.split(pat=",", expand=True).add_prefix('CalcMORSE_value_').fillna(0)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.join(Descriptores3D_Dataframe_subset3)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.drop('CalcMORSE', axis=1)

Descriptores3D_Dataframe_subset4 = Descriptores3D_Dataframe.astype('string')
Descriptores3D_Dataframe_subset4['CalcWHIM'] = Descriptores3D_Dataframe_subset4['CalcWHIM'].str.replace("[", "").str.replace("]", "")
Descriptores3D_Dataframe_subset4 = Descriptores3D_Dataframe_subset4['CalcWHIM'].str.split(pat=",", expand=True).add_prefix('CalcWHIM_value_').fillna(0)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.join(Descriptores3D_Dataframe_subset4)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.drop('CalcWHIM', axis=1)

Descriptores3D_Dataframe_subset5 = Descriptores3D_Dataframe.astype('string')
Descriptores3D_Dataframe_subset5['CalcGETAWAY'] = Descriptores3D_Dataframe_subset5['CalcGETAWAY'].str.replace("[", "").str.replace("]", "")
Descriptores3D_Dataframe_subset5 = Descriptores3D_Dataframe_subset5['CalcGETAWAY'].str.split(pat=",", expand=True).add_prefix('CalcGETAWAY_value_').fillna(0)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.join(Descriptores3D_Dataframe_subset5)
Descriptores3D_Dataframe = Descriptores3D_Dataframe.drop('CalcGETAWAY', axis=1)

Descriptores3D_Dataframe.to_csv("./DescriptorsCSV/Descriptores3D_Dataframe.csv")
