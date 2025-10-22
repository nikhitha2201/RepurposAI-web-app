#/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def similarity(smile1, smile2, sim_type = 'dice'):
    ms = [Chem.MolFromSmiles(smile1), Chem.MolFromSmiles(smile2)]
    fpgen = AllChem.GetRDKitFPGenerator()
    fps = [fpgen.GetFingerprint(x) for x in ms]
    if sim_type == 'tanimoto':
        sm = DataStructs.TanimotoSimilarity(fps[0],fps[1])
    elif sim_type == 'dice':
        sm = DataStructs.DiceSimilarity(fps[0],fps[1])  
    elif sim_type == 'cosine':
        sm = DataStructs.CosineSimilarity(fps[0],fps[1])
    elif sim_type == 'sokal':
        sm = DataStructs.SokalSimilarity(fps[0],fps[1])
    elif sim_type == 'russel':
        sm = DataStructs.RusselSimilarity(fps[0],fps[1])
    elif sim_type == 'kulczynski':
        sm = DataStructs.KulczynskiSimilarity(fps[0],fps[1])
    elif sim_type == 'mcconnaughey':
        sm = DataStructs.McConnaugheySimilarity(fps[0],fps[1])
    else:
        raise ValueError("Unknown similarity type: {}".format(sim_type))
    return [{"compound1": smile1,"compound2": smile2, "similarity": sm}]

if __name__ == "__main__":
    ## testing
    types = ['tanimoto', 'dice', 'cosine', 'sokal', 'russel', 'kulczynski', 'mcconnaughey']
    smile1 = "CCO"
    smile2 = "CCN"
    for type in types:
        print(similarity(smile1, smile2, sim_type=type))
 