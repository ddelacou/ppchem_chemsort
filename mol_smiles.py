from rdkit import Chem

Smiles = "CCO"
mol = Chem.MolFromSmiles(Smiles)
print(mol)