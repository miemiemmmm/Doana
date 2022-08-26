import rdkit.Chem as Chem
from IPython.display import Image
from rdkit.Chem import AllChem, Draw

def drawGridMols(results, legends=[], mpr=5):
  print(f"Drawing grid molecules {len(results)}")
  results = [i.replace("*","") for i in results]
  mols = [Chem.MolFromSmiles(m) for m in results]; 
  for m in mols: 
    tmp = AllChem.Compute2DCoords(m);
  img = Draw.MolsToGridImage(mols, molsPerRow=mpr, subImgSize=(400,400), maxMols=100,
                             legends=legends, returnPNG=False)
  img.save("/tmp/molgrid.png")
  return Image("/tmp/molgrid.png")



