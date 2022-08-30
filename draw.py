import rdkit.Chem as Chem
from IPython.display import Image
from rdkit.Chem import AllChem, Draw
import requests
import json
import base64
from PIL import Image as pilimage

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

def show_png(themol, savefile="/tmp/example_image.png"):
  json_data = {
    'structure': themol,
    'parameters': "png",
    'filterChain': [
      {
        'filter': 'standardizer',
        'parameters': {
          'standardizerDefinition': 'aromatize',
        },
      },
    ],
  }
  response = requests.post('https://marvinjs-demo.chemaxon.com/rest-v1/util/calculate/molExport', json=json_data)
  if response.status_code != 200: print(response.status_code, response.reason, response.url, response.text)
  imgstr = json.loads(response.text)["binaryStructure"];
  image = base64.b64decode(imgstr);
  with open(savefile, 'wb') as f:
    f.write(image)
  image = pilimage.open(savefile)
  image.show()

