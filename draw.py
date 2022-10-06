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

def MolDensityMap(mol, weights, colorMap="bwr", scale=-1, size=(400, 400), sigma=0.03,
  coordScale=1.5, step=0.01, colors='k', contourLines=10, alpha=0.5, **kwargs):
  """
  Map the weight to molecules and get the density map 
  """
  import math 
  import numpy as np 
  fig = Draw.MolToMPL(mol, coordScale=coordScale, size=size, **kwargs)
  text_kwargs = dict(ha='center', va='center', fontsize=28, color='C1')
  if sigma is None:
    bond = mol.GetBondWithIdx(0)
    idx1 = bond.GetBeginAtomIdx()
    idx2 = bond.GetEndAtomIdx()
    sigma = 0.3 * math.sqrt(sum([(mol._atomPs[idx1][i] - mol._atomPs[idx2][i])**2 for i in range(2)]))
    sigma = round(sigma, 2)

  # Calculate the coutour
  x, y, z = Draw.calcAtomGaussians(mol, sigma, weights=weights, step=step)

  # Scaling
  if scale <= 0.0:
    maxScale = max(abs(np.min(z)), math.fabs(np.max(z)))
  else:
    maxScale = scale
  fig.axes[0].imshow(z, cmap=colorMap, interpolation='bilinear', origin='lower',
                   extent=(0, 1, 0, 1), vmin=-maxScale, vmax=maxScale)
  # Contour lines
  # Only draw them when at least one weight is not zero
  if len([w for w in weights if w != 0.0]):
    contourset = fig.axes[0].contour(x, y, z, contourLines, colors=colors, alpha=alpha, **kwargs)
  for j, c in enumerate(contourset.collections):
    if contourset.levels[j] == 0.0:
      c.set_linewidth(0.0)
    elif contourset.levels[j] < 0:
      c.set_dashes([(0, (3.0, 3.0))])
  fig.axes[0].set_axis_off()
  return fig


