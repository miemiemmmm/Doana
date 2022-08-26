import requests
import json
import base64
from PIL import Image as pilimage

def toformat(themol, theformat):
  json_data = {
    'structure': themol,
    'parameters': theformat,
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
  if response.status_code != 200: 
    print(response.status_code, response.reason, response.url); 
    return
  else: 
    struct = json.loads(response.text)["structure"]
    return struct

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
  if response.status_code != 200: 
    print(response.status_code, response.reason, response.url, response.text)
    return
  else:
    imgstr = json.loads(response.text)["binaryStructure"];
    image = base64.b64decode(imgstr);
    with open(savefile, 'wb') as f:
      f.write(image)
    image = pilimage.open(savefile)
    image.show()

def hydrogenize(struct):
  json_data = {
    'structure': struct,
    'parameters': {
      'method': 'HYDROGENIZE',
      "ph":7,
    },
  }
  response = requests.post('https://marvinjs-demo.chemaxon.com/rest-v1/util/convert/hydrogenizer', json=json_data)
  if response.status_code == 200:
    return response.text
  else:
    return None

def gen3d(struct):
  json_data = {
    'structure': struct,
    "inputFormat": "smiles",
    "outputFormat": "mrv",
    'parameters': {
      'dim': 3,
    },
  }
  response = requests.post('https://marvinjs-demo.chemaxon.com/rest-v1/util/convert/clean', json=json_data)
  if response.status_code == 200:
    return response.text
  else:
    print(response.status_code, response.reason, response.url);




