import requests
import json
import base64
from PIL import Image as pilimage
import numpy as np 
import os

def clean_mol2(mol2file):
  # Will put the output file to the same directory of to mol2 file
  import subprocess;
  print("Running the clean mol2 file function"); 
  cleaned_mol2file = mol2file.replace(".mol2", "_clean.mol2"); 
  cleanMOL2LP(mol2file, cleaned_mol2file); 
  sdffile = mol2file.replace(".mol2", ".sdf"); 
  sp = subprocess.Popen(["obabel", cleaned_mol2file, "-O", sdffile]); 
  sp.wait(); 
  
  if os.path.isfile(sdffile):
    print("Found the sdf output;");
    return sdffile
  elif os.path.isfile(cleaned_mol2file):
    print("Found the clean mol2 output;"); 
    return cleaned_mol2file
  else :
    print(f"Failed to clean mol2 file {mol2file}, using the unchanged mol2 file ")
    return mol2file

def cleanMOL2LP(mol2input, mol2out):
  import re
  # Read all molecules from the mol2 file 
  with open(mol2input, 'r') as file1:
    # Only keep the MOLECULE, ATOM and BOND field 
    mols_raw = ['@<TRIPOS>MOLECULE\n'+i for i in file1.read().split('@<TRIPOS>MOLECULE\n') if len(i)>1 and i[0]!="#"]

  try: 
    mols_cleaned = ["@<TRIPOS>MOLECULE" + mols_raw[i].split("@<TRIPOS>MOLECULE")[1].split("@<TRIPOS>")[0]+"@<TRIPOS>ATOM"+ \
    mols_raw[i].split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0] + "@<TRIPOS>BOND" + \
    mols_raw[i].split("@<TRIPOS>BOND")[1].split("@<TRIPOS>")[0] \
    for i in range(len(mols_raw))]
  except: 
    mols_cleaned = ["@<TRIPOS>MOLECULE" + mols_raw[i].split("@<TRIPOS>MOLECULE")[1].split("@<TRIPOS>")[0]+"@<TRIPOS>ATOM"+ \
    mols_raw[i].split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0] + "@<TRIPOS>BOND" + \
    mols_raw[i].split("@<TRIPOS>BOND")[1] \
    for i in range(len(mols_raw))]

  mols = [i for i in mols_cleaned if len(i) > 1 ]
  # Find dummy atoms and index them
  with open(mol2out,'w') as file1:
    for molidx, mol in enumerate(mols):
      # Read mol2 parts
      molhead = mol.split("@<TRIPOS>MOLECULE")[1].split("@<TRIPOS>")[0].strip("\n").split("\n");
      molname = mol.split("\n")[1]
      atms    = np.array([i.split() for i in mol.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].split("\n") if len(i) >0])
      bds     = np.array([i.split() for i in mol.split("@<TRIPOS>BOND")[1].split("@<TRIPOS>")[0].split("\n") if len(i) >0])
      # Find lone pair / dummy atoms 
      atomtypes = np.array([i.lower() for i in atms[:,5]])
      idxs = np.where((atomtypes == 'lp') | (atomtypes == 'du') | (atomtypes == 'xx'))[0]
      if len(idxs) > 0:
        try:
          exp_idxs = [i+1 for i in idxs]
          charges = atms[:,8][idxs].astype(float)
          # find bond partner         
          partneridx = [np.where( ((bds[:,1]==str(i+1))|(bds[:,2]==str(i+1))) )[0] for i in idxs ]
          partneridx = np.array([i[0] for i in partneridx])
          partnerbds = bds[partneridx][:,1:3].astype(int)
          partners = [i[1] if (i[0] in exp_idxs) else i[0] for i in partnerbds]
          # add charge to its partner 
          for i in range(len(idxs)):
            print('Editing the molecule {} -> {}:{}:{} <= {}:{}:{}'.format( molname, atms[partners[i]-1,1], partners[i], atms[partners[i]-1,8], atms[idxs[i],1],idxs[i], charges[i]))
            charge_aft = float(atms[partners[i]-1,8]) + charges[i]
            atms[partners[i]-1,8] = str(np.round(charge_aft,3))
          # delete lone pair / dummy atoms 
          atms = np.delete(atms, idxs, axis=0);
          atommap = {}
          for i in range(len(atms)):
            atommap[atms[i,0]] = str(i+1)
            atms[i,0] = str(i+1)
          bds  = np.delete(bds , partneridx, axis=0);
          for i in range(len(bds)):
            bds[i,0] = str(i+1)
            bds[i,1] = atommap[bds[i,1]]
            bds[i,2] = atommap[bds[i,2]]
          molhead[1] = f'{len(atms)} {len(bds)}'
          molhead = "\n".join(molhead);
          stratms = '\n'.join([' '.join(i) for i in atms]);
          strbd   = '\n'.join([' '.join(i) for i in bds]);
          molstr = f'@<TRIPOS>MOLECULE\n{molhead}\n@<TRIPOS>ATOM\n{stratms}\n@<TRIPOS>BOND\n{strbd}\n';
        except:
          molstr = ""
          print(f"Failed to deal with compound {molidx+1}")
      else :
        stratms = '\n'.join([' '.join(i) for i in atms])
        strbd   = '\n'.join([' '.join(i) for i in bds])
        molhead = "\n".join(molhead);
        molstr = f'@<TRIPOS>MOLECULE\n{molhead}\n@<TRIPOS>ATOM\n{stratms}\n@<TRIPOS>BOND\n{strbd}\n';
      # Write the mol2 to a new mol2 file
      file1.write(molstr)

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



def AvailabilityCheck(libraryfile, libraryout, threshold=0, forceMCULE=False, forceCS=False):
  from rdkit import Chem
  from Doana import prices
  with open(libraryfile, "r") as file1: 
    lib = [f"@<TRIPOS>MOLECULE{i}" for i in file1.read().strip("@<TRIPOS>MOLECULE").split("@<TRIPOS>MOLECULE")]; 
  with open(libraryout, "w") as fileout: 
    for index, i in enumerate(lib): 
      try: 
        themol = Chem.MolFromMol2Block(i); 
        thesmi = Chem.MolToSmiles(themol)
        compounds = [thesmi]; 
        compounds = [i.replace("*", "") for i in compounds]
        try: 
          querier.querymols(compounds)
        except: 
          querier = prices.MolQuerier()
          querier.MCULE_TOKEN ="cd69ac4c24cccc1e99b5d6d8f0cb4c267ad84f47"
          querier.querymols(compounds)
        price_count = 0; 
        if isinstance(querier.pricedic[thesmi]["mcule_price"], str):
          print("not found mcule prices" )
          if forceMCULE == True: 
            continue
        elif isinstance(querier.pricedic[thesmi]["mcule_price"], list):
          price_count += len(querier.pricedic[thesmi]["mcule_price"]); 

        if isinstance(querier.pricedic[thesmi]["cs_price"], str):
          print("Not found ChemSpace prices"); 
          if forceCS == True: 
            continue
        elif isinstance(querier.pricedic[thesmi]["cs_price"], list):
          price_count += len(querier.pricedic[thesmi]["cs_price"]); 
        print(f"Compound {index} found {price_count} prices")
        if price_count > threshold: 
          fileout.write(i); 
      except: 
        continue 

def ConcVINAResults(basepath):
  with open(f"{basepath}/parms.json", "r") as file1:
    parms = json.load(file1)
    print(parms)
  sess_id = parms["session_id"]
  jobnr = int(parms['jobnr']);
  molout = open(f"{basepath}/all_mols.mol2", 'w')
  scoreout = open(f"{basepath}/all_scores.dat", 'w')
  c = 0;
  for taski in range(1, jobnr+1):
    if os.path.isfile(f"{basepath}/task{taski}_worked.mol2"):
      if os.path.isfile(f"{basepath}/task{taski}_outcome.json"):
        # Check pose number and score number.
        with open(f"{basepath}/task{taski}_worked.mol2", "r") as molfile, open(f"{basepath}/task{taski}_outcome.json", "r") as scorefile:
          mollines = molfile.read()
          mols = [f"@<TRIPOS>MOLECULE{i}" for i in mollines.strip("@<TRIPOS>MOLECULE").split("@<TRIPOS>MOLECULE")]
          molnr = len(mols)

          scores = json.load(scorefile)
          poses = list(scores["poses"].keys())
          posenr = len(poses)
          if molnr != posenr:
            print(f"The number of molecules and the number of score in task {taski} is not equil, skipping this task.")
          else:
            print(f"Processing task {taski}; ");
            molout.write("".join(mols))
            for p in poses:
              c+=1
              thesmi = scores["poses"][p]["smiles"];
              thescores = scores["poses"][p]["scores"];
              themolname = scores["poses"][p]["molname"];
              thedate = scores["poses"][p]["molfinishdate"].replace(" ", "-")
              scorestr = " ".join([str(i) for i in thescores])
              scoreout.write(f"{c} {themolname} {thesmi} {scorestr} {sess_id} {thedate}\n")
      else:
        print(f"Score file {taski} not found")
    else:
      print(f"MOL2 output file {taski} not found")
  molout.close()
  scoreout.close()
