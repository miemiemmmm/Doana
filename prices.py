import json
import requests
import urllib.parse
import pandas as pd 
from bs4 import BeautifulSoup 
    
class MolQuerier:
  """
  Use example: 
  compounds = ["CCOC(=O)c1c(C)[nH]c(c1C)C(=O)[C@@H](N1CCN(CC1)c1ccc(cn1)C(F)(F)F)C"]
  compounds = [i.replace("*", "") for i in compounds]
  try:
     ret_response = querier.querymols(compounds)
  except: 
     querier = MolQuerier()
     ret_response = querier.querymols(compounds)
  querier.formatprint()
  """
  def __init__(self):
    self.CSORIGIN = 'https://chem-space.com'; 
    self.CSREFERER = 'https://chem-space.com/search'; 

    self.CSCATMAP = {"CSCS": "Custom Request", "CSMS":"Make-On-Demand Screening Compounds", 
                     "CSSB":"In-Stock Building Blocks", "CSSS":"In-stock Screening Compounds" }
    self.pricedic = {}
    self.EmptyEntry = {
      "index":"", "can_smiles":"", 
      "mcule_found":"", "mcule_id":"", "mcule_price":"", "mcule_url":"", "mcule_smiles":"", "mcule_comment":"", 
      "cs_found":"", "cs_id":"", "cs_price":"", "cs_url":"", "cs_smiles":"", "cs_comment":"", 
    }
    self.csrf_token, self.cookie_csrf = self.gettoken();
    self.interval = 48; 

  def mculequery(self):
    headers = {'Authorization': f"Token {self.MCULE_TOKEN}"}
    for compi in range(self.compoundnr):
      thissmiles = self.compounds[compi];
      comp = urllib.parse.quote(self.compounds[compi], safe=''); 
      response = requests.get(f'https://mcule.com/api/v1/search/lookup/?query={comp}', headers=headers); 
      if response.status_code == 200: 
        # Query success - Stage 1: query the MCULE-ID
        retdic = json.loads(response.text)
        if len(retdic["results"]) == 0: 
          print(f"MCULE: Compound {compi+1} found no hits")
          self.pricedic[thissmiles]["mcule_comment"] = f"MCULE: Compound {compi+1} found no hits"; 
          self.pricedic[thissmiles]["mcule_found"] = False; 
          self.pricedic[thissmiles]["mcule_id"] = "Empty";
          self.pricedic[thissmiles]["mcule_price"] = "Empty"
          continue 
        elif len(retdic["results"]) > 1: 
          retnr = len(retdic["results"])
          print(f"MCULE: Compound {compi+1} found multiple ({retnr}) hits")
          self.pricedic[thissmiles]["mcule_comment"] = f"MCULE: Compound {compi+1} found multiple ({retnr}) hits"
          
        self.pricedic[thissmiles]["mcule_found"] = True; 
        self.pricedic[thissmiles]["mcule_smiles"] = retdic["results"][0]["smiles"]; 
        self.pricedic[thissmiles]["mcule_id"]  = molID = retdic["results"][0]["mcule_id"]; 
        self.pricedic[thissmiles]["mcule_url"] = retdic["results"][0]["url"]
        
        response2 = requests.get(f"https://mcule.com/api/v1/compound/{molID}/prices/", headers=headers)
        if response2.status_code == 200:
          # Query success - Stage 2: query the molecule price TODO: purity
          pricedic = json.loads(response2.text)
          if "best_prices" in pricedic.keys() and len(pricedic["best_prices"]) > 0:
            self.pricedic[thissmiles]["mcule_price"] = []
            for idx, i in enumerate(pricedic["best_prices"]):
              if "currency" in i.keys():
                thecurrency = i["currency"]; 
              else: 
                thecurrency = "USD"; 
              if "unit" in i.keys():
                theunit = i["unit"]; 
              else: 
                theunit = "mg"; 
              if "delivery_time_working_days" in i.keys():
                thetime = str(i["delivery_time_working_days"])+" days"; 
              else:
                thetime = "N/A"; 
              if "price" in i.keys():
                theprice = str(i["price"])+f" {thecurrency}"; 
              else: 
                theprice = "N/A"; 
              if "amount" in i.keys():
                theamount = str(i["amount"])+f" {theunit}"; 
              else:
                theamount = "N/A"; 
              self.pricedic[thissmiles]["mcule_price"].append({"price":theprice , "amount":theamount , "time":thetime, "purity": "", "supplier":"MCULE"})
          elif "best_prices" in pricedic.keys() and len(pricedic["best_prices"]) == 0: 
            self.pricedic[thissmiles]["mcule_price"] = "Empty"
          else: 
            self.pricedic[thissmiles]["mcule_comment"] += f"MCULE Query Stage 2 success, however, no price found"
        else: 
          print(f"MCULE: Price query failed ({molID}): ", response2.status_code, response2.reason, response2.url)
          self.pricedic[thissmiles]["mcule_comment"] += f"MCULE Query Stage 2 failed: StatusCode {response.status_code}; Reason: {response.reason};"
          
      else: 
        # Stage 1 query failed 
        print(f"MCULE: Compound {compi+1} query failed: StatusCode {response.status_code}, Reason: {response.reason}; URL: {response.url}")
        self.pricedic[thissmiles]["mcule_comment"] = f"MCULE Query Stage 1 Failed: StatusCode {response.status_code}; Reason: {response.reason};"
        self.pricedic[thissmiles]["mcule_found"] = False; 
        self.pricedic[thissmiles]["mcule_id"] = "Failed";
        self.pricedic[thissmiles]["mcule_price"] = "Empty"; 
  
  def csquery(self):
    """
        Update 13th June: the Chemspace do not reconize quoted structure in requests
    """
    cookies = { '_csrf': self.cookie_csrf }
    headers = {'X-CSRF-Token': self.csrf_token, 'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8', }
    csrf_quote = urllib.parse.quote(self.csrf_token, safe=''); 

    for compi in range(self.compoundnr):
      thissmiles = self.compounds[compi];
      comp_str = urllib.parse.quote(self.compounds[compi], safe=''); 
      data = f'_csrf={csrf_quote}&search={comp_str}'; 
      response = requests.post('https://chem-space.com/search/text', cookies=cookies, headers=headers, data=data); 
      if response.status_code == 200: 
        # Query success - Stage 1: query the ChemSpace ID or List 
        retdic = json.loads(response.text)
        if "message" in retdic.keys():
          retmsg = retdic['message']; 
          print(f"Chem-Space: Compound {compi+1} found no hits"); 
          self.pricedic[thissmiles]["cs_found"] = False; 
          self.pricedic[thissmiles]["cs_comment"] = f"Chem-Space: Compound {compi+1} found no hits; Return message:{retmsg}; "; 
          continue
        if "url" not in retdic.keys():
          retmsg = retdic['message'];
          print(f"Chem-Space: Compound {compi+1} found no hits");
          self.pricedic[thissmiles]["cs_found"] = False;
          self.pricedic[thissmiles]["cs_comment"] = f"Chem-Space: Compound {compi+1} found no hits; URL not returned; ";
          continue
        returl = retdic["url"]; 
        if "list" in returl:
          # Multiple entries returned / multiple pages
          finalurl    = f"{self.CSORIGIN}{returl}?per_page={self.interval}&view=list&page=1"; 
          result_page = requests.get(finalurl); 
          thepage = BeautifulSoup(result_page.content, 'html.parser'); 
          # self.savepage(result_page, f"/tmp/chemspace_list{compi+1}.html"); 
          
          self.CSCURRENCY = thepage.find("ul", {"id":"currency"}).find("strong").text.upper(); 
          try:
            foundNr = thepage.find("span", {"class":"found"}).find_all("b")[-1].getText(); 
            foundNr = int(foundNr)
            cs_ids = [i.find("a", {"data-role":"enquire-btn"})["data-cs_id"] for i in thepage.find_all("dl", {"class":"search-result-item"})]; 
          except:
            print("Chem-Space: Error: The following url returned a list however, no molecule entry found: ", finalurl); 
            continue; 
          
          if foundNr > 0:
            # Use the first ChemSpace-ID as the field entry
            self.pricedic[thissmiles]["cs_found"] = True; 
            self.pricedic[thissmiles]["cs_id"] = cs_ids[0];
          if int(foundNr) > 1:
            print(f"Chem-Space: Compound {compi+1} found {foundNr} entries:", ",".join(cs_ids)); 
            cs_ids_str = ",".join(cs_ids)
            self.pricedic[thissmiles]["cs_comment"] += f"Found {foundNr} entries in ChemSpace database including: {cs_ids_str}"

          mollist = thepage.find_all("dl", {"class":"search-result-item"}); 
          themolinfo = mollist[0].find("dt", {"class":"item-figure"}); 
          csid = self.pricedic[thissmiles]["cs_id"] = themolinfo.find("input", {"class":"idnumber"})["data-mark"]; 
          self.pricedic[thissmiles]["cs_url"] = f"{self.CSORIGIN}/{csid}"; 
          
          raw_smiles = [i["href"].split("=")[1] for i in themolinfo.find_all("a") if "smiles" in str(i)]; 
          thesmile = self.pricedic[thissmiles]["cs_smiles"] = urllib.parse.unquote(raw_smiles[0]); 
          
          self.pricedic[thissmiles]["cs_price"] = []
          for mol in mollist:
            table = mol.find("table", {"class":"offers-table"}); 
            if table: 
              columns = [i.getText() for i in table.find("thead").find_all("th")]; 
              for i in range(len(columns)):
                if "supplier" in columns[i]:
                  supplcol = i; 
                else:
                  supplcol = 0; 
                if "time"     in columns[i]:
                  timecol = i; 
                else: 
                  timecol = 1; 
                if "price"    in columns[i]:
                  pricecol = i; 
                else:
                  pricecol = 2;
                if "pack"     in columns[i]:
                  packcol = i; 
                else: 
                  packcol = 3; 

              #   {"time":thetime, "price":theprice, "amount":thepack, "purity": thepurity, "supplier": thesuppl}
              rows = table.find("tbody").find_all("tr", {"data-role": "cart-item-container"});
              for rowi in rows: 
                thesuppl  = rowi.find_all("td")[supplcol].find("div", {"class":"shorten-supplier-name"})["data-content"]; 
                try:
                  thetime = rowi.find_all("td")[timecol]["data-value"];
                  thetime = int(thetime);
                  if thetime >= 365:
                    raise
                  thetime = f"{thetime} days"; 
                except:
                  # "TBD"
                  thetime = rowi.find_all("td")[timecol].getText().strip(); 

                # Package size condition is either a selection element or a simple td element
                if rowi.find_all("td")[packcol].find("select", {"data-role":"pack-switcher"}):
                  # Select as price field 
                  for option in rowi.find_all("td")[packcol].find("select", {"data-role":"pack-switcher"}).find_all("option"):
                    thepack  = json.loads(option["data-pack_size"])["text"]; 
                    if "POA" in option["data-price"]:
                      theprice = "POA"
                    else: 
                      theprice = json.loads(option["data-price"])["value"].__str__()+f" {self.CSCURRENCY}"; 
                    self.pricedic[thissmiles]["cs_price"].append({"time":thetime, "price":theprice, "amount":thepack, "supplier": thesuppl})                  
                else:
                  # Plain td cell as price field 
                  thepack  = rowi.find_all("td")[packcol].getText().strip(); 
                  if "POA" in str(rowi.find_all("td")[pricecol]):
                    theprice = rowi.find_all("td")[pricecol].getText().strip(); 
                  else: 
                    theprice = rowi.find_all("td")[pricecol].getText().strip() + f" {self.CSCURRENCY}"; 
                  self.pricedic[thissmiles]["cs_price"].append({"time":thetime, "price":theprice, "amount":thepack, "supplier": thesuppl})
            else:
              pass
        elif "/CS" in returl:
          # Only One Molecule entry returned
          finalurl    = f"{self.CSORIGIN}{returl}"; 
          result_page = requests.get(finalurl); 
          bs = BeautifulSoup(result_page.content, 'html.parser')
          self.pricedic[thissmiles]["cs_found"] = True; 
          self.pricedic[thissmiles]["cs_id"]    = bs.find("div", {"class":"scheme"})["data-cscid"]
          self.pricedic[thissmiles]["cs_url"]   = finalurl; 
          self.pricedic[thissmiles]["cs_smiles"] = bs.find("li", {"class":"str-format", "data-format":"SMILES", "data-target-id":"smiles"})["data-value"]
          self.CSCURRENCY = bs.find("ul", {"id":"currency"}).find("strong").text.upper(); 
          #self.savepage(result_page, f"/tmp/chemspace_comp{compi+1}.html")

          if bs.find("div", {"class":"empty"}):
            # No price / provider available
            self.pricedic[thissmiles]["cs_price"] = "Empty"
          elif bs.find("table", {"class":"structureSuppliersList", "id":"suppliersList"}):
            # Available price / provider 
            columns = [i.find("a") for i in  bs.find("table", {"class":"structureSuppliersList"}).find_all("th")]
            columns = [i.getText().lower() for i in columns if i]
            self.pricedic[thissmiles]["cs_price"] = []
            for i in range(len(columns)): 
              if "supplier" in columns[i]:
                supplcol = i; 
              else:
                supplcol = 0; 
              if "purity"   in columns[i]:
                puritycol = i; 
              else:
                puritycol = 3; 
              if "pack"     in columns[i]:
                packcol = i; 
              else: 
                packcol = 4; 
              if "price"    in columns[i]:
                pricecol = i; 
              else:
                pricecol = 5;
              if "time"     in columns[i]:
                timecol = i; 
              else: 
                timecol = 1; 
            carttags = bs.find("tbody").find_all("tr", {"data-role": "cart-item-container"})
            for cart in carttags: 
              thecells = [i.text for i in cart.find_all("td")]; 
              thesuppl  = thecells[supplcol]
              thepurity = thecells[puritycol]
              thepack   = thecells[packcol]
              theprice  = thecells[pricecol]+f" {self.CSCURRENCY}"
              thetime   = thecells[timecol]
              self.pricedic[thissmiles]["cs_price"].append({"time":thetime, "price":theprice, "amount":thepack, "purity": thepurity, "supplier": thesuppl})
      else: 
        # Stage 1 query failed 
        print(f"ChemSpace: Compound {compi+1} query failed: StatusCode {response.status_code}, Reason: {response.reason}; URL: {response.url}")
        self.pricedic[thissmiles]["cs_found"] = False; 
        self.pricedic[thissmiles]["cs_comment"] = f"ChemSpace Query Stage 1 Failed: StatusCode {response.status_code}; Reason: {response.reason};"; 
      if (len(self.pricedic[thissmiles]["cs_price"]) == 0):
        # If no prices found, change it to stirng Empty
        self.pricedic[thissmiles]["cs_price"] = "Empty"


  def matchquery(self, retlist): 
    print(f"Matching {len(retlist)} return compounds with {len(self.compounds)} initial compound entries")
    values = []
    print(retlist)
    for i in retlist:
      found = 0
      if i == "None": 
        values.append("None")
        continue
      for j in self.compounds: 
        mol1 = Chem.MolFromSmarts(j)
        mol2 = Chem.MolFromSmarts(i)
        try: 
          have_match = mol1.HasSubstructMatch(mol2)
          samenha = (mol1.GetNumHeavyAtoms() == mol2.GetNumHeavyAtoms())
          if have_match and samenha:
            values.append(j)
            found =1
            break
        except: 
          continue
      if found == 0:
        values.append("None")
    if len(values) != len(retlist):
      raise Exception(f"Length of returned compound list is not aligned {len(values)}, {len(retlist)}")
    return values

  def querymols(self, compounds):
    self.compounds = compounds; 
    self.compoundnr = len(compounds);
    self.pricetable = pd.DataFrame(); 
    self.pricetable["index"] = range(len(self.compounds)); 
    self.pricetable["can_smile"] = self.compounds; 
    
    for idx,mol in zip(range(self.compoundnr), self.compounds):
      self.pricedic[mol] = json.loads(json.dumps(self.EmptyEntry))
      self.pricedic[mol]["index"] = idx+1; 
      self.pricedic[mol]["can_smiles"] = mol; 

    self.mculequery();  # Collect MCULE 
    self.csquery();     # Collect from ChemSpace 


  def getmol(self, struct, fomat):
    json_data = {
      'structure': struct,
      'parameters': fomat,
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
    if response.status_code == 200:
      if fomat == "inchi":
        return json.loads(response.text)["structure"].split('\n')[0]
      else :
        return json.loads(response.text)["structure"]
    else:
      print("Failed to interpret the molecule")
      print(response.status_code, response.reason, response.url)


  def savepage(self, response, path):
    with open(path, "w") as file1: file1.write(BeautifulSoup(response.content, 'html.parser').prettify())

  def gettoken(self):
    ret = requests.get("https://chem-space.com")
    cookie_csrf = ret.cookies.get_dict()["_csrf"]; 
    if ret.status_code != 200: print("Failed to get the csrf values:", ret.url, ret.status_code, ret.reason)
    bs = BeautifulSoup(ret.content, 'html.parser')
    metainfo = bs.find_all("meta")
    for i in metainfo:
      if i.has_attr("name"):
        if i["name"] == "csrf-token":
          csrf_token = i["content"]
    return csrf_token, cookie_csrf
  def gethtml(self):
    finestr=""
    headers=["Smiles", "MCULE-ID", "MCULE-prices", "ChemSpace-ID", "ChemSpace-prices"]
    data = []
    for i in range(self.compoundnr):
      smi_i = self.compounds[i]
      smiles  = self.pricedic[smi_i]["can_smiles"]
      mid     = self.pricedic[smi_i]["mcule_id"]
      if len(self.pricedic[smi_i]["mcule_price"]) == 0 or isinstance(self.pricedic[smi_i]["mcule_price"], str):
        mprice  = self.pricedic[smi_i]["mcule_price"]
      else:
        mprice  = "\n".join([f"{j['price']}/{j['amount']}/{j['time']}" for j in self.pricedic[smi_i]["mcule_price"]])
      csid    = self.pricedic[smi_i]["cs_id"]
      if len(self.pricedic[smi_i]["cs_price"]) == 0 or isinstance(self.pricedic[smi_i]["cs_price"], str):
        csprice  = self.pricedic[smi_i]["cs_price"]
      else:
        csprice = "\n".join([f"{j['price']}/{j['amount']}/{j['time']}" for j in self.pricedic[smi_i]["cs_price"]])
      data.append([smiles, mid, mprice,csid,csprice])
      finestr += f"{smiles},{mid},{mprice},{csid},{csprice}\n"; 
    pdtable = pd.DataFrame(data, columns=headers) 
    
    return pdtable.to_html().replace("\\n", "<br>")
  def getcsv(self):
    finestr="Smiles, MCULE-ID, MCULE-prices, ChemSpace-ID, ChemSpace-prices"; 
    for i in range(self.compoundnr):
      smi_i = self.compounds[i]
      smiles  = self.pricedic[smi_i]["can_smiles"]
      mid     = self.pricedic[smi_i]["mcule_id"]
      if len(self.pricedic[smi_i]["mcule_price"]) == 0 or isinstance(self.pricedic[smi_i]["mcule_price"], str):
        mprice  = self.pricedic[smi_i]["mcule_price"]
      else:
        mprice  = "\n".join([f"{j['price']}/{j['amount']}/{j['time']}" for j in self.pricedic[smi_i]["mcule_price"]])
      csid    = self.pricedic[smi_i]["cs_id"]
      if len(self.pricedic[smi_i]["cs_price"]) == 0 or isinstance(self.pricedic[smi_i]["cs_price"], str):
        csprice  = self.pricedic[smi_i]["cs_price"]
      else:
        csprice = "\n".join([f"{j['price']}/{j['amount']}/{j['time']}" for j in self.pricedic[smi_i]["cs_price"]])
      finestr += f"{smiles},{mid},{mprice},{csid},{csprice}\n"; 
    return finestr
  def writejson(self, fileout):
    with open(fileout,"w") as file1: 
      json.dump(self.pricedic, file1, indent=2)
