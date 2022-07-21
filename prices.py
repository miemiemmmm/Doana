import json
import requests
import urllib.parse
from bs4 import BeautifulSoup 
# from rdkit import Chem
import pandas as pd 

import sys

    
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
        self.csrf_token, self.cookie_csrf = self.gettoken();
        self.cspages  = []; 
        self.interval = 24
    
    def mculequery(self):
        MCULEfound = []; MCULEsmiles = []; MCULEid = []; MCULEurl = []; MCULEprices = [];
        headers = {'Authorization': f"Token {self.MCULE_TOKEN}"}
        for compi in range(self.compoundnr):
            comp = urllib.parse.quote(self.compounds[compi], safe='')
            response = requests.get(f'https://mcule.com/api/v1/search/lookup/?query={comp}', headers=headers) 
            if response.status_code == 200: 
                retdic = json.loads(response.text)
                if len(retdic["results"]) == 0: 
                    print(f"MCULE: Compound {compi+1} found no hits")
                    MCULEfound.append("Empty"); 
                    MCULEsmiles.append("None"); 
                    MCULEid.append("None"); 
                    MCULEurl.append("None"); 
                    MCULEprices.append("None")
                    continue 
                molsmi = retdic["results"][0]["smiles"]
                molID  = retdic["results"][0]["mcule_id"]
                molurl = retdic["results"][0]["url"]
                response2 = requests.get(f"https://mcule.com/api/v1/compound/{molID}/prices/", headers=headers)
                if response2.status_code == 200:
                    pricedic = json.loads(response2.text)
                    if "best_prices" in pricedic.keys() and len(pricedic["best_prices"]) > 0:
                        finalstr = f"Supplier MCULE: "
                        for idx, i in enumerate(pricedic["best_prices"]):
                            if "price" in i.keys() and "currency" in i.keys(): 
                                theprice = i["price"]
                                thecurrency = i["currency"]
                                finalstr += f"{theprice} {thecurrency}-"
                            if "amount" in i.keys() and "unit" in i.keys(): 
                                theamount = i["amount"]; 
                                theunit   = i["unit"];
                                finalstr += f"{theamount} {theunit}-"
                            if "delivery_time_working_days" in i.keys():
                                thetime = i["delivery_time_working_days"]
                                finalstr += f"{thetime}days"
                            if (idx < len(pricedic["best_prices"])-1):
                                finalstr += "_"
                        MCULEfound.append("Found"); 
                        MCULEsmiles.append(molsmi); 
                        MCULEid.append(molID); 
                        MCULEurl.append(molurl); 
                        MCULEprices.append(finalstr)
                    elif "best_prices" in pricedic.keys() and len(pricedic["best_prices"]) == 0: 
                        MCULEfound.append("Found"); 
                        MCULEsmiles.append(molsmi); 
                        MCULEid.append(molID); 
                        MCULEurl.append(molurl); 
                        MCULEprices.append("Empty")
                    else: 
                        MCULEfound.append("Found"); 
                        MCULEsmiles.append(molsmi); 
                        MCULEid.append(molID); 
                        MCULEurl.append(molurl); 
                        MCULEprices.append("Wrong price query")
                else: 
                    print(f"MCULE: Price query failed ({molID}): ", response2.status_code, response2.reason, response2.url)
                    MCULEfound.append("Found"); 
                    MCULEsmiles.append(molsmi); 
                    MCULEid.append(molID); 
                    MCULEurl.append(molurl); 
                    MCULEprices.append("Failed")
            else: 
                print(f"MCULE: Compound {compi} query failed: ", response.status_code, response.reason, response.url)
                MCULEfound.append("Failed"); 
                MCULEsmiles.append("None"); 
                MCULEid.append("None"); 
                MCULEurl.append("None"); 
                MCULEprices.append("None")

        self.pricetable["mcule_smile"] = MCULEsmiles
        self.pricetable["mcule_found"] = MCULEfound
        self.pricetable["mcule_id"]    = MCULEid
        self.pricetable["mcule_url"]   = MCULEurl
        self.pricetable["mcule_price"] = MCULEprices
        return None
            
    def csquery(self):
        """
            Update 13th June: the Chemspace do not reconize quoted structure in requests
        """
        cookies = { '_csrf': self.cookie_csrf }
        headers = {'X-CSRF-Token': self.csrf_token, 
                   'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8', }
        csrf_quote = urllib.parse.quote(self.csrf_token, safe=''); 
        
        CSfound  = []; 
        CSsmiles = []; 
        CSid     = []; 
        CSurl    = []; 
        CSprices = [];
        
        for compi in range(self.compoundnr):
            comp_str = urllib.parse.quote(self.compounds[compi], safe='')
            data = f'_csrf={csrf_quote}&search={comp_str}'
            retquery = requests.post('https://chem-space.com/search/text', cookies=cookies, headers=headers, data=data)
            
            if retquery.status_code == 200: 
                retdic = json.loads(retquery.text)
                
                if "message" in retdic.keys() and (retdic['message'] == "Chemicals not found"):
                    print(f"Chem-Space: Compound {compi+1} found no hits"); 
                    CSfound.append("Empty")
                    CSsmiles.append("None")
                    CSid.append("None")
                    CSurl.append("None")
                    CSprices.append("None")
                    continue
                returl = retdic["url"]; 
                if "list" in returl:
                    # Multiple entries returned / multiple pages
                    finalurl    = f"{self.CSORIGIN}{returl}?per_page={self.interval}&view=plate&page=1"; 
                    result_page = requests.get(finalurl); 
                    thepage = BeautifulSoup(result_page.content, 'html.parser')
                    self.CSCURRENCY = thepage.find("ul", {"id":"currency"}).find("strong").text.upper(); 
                    # self.savepage(result_page, f"/tmp/chemspace{compi+1}.html")
                    cspages = [finalurl]; 
                    pageinfo = thepage.find("ul", {"class":"port"})
                    foundNr = thepage.find("span", {"class":"found"}).find_all("b")[-1].getText()
                    smilst  = []
                    csidlst = []
                    urllst  = []
                    thissupp = {}
                    mollist = thepage.find_all("dl", {"class":"search-result-item"}); 
                    for i,j in zip(range(int(foundNr)), mollist):
                        # Different chemcial identities of the same molecule
                        # i.e. Different salt form
                        csid = j.find("input", {"class":"idnumber"})["data-mark"]
                        tmpatags = j.find_all("a")
                        thesmile = [urllib.parse.unquote(x["href"].split("=")[-1]) for x in tmpatags if (x.has_attr("href") and "smiles" in x["href"])][0]
                        category = j.find("div", {"class":"category-type"})["class"][-1]
                        carttags = j.find_all("tr", {"data-role": "cart-item-container"})
                        smilst.append(thesmile)
                        csidlst.append(csid)
                        urllst.append(f"{self.CSORIGIN}/{csid}")
                        thissupp[csid] = {}
                        for cartidx, cart in enumerate(carttags):
                            # multiple suppliers in one cart
                            # Only packages size and price vary 
                            leadtime = cart.find("td", {"data-role":"lead_time"}).text.strip()
                            supplier = cart.find("div", {"class": "shorten-supplier-name"}).text.strip()
                            price_default = cart.find("td", {"data-role":"price"}).text.strip()
                            psize_default = cart.find("td", {"data-role":"pack_size"}).text.strip()
                            thissupp[csid][f"supplier{cartidx+1}"] = supplier
                            psize_selector = cart.find("select", {"data-role":"pack-switcher"})
                            if psize_selector:
                                # Found package selector which means it has multiple prices. 
                                options  = psize_selector.find_all("option")
                                supplier_price = ""
                                for o in options: 
                                    if o.has_attr("data-pack_size"):
                                        supplier_price += (json.loads(o["data-pack_size"])["text"])
                                    if o.has_attr("data-price"):
                                        supplier_price += ("-" + json.loads(o['data-price'])["value"].__str__() + self.CSCURRENCY)
                                    if o.has_attr("data-lead_time"):
                                        supplier_price += ("-" + json.loads(o['data-lead_time'])["value"].__str__()+" days")
                                    if len(supplier_price) > 2 and (options.index(o)+1 < len(options)):
                                        supplier_price += "_"
                            else: 
                                # Not found package selector which means it only have one price. 
                                supplier_price = f"{psize_default}-{price_default}-{leadtime}"
                            if len(supplier_price) == 0:
                                supplier_price = "None"; 
                            thissupp[csid][f"prices{cartidx+1}"] = supplier_price
                    supplierstr = ""
                    for moli in thissupp.keys():
                        priceinds = [int(k.replace("supplier", '')) for k in thissupp[moli].keys() if "supplier" in k]
                        pricestr = ""
                        for ind in priceinds:
                            pricestr += thissupp[moli][f"supplier{ind}"]+": "+thissupp[moli][f"prices{ind}"]
                            if priceinds.index(ind) < len(priceinds)-1:
                                pricestr += "_"
                        supplierstr += pricestr
                        if [_ for _ in thissupp.keys()].index(moli) < len([_ for _ in thissupp.keys()])-1:
                            supplierstr += "_"
                    if len(csidlst) > 0: 
                        CSfound.append("Found")
                        CSsmiles.append("_".join(smilst))
                        CSid.append("_".join(csidlst))
                        CSurl.append("_".join(urllst))
                        CSprices.append(supplierstr)    #thissupp.__str__()
                    else: 
                        CSfound.append("Empty")
                        CSsmiles.append("None")
                        CSid.append("None")
                        CSurl.append("None")
                        CSprices.append("None")
                    
                elif "/CS" in returl:
                    # Only One Molecule entry returned
                    finalurl    = f"{self.CSORIGIN}{returl}"; 
                    result_page = requests.get(finalurl); 
                    bs = BeautifulSoup(result_page.content, 'html.parser')
                    self.CSCURRENCY = bs.find("ul", {"id":"currency"}).find("strong").text.upper(); 
                    csid = bs.find("div", {"class":"scheme"})["data-cscid"]
                    smival = bs.find("li", {"class":"str-format", "data-format":"SMILES", "data-target-id":"smiles"})["data-value"]
                    self.savepage(result_page, f"/tmp/comp1_chemspace{compi+1}.html")
                    
                    if bs.find("div", {"class":"empty"}):
                        # Returned Empty result AVAILABILITY
                        CSfound.append("Empty")
                        CSsmiles.append("None")
                        CSid.append("None")
                        CSurl.append("None")
                        CSprices.append("None")
                        
                    elif bs.find("table", {"class":"structureSuppliersList", "id":"suppliersList"}):
                        carttags = bs.find_all("tr", {"data-role": "cart-item-container"})
                        thissupp = ""
                        self.cssuppliermap ={}
                        suppliers = []
                        for cart in carttags:
                            supplier = [i for i in cart.find_all("span") if i.has_attr("data-supplier_id")][0]; 
                            supplierid = supplier["data-supplier_id"];
                            suppliertxt = supplier.text; 
                            
                            self.cssuppliermap[supplierid] = suppliertxt; 
                            if supplierid not in suppliers:
                                suppliers.append(supplierid); 
                                thissupp += (suppliertxt+": "); 
                            tdinfo = cart.find_all("td"); 
                            thetime  = tdinfo[1].text; 
                            thepack  = tdinfo[3].text; 
                            theprice = tdinfo[4].text; 
                            if len(thetime) > 0:
                                thissupp += f"{thepack}-{theprice}{self.CSCURRENCY} - {thetime}"; 
                            else :
                                thissupp += f"{thepack}-{theprice}{self.CSCURRENCY}"; 

                            if len(thissupp) > 2 and (carttags.index(cart)+1 < len(carttags)):
                                thissupp += "_"
                        CSfound.append("Found")
                        CSsmiles.append(smival)
                        CSid.append(csid)
                        CSurl.append(finalurl)
                        CSprices.append(thissupp)
            else: 
                CSfound.append("Failed")
                CSsmiles.append("None")
                CSid.append("None")
                CSurl.append("None")
                CSprices.append("None")
        self.pricetable["cs_price"] = CSprices
        self.pricetable["cs_smile"] = CSsmiles
        self.pricetable["cs_found"] = CSfound
        self.pricetable["cs_id"]    = CSid
        self.pricetable["cs_url"]   = CSurl
        return
                
            
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
    def formatprint(self):
        finestr=""
        for i in range(len(self.pricetable["can_smile"])):
            smiles  = self.pricetable["can_smile"][i]
            mid     = self.pricetable["mcule_id"][i]
            mprice  = self.pricetable["mcule_price"][i]
            csid    = self.pricetable["cs_id"][i]
            csprice = self.pricetable["cs_price"][i]
            print(f"{smiles}, {mid},{mprice},{csid},{csprice}")
            finestr += f"{smiles}, {mid},{mprice},{csid},{csprice}\n"; 
        return finestr
    


