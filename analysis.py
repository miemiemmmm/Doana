import prolif as plf
import numpy as np 
import matplotlib as plt 
import pytraj as pt 

import sys
import random
import pickle
import MDAnalysis as mda
import rdkit.Chem as Chem
from scipy.spatial import distance_matrix
import os

class PLIFGen:
    def __init__(self, parmdic):
        # keys # ['ligfile', 'profile', 'resultdat', 'resultmol2']
        self.parms = parmdic;

        # Clean the lone pair and dummy atoms in mol2 if necessary;
        if ("mol2" in parmdic["resultmols"]) and ("cleanscript" in parmdic.keys()):
            self.CLEANMOL2SCRIPT = parmdic["cleanscript"]
            self.ligandfile = self.clean_mol2(parmdic["resultmols"])
        elif ("mol2" in parmdic["resultmols"]) and ("cleanscript" not in parmdic.keys()):
            self.CLEANMOL2SCRIPT = "/home/miemie/Dropbox/Documents/PortableWork/myscripts/clean_SEEDmol2.sh"
            self.ligandfile = self.clean_mol2(parmdic["resultmols"])
        else:
            self.ligandfile = parmdic["resultmols"]

        self.FP_TYPE1 = ["HBDonor", "HBAcceptor", "PiStacking", "CationPi","Cationic"]
        self.FP_TYPE2 = ['Hydrophobic']
        self.VDWRADII = {'AA': 1.85, 'AG': 1.72, 'AL': 1.84, 'AR': 1.88, 'AT': 2.02, 'AU': 1.66, 'B': 1.92, 'BA': 2.68,
            'BE': 1.53, 'BI': 2.07, 'BR': 1.85, 'C': 1.7, 'CA': 2.31, 'CD': 1.58, 'CL': 1.75, 'Cl': 1.75,
            'CS': 3.43, 'CU': 1.4, 'F': 1.47, 'FR': 3.48, 'GA': 1.87, 'GE': 2.11, 'H': 1.1,
            'HE': 1.4, 'HH': 1.55, 'I': 1.98, 'IN': 1.93, 'K': 2.75, 'KR': 2.02, 'LI': 1.82, 'MG': 1.73,
            'N': 1.55, 'NA': 2.27, 'NE': 1.54, 'NI': 1.63, 'O': 1.52, 'P': 1.8, 'PB': 2.02, 'PD': 1.63,
            'PO': 1.97, 'PT': 1.75, 'RA': 2.83, 'RN': 2.2, 'RR': 3.03, 'S': 1.8, 'SB': 2.06, 'SE': 1.9,
            'SI': 2.1, 'SN': 2.17, 'SR': 2.49, 'TE': 2.06, 'TL': 1.96, 'U': 1.86, 'XE': 2.16, 'ZN': 1.39,
        }

        self.profile = parmdic["profile"]
        mda_prot = mda.Universe(self.profile, top=self.profile, guess_bonds=True, vdwradii=self.VDWRADII);
        elements = mda.topology.guessers.guess_types(mda_prot.atoms.names);
        mda_prot.add_TopologyAttr("elements", elements);
        prot = plf.Molecule.from_mda(mda_prot);

        self.failed_mol = []
        if '.sdf' in self.ligandfile:
            # SDF format
            self.lig_suppl, self.success_mol, self.failed_mol = self.sdf_supplier(self.ligandfile)
            if len(self.failed_mol) > 0: print(f"{len(self.failed_mol)} molecules failed in RDKit SDF loading:", *self.failed_mol)
        elif '.mol2' in self.ligandfile:
            # MOL2 format
            self.lig_suppl = [i for i in plf.mol2_supplier(self.ligandfile)]
            self.success_mol =[i for i in range(len(self.lig_suppl)) ]
            self.failed_mol = []
        print(f"{len(self.lig_suppl)} molecules put to the library supplier (from {self.ligandfile})")

        self.datfile = parmdic["resultdat"]
        # Generate protein-ligand interaction fingerprint
        fp_hpl = plf.Fingerprint(self.FP_TYPE1)
        fp_hpl.run_from_iterable(self.lig_suppl, prot)
        fp_data = fp_hpl.to_dataframe()
        self.fp_data = self.addinfo(fp_data)

        fp_hpb = plf.Fingerprint(self.FP_TYPE2)
        fp_hpb.run_from_iterable(self.lig_suppl, prot)
        fp_data2 = fp_hpb.to_dataframe()
        self.fp_data2 = self.addinfo(fp_data2)


    def clean_mol2(self, mol2file, verbose=False):
        # Will put the output file to the same directory of to mol2 file
        if verbose:
            osstr = f"bash {self.CLEANMOL2SCRIPT} {mol2file}";
        else :
            osstr = f"bash {self.CLEANMOL2SCRIPT} {mol2file} 2>&1 > /dev/null";
        os.system(osstr)
        sdffile = mol2file.replace(".mol2", ".sdf")
        cleaned_mol2file = mol2file.replace(".mol2", "_clean.mol2")
        if os.path.isfile(cleaned_mol2file):
            return cleaned_mol2file
        else:
            if os.path.isfile(sdffile):
                return sdffile
            else :
                print("Failed to clean mol2 file {mol2file}, using the unchanged mol2 file ")
                return mol2file

    def addinfo(self, df, topn=-1):
        # Add the Columns about the score and pos_id
        with open(self.datfile, "r") as file1: thelist = [i.strip('\n').split() for i in file1.read().split('\n')]
        score_table = [[i+1]+thelist[i] for i in self.success_mol]
        score_table    = np.array(score_table);
        col_pos_id     = score_table[:,1].astype(str);
        col_smiles     = score_table[:,2].astype(str);
        col_seed_total = score_table[:,3].astype(str);
        col_seed_vdw   = score_table[:,4].astype(str);
        col_ranking    = score_table[:,0].astype(int);
        col_nha        = score_table[:,9].astype(str);
        col_date       = score_table[:,10].astype(str);
        col_source_campaign   = score_table[:,11].astype(str);

        df["Rank"] = col_ranking
        df["seed_total"]   = col_seed_total;
        df["seed_vdw"] = col_seed_vdw;
        df["pos_id"]   = col_pos_id;
        df["nha"] = col_nha;
        df["date"] = col_date;
        df["source_campaign"] = col_source_campaign;
        df["can_smile"] = col_smiles;
        return df

    def sdf_supplier(self, path):
        suppl = Chem.SDMolSupplier(path, removeHs=False);
        indexes = [ i for i in range(suppl.__len__()) ]
        planned_indexes = []
        if ("onlymols" in self.parms.keys()) and (self.parms["onlymols"] != ""):
            if type(self.parms["onlymols"]) == int:
                planned_indexes = indexes[:self.parms["onlymols"]];
            else:
                planned_indexes = [i for i in self.parms["onlymols"]]
        else:
            planned_indexes = indexes;
        failed_mol = [];  mols = [];  success_mol = [];

        for i, mol in zip(indexes, suppl):
            if i in planned_indexes:
                try:
                    plfmol = plf.Molecule.from_rdkit(mol)
                    mols.append(plfmol)
                    success_mol.append(i)
                except:
                    failed_mol.append(i)
        return mols, success_mol, failed_mol;

    def docking_OI(self, refmol, dist_cutoff, printrecords=True):
        molfilename = refmol.__str__();
        try:
            if "sdf" in refmol:
                refmol, success, failed = self.sdf_supplier(refmol)
            elif "mol2" in refmol:
                refmol = [i for i in plf.mol2_supplier(refmol)]
                success = [1]
            refmol = Chem.RemoveHs(refmol)
            refxyz = [_.GetPositions() for _ in refmol.GetConformers()][0]
        except:
            print("Cannot Sanitize the molecule by RDKit. Only use the coordinates of heave atoms. \nTemporary file will be write to /tmp/testref.xyz. ")
            os.system(f"obabel {molfilename} -d -O /tmp/testref.xyz")
            with open("/tmp/testref.xyz", "r") as file1: refxyz = [i.split()[1:] for i in file1.read().strip("\n").split("\n")[2:]]
            refxyz = np.array(refxyz).astype(float)
            success = [1]

        if len(success) == 0:
            raise Exception(f"Error: Found {len(success)} molecules in the input {refmol}")
            return
        elif len(success) > 1:
            raise Exception(f"Error: Found {len(success)} molecules in the input {refmol}. Only the first molecule will be used.")
            refmol = refmol[0]
        else:
            refmol = refmol[0]

        if self.fp_data.__len__() != len(self.success_mol):
            raise Exception(f"Fingerprint dataframe 1 length ({self.fp_data.__len__()}) not equil to ligand supplier{len(self.success_mol)}")
        if self.fp_data2.__len__() != len(self.success_mol):
            raise Exception(f"Fingerprint dataframe 1 length ({self.fp_data.__len__()}) not equil to ligand supplier{len(self.success_mol)}")

        OI_results = [];  OI_ref  = [];  OI_avg  = [];  rnha = [];

        for i, resultmol in zip(range(len(self.lig_suppl)), self.lig_suppl):
            resultmol = Chem.RemoveHs(resultmol);
            testxyz = [_.GetPositions() for _ in resultmol.GetConformers()][0]

            nha_ratio = len(testxyz) / len(refxyz)

            distances = distance_matrix(refxyz, testxyz)
            distlt = distances < dist_cutoff
            count1 = np.count_nonzero(np.any(distlt, axis = 0))
            count2 = np.count_nonzero(np.any(distlt, axis = 1))
            rat1  = count1 / len(testxyz)
            rat2  = count2 / len(refxyz)
            r_avg = np.mean([rat1,rat2])

            rnha.append(nha_ratio)
            OI_results.append(rat1);
            OI_ref.append(rat2);
            OI_avg.append((rat1+rat2)/2)
            if printrecords: print(f"Rank: {self.fp_data.iloc[i].Rank.values[0]} | Ratio_t: {round(rat1,3)} | Ratio_r: {round(rat2,3)} | Ravg {round(r_avg,3)} \
| Rnha: {round(nha_ratio,3)} ({len(testxyz)}/{len(refxyz)})")
        self.fp_data["OI_avg"] = self.fp_data2["OI_avg"] = self.OI_avg = OI_avg
        self.fp_data["OI_ref"] = self.fp_data2["OI_ref"] = self.OI_ref = OI_ref
        self.fp_data["OI_results"] = self.fp_data2["OI_results"] = self.OI_results = OI_results

    def savepkl(self):
        dictosave={
            "failed_mol": self.failed_mol,
            "success_mol" : self.success_mol,
            "parms" : self.parms,
                  }
        try:
            dictosave['fp_data'] = self.fp_data;
            dictosave[self.FP_TYPE1] = self.FP_TYPE1;
        except:
            pass
        try:
            dictosave['fp_data2'] = self.fp_data2;
            dictosave[self.FP_TYPE2] = self.FP_TYPE2;
        except:
            pass

        with open(self.parms["outpkl"], "wb") as fileout:
            pickle.dump(dictosave, fileout)
            print("Saved the fingerprint to file:", self.parms["outpkl"])




# import pickle
# import numpy as np
# import hashlib
# import time
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import AllChem
from IPython.display import Image



class PLIFRead:
    def __init__(self, pklfile):
        with open(pklfile, "br") as file1: resultdic = pickle.load(file1);
        self.fpdata  = resultdic["fp_data"];
        self.fpdata2 = resultdic["fp_data2"];
        self.molidxs = np.array(resultdic["success_mol"]).astype(int);

    def show_cols(self):
        print("Data set 1")
        for i, j in enumerate(self.fpdata.columns):
            if "UNL1" in j: fptype = "_".join(j[1:])
            else: fptype = j[0]+"\t"
            if i%3-2 ==0: print(i,fptype, np.count_nonzero(self.fpdata[j]))
            else: print(i,fptype, np.count_nonzero(self.fpdata[j]), " | ", end="")
        print("Data set 2")
        for i, j in enumerate(self.fpdata2.columns):
            if "UNL1" in j: fptype = "_".join(j[1:])
            else: fptype = j[0]+"\t"
            if i%3-2 ==0: print(i,fptype, np.count_nonzero(self.fpdata2[j]))
            else : print(i,fptype, np.count_nonzero(self.fpdata2[j]) , " | ",end="")
        print("")

    def plifplot(self, theax, the_array, ticks=[], labelsize=12, rotation=-20):
        # Note: Vertically flip the array because top row in dataframe starts from the botton in figure
        flip_array = np.flip(the_array.to_numpy(), axis=0)
        theax.pcolormesh(flip_array, cmap='Greys')
        theax.axes.xaxis.set_ticks([])
        theax.axes.yaxis.set_ticks([])

        # Customize tick labels
        if len(ticks) >0:
            tick_pos = [i+0.5 for i in range(len(ticks))];
            theax.xaxis.set_major_formatter(ticker.NullFormatter());
            theax.xaxis.set_minor_locator(ticker.FixedLocator(tick_pos));
            theax.xaxis.set_minor_formatter(ticker.FixedFormatter(ticks));
            theax.tick_params(which='minor', axis="x", labelsize=labelsize, rotation=rotation);
        return theax

    def onlytopn(self, df, topn=8, summary=False):
        # Return the top N fingerprints
        true_count  = np.count_nonzero(df==True, axis = 0)  # convert to True/False status
        enum = np.array([[i,j] for i, j in enumerate(true_count)])
        enum_sort = enum[enum[:,1].argsort()[::-1]]
        ind = enum_sort[:topn,0]
        topnames = [df.columns[i] for i in ind]
        df = df.loc[:,topnames];
        percents = [ round(true_count[i]/len(self.molidxs),3) for i in ind]
        return df, ind, percents

    def sortby(self, fp):
        if isinstance(self.sel_cols, str) and self.sel_cols=="top":
            _, ind,_ = self.onlytopn(fp.UNL1)
            sort_columns = [fp.columns[i] for i in ind];
            fp = fp.sort_values(by=sort_columns , axis=0, ignore_index=True, ascending=False);
        elif len(self.sel_cols) > 0 and not isinstance(self.sel_cols, str):
            sort_columns = [fp.columns[i] for i in self.sel_cols];
            fp = fp.sort_values(by=sort_columns , axis=0, ignore_index=True, ascending=False);
            print("Sorting the table by the following fingerprint: \n", sort_columns);
        else:
            # sort by ranking by default (There is identical SEED value dueting rounded float)
            sort_columns = ["Rank"];
            fp = fp.sort_values(by=sort_columns , axis=0, ignore_index=True, ascending=True);
            print("Sorting the table by its ranking");
        fp = fp.reset_index(drop=True)
        return fp

    def summarize(self, ):
        print("Summary of the dataframe: ")
        print("    Fingerprint Number: {} ; Compound Number:{}".format(ColumsNr, RowNr))
        for i in range(len(topcounts)):
            print("    Column index: {}, Counts: {}, percent: {}, name: {}, type: {}".format(
                ret_df["index"][i], ret_df["Counts"][i], round(ret_df["Counts"][i]/RowNr,3), ret_df["resname"][i], ret_df["fptype"][i]))
        return ret_df

    def Docking_prolif(self, topn=8, sel_cols=[], labelsize = 13, rotation=0):
        # Get the DataFrame and its meta-information
        fig, ax = plt.subplots(2, 1, figsize=(12,12))
        self.sel_cols = sel_cols;
        # Sort values
        fpdata1 = self.fp1 = self.sortby(self.fpdata);
        FP1 = fpdata1.UNL1;
        # Display only top N fingerprints
        FP1, ind1, per1 = self.onlytopn(FP1, topn);
        # Get the tick labels
        FP1_ticks = [f"{j[0]}\n{j[1]}\n{per1[i]} fp {ind1[i]}" for i, j in enumerate(fpdata1.UNL1.columns[ind1])]
        # Plot the protein_ligand fingerprint
        ax[0] = self.plifplot(ax[0] , FP1, ticks=FP1_ticks, rotation=rotation);

        fpdata2 = self.fp2 = self.sortby(self.fpdata2);
        FP2 = fpdata2.UNL1;
        FP2, ind2, per2 = self.onlytopn(FP2, topn);
        FP2_ticks = [f"{j[0]}\n{j[1]}\n{per2[i]}  fp {ind2[i]}" for i,j in enumerate(fpdata2.UNL1.columns[ind2])]
        ax[1] = self.plifplot(ax[1] , FP2, ticks=FP2_ticks, rotation=rotation);
        return fig, ax
    def HBondFilter(self, operator, val):
        titles = [i for i in self.fpdata.columns]
        hbidxs = [i for i in range(len(titles)) if (('HBDonor' in titles[i]) or ('HBAcceptor' in titles[i])) ]
        # Convert the fingerprint to Numpy array and count the number of hydrogen bonds
        hbFingerprint = self.fpdata.iloc[:, hbidxs].to_numpy();
        hbCount = np.count_nonzero(hbFingerprint, axis=1);
        if operator == "gt": hbstatus = np.where(hbCount > val);
        elif operator == "ge": hbstatus = np.where(hbCount >= val);
        elif operator == "eq": hbstatus = np.where(hbCount == val);
        elif operator == "le": hbstatus = np.where(hbCount <= val);
        elif operator == "lt": hbstatus = np.where(hbCount < val);
        print(f"Hydrogen Bond Count Selector: {len(hbstatus[0])} poses are kept")
        return hbstatus[0]
    def colFilter(self, dataset, colNr, operator, val):
        thecol = [i for i in dataset.columns][colNr];
        colvalues = dataset[thecol].astype(float);
        if operator == "gt": status = np.where(colvalues > val);
        elif operator == "ge": status = np.where(colvalues >= val);
        elif operator == "eq": status = np.where(colvalues == val);
        elif operator == "le": status = np.where(colvalues <= val);
        elif operator == "lt": status = np.where(colvalues < val);
        print(f"Column Value Selector: {len(status[0])} poses are kept")
        return status[0]

    def FPSelectionFilter(self, fpdata, fp_sel, operator="any"):
        # Available operator, all and any,
        selcols = [fpdata.columns[i] for i in fp_sel]
        fp_status = (fpdata.iloc[:,fp_sel] == True).to_numpy()
        if operator == "any":
            fp_status = [np.any(i) for i in fp_status]
        if operator == "all":
            fp_status = [np.all(i) for i in fp_status]
        status = np.where( np.array(fp_status) == True )[0]
        print(f"Fingerprint Selector: {len(status)} poses are kept")
        return status

    def getIntersection(self, *args):
        intersec = set(args[0]);
        if len(args) > 1 :
            for i in args[1:]:
                intersec = intersec.intersection(set(i))
        print(f"Intersection Operator: Source from {len(args)} datasets, {len(intersec)} poses are kept")
        return [i for i in intersec]

    def getPOSIDByStatus(self, status):
        return self.fpdata.iloc[status, :].loc[:,"pos_id"].to_numpy().astype(str)
    def getPOSIDQuery(self, posids):
        pos_idstr = ",".join([str(i) for i in posids])
        return f"pos_id in ({pos_idstr})"
    def getPropByPOSID(self, pos_ids):
        result = np.array([ self.fpdata[self.fpdata['pos_id'] == i].to_numpy()[0][-11:] for i in pos_ids]).astype(str)
        propStr = "\n".join([", ".join(i) for i in result])
        return propStr
    def getSmiByPOSID(self, pos_ids):
        result = [ self.fpdata[self.fpdata['pos_id'] == i]["can_smile"].values.tolist()[0] for i in pos_ids]
        return [i for i in result]
    def getDetailsByPOSID(self,fpdata, pos_ids):
        rows = [fpdata[fpdata['pos_id'] == i] for i in pos_ids]
        for i in rows:
            heads = i.columns[(i == True).to_numpy().squeeze()].tolist()
            print(i.loc[:, heads])

    def drawMOLSByPOSID(self,pos_ids, mpr=5):
        result = np.array([ self.fpdata[self.fpdata['pos_id'] == i]["can_smile"].to_numpy()[0] for i in pos_ids]).astype(str)
        result = [i.replace("*","") for i in result]
        mols = [Chem.MolFromSmiles(m) for m in result]
        for m in mols: tmp = AllChem.Compute2DCoords(m)
        img = Draw.MolsToGridImage(mols, molsPerRow=mpr, subImgSize=(400,400), maxMols=100,
                                   legends=[f"pos_id: {i}" for i in pos_ids],returnPNG=False)
        img.save("/tmp/molgrid.png")
        return Image("/tmp/molgrid.png")
    def drawHist(self, dataframe, col, n_bins=10):
        dataset = dataframe[dataframe.columns[col]].to_numpy().astype(float);
        fig, ax = plt.subplots(1, 1,figsize=(10, 5), sharey=True, tight_layout=True)
        N, bins, patches = ax.hist(dataset, bins=n_bins)
        return fig, ax, N, bins, patches
    def matchSubStructure(self, poslst, structstr):
        substruct = Chem.MolFromSmarts(structstr); 
        smilst = self.getSmiByPOSID(poslst);  
        smilst = [i.replace("*", "") for i in smilst]; 
        retmols = []
        for i,j in zip(poslst,smilst):
            try:
                tmpmol = Chem.MolFromSmiles(j)
                if tmpmol.HasSubstructMatch(substruct):
                    retmols.append(i)
            except:
                pass
        print(f"Substructure Match: {len(retmols)} poses kept from {len(poslst)}")
        return retmols; 
        

