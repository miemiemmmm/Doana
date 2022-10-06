import os
import json
import pickle
import numpy as np
import rdkit.Chem as Chem
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import pytraj as pt
import matplotlib.ticker as ticker
from IPython.display import Image

VDWRADII = {'AA': 1.85, 'AG': 1.72, 'AL': 1.84, 'AR': 1.88, 'AT': 2.02, 'AU': 1.66, 'B': 1.92, 'BA': 2.68,
  'BE': 1.53, 'BI': 2.07, 'BR': 1.85, 'C': 1.7, 'CA': 2.31, 'CD': 1.58, 'CL': 1.75, 'Cl': 1.75,
  'CS': 3.43, 'CU': 1.4, 'F': 1.47, 'FR': 3.48, 'GA': 1.87, 'GE': 2.11, 'H': 1.1,
  'HE': 1.4, 'HH': 1.55, 'I': 1.98, 'IN': 1.93, 'K': 2.75, 'KR': 2.02, 'LI': 1.82, 'MG': 1.73,
  'N': 1.55, 'NA': 2.27, 'NE': 1.54, 'NI': 1.63, 'O': 1.52, 'P': 1.8, 'PB': 2.02, 'PD': 1.63,
  'PO': 1.97, 'PT': 1.75, 'RA': 2.83, 'RN': 2.2, 'RR': 3.03, 'S': 1.8, 'SB': 2.06, 'SE': 1.9,
  'SI': 2.1, 'SN': 2.17, 'SR': 2.49, 'TE': 2.06, 'TL': 1.96, 'U': 1.86, 'XE': 2.16, 'ZN': 1.39,
}


class PLIFGen_Dock:
  """
    Both SDF and MOL2 formats are supported as input
    Example:
    >>> from Doana import analysis
    >>> parmdic = {
      'reflig'     : "/home/yzhang/Documents/Teachings/Tim_Ruth/PaIVKLkS/tmp_Sampling_target_clean.mol2",
      'profile'    : "/home/yzhang/Documents/Teachings/Tim_Ruth/PaIVKLkS/tmp_Sampling_target.pdb",
      'resultmols' : "/home/yzhang/Documents/Teachings/Tim_Ruth/PaIVKLkS/PaIVKLkS_SEEDdock.mol2",
      'resultdat'  : "/home/yzhang/Documents/Teachings/Tim_Ruth/PaIVKLkS/PaIVKLkS_SEEDdock.dat",
      'outpkl'     : "/home/yzhang/Documents/Teachings/Tim_Ruth/PaIVKLkS/PaIVKLkS_SEEDdock.pkl",
      'onlymols'   : list(range(1,30))+list(range(50,69)),
    }
    >>> wrapper = analysis.PLIFGen_Dock(parmdic)
    >>> wrapper.gen()
    >>> OIdic = wrapper.calc_OI(1.5, printrecords=False)
    >>> wrapper.savedata()
  """
  def __init__(self, parmdic):
    import prolif as plf
    import MDAnalysis as mda
    self.parms = parmdic;
    self.refmol = [];
    self.lig_suppl = [];
    self.ligandfile = parmdic["resultmols"];
    self.type = (self.parms["type"] if "type" in self.parms.keys() and len(self.parms["type"]) > 0 else "seed")
    # Default parameters
    self.FP_TYPE1 = ["HBDonor", "HBAcceptor", "PiStacking", "CationPi","Cationic"]
    self.FP_TYPE2 = ['Hydrophobic']

    # Initial receptor PDB structure
    try:
      self.profile = self.parms["profile"]
      print(f"Loading receptor PDB structure: {self.profile}")
      mda_prot = mda.Universe(self.profile, top=self.profile, guess_bonds=True, vdwradii=VDWRADII);
      elements = mda.topology.guessers.guess_types(mda_prot.atoms.names);
      mda_prot.add_TopologyAttr("elements", elements);
      self.prot = plf.Molecule.from_mda(mda_prot);
    except:
      print("Failed to load the PDB into prolif");
      print("Please carefully read the error message and modify the PDB file.");
      raise

    # Initial pose library
    if "reflig" in self.parms.keys() and len(self.parms["reflig"]) > 0:
      # Initialize reference ligand while initializing the PLIF generator
      self.setreflig(self.parms["reflig"])

    # Read ligand molecules
    self.failed_mol = []
    if '.sdf' in self.ligandfile:
      # SDF format
      if ("onlymols" in self.parms.keys()) and (self.parms["onlymols"] != ""):
        self.parms["onlymols"] = planned_indexes = [i for i in self.parms["onlymols"]];
        print(f"Will only input the following molecules: {planned_indexes}");
        self.lig_suppl, self.success_mol, self.failed_mol = self.sdf_supplier(self.ligandfile, inputmols=planned_indexes)
      else:
        print("Will input all molecules")
        self.lig_suppl, self.success_mol, self.failed_mol = self.sdf_supplier(self.ligandfile)
    elif '.mol2' in self.ligandfile:
      # MOL2 format
      if ("onlymols" in self.parms.keys()) and (self.parms["onlymols"] != ""):
        self.parms["onlymols"] = planned_indexes = [i for i in self.parms["onlymols"]];
        print(f"Will only input the following molecules: {planned_indexes}");
        self.lig_suppl, self.success_mol, self.failed_mol = self.mol2_supplier(self.ligandfile, inputmols=planned_indexes);
      else:
        print("Will input all molecules")
        self.lig_suppl, self.success_mol, self.failed_mol = self.mol2_supplier(self.ligandfile)
    if len(self.failed_mol) > 0:
        print(f"{len(self.failed_mol)} molecules failed in RDKit SDF loading:", *self.failed_mol)
    print(f"{len(self.lig_suppl)} molecules put to the library supplier (from {self.ligandfile})")
    self.datfile = parmdic["resultdat"]

  def sdf_supplier(self, path, inputmols=False):
    """
      Apart from output molecule, successed and failed molecules should be recorded
    """
    suppl = Chem.SDMolSupplier(path, removeHs=False);
    indexes = [ i for i in range(suppl.__len__()) ];
    planned_indexes = [];
    if isinstance(inputmols, list):
      planned_indexes = inputmols;
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
      if i > max(planned_indexes):
        break
    return mols, success_mol, failed_mol;

  def mol2_supplier(self, path, inputmols=False):
    """
      Apart from output molecule, successed and failed molecules should be recorded
    """
    with open(path, "r") as file1:
      mol2blocks = file1.read().strip("@<TRIPOS>MOLECULE").split("@<TRIPOS>MOLECULE");
    indexes = [ i for i in range(len(mol2blocks))];
    planned_indexes = [];
    if isinstance(inputmols, list):
      planned_indexes = inputmols;
    else:
      planned_indexes = indexes;
    failed_mol = [];  mols = [];  success_mol = [];
    for i, moli in zip(indexes, mol2blocks):
      mol2str = "@<TRIPOS>MOLECULE"+moli;
      if i in planned_indexes:
        try:
          mol = Chem.MolFromMol2Block(mol2str, removeHs=False)
          mols.append(plf.Molecule.from_rdkit(mol));
          success_mol.append(i)
        except:
          failed_mol.append(i)
    return mols, success_mol, failed_mol;

  def setreflig(self, ligfile):
    """
      Not fixed to self.ligandfile because user might want a second reference ligand to use.
      Only needs the XYZ for now
    """
    try:
      if "sdf" in ligfile:
        refmol, success, failed = self.sdf_supplier(ligfile);
      elif "mol2" in ligfile:
        refmol, _s, _f = [i for i in self.mol2_supplier(ligfile)];
      fp_hpl = plf.Fingerprint(self.FP_TYPE1+self.FP_TYPE2)
      fp_hpl.run_from_iterable(refmol, self.prot)
      fp_data = fp_hpl.to_dataframe()

      refmol = Chem.RemoveHs(refmol[0]);
      self.refxyz = [_.GetPositions() for _ in refmol.GetConformers()][0];
      print(f"Reference ligand interactions: {fp_data}")

    except:
      print("Cannot Sanitize the molecule by RDKit. Trying to use pytraj")
      import pytraj as pt;
      tmptraj = pt.load(ligfile, top=ligfile, mask="!@H*");
      self.refxyz = tmptraj.xyz[0].astype(float);
    if len(self.refxyz) > 0:
      print("Successfully loaded the reference ligand ")
  def addinfo(self, df, topn=-1):
    # Add the Columns about the score and pos_id
    with open(self.datfile, "r") as file1:
      thelist = [i.strip('\n').split() for i in file1.read().split('\n')]
    score_table = [[i + 1] + thelist[i] for i in self.success_mol]
    score_table = np.array(score_table);
    if self.type == "seed":
      col_pos_id     = score_table[:,1].astype(str);
      col_smiles     = score_table[:,2].astype(str);
      col_seed_total = score_table[:,3].astype(str);
      col_seed_vdw   = score_table[:,4].astype(str);
      col_seed_delec = score_table[:,5].astype(str);
      col_seed_elinw = score_table[:,6].astype(str);
      col_seed_resdes  = score_table[:,7].astype(str);
      col_seed_fragdes = score_table[:,8].astype(str);
      col_ranking    = score_table[:,0].astype(int);
      col_nha        = score_table[:,9].astype(str);
      col_date       = score_table[:,10].astype(str);
      col_source_campaign   = score_table[:,11].astype(str);

      df["Rank"] = col_ranking
      df["seed_total"]   = col_seed_total;
      df["seed_vdw"] = col_seed_vdw;
      df["seed_delec"] = col_seed_delec;
      df["seed_elinw"] = col_seed_elinw;
      df["seed_resdes"] = col_seed_resdes;
      df["seed_fragdes"] = col_seed_fragdes;
      df["pos_id"]   = col_pos_id;
      df["nha"] = col_nha;
      df["date"] = col_date;
      df["source_campaign"] = col_source_campaign;
      df["can_smile"] = col_smiles;
      return df
    elif self.type == "vina":
      col_pos_id = score_table[:, 1].astype(str);
      col_smiles = score_table[:, 3].astype(str);
      col_vina_total = score_table[:, 4].astype(str);
      col_inter = score_table[:, 5].astype(str);
      col_intra = score_table[:, 6].astype(str);
      col_torsion = score_table[:, 7].astype(str);
      col_ranking = score_table[:, 1].astype(int);
      col_date = score_table[:, 10].astype(str);
      col_source_campaign = score_table[:, 9].astype(str);

      df["Rank"] = col_ranking
      df["vina_total"] = col_vina_total;
      df["vina_inter"] = col_inter;
      df["vina_intra"] = col_intra;
      df["vina_torsion"] = col_torsion;
      df["pos_id"] = col_pos_id;
      df["date"] = col_date;
      df["source_campaign"] = col_source_campaign;
      df["can_smile"] = col_smiles;
      return df

  def gen(self):
    # Generate protein-ligand interaction fingerprint
    fp_hpl = plf.Fingerprint(self.FP_TYPE1)
    fp_hpl.run_from_iterable(self.lig_suppl, self.prot)
    fp_data = fp_hpl.to_dataframe()
    self.fp_data = self.addinfo(fp_data)

    fp_hpb = plf.Fingerprint(self.FP_TYPE2)
    fp_hpb.run_from_iterable(self.lig_suppl, self.prot)
    fp_data2 = fp_hpb.to_dataframe()
    self.fp_data2 = self.addinfo(fp_data2)

  def savedata(self):
    """
      Serialize the PLIF data
    """
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

  def calc_OI(self, dist_cutoff, printrecords=True):
    """
      Additional computations: Overlapping Index
    """
    if self.fp_data.__len__() != len(self.success_mol):
      raise Exception(f"Fingerprint dataframe 1 length ({self.fp_data.__len__()}) not equil to ligand supplier{len(self.success_mol)}")
    if self.fp_data2.__len__() != len(self.success_mol):
      raise Exception(f"Fingerprint dataframe 1 length ({self.fp_data.__len__()}) not equil to ligand supplier{len(self.success_mol)}")
    OI_results = [];  OI_ref  = [];  OI_avg  = [];  rnha = [];
    for i, resultmol in zip(range(len(self.lig_suppl)), self.lig_suppl):
      resultmol = Chem.RemoveHs(resultmol);
      testxyz = [_.GetPositions() for _ in resultmol.GetConformers()][0]

      nha_ratio = len(testxyz) / len(self.refxyz)

      distances = distance_matrix(self.refxyz, testxyz)
      distlt = distances < dist_cutoff
      count1 = np.count_nonzero(np.any(distlt, axis = 0))
      count2 = np.count_nonzero(np.any(distlt, axis = 1))
      rat1  = count1 / len(testxyz)
      rat2  = count2 / len(self.refxyz)
      r_avg = np.mean([rat1,rat2])

      rnha.append(nha_ratio)
      OI_results.append(rat1);
      OI_ref.append(rat2);
      OI_avg.append((rat1+rat2)/2)
      if printrecords: print(f"Rank: {self.fp_data.iloc[i].Rank.values[0]} | Ratio_t: {round(rat1,3)} | Ratio_r: {round(rat2,3)} | Ravg {round(r_avg,3)} \
| Rnha: {round(nha_ratio,3)} ({len(testxyz)}/{len(self.refxyz)})")
    self.fp_data["OI_avg"] = self.fp_data2["OI_avg"] = self.OI_avg = OI_avg
    self.fp_data["OI_ref"] = self.fp_data2["OI_ref"] = self.OI_ref = OI_ref
    self.fp_data["OI_results"] = self.fp_data2["OI_results"] = self.OI_results = OI_results


class PLIFGen_MD:
  """
  >>> thedic = {
        "masklig" : ":LIG",
        # The mask cannot contain separate molecules
        # NOTE ESPECIALLY: ligand, water, ions, cofactors
        # If maskpro is not defined, the protein will be chosen
        "maskpro" : ":LIG<:7&!:T3P,LIG,K+,CL-",
        "stride" : 1,
        "topfile" : "/home/miemie/Downloads/C0084GULFky4/C0084GULFky4_PDB.pdb", 
        "trajfile" : "/home/miemie/Downloads/C0084GULFky4/C0084GULFky4_TRJ.nc",
        "outpkl":"/tmp/test.pkl", 
      }
  >>> mdplif = PLIFGen_MD(thedic)
  >>> mdplif.gen()
  >>> mdplif.calc_OI(dist_cutoff=1.5, use_mean=True)
  >>> plt.plot(mdplif.fp_data["OI_avg"])
  >>> mdplif.savedata()
  """
  def __init__(self, parmdic):
    import MDAnalysis as mda; 
    from Doana import utils; 
    self.parms = parmdic;
    self.FP_TYPE1 = ['HBDonor', 'HBAcceptor', 'PiStacking', 'CationPi', 'Cationic']
    self.FP_TYPE2 = ['Hydrophobic']
    
    self.topfile = parmdic["topfile"]; 
    self.trajfile = parmdic["trajfile"]; 
    self.stride = parmdic["stride"]; 
    
    self.mda_prot = mda.Universe(self.topfile, self.trajfile, guess_bonds=True, vdwradii=VDWRADII)
    elements = mda.topology.guessers.guess_types(self.mda_prot.atoms.names)
    self.mda_prot.add_TopologyAttr('elements', elements)
    
    #     mda_ref = mda.Universe(self.topfile, self.trajfile, guess_bonds=True, vdwradii=VDWRADII)
    #     self.mda_prot.trajectory[-1]
    #     mda_ref.trajectory[0]
    #     print(f"will run alignment")
    #     alignment = align.AlignTraj(self.mda_prot, mda_ref, select="name CA", weights="mass")
    #     alignment.run()
    #     print(f"RMSD new is {rmsd_new}")
    
    traj_pt = pt.load(self.topfile); 
    traj_pt.top.set_reference(traj_pt[0]); 
    traj_pt.superpose("@CA"); 
    
    # Select the ligand and the protein moiety to calculate. 
    try: 
      self.lig, self.ligindex = utils.PTToMDASelect(traj_pt, self.mda_prot, self.parms["masklig"]); 
    except: 
      print("Failed to select, selecting the devault value <resname LIG or resname MDL>"); 
      self.lig = self.mda_prot.atoms.select_atoms('resname LIG or resname MDL');
      self.ligindex = np.array([i.index for i in self.lig]); 
    self.ligmask = "@"
    for x in self.ligindex: 
      self.ligmask += f"{x},"
    
    try: 
      self.prot, self.protindex = utils.PTToMDASelect(traj_pt, self.mda_prot, self.parms["maskpro"]); 
    except: 
      print("Failed to select, selecting the default values: <protein>"); 
      self.prot = self.mda_prot.atoms.select_atoms('protein'); 
      self.protindex = np.array([i.index for i in self.prot]); 
    self.promask = "@"
    for x in self.protindex: 
      self.promask += f"{x},"
    
  def calc_OI(self, dist_cutoff=1, printrecords=True, use_mean=False):
    """
      Additional computations: Overlapping Index
      It is strange that the alignment of MDAnalysis is not working
    """
    traj_pt = pt.load(self.trajfile, top=self.topfile, stride=self.stride); 
    traj_pt.top.set_reference(traj_pt[0]); 
    traj_pt.superpose("@CA"); 
    OI_avg, OI_details = OverlappingIndex(traj_pt, self.ligmask, dist_cutoff=dist_cutoff, use_mean=use_mean)
    self.fp_data["OI_avg"] = self.fp_data2["OI_avg"] = OI_avg
    self.fp_data["OI_ref"] = self.fp_data2["OI_ref"] = OI_details['OI_reference']
    self.fp_data["OI_results"] = self.fp_data2["OI_results"] = OI_details['OI_ligand']
    
  
  
  def calc_OI2(self, dist_cutoff=1, printrecords=True, use_mean=False):
    """
      Additional computations: Overlapping Index
    """
    OI_results = [];  
    OI_ref  = [];  
    OI_avg  = [];  
    rnha = [];
    
    if use_mean==True:
      self.refxyz = [i.positions[self.ligindex] for i in self.mda_prot.trajectory]
      self.refxyz = np.mean(self.refxyz, axis=0)
    else: 
      self.refxyz = self.mda_prot.trajectory[0].positions[self.ligindex]
    
    print("==>", self.refxyz)
    print("==>", self.refxyz.shape)
    
    
    for idx, theframe in enumerate(self.mda_prot.trajectory[::self.stride]): 
      testxyz = theframe.positions[self.ligindex]; 
      nha_ratio = len(testxyz) / len(self.refxyz)

      distances = distance_matrix(self.refxyz, testxyz)
      distlt = distances < dist_cutoff
      count1 = np.count_nonzero(np.any(distlt, axis = 0))
      count2 = np.count_nonzero(np.any(distlt, axis = 1))
      rat1  = count1 / len(testxyz)
      rat2  = count2 / len(self.refxyz)
      r_avg = np.mean([rat1,rat2])

      rnha.append(nha_ratio)
      OI_results.append(rat1);
      OI_ref.append(rat2);
      OI_avg.append((rat1+rat2)/2)
    self.fp_data["OI_avg"] = self.fp_data2["OI_avg"] = OI_avg
    self.fp_data["OI_ref"] = self.fp_data2["OI_ref"] = OI_ref
    self.fp_data["OI_results"] = self.fp_data2["OI_results"] = OI_results

  def gen(self):
    import prolif as plf
    # use default interactions
    fp = plf.Fingerprint(self.FP_TYPE1)
    fp.run(self.mda_prot.trajectory[::self.stride], self.lig, self.prot, residues='all')
    self.fp_data = fp.to_dataframe()

    fp2 = plf.Fingerprint(self.FP_TYPE2);
    fp2.run(self.mda_prot.trajectory[::self.stride], self.lig, self.prot, residues='all')
    self.fp_data2 = fp2.to_dataframe()
  
  def savedata(self):
    """
      Serialize the PLIF data
    """
    dictosave={
      "parms" : self.parms,
    }
    try:
      dictosave['fp_data'] = self.fp_data;
      dictosave["fptypes"] = self.FP_TYPE1;
    except:
      pass
    try:
      dictosave['fp_data2'] = self.fp_data2;
      dictosave["fptypes2"] = self.FP_TYPE2;
    except:
      pass

    with open(self.parms["outpkl"], "wb") as fileout:
      pickle.dump(dictosave, fileout)
      print("Saved the fingerprint to file:", self.parms["outpkl"])

class PLIFRead_Dock:
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
    # Available operator, all and any
    selcols = [fpdata.columns[i] for i in fp_sel]
    fp_status = (fpdata.iloc[:,fp_sel] == True).to_numpy();
    if operator == "or":
      print("Fingerprint Selector: Using and operator, the poses are kept if they have any of defined fingerprint.");
      fp_status = [np.any(i) for i in fp_status];
    if operator == "and":
      print("Fingerprint Selector: Using and operator, the poses are kept if they have all of defined fingerprint.");
      fp_status = [np.all(i) for i in fp_status];
    status = np.where( np.array(fp_status) == True )[0];
    status = np.array(list(set(status)));
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
    from rdkit.Chem import Draw, AllChem
    result = np.array([ self.fpdata[self.fpdata['pos_id'] == i]["can_smile"].to_numpy()[0] for i in pos_ids]).astype(str)
    result = [i.replace("*","") for i in result]
    mols = [Chem.MolFromSmiles(m) for m in result]
    for m in mols: tmp = AllChem.Compute2DCoords(m)
    img = Draw.MolsToGridImage(mols, molsPerRow=mpr, subImgSize=(400,400), maxMols=100,
                               legends=[f"pos_id: {i}" for i in pos_ids],returnPNG=False)
    img.save("/tmp/molgrid.png")
    return Image("/tmp/molgrid.png")

  def drawHist(self, dataframe, n_bins=10):
    dataset = dataframe.to_numpy().astype(float);
    fig, ax = plt.subplots(1, 1,figsize=(10, 5), sharey=True, tight_layout=True)
    N, bins, patches = ax.hist(dataset, bins=n_bins)
    return fig, ax, N, bins, patches

  def drawHistByCol(self, dataframe, col, n_bins=10):
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

class PLIFRead_MD:
  """
  >>> parms = {
        "picklefile"  : "/tmp/test.pkl",
        "outputdir"   : "/tmp/",
        "file_prefix" : "MDPLIFREADER",
        "TOPN"        : 7,
      }
  >>> MDreader = da.analysis.PLIFRead_MD(parms)
  >>> MDreader.savefig()
  """
  def __init__(self, parms):
    self.parms = parms
    with open(self.parms['picklefile'], "rb") as file1:
      self.fp_dic = pickle.load(file1);

    self.fp_data = self.fp_dic["fp_data"];
    self.fp_data2 = self.fp_dic["fp_data2"];
    self.rownr = len(self.fp_data)

    try:
      self.outputdir = self.parms["outputdir"]
    except:
      print("The outputdir is not defined, using the current working directory. ")
      self.outputdir = "./";
    self.file_prefix = self.parms["file_prefix"]
    self.fig, self.ax = self.MD_prolif(topn=self.parms['TOPN'], summary=True);

  def MD_prolif(self, topn=8, sel_fps=[], labelsize = 13, rotation=-20, summary=True, sortbyfp=False):
    # Get the DataFrame and its meta-information
    fig, ax = plt.subplots(2, 1, figsize=(12,12))
    fig.tight_layout(pad=7.0)
    resname="LIG"
    for i in self.fp_data:
      resnames = [j for j in i if "LIG" in j]
      if len(resnames) > 0:
        resname = resnames[0]

    FP1 = self.fp_data[resname];
    topn1 = self.TopFPNumber(FP1, topn);
    topdf, ind, percents = self.onlytopn(FP1, topn=topn1, summary=summary);
    ticks1 = ["\n".join(i)+f"\n{j*100:.1f}%" for i,j in zip(FP1.columns, percents)];
    self.plifplot(ax[0], topdf, rotation=rotation, ticks=ticks1)
    ax[0].set_title(f"Non hydrophobic interactions ({self.file_prefix})");

    FP2 = self.fp_data2[resname];
    topn2 = self.TopFPNumber(FP2, topn);
    topdf2, ind2, percents2 = self.onlytopn(FP2, topn=topn2, summary=summary);
    ticks2 = ["\n".join(i)+f"\n{j*100:.1f}%" for i,j in zip(FP2.columns, percents2)];
    self.plifplot(ax[1], topdf2, rotation=rotation, ticks=ticks2)
    ax[1].set_title(f"Hydrophobic interactions ({self.file_prefix})");
    return fig, ax

  def plifplot(self, theax, the_array, ticks=[], labelsize=12, **kwarg):
    # Note: Vertically flip the array because top row in dataframe starts from the botton in figure
    if "rotation" in kwarg.keys():
      rotation = kwarg["rotation"]
    else:
      rotation = -20
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
    percents = [ round(true_count[i]/self.rownr,3) for i in ind]
    return df, ind, percents
  def TopFPNumber(self, fingerprint, topn):
    if fingerprint.shape[1] < topn:
      TOPN1 = fingerprint.shape[1];
      print(f"FP number {TOPN1} less than desired number {topn} of fingerprint, setting top Number to {TOPN1}")
    else:
      TOPN1 = topn;
    return TOPN1

  def savefig(self):
    self.fig.savefig(self.outputdir+"/"+self.file_prefix+'_PLIFs.png', dpi=300, transparent=False, bbox_inches='tight')



def OverlappingIndex(traj, ligand_mask, dist_cutoff=1, use_mean=True):
  OI_lig = []; OI_ref  = []; OI_avg  = []; 
  sel = traj.top.select(ligand_mask)
  if use_mean==True:
    ligarr = np.array([traj.xyz[i][sel] for i in range(len(traj.xyz))]); 
    refxyz = np.mean(ligarr, axis=0); 
  else: 
    refxyz = traj[0].xyz[sel]
  for i in range(len(traj)):
    theframe = traj[i];
    this_ligxyz = theframe.xyz[sel]
    distances = distance_matrix(refxyz, this_ligxyz)
    distlt = distances < dist_cutoff
    dist_atoms_test = np.any(distlt, axis = 0)
    dist_atoms_ref = np.any(distlt, axis = 1)
    thisligOI = np.count_nonzero(dist_atoms_test)/len(this_ligxyz); 
    thisrefOI = np.count_nonzero(dist_atoms_ref)/len(refxyz); 
    OI_lig.append(thisligOI);
    OI_ref.append(thisrefOI);
    OI_avg.append((thisligOI+thisrefOI)/2); 
  OI_lig = np.array(OI_lig); 
  OI_ref = np.array(OI_ref); 
  OI_avg = np.array(OI_avg); 
  return OI_avg, {'OI_ligand': OI_lig, 'OI_reference': OI_ref, 'OI_average':OI_avg, "selection":sel}

def PairwiseDist(traj, mask1, mask2, use_mean=False):
  selmask1 = traj.top.select(mask1);
  selmask2 = traj.top.select(mask2);
  if use_mean == True:
    frame_mean = np.mean(traj.xyz, axis=0); 
    this_ligxyz = frame_mean[selmask1];
    this_proxyz = frame_mean[selmask2];
    ref_frame = distance_matrix(this_ligxyz, this_proxyz);
  else: 
    pdist = pt.pairwise_distance(traj, mask_1=mask1, mask_2=mask2); 
    ref_frame = pdist[0][0]; 

  minindex = [np.where(ref_frame[i] == np.min(ref_frame[i]))[0][0] for i in range(len(ref_frame))]
  absolute_index = [selmask2[i] for i in minindex]
  min_dists = np.min(ref_frame, axis=1)

  # For mask selection, remember to add 1 because it is the read index number
  # print(f"processing atom {resnames[i]}:{resids[i]+1}_{atoms[i]}@{indexs[i]+1} and atom {resnames[j]}:{resids[j]+1}_{atoms[j]}@{indexs[j]+1}: mean {np.mean(dist_tmp).round(2)}, std {np.std(dist_tmp).round(2)}")
  distlist = []
  for i,j in zip(selmask1, absolute_index):
    dist_tmp = pt.distance(traj, f"@{i+1} @{j+1}")
    distlist.append(dist_tmp)

  distarr = np.array(distlist)
  stds = np.std(distarr, axis=1)
  means = np.mean(distarr, axis=1)
  overall_std = np.mean(stds)
  ranking = stds.argsort().argsort(); 

  atom_names = np.array([i.name for i in traj.top.atoms])
  atom_ids = np.array([i.index for i in traj.top.atoms])
  resids = np.array([i.resid for i in traj.top.atoms])
  resnames = np.array([i.resname for i in traj.top.atoms])

  lig_atom_names = atom_names[selmask1]
  lig_atomids  = atom_ids[selmask1]
  lig_resids   = resids[selmask1]
  lig_resnames = resnames[selmask1]
  pro_atom_names = atom_names[absolute_index]
  pro_atomids  = atom_ids[absolute_index]
  pro_resids   = resids[absolute_index]
  pro_resnames = resnames[absolute_index]

  return distarr, ranking, {
    "stds":stds, "means":means, "overall_std":overall_std, 
    "lig_atom_names":lig_atom_names,"lig_atomids":lig_atomids,"lig_resids":lig_resids,"lig_resnames":lig_resnames,
    "pro_atom_names":pro_atom_names,"pro_atomids":pro_atomids,"pro_resids":pro_resids,"pro_resnames":pro_resnames,
  }



