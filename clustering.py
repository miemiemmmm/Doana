import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from rdkit import Chem
from rdkit.Chem import AllChem

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


def GetOneCluster(raw_data:pd.DataFrame, key:str, idx:int):
    if key not in raw_data.keys():
        raise("Fatal: not found the group by keys")
    cluster_rows = []
    for k,v in raw_data.groupby(key):
        if k == idx:
            for r, row in v.iterrows():
                cluster_rows.append([row.SMILES,  k, len(v), row.ID])
    ret_df = pd.DataFrame(cluster_rows, columns=["SMILES", "Cluster", "Num", "ID"])
    return ret_df


def GetEachCluster(raw_data:pd.DataFrame, key:str, idx_method):
    """
    There have to be SMILES and ID column in the input dataframe
    """
    if key not in raw_data.keys():
        raise("Fatal: not found the group by keys")
    cluster_rows = []
    for k,v in raw_data.groupby(key):
        _v = v.reset_index(drop=True)
        idx = idx_method(_v)
        print("Min distance idx: ", idx)
        cluster_rows.append([
            v.SMILES.values[idx],
            k,
            len(v),
            v.ID.values[idx],
        ])
    ret_df = pd.DataFrame(cluster_rows, columns=["SMILES", "Cluster", "Num", "ID"])
    return ret_df


# Convert SMILES to RDKit molecules and then to ECFP fingerprints
raw = pd.read_csv("/home/yzhang/Downloads/BFL-1_compounds.csv")
molecules = [Chem.MolFromSmiles(smiles) for smiles in raw.Smiles]
fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for mol in molecules]

# Dimensionality reduction and clustering
pca = PCA(n_components=10)
reduced_data_pca = pca.fit_transform(fingerprints)

kmeans = KMeans(n_clusters=20)
clusters = kmeans.fit_predict(reduced_data_pca)
centroids = kmeans.cluster_centers_

uniq = np.unique(clusters, return_counts=True)

final_data = []
for i in range(len(clusters)):
    final_data.append([
        reduced_data_pca[i][0], reduced_data_pca[i][1],
        clusters[i],
        raw.Smiles.values[i],
        raw.ID.values[i]
    ])
data_final = pd.DataFrame(final_data, columns = ["pc1", "pc2", "cluster", "SMILES", "ID"])

# Distance to centroid
distances = []
for i in range(len(clusters)):
    distances.append(np.linalg.norm(reduced_data_pca[i] - centroids[clusters[i]]))
data_final["distance"] = distances

plt.figure(figsize=(6,6))
sns.scatterplot(data=data_final, x="pc1", y="pc2", hue="cluster", palette='rainbow')
legends = []
for i, centroid in enumerate(centroids):
    plt.scatter(centroid[0], centroid[1], c='black', marker='x')
    plt.text(centroid[0]+0.05, centroid[1]-0.04, f"c{i+1}", c='black')
    legends.append(f"c{i+1}")
plt.legend([], frameon=False)
plt.title('2D PCA of Molecular Fingerprints')
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')
plt.show()


def GetMinDistIdx(df:pd.DataFrame):
    return df.distance.idxmin()
print(GetMinDistIdx(data_final))


data2grid = GetEachCluster(data_final, "cluster", GetMinDistIdx)
# data2grid = GetOneCluster(data_final, "cluster", 14)
print(data2grid)
# mols2grid.display(data2grid,subset=["img","Num"])
