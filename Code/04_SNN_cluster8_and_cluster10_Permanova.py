# -*- coding: utf-8 -*-
"""
@author: Wen-Hsuan Yu
"""

from pathlib import Path
import pandas as pd, numpy as np
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.stats.distance import permanova

script_folder = Path.cwd()
outputs_folder = script_folder.parent / 'Outputs'

fname = outputs_folder / 'Seurat_integration_PCA_cell_embeddings.txt'
pca = pd.read_csv(fname, sep='\t', header=0, index_col=0, encoding='utf-8')
pca = pca.iloc[:, :18]

fname = outputs_folder / 'Seurat_integration_SNN_clusters.txt'
clusters = pd.read_csv(fname, sep='\t', header=0, index_col=0, encoding='utf-8')

fname = outputs_folder / 'WT_and_KO_cells_celltypes.txt'
celltypes = pd.read_csv(fname, sep='\t', header=0, index_col=0, encoding='utf-8')

cluster2celltype = {}
for C in [8, 10]:
    cells = clusters.index[clusters['Cluster'] == C].tolist()
    temp = {}
    for barcode in cells:
        celltype = celltypes.loc[barcode, 'Maintype']
        if celltype not in temp:
            temp.update({celltype:[]})
        temp[celltype].append(barcode)
    cluster2celltype.update({C:temp})

fname = outputs_folder / 'Permanova_results.xlsx'
with pd.ExcelWriter(fname) as writer:
    for C in [8, 10]:
        frames = []
        for (celltype_A, celltype_B) in list(combinations(sorted(cluster2celltype[C].keys()), 2)):
            cells = cluster2celltype[C][celltype_A] + cluster2celltype[C][celltype_B]
            grouping = [celltype_A]*len(cluster2celltype[C][celltype_A]) + [celltype_B]*len(cluster2celltype[C][celltype_B])
            
            X = pca.loc[cells, :].copy()
            dm = squareform(pdist(X, metric='euclidean'))
            dist_mat = DistanceMatrix(dm, cells)
            np.random.seed(0)
            result = permanova(dist_mat, grouping, permutations=1000)
            result.name = ('%s vs %s'%(celltype_A, celltype_B))
            frames.append(result)
        result = pd.concat(frames, axis='columns')
        sheet_name = 'cluster %d'%(C)
        result.T.to_excel(writer, sheet_name=sheet_name, header=True, index=True, 
                          index_label='', encoding='utf-8')
