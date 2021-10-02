# -*- coding: utf-8 -*-
"""
@author: Wen-Hsuan Yu
"""

from pathlib import Path
import pandas as pd, numpy as np

script_folder = Path.cwd()
outputs_folder = script_folder.parent / 'Outputs'
input_data_folder = outputs_folder / 'KOvWT_DEG_inputs'

subtypes2newcelltype = {'ILC (ILC1.CD127+)':'ILC1', 'ILC (LIV.ILC1.DX5-)':'ILC1', 
                        'ILC (LPL.NCR+ILC1)':'ILC1', 'ILC (ILC2)':'ILC2', 
                        'ILC (LPL.NCR+ILC3)':'ILC3', 'ILC (ILC3.LTI.CD4+)':'ILC3', 
                        'ILC (ILC3.LTI.CD4-)':'ILC3', 'ILC (ILC3.LTI.4+)':'ILC3', 
                        'ILC (LIV.NK.DX5+)':'NK cells', 'ILC (LPL.NCR+CNK)':'NK cells', 
                        'T cells (T.4FP3+25+)':'Treg', 'T cells (T.Tregs)':'Treg', 
                        'T cells (T.CD4.1H)':'CD4+T cells', 
                        'T cells (T.CD4.24H)':'CD4+T cells', 
                        'T cells (T.CD4.48H)':'CD4+T cells', 
                        'T cells (T.CD4.5H)':'CD4+T cells', 
                        'T cells (T.CD4.96H)':'CD4+T cells', 
                        'T cells (T.CD4.CTR)':'CD4+T cells', 
                        'T cells (T.4EFF49D+11A+.D8.LCMV)':'CD4+T cells', 
                        'T cells (T.4MEM49D+11A+.D30.LCMV)':'CD4+T cells', 
                        'T cells (T.4NVE44-49D-11A-)':'CD4+T cells', 
                        'T cells (T.4SP24-)':'CD4+T cells', 
                        'T cells (T.4SP24int)':'CD4+T cells', 
                        'T cells (T.4SP69+)':'CD4+T cells', 
                        'T cells (T.CD4+TESTNA)':'CD4+T cells', 
                        'T cells (T.CD4+TESTDB)':'CD4+T cells', 
                        'T cells (T.CD4CONTROL)':'CD4+T cells', 
                        'T cells (T.CD4TESTJS)':'CD4+T cells', 
                        'T cells (T.CD4TESTCJ)':'CD4+T cells', 
                        'T cells (T.4MEM)':'CD4+T cells', 
                        'T cells (T.4Mem)':'CD4+T cells', 
                        'T cells (T.4MEM44H62L)':'CD4+T cells', 
                        'T cells (T.4Nve)':'CD4+T cells', 
                        'T cells (T.4NVE)':'CD4+T cells', 
                        'T cells (T.4)':'CD4+T cells', 
                        'T cells (T.4.Pa)':'CD4+T cells', 
                        'T cells (T.4.PLN)':'CD4+T cells', 
                        'T cells (T.4FP3-)':'CD4+T cells', 
                        'T cells (T.CD8.1H)':'CD8+T cells', 
                        'T cells (T.CD8.24H)':'CD8+T cells', 
                        'T cells (T.CD8.48H)':'CD8+T cells', 
                        'T cells (T.CD8.5H)':'CD8+T cells', 
                        'T cells (T.CD8.96H)':'CD8+T cells', 
                        'T cells (T.CD8.CTR)':'CD8+T cells', 
                        'T cells (T.8EFF.TBET+.OT1LISOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.TBET-.OT1LISOVA)':'CD8+T cells', 
                        'T cells (T.8EFFKLRG1+CD127-.D8.LISOVA)':'CD8+T cells', 
                        'T cells (T.8MEMKLRG1-CD127+.D8.LISOVA)':'CD8+T cells', 
                        'T cells (T.8SP24-)':'CD8+T cells', 
                        'T cells (T.8SP24int)':'CD8+T cells', 
                        'T cells (T.8SP69+)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D15.VSVOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D5.VSVOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.VSVOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D8.VSVOVA)':'CD8+T cells', 
                        'T cells (T.8MEM)':'CD8+T cells', 
                        'T cells (T.8Mem)':'CD8+T cells', 
                        'T cells (T.8MEM.OT1.D106.VSVOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D45VSV)':'CD8+T cells', 
                        'T cells (T.8Nve)':'CD8+T cells', 
                        'T cells (T.8NVE)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D10LIS)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D10.LISOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D15LIS)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D15.LISOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1LISO)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.LISOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D8LISO)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.D8.LISOVA)':'CD8+T cells', 
                        'T cells (T.8MEM.OT1.D100.LISOVA)':'CD8+T cells', 
                        'T cells (T.8MEM.OT1.D45.LISOVA)':'CD8+T cells', 
                        'T cells (T.8NVE.OT1)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.12HR.LISOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.24HR.LISOVA)':'CD8+T cells', 
                        'T cells (T.8EFF.OT1.48HR.LISOVA)':'CD8+T cells'}
celltype_rename = {'CD8+T cells':'CD8+ T cell', 'Macrophages':'Macrophage', 
                   'B cells':'B cell', 'NK cells':'NK cell', 
                   'CD4+T cells':'CD4+ T cell', 'Monocytes':'Monocyte', 
                   'Other T cells':'Other T cell', 'Neutrophils':'Neutrophil', 
                   'Mast cells':'Mast cell', 'Basophils':'Basophil'}
subtype_rename = {'Tgd.vg5+.act':'Tgd.VG5+.ACT', 'B1a':'B1A', 
                  'B.Fo':'B.FO', 'T.8Nve':'T.8NVE'}

fname = outputs_folder / 'SingleR_ImmGen_fine_single_labels.txt'
ImmGen = pd.read_csv(fname, sep='\t', header=0, index_col=0, encoding='utf-8')

temp = {}
for barcode in ImmGen.index:
    label = ImmGen.loc[barcode, 'labels']
    if label in subtypes2newcelltype:
        maintype = subtypes2newcelltype[label]
    elif 'T cells' == label[:7]:
        maintype = 'Other T cells'
    else:
        maintype = label.split(' (')[0]
    if maintype in celltype_rename:
        maintype = celltype_rename[maintype]
    subtype = label.split(' (')[-1][:-1]
    if maintype in celltype_rename:
        maintype = celltype_rename[maintype]
    if subtype in subtype_rename:
        subtype = subtype_rename[subtype]
    temp.update({barcode:{'Maintype':maintype, 'Subtype':subtype}})
ImmGen = pd.DataFrame.from_dict(temp, orient='index')
fname = outputs_folder / 'WT_and_KO_cells_celltypes.txt'
ImmGen.to_csv(fname, sep='\t', header=True, index=True, index_label='Cell barcode', encoding='utf-8')

fname = outputs_folder / 'Seurat_integration_RNA_counts.txt'
data = pd.read_csv(fname, sep='\t', header=0, index_col=0, dtype=np.str_, encoding='utf-8')

celltype_list = sorted(set(ImmGen['Maintype'].values.tolist()))
for celltype in celltype_list:
    barcodes = sorted(ImmGen.index[ImmGen['Maintype']==celltype].tolist())
    cond = []
    for barcode in barcodes:
        if 'WT' in barcode:
            cond.append('0')
        else:
            cond.append('1')
    n_comb_0, n_comb_1 = cond.count('0'), cond.count('1')
    if (n_comb_0 > 2) and (n_comb_1 > 2):
        fname = input_data_folder / ('Differential_expression_analysis_%s_group.txt'%(celltype))
        with fname.open(mode='w', encoding='utf-8') as fw:
            fw.write('%s\n'%('\t'.join(cond)))
        fname = input_data_folder / ('Differential_expression_analysis_%s_matrix.txt'%(celltype))
        temp = data.loc[:, barcodes]
        temp.to_csv(fname, sep='\t', header=True, index=True, index_label='', encoding='utf-8')
        