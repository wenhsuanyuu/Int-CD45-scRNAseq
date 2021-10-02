# -*- coding: utf-8 -*-
"""
@author: Wen-Hsuan Yu
"""

from pathlib import Path
import pandas as pd

script_folder = Path.cwd()
# Xu et al. 2019 (GEO accession:GSE124880)
data_folder = script_folder.parent / 'Gene_Barcode_Matrices/GSE124880_Xu'
outputs_folder = script_folder.parent / 'Outputs'

abbr_dict = {'LP':'Lamina propria', 'PP':'Peyer\'s patche'}

fname = data_folder / 'GSE124880_PP_LP_mm10_count_barcode.tsv.gz'
GSE124880_barcode = pd.read_csv(fname, sep='\t', header=None, index_col=None, encoding='utf-8')

fname = data_folder / 'Food-Allergy-PP-LP_cluster_v1.txt'
Food_Allergy_cluster = pd.read_csv(fname, sep='\t', header=[0, 1], index_col=0, encoding='utf-8')
Food_Allergy_cluster.columns = ['Cluster']

sample_type_dict = {'IL2_IL7_rep1':'Healthy', 
                    'IL2_IL7_rep2':'Healthy', 
                    'IL2_IL7_rep3':'Healthy', 
                    'IL2_IL7_IL25_rep1':'Healthy', 
                    'IL2_IL7_IL25_rep2':'Healthy', 
                    'IL2_IL7_IL25_rep3':'Healthy', 
                    'IL2_IL7_IL25_CGRP0.4ug_rep1':'Healthy', 
                    'IL2_IL7_IL25_CGRP0.4ug_rep2':'Healthy', 
                    'IL2_IL7_IL25_CGRP0.4ug_rep3':'Healthy', 
                    'IL2_IL7_CGRP0.4ug_rep1':'Healthy', 
                    'IL2_IL7_CGRP0.4ug_rep2':'Healthy', 
                    'IL2_IL7_CGRP0.4ug_rep3':'Healthy', 
                    'IL2_IL7_CGRP0.4ug_rep4':'Healthy', 
                    '04242017_PP_location_WT_Allergy_v2_AD':'Inflammation', 
                    '04242017_PP_location_WT_Allergy_v2_AI':'Inflammation', 
                    '04242017_PP_location_WT_Allergy_v2_AJ':'Inflammation', 
                    '04242017_PP_location_WT_Allergy_v2_CD':'Healthy', 
                    '04242017_PP_location_WT_Allergy_v2_CI':'Healthy', 
                    '04242017_PP_location_WT_Allergy_v2_CJ':'Healthy', 
                    '05152017_PP_location_WTAllergy_v2_AD':'Inflammation', 
                    '05152017_PP_location_WTAllergy_v2_AI':'Inflammation', 
                    '05152017_PP_location_WTAllergy_v2_AJ':'Inflammation', 
                    '05152017_PP_location_WTAllergy_v2_CD':'Healthy', 
                    '05152017_PP_location_WTAllergy_v2_CI':'Healthy', 
                    '05152017_PP_location_WTAllergy_v2_CJ':'Healthy', 
                    '06132017_PP_WT_Allergy_v2_allergy':'Inflammation', 
                    '06132017_PP_WT_Allergy_v2_ctrl':'Healthy', 
                    '07142017_PP_WT_Allergy_v2_Allergy':'Inflammation', 
                    '07142017_PP_WT_Allergy_v2_Control':'Healthy', 
                    '09012017_PP_WT_Allergy_v2_Allergy_1':'Inflammation', 
                    '09012017_PP_WT_Allergy_v2_Allergy_2':'Inflammation', 
                    '09012017_PP_WT_Allergy_v2_Control_1':'Healthy', 
                    '09012017_PP_WT_Allergy_v2_Control_2':'Healthy', 
                    '10172017_ctrl_allergy_nonPP_nonTB_v2_Allergy_1':'Inflammation', 
                    '10172017_ctrl_allergy_nonPP_nonTB_v2_Allergy_2':'Inflammation', 
                    '10172017_ctrl_allergy_nonPP_nonTB_v2_Ctrl_1':'Healthy', 
                    '10172017_ctrl_allergy_nonPP_nonTB_v2_Ctrl_2':'Healthy', 
                    '11292017_IgDnegSI_nonTBSI_v2_Allergy_IgDLow':'Inflammation', 
                    '11292017_IgDnegSI_nonTBSI_v2_Allergy_nonTB':'Inflammation', 
                    '11292017_IgDnegSI_nonTBSI_v2_Ctrl_IgDLow':'Healthy', 
                    '11292017_IgDnegSI_nonTBSI_v2_Ctrl_nonTB':'Healthy'}
cellsubset_label = ['Resting CD4+ T cell', 'Resting B cell', 'ILC3', 
                    'LTi cell', 'CD8+ T cell', 'pDC', 'Activated CD4+ T cell', 
                    'NK cell', 'GC B cell (DZ)', 'GC B cell (LZ)', 'ILC1', 
                    'ILC2', 'CD4+ T cell (low UMI count)', 'DC (CD103+CD11b+)', 
                    'ILC3 (low UMI count)', 'LTi cell (low UMI count)', 
                    'Gamma delta T cell (Xcl1+)', 'Unresolved', 'NKT cell', 
                    'Plasma cell', 'DC (CD103+CD11b-)', 'Macrophage', 
                    'Resting B cell (low UMI count)', 'DC (CD103-C1)', 
                    'Endothelial cell', 'Mast cell', 'NK cell (low UMI count)', 
                    'GC B cell (low UMI count)', 'Doublets', 'Epithelial cell C1', 
                    'DC (CD103-C2)', 'Stromal cell (DN)', 'Fibroblast', 'Unresolved', 
                    'Epithelial cell C2', 'Neutrophil', 'Unresolved', 'Gamma delta T cell (Gzma+)', 
                    'Epithelial cell (low UMI count)', 'Macrophage (low UMI count)', 
                    'Doublets', 'Doublets', 'Doublets', 'Basophil', 
                    'T precursor-like cell', 'Lymphatic endothelial-like cell']

barcode_dict = {}
for barcode in Food_Allergy_cluster.index:
    label = barcode.split('_')[0].split('-')[-1]
    cond = sample_type_dict['_'.join(barcode.split('_')[1:-1])]
    region = abbr_dict[barcode.split('_')[0].split('-')[-1]]
    subsetID = Food_Allergy_cluster.loc[barcode, 'Cluster']
    celltype = cellsubset_label[subsetID-1]
    barcode = '_'.join(barcode.split('_')[1:])
    barcode_new = '%s_%s'%(label, barcode)
    barcode_dict.update({barcode:{'Barcode':barcode_new, 'Condition':cond, 
                                  'Region':region, 'SubsetID':subsetID, 
                                  'Celltype':celltype}})
barcode_info = pd.DataFrame.from_dict(barcode_dict, orient='index')

print('Healthy', sum(barcode_info['Condition']=='Healthy'))
print('Inflammation', sum(barcode_info['Condition']=='Inflammation'))
print('Lamina propria', sum(barcode_info['Region']=='Lamina propria'))
print('Peyer\'s patche', sum(barcode_info['Region']=='Peyer\'s patche'))

fname = outputs_folder / 'GSE124880_Xu_barcode_metadata.txt'
barcode_info.to_csv(fname, sep='\t', header=True, index=False, index_label='Cell barcode', encoding='utf-8')

temp = barcode_info.loc[GSE124880_barcode[0].values, ['Barcode', 'Region']].copy()
fname = data_folder / 'GSE124880_PP_LP_mm10_count_barcode_v2.tsv.gz'
temp.to_csv(fname, sep='\t', header=False, index=False, encoding='utf-8')
