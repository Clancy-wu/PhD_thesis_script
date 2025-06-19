#/bin/python
import os
import pandas as pd
from scipy.stats import permutation_test
import numpy as np

clinic_info_73 = pd.read_csv('subs_146_info.csv')
brainnetome_surface = pd.read_csv('brainnetome_surface.csv')

func_vertex = pd.read_csv('BrainGraphRDS/CFS_thickness_scn_vertex.csv')

select_col = 'dist'
net_df = func_vertex

def mean_diff(x, y):
    return np.mean(x) - np.mean(y)

def perm_test(group1, group2):
    try:
        res = permutation_test((group1, group2), mean_diff, permutation_type='independent', n_resamples=1000)
        res_out = np.around(res.pvalue,8)
    except:
        res_out = 1
    return res_out

def corr_between_net_behavior_vertex(net_df, select_col):
    net_select_df = net_df.loc[:, [select_col,'region','all_group','threshold']]
    net_select_df.columns = ['net_value', 'region','all_group','threshold']
    merge_df = pd.merge(brainnetome_surface.loc[:,['name','index']], net_select_df, left_on='name', right_on='region')
    merge_df = merge_df.sort_values(['index', 'threshold'])
    patient_before = merge_df[merge_df['all_group']=='patient_before']
    patient_after = merge_df[merge_df['all_group']=='patient_after']
    health_before = merge_df[merge_df['all_group']=='health_before'] 
    health_after = merge_df[merge_df['all_group']=='health_after']

    ## all regions
    all_regions = brainnetome_surface['name'].values
    p_before_h_before_p = []
    p_before_p_after_p = []
    h_before_h_after_p = []

    for i in range(210):
        i_name = all_regions[i]

        i_group1 = patient_before[patient_before['name']==i_name]['net_value'].values
        i_group2 = health_before[health_before['name']==i_name]['net_value'].values
        p_before_h_before_p.append(perm_test(i_group1, i_group2))

        i_group1 = patient_before[patient_before['name']==i_name]['net_value'].values
        i_group2 = patient_after[patient_after['name']==i_name]['net_value'].values
        p_before_p_after_p.append(perm_test(i_group1, i_group2))

        i_group1 = health_before[health_before['name']==i_name]['net_value'].values
        i_group2 = health_after[health_after['name']==i_name]['net_value'].values
        h_before_h_after_p.append(perm_test(i_group1, i_group2))

    out_df = pd.DataFrame({
        'p_before_h_before_p': p_before_h_before_p,
        'p_before_p_after_p': p_before_p_after_p,
        'h_before_h_after_p': h_before_h_after_p
    })
    return(out_df)