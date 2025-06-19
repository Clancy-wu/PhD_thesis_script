import os
from shutil import copy
import pandas as pd

org_dir = '/home/clancy/Desktop/cfs_func_fiber_network/NetResults_invnodal'
new_dir = 'NetResults_invnodal'
# 152
file_152 = pd.read_csv('all_subs_152_info.csv')
subs = file_152['participant_id'].values

def generate_mat_from_name(sub_name):
    sub_org_dir = os.path.join(org_dir, 'CFS', sub_name)
    sub_org_func = os.path.join(sub_org_dir, f'raw_func_{sub_name}_invnodal.txt')
    sub_org_fiber = os.path.join(sub_org_dir, f'raw_fiber_{sub_name}_invnodal.txt')
    ## only surface, 210
    # func
    sub_func_df = pd.read_csv(sub_org_func, sep=' ', header=None)
    sub_func = sub_func_df.iloc[0:210, 0:210]
    sub_func_name = os.path.join(new_dir, sub_name, f'raw_func_{sub_name}_invnodal.txt')
    os.makedirs(os.path.dirname(sub_func_name), exist_ok=True)
    sub_func.to_csv(sub_func_name, sep=' ', index=None, header=None)
    # fiber
    sub_fiber_df = pd.read_csv(sub_org_fiber, sep=' ', header=None)
    sub_fiber = sub_fiber_df.iloc[0:210, 0:210]
    sub_fiber_name = os.path.join(new_dir, sub_name, f'raw_fiber_{sub_name}_invnodal.txt')
    sub_fiber.to_csv(sub_fiber_name, sep=' ', index=None, header=None)
    return 0

for sub in subs:
    a = generate_mat_from_name(sub)
    print(a)


##################################
import os
from shutil import copy
import pandas as pd
import nibabel as nib
import numpy as np

org_dir = '/home/clancy/Desktop/doctoral_thesis/NetResults_invnodal'
new_dir = 'BrainGraphRDS'
# 146
file_146 = pd.read_csv('subs_146_info.csv')
file_146 = file_146[file_146['subject']!='Subject_019']
## sub-sub059 exlucded, sub-sub059 = Subject_019, so sub-sub019 excluded
subs = file_146['participant_id'].values

def create_thick_coupling_mat(sub_name):
    ## lh thickness
    sub_lh_file = f'/home/clancy/data/BigData/chronic_fatigue/fsaverage_thickness_6/{sub_name}_space-fsaverage_hemi-L.thickness.gii'
    sub_lh = nib.load(sub_lh_file) # sub_lh.darrays[0].data # 163842
    sub_lh_data = sub_lh.darrays[0].data
    lh_atlas_file = '/home/clancy/TemplateFlow/BN_Atlas_freesurfer/fsaverage/label/lh.BN_Atlas.annot'
    lh_atlas = nib.freesurfer.read_annot(lh_atlas_file) # lh_atlas[0] # 1,3...209 # 163842
    lh_atlas_data = lh_atlas[0]
    roi_results = []
    roi_results.append(('subject', sub_name))
    for i in range(1,211):
        i_roi = 'ROI-'+str(i)
        i_mean = np.mean(sub_lh_data[lh_atlas_data==i])
        roi_results.append((i_roi, str(i_mean)))
    roi_lh = pd.DataFrame(roi_results)
    roi_lh = roi_lh[roi_lh[1]!='nan']
    ## rh thickness
    sub_rh_file = f'/home/clancy/data/BigData/chronic_fatigue/fsaverage_thickness_6/{sub_name}_space-fsaverage_hemi-R.thickness.gii'
    sub_rh = nib.load(sub_rh_file) # sub_rh.darrays[0].data # 163842
    sub_rh_data = sub_rh.darrays[0].data
    rh_atlas_file = '/home/clancy/TemplateFlow/BN_Atlas_freesurfer/fsaverage/label/rh.BN_Atlas.annot'
    rh_atlas = nib.freesurfer.read_annot(rh_atlas_file) # rh_atlas[0] # 163842
    rh_atlas_data = rh_atlas[0]
    roi_results = []
    for i in range(1,211):
        i_roi = 'ROI-'+str(i)
        i_mean = np.mean(sub_rh_data[rh_atlas_data==i])
        roi_results.append((i_roi, str(i_mean)))
    roi_rh = pd.DataFrame(roi_results)
    roi_rh = roi_rh[roi_rh[1]!='nan']
    roi_both = pd.concat((roi_lh, roi_rh), axis=0)
    return roi_both

df_example = create_thick_coupling_mat(subs[0])
df_out = df_example[0]
for sub_name in subs:
    df = create_thick_coupling_mat(sub_name)
    df_out = pd.concat((df_out, df[1]), axis=1)

desire_order = ['subject']+['ROI-'+str(i) for i in range(1,211)]
df_sorted = df_out.set_index(0).reindex(desire_order).reset_index()
df_sorted.T.to_csv('BrainGraphRDS/CFS_thickness_mat_exclude_Subject_019.csv', header=False, index=False)

##################################
from scipy.stats import spearmanr

org_dir = '/home/clancy/Desktop/doctoral_thesis/NetResults_invnodal'
new_dir = 'BrainGraphRDS'
# 146
file_146 = pd.read_csv('subs_146_info.csv')
file_146 = file_146[file_146['subject']!='Subject_019']
## sub-sub059 exlucded, sub-sub059 = Subject_019, so sub-sub019 excluded
subs = file_146['participant_id'].values

def create_thick_coupling_mat(sub_name):
    sub_func_file = f'/home/clancy/Desktop/doctoral_thesis/NetResults_invnodal/{sub_name}/raw_func_{sub_name}_invnodal.txt'
    sub_func = pd.read_csv(sub_func_file, sep=' ', header=None)
    sub_fiber_file = f'/home/clancy/Desktop/doctoral_thesis/NetResults_invnodal/{sub_name}/raw_fiber_{sub_name}_invnodal.txt'
    sub_fiber = pd.read_csv(sub_fiber_file, sep=' ', header=None)
    roi_results = []
    for i in range(210):
        i_cor = spearmanr(sub_func.iloc[i,:], sub_fiber.iloc[i,:])
        i_z = np.arctanh(i_cor[0])
        roi_results.append(i_z)
    sub_out = np.append(sub_name, roi_results)
    return sub_out

row_name = ['subject']+['ROI-'+str(x) for x in range(1,211)]
df_out = pd.DataFrame(row_name)
for sub_name in subs:
    sub_out = pd.DataFrame(create_thick_coupling_mat(sub_name))
    df_out = pd.concat((df_out, sub_out), axis=1)

df_out.T.to_csv('BrainGraphRDS/CFS_coupling_mat_exclude_Subject_019.csv', header=False, index=False)
