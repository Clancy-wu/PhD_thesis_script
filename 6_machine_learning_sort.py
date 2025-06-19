import os
import re
import pandas as pd
import numpy as np
atlas_info = pd.read_csv('brainnetome_surface.csv')
def get_roi_name(roi_num):
    #roi_num = 'ROI-13'
    num = re.findall('\d+', roi_num)[0]
    roi_name = atlas_info[atlas_info['index']==int(num)]['name'].values
    return str(roi_name[0])
#############
func_file = 'BrainGraphRDS/CFS_func_weighted_subject_vertex.csv'
func_df = pd.read_csv(func_file)
func_df = func_df[func_df['threshold']==0.25]
df = func_df[func_df['region']==get_roi_name('ROI-193')].loc[:, ['participant_id','E.nodal']]
df.to_csv('mechine_learning/fun_ROI-193_Enodal.csv', index=None)
#############
fiber_file = 'BrainGraphRDS/CFS_fiber_weighted_subject_vertex.csv'
fiber_df = pd.read_csv(fiber_file)
fiber_df = fiber_df[fiber_df['threshold']==0.25]
df = fiber_df[fiber_df['region']==get_roi_name('ROI-29')].loc[:, ['participant_id','E.nodal.wt']]
df.to_csv('mechine_learning/fiber_ROI-29_Enodalwt.csv', index=None)
#############
thickness_file = 'Measurements/CFS_thickness_mat_exclude_Subject_019.csv'
thick_df = pd.read_csv(thickness_file)
df = thick_df.loc[:, ['subject', 'ROI-189', 'ROI-190', 'ROI-192', 'ROI-193', 'ROI-198']]
df.columns = ['participant_id', 'ROI-189', 'ROI-190', 'ROI-192', 'ROI-193', 'ROI-198']
df.to_csv('mechine_learning/thickness_exclude_Subject_019.csv', index=None)

