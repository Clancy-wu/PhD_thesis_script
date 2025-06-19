import os
from re import findall
from nilearn import image
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from itertools import product
from glob import glob
cfs_dir = '/home/clancy/data/BigData/chronic_fatigue'

def run(f, this_iter):
    with ProcessPoolExecutor(max_workers=14) as executor:
        results = list(tqdm(executor.map(f, this_iter), total=len(this_iter)))
    return results

def get_roi_alff(sub_file):
    # ROI-29, ROI-115, ROI-189, ROI-190, ROI-192, ROI-193, ROI-198
    roi_index = [29, 115, 189, 190, 192, 193, 198]
    sub_name = findall(r'(sub-sub\d+)_task-rest', os.path.basename(sub_file))[0]
    atlas_file = f'{cfs_dir}/tpl-MNI152NLin2009cAsym_res-02_atlas-BN_desc-246_dseg.nii.gz'
    atlas_img = image.get_data(atlas_file)
    sub_img = image.get_data(sub_file)
    roi_results = []
    for i in roi_index:
        alff = np.mean(sub_img[atlas_img == i])
        roi_results.append(alff)
    return np.append(sub_name, np.array(roi_results))    

all_files = glob(f'{cfs_dir}/xcp_d_155/sub-sub*/func/sub-sub*_task-rest_space-MNI152NLin2009cAsym_res-2_desc-denoisedSmoothed_bold.nii.gz')
future_results = run(get_roi_alff, all_files)
alff_results = pd.DataFrame(future_results)
alff_results.columns = ['participant_id', 'ROI-29', 'ROI-115', 'ROI-189', 'ROI-190', 'ROI-192', 'ROI-193', 'ROI-198']
clinic_info = pd.read_csv('subs_146_info.csv')
roi_df = pd.merge(clinic_info, alff_results, on='participant_id', how='left')
roi_df.to_csv('alff_results.csv', index=None)

## compute brain size
def get_sub_brain_size(sub_file):
    mask_img = image.load_img(sub_file)
    voxel_size = mask_img.header.get_zooms()
    brain_size = np.around(np.sum(mask_img.get_fdata()==1) * voxel_size[0] * voxel_size[1] * voxel_size[2], 3)
    sub_name = findall(r'(sub-sub\d+)_desc-brain', os.path.basename(sub_file))[0]
    return np.append(sub_name, brain_size)

all_files = [i for i in glob(f'{cfs_dir}/fmriprep_155/sub-sub*/anat/sub-sub*_desc-brain_mask.nii.gz') if 'MNI152NLin2009cAsym' not in i]
future_results = run(get_sub_brain_size, all_files)
df = pd.DataFrame(future_results, columns=['participant_id', 'brain_size'])
df.to_csv('Measurements/whole_brain.csv', index=None)


## compute gray matter volume
import ants
def get_sub_brain_size(sub_name):
    # sub_file
    bn_file = f'{cfs_dir}/tpl-MNI152NLin2009cAsym_res-02_atlas-BN_desc-246_dseg.nii.gz'
    bn_img = ants.image_read(bn_file)
    sub_trans_file = f'{cfs_dir}/fmriprep_155/{sub_name}/anat/{sub_name}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5'
    sub_t1_file = f'{cfs_dir}/fmriprep_155/{sub_name}/anat/{sub_name}_desc-preproc_T1w.nii.gz'
    sub_t1 = ants.image_read(sub_t1_file)
    sub_trans = ants.read_transform(sub_trans_file)
    bn_native = sub_trans.apply_to_image(bn_img, reference=sub_t1, interpolation='multilabel')
    bn_spacing = bn_native.spacing
    voxel_size = bn_spacing[0] * bn_spacing[1] * bn_spacing[2]
    gm_size = np.around(np.sum(bn_native.numpy() > 0) * voxel_size, 3)
    return np.append(sub_name, gm_size)

all_subs = [i for i in os.listdir(f'{cfs_dir}/BIDS_155') if 'sub-sub' in i]
future_results = []
i=1
for sub_name in all_subs:
    print(i)
    future_results.append(get_sub_brain_size(sub_name))
    i += 1
df = pd.DataFrame(future_results, columns=['participant_id', 'gm_size'])
df.to_csv('Measurements/gray_matter_size.csv', index=None)
