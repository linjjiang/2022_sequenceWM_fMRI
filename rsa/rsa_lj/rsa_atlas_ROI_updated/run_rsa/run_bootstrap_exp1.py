#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# relevant imports
import numpy as np
from scipy import io
import matplotlib.pyplot as plt
from matplotlib import rcParams
import rsatoolbox
import rsatoolbox.data as rsd # abbreviation to deal with dataset
import rsatoolbox.rdm as rsr
import os
import seaborn as sns
import sklearn as sk
import math
import pandas as pd
import pickle
import copy
import multiprocessing as mp


# In[ ]:


# change the working directory to be the timecourse data
#os.chdir('/mnt/Data1/linjdata1/vswmda/scan_data/rsa/full_GLM_mgs_0.05_50/')
os.chdir('/gpfs/scratch/linjjiang/scan_data/rsa/full_GLM_mgs_0.05_50/')


# In[ ]:


subjects = ['f09','f10','f11','f12','f15','f16','f17','f18','f19']
epochs = ['delay','response','stimulus']
exp1_subjects = ['f09','f10','f11','f12','f15','f16']
exp2_subjects = ['f17','f18','f19']


# In[ ]:


# get ROI name
# if full_GLM_mgs_xxx
order = ['area4', 'v1', 'v2', 'ips', 'fef', 'sfg', 'mfg', 'ifg',
         'ips0', 'ips1', 'ips2', 'ips3', 'ips4', 'ips5', 'spl1']

# # if full_GLM_atlas_roi
# order = ['area4-ju50',
#          'area8-hcp','area9-hcp','area9|46-hcp','area44|45|47l-hcp','fef-hcp',
         
#          'v1-wang25','v2-wang25','ips-wang15','ips0-wang15','ips1-wang15',
#          'ips2-wang15','ips3-wang15','ips4-wang15','ips5-wang15','spl1-wang15','fef-wang25',
         
#          'ipcs-md','spcs-md','amfg-md','pmfg-md','ifg-md'
#          ]


# In[ ]:


# change the following as needed
roi_labels = {'area4': 'M1', 'v1': 'V1', 'v2': 'V2', 
              
              'ips0': 'IPS0', 'ips1': 'IPS1', 'ips2': 'IPS2', 
              'ips3': 'IPS3', 'ips4': 'IPS4', 'ips5': 'IPS5', 'spl1': 'SPL1',
              
          'ips': 'IPS', 'fef': 'FEF', 'sfg': 'SFG', 
          'mfg': 'MFG', 'ifg': 'IFG'}
# roi_order = list(roi_labels.values())
subj_labels = {index: subject for index, subject in enumerate(subjects, start=0)}
epoch_labels = {0: 'stimulus', 1: 'delay', 2: 'response'}


# In[ ]:


with open('rdm_crossnobis_split_by_sess_all.pkg','rb') as f:
    RDM = pickle.load(f)


# In[ ]:


#store the subject id and roi name in dictionary, then use the dictionary to retrieve rdms to plot
subj_dict = dict()
roi_dict = dict()
epoch_dict = dict()
sess_dict = dict()

for idx, rdm in enumerate(RDM):
    subj = rdm.rdm_descriptors['subj'][0]
    roi = rdm.rdm_descriptors['roi'][0]
    epoch = rdm.rdm_descriptors['epoch']
    sess = rdm.rdm_descriptors['session']
    
    if subj not in subj_dict:
        subj_dict[subj] = []
    subj_dict[subj].append(idx)
    
    if roi not in roi_dict:
        roi_dict[roi] = []
    roi_dict[roi].append(idx)

    if epoch not in epoch_dict:
        epoch_dict[epoch] = []
    epoch_dict[epoch].append(idx)

    if sess not in sess_dict:
        sess_dict[sess] = []
    sess_dict[sess].append(idx)


# # Experiment 1

# In[ ]:


import numpy as np

# Define conditions
conditions = ['3WL', '3TL', '3WR', '3TR', '5WL', '5TL', '5WR', '5TR']

def create_orthogonal_rdm():
    size = len(conditions)
    rdms = {key: np.zeros((size, size)) for key in ['LR', 'WT', 'Ecc', 'LR-WT', 'Ecc-LR', 'Ecc-WT', 'All']}

    for i, cond1 in enumerate(conditions):
        for j, cond2 in enumerate(conditions):
            ecc1, shape1, pattern1 = cond1[0], cond1[1], cond1[2]
            ecc2, shape2, pattern2 = cond2[0], cond2[1], cond2[2]

            # Individual differences
            if pattern1 != pattern2 and shape1 == shape2 and ecc1 == ecc2:
                rdms['LR'][i, j] = 1
            if shape1 != shape2 and pattern1 == pattern2 and ecc1 == ecc2:
                rdms['WT'][i, j] = 1
            if ecc1 != ecc2 and shape1 == shape2 and pattern1 == pattern2:
                rdms['Ecc'][i, j] = 1

            # Two-way interaction not covered by main effects
            if pattern1 != pattern2 and shape1 != shape2 and ecc1 == ecc2:
                rdms['LR-WT'][i, j] = 1
            if ecc1 != ecc2 and pattern1 != pattern2 and shape1 == shape2:
                rdms['Ecc-LR'][i, j] = 1
            if ecc1 != ecc2 and shape1 != shape2 and pattern1 == pattern2:
                rdms['Ecc-WT'][i, j] = 1

            # Three-way interaction (All)
            if ecc1 != ecc2 and shape1 != shape2 and pattern1 != pattern2:
                rdms['All'][i, j] = 1

    return rdms

model_rdms = create_orthogonal_rdm()

# # Example of how to check orthogonality
# for key1 in model_rdms:
#     for key2 in model_rdms:
#         if key1 != key2:
#             dot_product = np.dot(model_rdms[key1].flatten(), model_rdms[key2].flatten())
#             print(f"Dot product between {key1} and {key2}: {dot_product}")


# In[ ]:


model_names=[]
rdm_arrays=[]
for key,value in model_rdms.items():
    model_names.append(key)
    rdm_arrays.append(value)
result = np.stack(rdm_arrays, axis=0)
print(result.shape)


# In[ ]:


model_rdms = rsatoolbox.rdm.RDMs(result,
                            rdm_descriptors={'model':model_names},
                            pattern_descriptors={'conds':conditions},
                            dissimilarity_measure='Crossnobis')


# In[ ]:


models_fixed = []
for i in range(len(model_names)):
    models_fixed.append(rsatoolbox.model.ModelFixed(model_names[i], model_rdms[i]))
print(models_fixed)


# In[ ]:





# In[ ]:


def run_bootstrap(args):
    epoch, roi, sess = args
    print(epoch,roi,sess)
    rdms = []
    for subj in exp1_subjects:
        idx_set = set(epoch_dict[epoch]) & set(roi_dict[roi]) & set(sess_dict[sess]) & set(subj_dict[subj])
        idx = list(idx_set) #set(subj_dict[subj]) &
        #print(subj,epoch,roi,sess)
        #print(len(idx))
        #print(idx)
        if len(idx) > 0:
            rdm = RDM[idx[0]]
            rdm.rdm_descriptors.pop('name',None)
            rdm.rdm_descriptors.pop('subj',None)

            temp = rdm.rdm_descriptors['session']
            if isinstance(temp, list)==False:
                rdm.rdm_descriptors['session'] = [rdm.rdm_descriptors['session']]

            temp = rdm.rdm_descriptors['epoch']
            if isinstance(temp, list)==False:
                rdm.rdm_descriptors['epoch'] = [rdm.rdm_descriptors['epoch']]

            temp = rdm.rdm_descriptors['session_name']
            if isinstance(temp, list)==False:
                rdm.rdm_descriptors['session_name'] = [rdm.rdm_descriptors['session_name']]

            rdms.append(rdm)

    rdm_test = rsatoolbox.rdm.rdms.concat(rdms)

#             results = rsatoolbox.inference.eval_bootstrap_rdm(models_fixed, 
#                                             rdm_test, method='cosine',
#                                             N=1000,boot_noise_ceil=True) #rsatoolbox.inference.eval_bootstrap_rdm(models, rdms_data, method='corr')
    results = rsatoolbox.inference.eval_dual_bootstrap(
            models_fixed, rdm_test, method='cosine', fitter=None,
            k_pattern=2, k_rdm=2, N=10000, n_cv=2,
            pattern_descriptor='index', rdm_descriptor='index',
            use_correction=True)

#     df_unique = pd.concat([df_unique,pd.DataFrame({
#                 'results': [results], 
#                 'roi': [roi], 
#                 'session': [sess],
#                 'epoch': [epoch]
#             })], axis=0) #

    return pd.DataFrame({
        'results': [results],
        'roi': [roi],
        'session': [sess],
        'epoch': [epoch]
    })
            


# In[ ]:


# spawn a pool of parallel workers
p = mp.Pool(processes=96)

import itertools
combinations = itertools.product(epoch_dict.keys(), roi_dict.keys(), sess_dict.keys())

results = p.map(run_bootstrap, combinations)


# In[ ]:


# close the pool of workers
p.close()


# In[ ]:


# Concatenate the results into a single dataframe
df_unique = pd.concat(results, axis=0)


# In[ ]:


with open('crossnobis_dual_bootstrap_10000_exp1.pkg','wb') as f:
    pickle.dump(df_unique,f)   


# In[ ]:


import pandas as pd

# Initialize an empty list to store the extracted data
data = []

# Model names
models = ['LR', 'WT', 'Ecc', 'LR-WT', 'Ecc-LR', 'Ecc-WT', 'All']

# Iterate through each row in df_unique
for index, row in df_unique.iterrows():
    results_object = row['results']
    roi = row['roi']
    session = row['session']
    epoch = row['epoch']
    
    # Extract the test results
    pairwise_results = results_object.test_all()[0]
    test_against_zero = results_object.test_all()[1]
    test_against_noise_ceiling = results_object.test_all()[2]
    eval_mean = results_object.get_means()
    eval_sem = results_object.get_sem()
    noise_ceil = results_object.get_noise_ceil()
    
    # Store pairwise test results
    for i, model1 in enumerate(models):
        for j, model2 in enumerate(models):
            data.append({
                'roi': roi,
                'session': session,
                'epoch': epoch,
                'model': f'{model1} vs {model2}',
                'comparison': 'pairwise',
                'p_value': pairwise_results[i, j],
                'mean': [eval_mean[i],eval_mean[j]],
                'sem': [eval_sem[i],eval_sem[j]],
                'noise_ceil': noise_ceil
                
            })
    
    # Store test against zero results
    for i, model in enumerate(models):
        data.append({
            'roi': roi,
            'session': session,
            'epoch': epoch,
            'model': model,
            'comparison': 'against_zero',
            'p_value': test_against_zero[i],
                'mean': eval_mean[i],
                'sem': eval_sem[i],
                'noise_ceil': noise_ceil
        })
    
    # Store test against noise ceiling results
    for i, model in enumerate(models):
        data.append({
            'roi': roi,
            'session': session,
            'epoch': epoch,
            'model': model,
            'comparison': 'against_noise_ceiling',
            'p_value': test_against_noise_ceiling[i],
                'mean': eval_mean[i],
                'sem': eval_sem[i],
                'noise_ceil': noise_ceil
        })

# Convert the list of dictionaries to a DataFrame
df_results = pd.DataFrame(data)

df_results[df_results['p_value'] <= 0.05].to_csv('crossnobis_dual_bootstrap_10000_exp1.csv')

