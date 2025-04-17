#!/usr/bin/env python
# coding: utf-8

# In[49]:


import numpy as np
import multiprocessing as mp
import os
import re
import nibabel as nib
import pickle


# In[50]:


# where are the residuals?
residual_path = '/gpfs/scratch/linjjiang/scan_data/spm/output/spatiotemporal_res'


# In[51]:


def list_subfolders(directory):
    subfolders = [f.name for f in os.scandir(directory) if f.is_dir()]
    paths =  [f.path for f in os.scandir(directory) if f.is_dir()]
    return subfolders,paths


# In[52]:


def extract_folder_name(path, target_folder):
    # Split the path into parts
    path_parts = path.split(os.sep)
    
    # Find the index of the target folder
    if target_folder in path_parts:
        index = path_parts.index(target_folder)
        
        # Return the next part of the path if it exists
        if index + 1 < len(path_parts):
            return path_parts[index + 1]
    
    return None


# In[53]:


def extract_subject_session(path):
    # Split the path into parts
    path_parts = path.split(os.sep)
    
    # Find the subject and session parts
    subject = None
    session = None
    
    if 'spatiotemporal_res' in path_parts:
        idx = path_parts.index('spatiotemporal_res')
        if idx + 2 < len(path_parts):
            subject = path_parts[idx + 1]
            session = path_parts[idx + 2]
    
    return subject, session


# In[54]:


def list_matching_files(directory):
    # Define the regular expression pattern
    pattern = re.compile(r'^Res_\d{4}\.nii$')
    
    # List all files in the directory
    all_files = os.listdir(directory)
    
    # Filter files that match the pattern
    matching_files = [os.path.join(directory, f) for f in all_files if pattern.match(f)]
    
    # Sort the files by name
    matching_files = sorted(matching_files)
    
    return matching_files


# In[55]:


subjects = ['f18','f19'] #
#

# In[56]:


final_paths = []
for subject in subjects:
    folder,path = list_subfolders(os.path.join(residual_path,subject))
    for pt in path:
        specific_file_path = os.path.join(pt, 'smgs', 'SPM.mat')
        if os.path.exists(specific_file_path):
            final_paths.append(os.path.join(pt, 'smgs'))


# In[57]:


print(len(final_paths))
print(final_paths[0])


# In[58]:


# spawn a pool of parallel workers
p = mp.Pool(processes=28)


# In[59]:


for final_path in final_paths:
    files = list_matching_files(final_path)
    subject, session = extract_subject_session(final_path)

    print(subject,session)

    # Initialize a list to hold the residuals
    residuals = []
    # Load each residual image and add to the residuals list
    for file in files:  # assuming you have 2000 residual images
        residual_img = nib.load(file)
        residual_data = residual_img.get_fdata()
        residuals.append(residual_data)  # flatten to a 1D array .flatten()
    # Convert the list of residuals to a 2D numpy array (T x V)
    residuals = np.array(residuals)

    # save residual file
    #/gpfs/scratch/linjjiang/scan_data/rsa/residual_from_spm
    with open(os.path.join('/gpfs/scratch/linjjiang/scan_data/rsa/residual_from_spm',
                   'res_'+subject+'_'+session+'.pkg'),'wb') as f:
        pickle.dump(residuals,f)
            


# In[ ]:


# close the pool of workers
p.close()

