This folder contains RSA scripts for Linjing's fMRI dissertation project, using the atlas ROI
Date: 09/10/2024
Contact: Linjing Jiang, linjjiang07@gmail.com

1. make_rsa_figures: make RSA figures - naming should be straightforward, e.g., make boostrap tables, make crossnibis figures (for combined session or each session), make correlation figures etc.

2. run_rsa: run RSA with the following order
- Step 0 and 0.1: calculate noise from GLM residual
- Step 1: Make RSA dataset
- Step 2: Calculate simple correlation
- Step 3: Calculate crossnobis distance and RDM
- Step 4: Model inference
- Step 6: rsa for eye movements (exploratory analysis, not included in my dissertation)
- run_bootstrap_exp1/2.slurm/py/ipynb: run bootstrapping on HPC using parallel computing. 

3. not_used: not used




