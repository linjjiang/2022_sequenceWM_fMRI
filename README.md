# Representational Geometry of Visuospatial Sequences: An fMRI Study

This is an fMRI Study I conducted as part of my Ph.D. dissertation: [Behavioral and Neural Correlates of Serial-Order Spatial Working Memory](https://www.proquest.com/docview/3123663226?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses)

We are preparing a manuscript to be submitted to a peer reviewed journal.

## Description

This study explored how the implicit structure of visuospatial sequences is encoded and maintained in the dorsal visual pathway.

In two fMRI experiments, we collected fMRI data from nine healthy young adults during a serial-order memory-guided saccade task with varying structures of visuospatial sequences. We applied representational similarity analysis on GLM-modeled fMRI data, integrated bootstrapped cross-validation and noise covariance modeling for robust pattern estimation and generalization. 

This repo contains the scripts used to perform RSA on fMRI timeseries, after preprocessing and GLM:
* make_roi: make different sets of ROI for RSA
* run_rsa: main scripts for running RSA.
* glm_example: example scripts for running General Linear Model before RSA.

## Getting Started

### Dependencies

* ubuntu 24.04
* Python 3.10 with:
	* NumPy, Pandas, Scipy
	* rsatoolbox: https://rsatoolbox.readthedocs.io/en/stable/
	* NiBabel: https://nipy.org/nibabel/
	* Nilearn: https://nilearn.github.io/dev/index.html#

### Installing

To be updated...

### Executing program

To be updated...

## Help & Potential Issues


## Authors

Linjing Jiang
[@linjjiang](https://github.com/linjjiang)

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MTI License - see the LICENSE.md file for details

## Acknowledgments

* [rsatoolbox](https://rsatoolbox.readthedocs.io/en/stable/)
* [NiBabel](https://nipy.org/nibabel/)
* [Nilearn](https://nilearn.github.io/dev/index.html#)
