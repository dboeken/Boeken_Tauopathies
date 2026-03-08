[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15635269/.svg)](https://doi.org/10.5281/zenodo.15635269)

# Boeken_Tauopathies

This repository contains the analysis code associated with the Tauopathies single-molecule characterisation project, led by Dorothea Böken. 

The associated manuscript is now published [here](https://www.sciencedirect.com/science/article/pii/S2211124726000124) as:

**Böken, D. et al. (2026). _Nanoscopic tau aggregates are not shared intermediates but disease-specific entities across tauopathies._ Cell Reports, 45(2), 116934.**  

The original preprint is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.10.658934v1).


## Prerequisites

This analysis assumes a standard installation of Python 3 (v 3.10.5). For specific package requirements, see the environment.yml file, or  create a new conda environment containing all packages by running ```conda create -f environment.yml```. In addition to the analysis contained here, some simple statistical tests were performed using [GraphPad Prism v 10.4](https://www.graphpad.com/scientific-software/prism/).

## Raw data

For convenience, example image files are provided here under the ```data``` folder. These data may be used to explore the workflows presented here as described below.


## Workflow

Example raw images are provided here within the ```data``` folder to test the included analysis scripts.

Individual analyses are presented within the ```src``` folder. Where processing order is important for individual analyses, scripts have been numbered and should be run in order before unnumbered counterparts.

