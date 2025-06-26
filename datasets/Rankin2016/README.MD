# Rankin2016

*Avian* / *Gait* / *Overground*

The data in this folder comes from a [publicly available dataset](https://datadryad.org/stash/dataset/doi:10.5061/dryad.fh3h6) released with the paper:

**Rankin JW, Rubenson J & Hutchinson JR (2016). Inferring muscle functional roles of the ostrich pelvic limb during walking and running using computer optimization. *J R Soc Interface*, 13: 20160035. doi: [10.1098/rsif.2016.0035](https://doi.org/10.1098/rsif.2016.0035)** 

Any additional use must be in line with that described in the original paper.

## Code

The Python scripts in this folder were developed using OpenSim v4.4 alongside Python 3.10 on Windows. The code may work with other versions or operating systems, but has not been tested for this. The script uses relative file paths to extract the data, so it must be run from this base directory.

## Steps for Extracting Data

The following steps describe the processes and code used to generate the data used for test cases related to the overall aim of this work (i.e. integrating ZMP estimates of GRFs into OpenSim). The final data is already available in this folder (see **Data** section below), however the processes to produce these data can be replicated using these steps.

1. **<u>Not recommended</u>**, but you can re-download the original data from [this link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.fh3h6) and separate the `.mot` files into a **data** folder and the `.osim` files into a *models* folder. For subsequent steps and scripts to work, the compliant tendon running ostrich model must be renamed to `Rankin2016_Ostrich.osim` and the illegal path points (those that contain a + symbol) be renamed — and be placed in the **model** folder. The **data** folder should contain four files, those being the two `.mot` files for each of the walking and running trials.
2. Run the **`1_cleanData.py`** script. This will create some new files for use in later scripts, including a model with muscles removed, states kinematics files, and GRF files that re-express the centre of pressure in the ground instead of body frames and enforces the ground plane to be at y = 0.
3. Run the **`2_processData.py`** script. ***TODO: add details on what this does...***

## Data

This dataset includes walking and running data from a single ostrich (*Stuthio camelus*) specimen, processed using an OpenSim Moco state tracking approach. The idea behind this dataset is the demonstrate how ZMP estimates of GRFs work in: (1) a non-human context; and (2) when an incomplete multi-body model is used (i.e. head and neck missing). 

Note that the original dataset included multiple versions of the ostrich model with different muscle-tendon unit parameters (e.g. rigid vs. non-rigid tendon for walking and running). As muscle-driven simulations were not a focus of this work, only the compliant tendon running model was retained and used.

***TODO: add any final descriptors of data folders...***

