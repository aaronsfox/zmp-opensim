# SchreiberMoissenet2019

*Human* / *Walking* / *Overground*

The data in this folder comes from a [publicly available dataset](https://doi.org/10.6084/m9.figshare.7734767) released with the paper:

**Schreiber C & Moissenet F (2019). A multimodal dataset of human gait at different walking speeds established on injury-free adult participants. *Sci Data*, 6: 111. doi: [10.1038/s41597-019-0124-4](https://doi.org/10.1038/s41597-019-0124-4)**

Any additional use must be in line with that described in the original paper.

## Code

The Python scripts in this folder were developed using OpenSim v4.4 alongside Python 3.10 on Windows. The code may work with other versions or operating systems, but has not been tested for this. The script uses relative file paths to extract the data, so it must be run from this base directory.

## Model

The model included in this folder and used to process data is a modified version from the paper:

**Rajagopal A et al. (2016). Full-Body Musculoskeletal Model for Muscle-Driven Simulation of Human Gait. *IEEE Trans Biomed Eng*, 63: 2068-2079. doi: [10.1109/tbme.2016.2586891](https://doi.org/10.1109/tbme.2016.2586891)**

The modifications to the model were:

- Removing the muscles as all simulations of these data will be torque-driven
- Adapting the marker set to match that of the experimental data being processed, with some markers on the scapulae out due to not being tracked in simulations.

## Steps for Extracting Data

The following steps describe the processes and code used to generate the data used for test cases related to the overall aim of this work (i.e. integrating ZMP estimates of GRFs into OpenSim). The final data is already available in this folder (see **Data** section below), however the processes to produce these data can be replicated using these steps.

1. Download the archive of participant data from [this link](https://doi.org/10.6084/m9.figshare.7734767). 
2. Create a `raw` folder in this directory and unzip each participants data into this folder. Each participant should therefore have their own folder within the `raw` directory (e.g. `raw\\2014001`) that contains the `.c3d` files from Schreiber & Moissenet (2019). The participant folder should simply contain all of the individual `.c3d` files for that participant to be appropriately extracted (e.g. `raw\\2014001\\2014001_C1_01.c3d`).
3. Run the **`1_extractData.py`** script. This will extract and convert the necessary experimental `.c3d` files from the raw folder into individual participant folders within this directory. Note that only the 0—0.4m·s<sup>-1</sup> and 0.8—1.2m·s<sup>-1</sup> speeds are extracted. After this you may delete the `raw` folder as it will no longer be used.
4. Run the **`2_processData.py`** script. ***TODO: add details on this processes...***. **Note that this processing can take a significant amount of time, so it is not necessarily recommended to repeat**.

## Data

This dataset includes overground walking at various speeds, processed using an OpenSim Moco marker tracking approach. The idea behind this dataset is to demonstrate how ZMP estimates of GRFs are impacted by: (1) walking at different speeds (i.e. 0—0.4m·s<sup>-1</sup> and 0.8—1.2m·s<sup>-1</sup>); and (2) a double support movement.

***TODO: add descriptors of data folders...***

***TODO Notes:***

- *Only participants with 3 'good' trials for each speed selected*
  - *A 'good' trial indicated by the foot contacts on both plates*