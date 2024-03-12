# HamnerDelp2013

*Human* / *Running* / *Treadmill*

The data in this folder comes from a [publicly available dataset](https://simtk.org/projects/nmbl_running) released with the paper:

**Hamner SR & Delp SL (2013). Muscle contributions to fore-aft and vertical body mass center accelerations over a range of running speeds. *J Biomech*, 46(4): 780-787. doi: [10.1016/j.jbiomech.2012.11.024](https://doi.org/10.1016/j.jbiomech.2012.11.024)** 

Any additional use must be in line with that described in the original paper.

## Code

The Python scripts in this folder were developed using OpenSim v4.4 alongside Python 3.10 on Windows. The code may work with other versions or operating systems, but has not been tested for this. The script uses relative file paths to extract the data, so it must be run from this base directory.

## Steps for Extracting Data

The following steps describe the processes and code used to generate the data used for test cases related to the overall aim of this work (i.e. integrating ZMP estimates of GRFs into OpenSim). The final data is already available in this folder (see **Data** section below), however the processes to produce these data can be replicated using these steps.

1. Download the individual subject datasets from [this link](https://simtk.org/projects/nmbl_running). This will result in separate `.zip` folders for `subject01`, `subject02`, `subject03`, `subject04`, `subject08`, `subject11`, `subject17`, `subject19` and `subject20`.
2. Create a `raw` folder in this directory and unzip each participants data into this folder. Each participant should therefore have their own folder within the `raw` directory (e.g. `raw\\subject01`) that contains the experimental data from Hamner & Delp (2013). The participant folder should have the data folders directly underneath it to be appropriately extracted (e.g. `raw\\subject01\\ExportedData`).
3. Run the **`1_extractData.py`** script. This will extract the necessary experimental data files from the raw folder into individual participant folders within this directory. Note that only the 3.0m·s<sup>-1</sup> and 5.0m·s<sup>-1</sup> datasets are extracted. After this you may delete the `raw` folder as it will no longer be used.
4. Run the **`2_processData.py`** script. This will generate torque-driven simulations of the running gait cycles at various speeds using two different approaches (i.e. *RRA* and *MocoTrack*) designed to leave varying degrees of residual forces and moments present in the data. **Note that this processing can take a decent amount of time, so it is not necessarily recommended to repeat**.

## Data

This dataset includes treadmill running at various speeds, processed with two different techniques (i.e. *RRA* and *MocoTrack*) designed to leave variable remaining levels of residual forces and moments in the dataset. The idea here is to demonstrate how ZMP estimates of GRFs may be impacted by: (1) treadmill running at different speeds (i.e. 3.0m·s<sup>-1</sup> vs. 5.0m·s<sup>-1</sup>); and (2) higher levels of residual forces and moments in the data. Each participant's processed data will end up in a folder labelled with the subject ID. Within these folders, data is separated into:

- **expData:** original `.trc` marker files, `.mot` and `.xml` files from running trials. This folder also includes a `.pkl` file created to indicate the timings of gait cycles in the data.
- **ik:** data from the original [Hamner & Delp paper](https://doi.org/10.1016/j.jbiomech.2012.11.024) processed using OpenSim's inverse kinematics tool.
- **moco:** set-up files and outputs from the Moco tracking approach. Some files from the expData folder are copied across for ease of use with the Moco tracking tool. Outputs from different running speeds and gait cycles are separated into folders. The folder also contains a `.pkl` file that records the computational time for each Moco track iteration.
- **model:** static trial data and scaled model from the original [Hamner & Delp paper](https://doi.org/10.1016/j.jbiomech.2012.11.024).
- **rra:** set-up files and outputs from the residual reduction algorithm (RRA) approach. Outputs from different running speeds and gait cycles are separated into folders. The folder also contains `.pkl` files that record the computational time and recommended model mass adjustments for each RRA iteration.