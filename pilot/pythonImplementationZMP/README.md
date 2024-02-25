TODO: add any details

Aim of this code is to generate ZMP GRF predictions from a kinematic solution using core OpenSim processes. The idea here should be to use as similar code to opensim C++ implementation of these calculations as possible to make it easy to transfer down the track.

For the run, uses processed data from dynamic consistency quest paper from Moco Tracking, given this was the most dynamically consistent of the bunch.

For the cut, uses RRA processed data from the anticipated vs. unanticipated project (PW10, A11) — which may contain greater residuals than the processed running data.

Uses OpenSim 4.4 with Python 3.10