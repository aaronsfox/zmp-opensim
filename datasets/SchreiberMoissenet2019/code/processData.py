# -*- coding: utf-8 -*-
"""

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code processes the data from the SchreiberMoissenet2019 dataset through
    an OpenSim Moco coordinate tracking framework. A torque-driven model is used to
    produce kinematics that track kinematics consistent with the applied external
    loads.
    
"""

# =========================================================================
# Import packages
# =========================================================================

import opensim as osim
import os
import numpy as np
import pandas as pd
import shutil
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pickle
import argparse

# =========================================================================
# Flags for running analyses
# =========================================================================

# # Set participant Id to run
# participant = '2014001'
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--participant', action = 'store', type = str, help = 'Enter the participant ID number')
args = parser.parse_args()
participant = args.participant

# =========================================================================
# Set-up
# =========================================================================

# Set matplotlib parameters
# rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['font.weight'] = 'bold'
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 16
rcParams['axes.linewidth'] = 1.5
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False
rcParams['legend.fontsize'] = 10
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.width'] = 1.5
rcParams['legend.framealpha'] = 0.0
rcParams['savefig.dpi'] = 300
rcParams['savefig.format'] = 'pdf'

# Add the utility geometry path for model visualisation
osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.join(os.getcwd(), '..', 'model', 'Geometry'))

# Set the participant list based on those in dataset folder
participantList = [ii for ii in os.listdir(os.path.join('..', 'data',)) if
                   os.path.isdir(os.path.join('..', 'data' ,ii))]

# Check if input participant is in list
if participant not in participantList:
    raise ValueError(f'No data found for participant ID {participant}. Check input for error...')

# Read in participant anthropometrics
anthropometrics = pd.read_csv(os.path.join('..', 'data', 'anthropometrics.csv'))

# Set the list of speed conditions to process
# Modify this if you want to extract different walking speeds
conditionList = [
    # 'C1',   # 0-0.4 m/s
    # 'C2',   # 0.4-0.8 m/s
    # 'C3',   # 0.8-1.2 m/s
    'C4',   # self-selected spontaneuous speed
    'C5',   # self-selected fast speed
    ]

# Create dictionary that sets optimal forces for torque model actuators
actForces = {'pelvis_tx': {'actuatorType': 'residual', 'optForce': 1},
             'pelvis_ty': {'actuatorType': 'residual', 'optForce': 1},
             'pelvis_tz': {'actuatorType': 'residual', 'optForce': 1},
             'pelvis_tilt': {'actuatorType': 'residual', 'optForce': 1},
             'pelvis_list': {'actuatorType': 'residual', 'optForce': 1},
             'pelvis_rotation': {'actuatorType': 'residual', 'optForce': 1},
             'hip_flexion_r': {'actuatorType': 'actuator', 'optForce': 300},
             'hip_adduction_r': {'actuatorType': 'actuator', 'optForce': 200},
             'hip_rotation_r': {'actuatorType': 'actuator', 'optForce': 100},
             'knee_angle_r': {'actuatorType': 'actuator', 'optForce': 300},
             'ankle_angle_r': {'actuatorType': 'actuator', 'optForce': 300},
             'hip_flexion_l': {'actuatorType': 'actuator', 'optForce': 300},
             'hip_adduction_l': {'actuatorType': 'actuator', 'optForce': 200},
             'hip_rotation_l': {'actuatorType': 'actuator', 'optForce': 100},
             'knee_angle_l': {'actuatorType': 'actuator', 'optForce': 300},
             'ankle_angle_l': {'actuatorType': 'actuator', 'optForce': 300},
             'lumbar_extension': {'actuatorType': 'actuator', 'optForce': 300},
             'lumbar_bending': {'actuatorType': 'actuator', 'optForce': 200},
             'lumbar_rotation': {'actuatorType': 'actuator', 'optForce': 100},
             'arm_flex_r': {'actuatorType': 'actuator', 'optForce': 300},
             'arm_add_r': {'actuatorType': 'actuator', 'optForce': 200},
             'arm_rot_r': {'actuatorType': 'actuator', 'optForce': 100},
             'elbow_flex_r': {'actuatorType': 'actuator', 'optForce': 300},
             'pro_sup_r': {'actuatorType': 'actuator', 'optForce': 100},
             'arm_flex_l': {'actuatorType': 'actuator', 'optForce': 300},
             'arm_add_l': {'actuatorType': 'actuator', 'optForce': 200},
             'arm_rot_l': {'actuatorType': 'actuator', 'optForce': 100},
             'elbow_flex_l': {'actuatorType': 'actuator', 'optForce': 300},
             'pro_sup_l': {'actuatorType': 'actuator', 'optForce': 100},
             }

# Set the list of markers designating floor level on the models feet
floorMarkers = ['R_FM1_ground', 'R_FM2_ground', 'R_FM5_ground',
                'R_FM1_mid_ground', 'R_FM2_mid_ground', 'R_FM5_mid_ground', 'R_FCC_ground',
                'L_FM1_ground', 'L_FM2_ground', 'L_FM5_ground',
                'L_FM1_mid_ground', 'L_FM2_mid_ground', 'L_FM5_mid_ground', 'L_FCC_ground']

# Set weights for optimisations
globalStatesTrackingWeight = 1e0  # TODO: could make this 10 given different problem make-up?
globalControlEffortGoal = 1e-3

# Set mesh interval
meshInterval = 25

# =========================================================================
# Run tracking simulation
# =========================================================================

# Print out to console as a bookmark in any log file
print(f'***** ----- RUNNING TRACKING SIM FOR PARTICIPANT {participant} ----- *****')

# Make directory for participant
os.makedirs(os.path.join('..', 'data', participant, 'tracking'), exist_ok=True)

# Loop through the conditions
for condition in conditionList:

    # Print out to console as a bookmark in any log file
    print(f'***** ----- RUNNING TRACKING SIM FOR CONDITION {condition} ----- *****')

    # Set the trial folder
    # Note assumes only one/first trial for each condition to be run
    trialFolder = [ii for ii in os.listdir(os.path.join('..', 'data', participant, 'ik')) if
                   os.path.isdir(os.path.join('..', 'data', participant, 'ik', ii)) and
                   ii.startswith(condition)][0]

    # Create a folder for the condition and trial
    os.makedirs(os.path.join('..', 'data', participant, 'tracking', trialFolder), exist_ok=True)

    # Load in event data for trial
    with open(os.path.join('..', 'data', participant, 'events',
                           f'{participant}_{trialFolder}.pkl'), 'rb') as pklFile:
        eventData = pickle.load(pklFile)

    # Set-up model and files
    # -------------------------------------------------------------------------

    # Navigate to current folder to avoid issues with file paths in external loads an elsewhere
    homeDir = os.getcwd()
    os.chdir(os.path.join('..', 'data', participant, 'tracking', trialFolder))

    # Copy external loads files to directory
    shutil.copyfile(os.path.join('..', '..', 'dynamic', f'{participant}_{trialFolder}_grf.xml'),
                    f'{participant}_{trialFolder}_grf.xml')
    shutil.copyfile(os.path.join('..', '..', 'dynamic', f'{participant}_{trialFolder}_grf.mot'),
                    f'{participant}_{trialFolder}_grf.mot')

    # Construct a model processor to use with the tool
    modelProc = osim.ModelProcessor(os.path.join('..', '..', 'scaling', f'{participant}_scaledModelAdjusted.osim'))

    # Append external loads
    modelProc.append(osim.ModOpAddExternalLoads(f'{participant}_{trialFolder}_grf.xml'))

    # Weld desired locked joints
    # Create vector string object
    weldVectorStr = osim.StdVectorString()
    [weldVectorStr.append(joint) for joint in [
        'radius_hand_r', 'radius_hand_l', 'subtalar_r', 'subtalar_l', 'mtp_r', 'mtp_l']]
    # Append to model processor
    modelProc.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))

    # Process the model to edit
    trackingModel = modelProc.process()

    # Add coordinate actuators to model
    for coordinate in actForces:
        # Create actuator
        actu = osim.CoordinateActuator()
        # Set name
        actu.setName(f'{coordinate}_{actForces[coordinate]["actuatorType"]}')
        # Set coordinate
        actu.setCoordinate(trackingModel.updCoordinateSet().get(coordinate))
        # Set optimal force
        actu.setOptimalForce(actForces[coordinate]['optForce'])
        # Set min and max control
        actu.setMinControl(np.inf * -1)
        actu.setMaxControl(np.inf * 1)
        # Append to model force set
        trackingModel.updForceSet().cloneAndAppend(actu)

    # Finalise model connections
    trackingModel.finalizeConnections()

    # Print model to file in tracking directory
    trackingModel.printToXML(f'{participant}_{trialFolder}_trackingModel.osim')
    trackingModel.initSystem()

    # Set-up the Moco tracking problem
    # -------------------------------------------------------------------------

    # Create tracking tool
    track = osim.MocoTrack()
    track.setName(f'{participant}_{trialFolder}')

    # Set model
    trackModelProc = osim.ModelProcessor(f'{participant}_{trialFolder}_trackingModel.osim')
    track.setModel(trackModelProc)

    # Set coordinates from IK in tracking tool

    # Read in IK data
    ikTable = osim.TimeSeriesTable(os.path.join('..', '..', 'ik', trialFolder, f'{trialFolder}_ik.mot'))

    # Convert IK columns to state names
    stateColLabels = osim.StdVectorString()
    for currCol in ikTable.getColumnLabels():
        try:
            stateName = trackingModel.getCoordinateSet().get(currCol).getAbsolutePathString() + '/value'
            stateColLabels.append(stateName)
        except:
            ikTable.removeColumn(currCol)
    ikTable.setColumnLabels(stateColLabels)

    # Convert IK data to radians
    ikTracking = osim.TableProcessor(ikTable).processAndConvertToRadians(trackingModel)

    # Print to file for reference
    osim.STOFileAdapter.write(ikTracking, 'ikTrackingData.sto')

    # Set kinematics and markers to track in tool
    track.setStatesReference(osim.TableProcessor(ikTracking))

    # Set general parameters
    track.set_allow_unused_references(True)
    track.set_track_reference_position_derivatives(True)
    track.set_apply_tracked_states_to_guess(True)

    # Set global tracking weights
    track.set_states_global_tracking_weight(globalStatesTrackingWeight)

    # Set times
    track.set_initial_time(ikTracking.getIndependentColumn()[0])
    track.set_final_time(ikTracking.getIndependentColumn()[-1])

    # Initialise to a Moco study and problem to finalise
    # -------------------------------------------------------------------------

    # Get study and problem
    study = track.initialize()
    problem = study.updProblem()

    # Update control effort goal
    # -------------------------------------------------------------------------

    # Get a reference to the MocoControlCost goal and set parameters
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(globalControlEffortGoal)
    effort.setExponent(2)

    # Update individual weights in control effort goal
    # Put high weight on residual use
    effort.setWeightForControlPattern('/forceset/.*_residual', 10.0)
    # Put standard weight on the torque actuators
    effort.setWeightForControlPattern('/forceset/.*_actuator', 1.0)

    # Update the states tracking goal
    # -------------------------------------------------------------------------

    # Get a reference to the MocoState tracking goal
    tracking = osim.MocoStateTrackingGoal.safeDownCast(problem.updGoal('state_tracking'))

    # Scale state weights based on coordinate values
    tracking.setScaleWeightsWithRange(True)

    # Add constraints to the problem
    # -------------------------------------------------------------------------

    # Constrain contact markers to at minimum be at ground level
    # Loop through markers to create path constraints
    for marker in floorMarkers:

        # Create constraint
        markerConstraint = osim.MocoOutputConstraint()
        markerConstraint.setName(f'{marker}_constraint')

        # Set path to marker location
        markerConstraint.setOutputPath(f'/markerset/{marker}|location')

        # Set output index to y-axis
        markerConstraint.setOutputIndex(1)

        # Create and set the bounds to ground level and a reasonable height
        markerBounds = osim.StdVectorMocoBounds()
        markerBounds.append(osim.MocoBounds(0.0, 0.25))
        markerConstraint.updConstraintInfo().setBounds(markerBounds)

        # Add to problem
        problem.addPathConstraint(markerConstraint)

    # Set bounds in problem
    # -------------------------------------------------------------------------

    # Set bounds on joint coordinates using the tracking data as a reference
    # Initial and final bounds are set to be within a percentage of the range their tracking values
    # Max and minimum bounds are set to be within a percentage of the range for the coordinate
    perRange = 0.2

    # Loop through states
    for stateName in list(ikTracking.getColumnLabels()):

        # Check for joint coordinate state
        if stateName.endswith('/value'):

            # Calculate the range, initial and final values for coordinate
            valRange = np.ptp(ikTracking.getDependentColumn(stateName).to_numpy())
            initialVal = ikTracking.getDependentColumn(stateName).to_numpy()[0]
            finalVal = ikTracking.getDependentColumn(stateName).to_numpy()[-1]
            maxVal = ikTracking.getDependentColumn(stateName).to_numpy().max()
            minVal = ikTracking.getDependentColumn(stateName).to_numpy().min()

            # Set the bounds
            problem.setStateInfo(stateName,
                                 [minVal - (valRange * perRange), maxVal + (valRange * perRange)],
                                 [initialVal - (valRange * perRange), initialVal + (valRange * perRange)],
                                 [finalVal - (valRange * perRange), finalVal + (valRange * perRange)])

    # Define and configure the solver
    # -------------------------------------------------------------------------
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())

    # Set solver options
    solver.set_optim_max_iterations(2000)
    solver.set_num_mesh_intervals(meshInterval)
    solver.set_optim_constraint_tolerance(1e-3)
    solver.set_optim_convergence_tolerance(1e-3)
    solver.resetProblem(problem)

    # Solve the problem
    # -------------------------------------------------------------------------
    trackingSolution = study.solve()

    # # Option to visualise solution
    # study.visualize(trackingSolution)

    # Save files and finalize
    # -------------------------------------------------------------------------

    # Write solution to file
    if trackingSolution.isSealed():
        trackingSolution.unseal()
    trackingSolution.write(f'{participant}_{trialFolder}_trackingSolution.sto')

    # Remove initial tracked states and markers file
    os.remove(f'{participant}_{trialFolder}_tracked_states.sto')

    # Return to home directory
    os.chdir(homeDir)

    # Print out to console as a bookmark in any log file
    print(f'***** ----- FINISHED TRACKING SIM FOR CONDITION {condition} ----- *****')

# Print out to console as a bookmark in any log file
print(f'***** ----- FINISHED TRACKING SIM FOR PARTICIPANT {participant} ----- *****')

# =========================================================================
# Finalise and exit
# =========================================================================

# Exit console to avoid exit code error
os._exit(00)

# %% ----- end of processData.py ----- %% #