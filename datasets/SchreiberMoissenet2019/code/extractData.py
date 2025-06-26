# -*- coding: utf-8 -*-
"""

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code is a preliminary step in extracting just the relevant data needed
    from the SchreiberMoissenet2019 dataset. There are a number of C3D files in
    each participants folder, and this code grabs the appropriate files for the
    desired speeds. Participant data is checked for trials that meet the conditions
    for analysis (see README for details). The c3d files are converted to OpenSim
    friendly formats, those being TRC and MOT file types. Each participant needs
    to have their data within a 'raw' folder and then in a folder named with the
    participant code (e.g. '2014001') for this script to work (see README in
    folder for details).

    This code also scales a musculoskeletal model and extracts the relevant event
    data for subsequent planned simulations with the dataset. The events are extracted
    as a toe-off to the subsequent toe-off while the participant is in contact with
    the force plates.

"""

# =========================================================================
# Import packages
# =========================================================================

import opensim as osim
import os
import glob
import shutil
import numpy as np
import btk
import pandas as pd
import random
import pickle

# =========================================================================
# Set-up
# =========================================================================

# Set the participant list based on the codes in the raw folder
# See README.MD in dataset folder for details on setting up raw data
participantList = [ii for ii in os.listdir(os.path.join('..', 'data', 'raw')) if os.path.isdir(os.path.join('..', 'data', 'raw', ii))]

# Read in participant anthropometrics
anthropometrics = pd.read_csv(os.path.join('..', 'data', 'anthropometrics.csv'))

# Set the list of speed conditions to extract
# Modify this if you want to extract different walking speeds
conditionList = [
    # 'C1',   # 0-0.4 m/s
    # 'C2',   # 0.4-0.8 m/s
    # 'C3',   # 0.8-1.2 m/s
    'C4',   # self-selected spontaneuous speed
    'C5',   # self-selected fast speed
    ]

# Create the two rotations needed for data
rot1 = osim.Rotation(np.deg2rad(-90), osim.Vec3(1,0,0))
rot2 = osim.Rotation(np.deg2rad(180), osim.Vec3(0,1,0))

# Read in trial data file
trials = pd.read_csv(os.path.join('..', 'data', 'trials.csv'))

# Add the utility geometry path for model visualisation
osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.join(os.getcwd(), '..', 'model', 'Geometry'))

# Create dictionaries for tools to avoid over-writing
scaleTool = {participantId: osim.ScaleTool() for participantId in participantList}
ikTool = {participantId: {condition: osim.InverseKinematicsTool() for condition in conditionList} for participantId in participantList}

# =========================================================================
# Create the relevant set-up files for OpenSim tools
# =========================================================================

# Create measurement set for scaling
# -------------------------------------------------------------------------

# Set the parameters for the measurement sets
measurementSetParams = {
    # Torso
    'torso_depth': {'markerPairs': [['CV7', 'SJN'], ['TV10', 'SXS'], ], 'bodyScale': ['torso'], 'axes': 'X'},
    'torso_height': {'markerPairs': [['CV7', 'TV10'], ['SJN', 'SXS'], ], 'bodyScale': ['torso'], 'axes': 'Y'},
    'torso_width': {'markerPairs': [['R_SAE', 'L_SAE'],], 'bodyScale': ['torso'], 'axes': 'Z'},
    # Pelvis
    'pelvis_depth': {'markerPairs': [['R_IAS', 'R_IPS'], ['L_IAS', 'L_IPS'], ], 'bodyScale': ['pelvis'], 'axes': 'X'},
    'pelvis_height': {'markerPairs': [['R_IAS', 'R_FTC'], ['L_IAS', 'L_FTC'], ], 'bodyScale': ['pelvis'], 'axes': 'Y'},
    'pelvis_width': {'markerPairs': [['R_IAS', 'L_IAS'], ['R_IPS', 'L_IPS'], ], 'bodyScale': ['pelvis'], 'axes': 'Z'},
    # Right thigh
    'r_thigh_length': {'markerPairs': [['R_FTC', 'R_FLE'], ], 'bodyScale': ['femur_r'], 'axes': 'Y'},
    'r_thigh_breadth': {'markerPairs': [['R_FLE', 'R_FME'], ], 'bodyScale': ['femur_r'], 'axes': 'XZ'},
    # Right patella
    'r_patella': {'markerPairs': [['R_FTC', 'R_FLE'], ], 'bodyScale': ['patella_r'], 'axes': 'XYZ'},
    # Right shank
    'r_tibia_length': {'markerPairs': [['R_TTC', 'R_TAM'], ], 'bodyScale': ['tibia_r', 'talus_r'], 'axes': 'Y'},
    'r_tibia_breadth': {'markerPairs': [['R_TAM', 'R_FAL'], ], 'bodyScale': ['tibia_r', 'talus_r'], 'axes': 'XZ'},
    # Right foot
    'r_foot_length': {'markerPairs': [['R_FCC', 'R_FM2'], ], 'bodyScale': ['calcn_r', 'toes_r'], 'axes': 'X'},
    'r_foot_width': {'markerPairs': [['R_FM1', 'R_FM5'], ], 'bodyScale': ['calcn_r', 'toes_r'], 'axes': 'YZ'},
    # Left thigh
    'l_thigh_length': {'markerPairs': [['L_FTC', 'L_FLE'], ], 'bodyScale': ['femur_l'], 'axes': 'Y'},
    'l_thigh_breadth': {'markerPairs': [['L_FLE', 'L_FME'], ], 'bodyScale': ['femur_l'], 'axes': 'XZ'},
    # Left patella
    'l_patella': {'markerPairs': [['L_FTC', 'L_FLE'], ], 'bodyScale': ['patella_l'], 'axes': 'XYZ'},
    # Left shank
    'l_tibia_length': {'markerPairs': [['L_TTC', 'L_TAM'], ], 'bodyScale': ['tibia_l', 'talus_l'], 'axes': 'Y'},
    'l_tibia_breadth': {'markerPairs': [['L_TAM', 'L_FAL'], ], 'bodyScale': ['tibia_l', 'talus_l'], 'axes': 'XZ'},
    # Left foot
    'l_foot_length': {'markerPairs': [['L_FCC', 'L_FM2'], ], 'bodyScale': ['calcn_l', 'toes_l'], 'axes': 'X'},
    'l_foot_width': {'markerPairs': [['L_FM1', 'L_FM5'], ], 'bodyScale': ['calcn_l', 'toes_l'], 'axes': 'YZ'},
    # Right upper arm
    'r_upperarm_length': {'markerPairs': [['R_SAA', 'R_UOA'], ], 'bodyScale': ['humerus_r'], 'axes': 'Y'},
    'r_upperarm_breadth': {'markerPairs': [['R_HLE', 'R_HME'], ], 'bodyScale': ['humerus_r'], 'axes': 'XZ'},
    # Right forearm
    'r_forearm_length': {'markerPairs': [['R_HLE', 'R_RSP'], ['R_HME', 'R_UHE'], ], 'bodyScale': ['radius_r', 'ulna_r'], 'axes': 'Y'},
    'r_forearm_breadth': {'markerPairs': [['R_RSP', 'R_UHE'], ], 'bodyScale': ['radius_r', 'ulna_r'], 'axes': 'XZ'},
    # Right hand
    'r_hand_length': {'markerPairs': [['R_RSP', 'R_HM2'], ['R_UHE', 'R_HM5'], ], 'bodyScale': ['hand_r'], 'axes': 'Y'},
    'r_hand_breadth': {'markerPairs': [['R_HM2', 'R_HM5'], ], 'bodyScale': ['hand_r'], 'axes': 'XZ'},
    # Left upper arm
    'l_upperarm_length': {'markerPairs': [['L_SAA', 'L_UOA'], ], 'bodyScale': ['humerus_l'], 'axes': 'Y'},
    'l_upperarm_breadth': {'markerPairs': [['L_HLE', 'L_HME'], ], 'bodyScale': ['humerus_l'], 'axes': 'XZ'},
    # Left forearm
    'l_forearm_length': {'markerPairs': [['L_HLE', 'L_RSP'], ['L_HME', 'L_UHE'], ], 'bodyScale': ['radius_l', 'ulna_l'], 'axes': 'Y'},
    'l_forearm_breadth': {'markerPairs': [['L_RSP', 'L_UHE'], ], 'bodyScale': ['radius_l', 'ulna_l'], 'axes': 'XZ'},
    # Left hand
    'l_hand_length': {'markerPairs': [['L_RSP', 'L_HM2'], ['L_UHE', 'L_HM5'], ], 'bodyScale': ['hand_l'], 'axes': 'Y'},
    'l_hand_breadth': {'markerPairs': [['L_HM2', 'L_HM5'], ], 'bodyScale': ['hand_l'], 'axes': 'XZ'},
}

# Create the measurement set
scaleMeasurementSet = osim.MeasurementSet()

# Append the measurements from parameters
for measureName in measurementSetParams.keys():
    # Create the measurement
    measurement = osim.Measurement()
    measurement.setName(measureName)
    # Append the marker pairs
    for ii in range(len(measurementSetParams[measureName]['markerPairs'])):
        measurement.getMarkerPairSet().cloneAndAppend(
            osim.MarkerPair(measurementSetParams[measureName]['markerPairs'][ii][0],
                            measurementSetParams[measureName]['markerPairs'][ii][1]))
    # Append the body scales
    for ii in range(len(measurementSetParams[measureName]['bodyScale'])):
        # Create body scale
        bodyScale = osim.BodyScale()
        bodyScale.setName(measurementSetParams[measureName]['bodyScale'][ii])
        # Create and set axis names
        axes = osim.ArrayStr()
        for jj in range(len(measurementSetParams[measureName]['axes'])):
            axes.append(measurementSetParams[measureName]['axes'][jj])
        bodyScale.setAxisNames(axes)
        # Apppend to body scale set
        measurement.getBodyScaleSet().cloneAndAppend(bodyScale)
    # Append the measurement to the set
    scaleMeasurementSet.cloneAndAppend(measurement)

# Create scale task set
# -------------------------------------------------------------------------

# Set the parameters for the scale marker and joint task sets
markerParams = {
    # Torso
    'SJN': {'weight': 5.0}, 'SXS': {'weight': 5.0}, 'CV7': {'weight': 5.0}, 'TV10': {'weight': 5.0},
    'R_SAE': {'weight': 0.0}, 'L_SAE': {'weight': 0.0},
    # Pelvis
    'L_IAS': {'weight': 10.0}, 'L_IPS': {'weight': 10.0}, 'R_IAS': {'weight': 5.0}, 'R_IPS': {'weight': 5.0},
    # Right thigh
    'R_FTC': {'weight': 2.5}, 'R_FLE': {'weight': 5.0}, 'R_FME': {'weight': 5.0},
    # Right shank
    'R_TTC': {'weight': 2.5}, 'R_FAX': {'weight': 2.5}, 'R_TAM': {'weight': 5.0}, 'R_FAL': {'weight': 5.0},
    # Right foot
    'R_FCC': {'weight': 5.0}, 'R_FM1': {'weight': 2.5}, 'R_FM2': {'weight': 2.5}, 'R_FM5': {'weight': 2.5},
    # Left thigh
    'L_FTC': {'weight': 2.5}, 'L_FLE': {'weight': 5.0}, 'L_FME': {'weight': 5.0},
    # Left shank
    'L_TTC': {'weight': 2.5}, 'L_FAX': {'weight': 2.5}, 'L_TAM': {'weight': 5.0}, 'L_FAL': {'weight': 5.0},
    # Left foot
    'L_FCC': {'weight': 5.0}, 'L_FM1': {'weight': 2.5}, 'L_FM2': {'weight': 2.5}, 'L_FM5': {'weight': 2.5},
    # Right upper arm
    'R_SAA': {'weight': 0.0}, 'R_HLE': {'weight': 1.0}, 'R_HME': {'weight': 1.0},
    # Right forearm
    'R_UOA': {'weight': 1.0}, 'R_RSP': {'weight': 1.0}, 'R_UHE': {'weight': 1.0},
    # Right hand
    'R_HM2': {'weight': 1.0}, 'R_HM5': {'weight': .0},
    # Left upper arm
    'L_SAA': {'weight': 0.0}, 'L_HLE': {'weight': 1.0}, 'L_HME': {'weight': 1.0},
    # Left forearm
    'L_UOA': {'weight': 1.0}, 'L_RSP': {'weight': 1.0}, 'L_UHE': {'weight': 1.0},
    # Left hand
    'L_HM2': {'weight': 1.0}, 'L_HM5': {'weight': 1.0},
}

# Create the task set
scaleTaskSet = osim.IKTaskSet()

# Append the tasks from the parameters
for taskName in markerParams.keys():
    # Create the task and add details
    task = osim.IKMarkerTask()
    task.setName(taskName)
    task.setWeight(markerParams[taskName]['weight'])
    if markerParams[taskName]['weight'] == 0.0:
        task.setApply(False)
    # Append to task set
    scaleTaskSet.cloneAndAppend(task)

# Create the IK task set for tracking
# -------------------------------------------------------------------------
""" NOTE: Currently just defaulting to all 1's... """

# Set the parameters for the IK task sets
ikTaskSetParams = {
    # Torso
    'SJN': {'weight': 1}, 'SXS': {'weight': 1}, 'CV7': {'weight': 1}, 'TV10': {'weight': 1}, 'R_SAE': {'weight': 0}, 'L_SAE': {'weight': 0},
    # Pelvis
    'L_IAS': {'weight': 1}, 'L_IPS': {'weight': 1}, 'R_IAS': {'weight': 1}, 'R_IPS': {'weight': 1},
    # Right thigh
    'R_FTC': {'weight': 1}, 'R_FLE': {'weight': 1}, 'R_FME': {'weight': 1},
    # Right shank
    'R_TTC': {'weight': 1}, 'R_FAX': {'weight': 1}, 'R_TAM': {'weight': 1}, 'R_FAL': {'weight': 1},
    # Right foot
    'R_FCC': {'weight': 1}, 'R_FM1': {'weight': 1}, 'R_FM2': {'weight': 1}, 'R_FM5': {'weight': 1},
    # Left thigh
    'L_FTC': {'weight': 1}, 'L_FLE': {'weight': 1}, 'L_FME': {'weight': 1},
    # Left shank
    'L_TTC': {'weight': 1}, 'L_FAX': {'weight': 1}, 'L_TAM': {'weight': 1}, 'L_FAL': {'weight': 1},
    # Left foot
    'L_FCC': {'weight': 1}, 'L_FM1': {'weight': 1}, 'L_FM2': {'weight': 1}, 'L_FM5': {'weight': 1},
    # Right upper arm
    'R_SAA': {'weight': 0}, 'R_HLE': {'weight': 1}, 'R_HME': {'weight': 1},
    # Right forearm
    'R_UOA': {'weight': 1}, 'R_RSP': {'weight': 1}, 'R_UHE': {'weight': 1},
    # Right hand
    'R_HM2': {'weight': 1}, 'R_HM5': {'weight': 1},
    # Left upper arm
    'L_SAA': {'weight': 0}, 'L_HLE': {'weight': 1}, 'L_HME': {'weight': 1},
    # Left forearm
    'L_UOA': {'weight': 1}, 'L_RSP': {'weight': 1}, 'L_UHE': {'weight': 1},
    # Left hand
    'L_HM2': {'weight': 1}, 'L_HM5': {'weight': 1},
}

# Create the task set
ikTaskSet = osim.IKTaskSet()

# Append the tasks from the parameters
for taskName in ikTaskSetParams.keys():
    # Create the task and add details
    task = osim.IKMarkerTask()
    task.setName(taskName)
    task.setWeight(ikTaskSetParams[taskName]['weight'])
    # Append to task set
    ikTaskSet.cloneAndAppend(task)

# =========================================================================
# Loop through participants to extract data
# =========================================================================

# Participant loop
for participant in participantList:

    # Determine if participant has the necessary good trials for the selected speeds
    # -------------------------------------------------------------------------

    # Create a dictionary to store boolean status for participant conditions
    # Defaults to False and will be changed if appropriate
    conditionStatus = {condition: False for condition in conditionList}
    conditionTrial = {condition: None for condition in conditionList}

    # An initial check is required to see if the force plate info is there for the participant
    if int(participant) in trials['subjectId'].to_list():

        # Check for trials from condition
        for condition in conditionList:

            # Get condition labels from dataframe
            conditionLabels = [label for label in trials.columns if label.startswith(condition + '_')]

            # Extract the trial FP contacts for current condition
            fpContacts = [trials.loc[(trials['subjectId'] == int(participant))][label].values[0] for label in conditionLabels]

            # Find where trial has the two FP contacts (i.e. 12 or 21)
            validTrials = [contact == 12 or contact == 21 for contact in fpContacts]

            # Check if condition is useable
            if np.sum(validTrials) > 0:

                # Convert condition status if appropriate
                conditionStatus[condition] = True

                # Make the random selection of the trial to use from the viable trials
                random.seed(int(participant))
                trialSelect = random.choice(list(np.where(validTrials)[0]))
                conditionTrial[condition] = conditionLabels[trialSelect]

        # Create variable whether to process participant based on all conditions being valid
        processParticipant = all(conditionStatus[condition] for condition in conditionList)

    # Don't process if participant isn't there
    else:

        processParticipant = False

    # Process participant if appropriate
    # -------------------------------------------------------------------------
    if processParticipant:

        # Create folders for participant
        # -------------------------------------------------------------------------

        # Starting directory
        os.makedirs(os.path.join('..', 'data', participant), exist_ok = True)

        # Static files directory
        os.makedirs(os.path.join('..', 'data', participant, 'static'), exist_ok = True)

        # Dynamic trials file directory
        os.makedirs(os.path.join('..', 'data', participant, 'dynamic'), exist_ok = True)

        # Event data directory
        os.makedirs(os.path.join('..', 'data', participant, 'events'), exist_ok=True)

        # Scaling directory
        os.makedirs(os.path.join('..', 'data', participant, 'scaling'), exist_ok=True)

        # Inverse kinematics directory
        os.makedirs(os.path.join('..', 'data', participant, 'ik'), exist_ok=True)

        # Get the static file
        # -------------------------------------------------------------------------

        # Get the path to the static file
        staticFile = glob.glob(os.path.join('..', 'data', 'raw', participant, '*_ST.c3d*'))[0]

        # Copy c3d across to static directory
        shutil.copyfile(staticFile, os.path.join('..', 'data', participant, 'static', os.path.split(staticFile)[-1]))

        # Construct opensim 3d object
        c3dFile = osim.C3DFileAdapter()
        c3dFile.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)

        # Read in the static trial
        staticC3D = c3dFile.read(staticFile)

        # Get markers table
        staticMarkers = c3dFile.getMarkersTable(staticC3D)

        # Rotate marker data
        for iRow in range(staticMarkers.getNumRows()):
            #Apply the two rotations
            staticMarkers.setRowAtIndex(iRow, rot1.multiply(staticMarkers.getRowAtIndex(iRow)))

        # There's a potentially annoying bug to deal with in scaling when only one row of marker data exists
        # Add a pseudo row of identical marker data to avoid this
        staticMarkers.appendRow(1, staticMarkers.getRow(0))

        # Write static markers to TRC file
        osim.TRCFileAdapter().write(staticMarkers, os.path.join('..', 'data', participant, 'static',
                                                                os.path.split(staticFile)[-1].split('.')[0]+'.trc'))

        # Scale musculoskeletal model for participant
        # -------------------------------------------------------------------------

        # Get static trial file
        staticTRC = glob.glob(os.path.join('..', 'data', participant, 'static', '*_ST.trc'))[0]

        # Set-up and run the scale tool
        # -------------------------------------------------------------------------

        # Set participant mass
        scaleTool[participant].setSubjectMass(
            anthropometrics[anthropometrics['subjectID'] == int(participant)]['mass'].values[0]
        )

        # Set generic model file
        scaleTool[participant].getGenericModelMaker().setModelFileName(os.path.join('..', 'model',
                                                                                    'Rajagopal2015_SchreiberMoissenet2019.osim'))

        # Set measurement set in model scaler
        scaleTool[participant].getModelScaler().setMeasurementSet(scaleMeasurementSet)

        # Set scale tasks in tool
        for ii in range(scaleTaskSet.getSize()):
            scaleTool[participant].getMarkerPlacer().getIKTaskSet().cloneAndAppend(scaleTaskSet.get(ii))

        # Set marker file
        scaleTool[participant].getMarkerPlacer().setMarkerFileName(staticTRC)
        scaleTool[participant].getModelScaler().setMarkerFileName(staticTRC)

        # Set options
        scaleTool[participant].getModelScaler().setPreserveMassDist(True)
        scaleOrder = osim.ArrayStr()
        scaleOrder.set(0, 'measurements')
        scaleTool[participant].getModelScaler().setScalingOrder(scaleOrder)

        # Set time ranges
        timeRange = osim.ArrayDouble()
        timeRange.set(0, 0)  # initial time
        timeRange.set(1, 1)  # final time
        scaleTool[participant].getMarkerPlacer().setTimeRange(timeRange)
        scaleTool[participant].getModelScaler().setTimeRange(timeRange)

        # Set output files
        scaleTool[participant].getModelScaler().setOutputModelFileName(
            os.path.join('..', 'data', participant, 'scaling', f'{participant}_scaledModel.osim'))
        scaleTool[participant].getModelScaler().setOutputScaleFileName(
            os.path.join('..', 'data', participant, 'scaling', f'{participant}_scaleSet.xml'))

        # Set marker adjustment parameters
        scaleTool[participant].getMarkerPlacer().setOutputMotionFileName(
            os.path.join('..', 'data', participant, 'scaling', f'{participant}_staticMotion.mot'))
        scaleTool[participant].getMarkerPlacer().setOutputModelFileName(
            os.path.join('..', 'data', participant, 'scaling', f'{participant}_scaledModelAdjusted.osim'))

        # Save and run scale tool
        scaleTool[participant].printToXML(
            os.path.join('..', 'data', participant, 'scaling', f'{participant}_scaleSetup.xml'))
        scaleTool[participant].run()

        # Adjust the model
        # -------------------------------------------------------------------------

        # Add marker locations underneath foot markers at floor level based on static motion
        # These may be useful later in determining foot-ground contact points

        # Set the list of markers to project to the floot
        floorMarkers = ['R_FM1', 'R_FM2', 'R_FM5', 'R_FCC',
                        'L_FM1', 'L_FM2', 'L_FM5', 'L_FCC']

        # Load the model
        scaledModel = osim.Model(
            os.path.join('..', 'data', participant, 'scaling', f'{participant}_scaledModelAdjusted.osim'))

        # Set model name
        scaledModel.setName(participant)

        # Load the static motion
        staticMotion = osim.TimeSeriesTable(
            os.path.join('..', 'data', participant, 'scaling', f'{participant}_staticMotion.mot'))

        # Initialise the model state
        modelState = scaledModel.initSystem()

        # Set model states from static motion
        for stateName in staticMotion.getColumnLabels():
            scaledModel.setStateVariableValue(modelState, stateName, staticMotion.getDependentColumn(stateName).to_numpy()[0])

        # Realize model to position stage
        scaledModel.realizePosition(modelState)

        # Loop through floor markers
        # Identify position in ground at floor level and translate back to associated body
        for marker in floorMarkers:

            # Get marker location in ground
            markerInGround = scaledModel.updMarkerSet().get(marker).getLocationInGround(modelState)

            # Set the y-axis to be zero for the marker in ground position (i.e. floor level)
            markerInGround.set(1, 0)

            # Translate this new position back to the markers associated body
            floorMarkerLoc = scaledModel.getGround().findStationLocationInAnotherFrame(
                modelState, markerInGround, scaledModel.updMarkerSet().get(marker).getParentFrame()
            )

            # Create new marker and append to model
            newMarker = osim.Marker()
            newMarker.setName(marker + '_ground')
            newMarker.setParentFrameName(scaledModel.updMarkerSet().get(marker).getParentFrameName())
            newMarker.set_location(floorMarkerLoc)
            scaledModel.getMarkerSet().cloneAndAppend(newMarker)

        # Create mid-foot markers between FM and heel points
        for marker in floorMarkers:

            # Add additional markers at the mid-foot with the FM markers
            # These sit at the mid distance to the heel marker along the X, Y and Z axes
            if '_FM' in marker:

                # Get marker location in ground
                markerInGround = scaledModel.updMarkerSet().get(marker).getLocationInGround(modelState)

                # Set the y-axis to be zero for the marker in ground position (i.e. floor level)
                markerInGround.set(1, 0)

                # Get the appropriate heel marker position
                heelMarkerInGround = scaledModel.updMarkerSet().get(marker[0] + '_FCC').getLocationInGround(modelState)

                # Find mid-point of the points in the ground along the x-axis and y-axis
                midPointX = markerInGround.get(0) - ((markerInGround.get(0) - heelMarkerInGround.get(0)) / 2)
                midPointY = markerInGround.get(1) - ((markerInGround.get(1) - heelMarkerInGround.get(1)) / 2)
                midPointZ = markerInGround.get(2) - ((markerInGround.get(2) - heelMarkerInGround.get(2)) / 2)

                # Get the mid-point marker location in the body frame
                midMarkerLoc = scaledModel.getGround().findStationLocationInAnotherFrame(
                    modelState, osim.Vec3(midPointX, midPointY, midPointZ),
                    scaledModel.updMarkerSet().get(marker).getParentFrame())

                # Create new marker and append to model
                newMarker = osim.Marker()
                newMarker.setName(marker + '_mid_ground')
                newMarker.setParentFrameName(scaledModel.updMarkerSet().get(marker).getParentFrameName())
                newMarker.set_location(midMarkerLoc)
                scaledModel.getMarkerSet().cloneAndAppend(newMarker)

        # Finalise model connections
        scaledModel.finalizeConnections()

        # Print to file (overwrites original adjusted model)
        scaledModel.printToXML(os.path.join('..', 'data', participant, 'scaling', f'{participant}_scaledModelAdjusted.osim'))

        # Get the dynamic file
        # -------------------------------------------------------------------------

        # Identify dynamic files to extract from those selected
        extractFiles = [os.path.join('..', 'data', 'raw', participant,
                                     f'{participant}_{conditionTrial[condition]}.c3d') for condition in conditionList]

        # Loop through files to copy and convert
        for dynaFile in extractFiles:

            # Copy c3d across to dynamic directory
            shutil.copyfile(dynaFile, os.path.join('..', 'data', participant, 'dynamic', os.path.split(dynaFile)[-1]))

            # Construct opensim 3d object
            c3dFile = osim.C3DFileAdapter()
            c3dFile.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)

            # Read in the dynamic trial
            dynaC3D = c3dFile.read(dynaFile)

            # Get markers table
            dynaMarkers = c3dFile.getMarkersTable(dynaC3D)

            # Rotate marker data
            for iRow in range(dynaMarkers.getNumRows()):
                # Apply the two rotations
                dynaMarkers.setRowAtIndex(iRow, rot1.multiply(dynaMarkers.getRowAtIndex(iRow)))

            # Check if moving in negative x-direction and data needs to be flipped for consistency
            # Check by examining average difference in pelvis marker (negativ difference needs to be flipped)
            if np.diff(dynaMarkers.flatten().getDependentColumn('L_IAS_1').to_numpy()).mean() < 0:
                flipDirection = True
            else:
                flipDirection = False

            # Rotate 180 degrees abot y-axis to flip direction if necessary
            if flipDirection:
                # Rotate marker data
                for iRow in range(dynaMarkers.getNumRows()):
                    # Apply the two rotations
                    dynaMarkers.setRowAtIndex(iRow, rot2.multiply(dynaMarkers.getRowAtIndex(iRow)))

            # Write dynamic markers to TRC file
            osim.TRCFileAdapter().write(dynaMarkers, os.path.join('..', 'data', participant, 'dynamic',
                                                                  os.path.split(dynaFile)[-1].split('.')[0]+'.trc'))

            # Get forces table
            dynaForces = c3dFile.getForcesTable(dynaC3D)

            # Rotate forces data
            for iRow in range(dynaForces.getNumRows()):
                dynaForces.setRowAtIndex(iRow, rot1.multiply(dynaForces.getRowAtIndex(iRow)))

            # Rotate 180 degrees abot y-axis to flip direction if necessary
            if flipDirection:
                # Rotate forces data
                for iRow in range(dynaForces.getNumRows()):
                    dynaForces.setRowAtIndex(iRow, rot2.multiply(dynaForces.getRowAtIndex(iRow)))

            # Flatten forces data
            forcesFlat = dynaForces.flatten()

            # Convert to numpy array
            # Pre-allocate numpy array based on data size
            dataArray = np.zeros((forcesFlat.getNumRows(),
                                  forcesFlat.getNumColumns()))

            # Extract data
            for forceInd in range(forcesFlat.getNumColumns()):
                dataArray[:,forceInd] = forcesFlat.getDependentColumn(forcesFlat.getColumnLabels()[forceInd]).to_numpy()

            # Replace nan's for COP and moment data with zeros
            np.nan_to_num(dataArray, copy = False, nan = 0.0)

            # Convert force point data from mm to m
            for forceName in list(forcesFlat.getColumnLabels()):
                if forceName.startswith('p') or forceName.startswith('m'):
                    #Get force index
                    forceInd = list(forcesFlat.getColumnLabels()).index(forceName)
                    #Convert to m units in data array
                    dataArray[:,forceInd] = dataArray[:,forceInd] / 1000

            # Build the new time series table
            forcesStorage = osim.Storage()

            # Get the time data
            time = forcesFlat.getIndependentColumn()

            # Create maps to replace text from force labels with
            # Force plate and type identifiers
            forceType = {}
            for ii in range(1,int(len(forcesFlat.getColumnLabels()) / 9)+1):
                forceType[f'f{ii}'] = f'ground_force_{ii}_v'
                forceType[f'p{ii}'] = f'ground_force_{ii}_p'
                forceType[f'm{ii}'] = f'ground_force_{ii}_m'
            # Axis identifiers
            forceAxis = {'1': 'x',
                         '2': 'y',
                         '3': 'z'}

            # Set labels in table
            newLabels = osim.ArrayStr()
            newLabels.append('time')
            for forceLabel in forcesFlat.getColumnLabels():
                # Split the label to get parts
                labelSplit = forceLabel.split('_')
                # Create new label
                forceLabel = f'{forceType[labelSplit[0]]}{forceAxis[labelSplit[1]]}'
                # Append to labels vector
                newLabels.append(forceLabel)
            forcesStorage.setColumnLabels(newLabels)

            # Add data
            for iRow in range(dataArray.shape[0]):
                row = osim.ArrayDouble()
                for iCol in range(dataArray.shape[1]):
                    row.append(dataArray[iRow,iCol])
                # Add data to storage
                forcesStorage.append(time[iRow], row)

            # Set name for storage object
            forcesStorage.setName(os.path.split(dynaFile)[-1].split('.')[0]+'_grf')

            # Write to file
            forcesStorage.printResult(forcesStorage,
                                      os.path.split(dynaFile)[-1].split('.')[0]+'_grf',
                                      os.path.join('..', 'data', participant, 'dynamic'),
                                      0.001, '.mot')

            # Create external loads files
            # Note that external loads are applied to feet based on timing of force plate contact and event name

            # Create the external loads
            forceXML = osim.ExternalLoads()

            # Convert forces to time-series table for easier use
            forcesTable = osim.TimeSeriesTable(os.path.join('..', 'data', participant,
                                                            'dynamic', os.path.split(dynaFile)[-1].split('.')[0]+'_grf.mot'))

            # Extract the vertical force data from the two plates
            vForce1 = forcesTable.getDependentColumn('ground_force_1_vy').to_numpy()
            vForce2 = forcesTable.getDependentColumn('ground_force_2_vy').to_numpy()

            # Identify contact indices and times on plates based on force threshold
            forceThreshold = 20

            # Get the c3d data with events using btk
            reader = btk.btkAcquisitionFileReader()
            reader.SetFilename(dynaFile)
            reader.Update()
            c3dData = reader.GetOutput()

            # Get events
            eventData = {'eventName': [], 'eventFrame': [], 'eventTime': []}
            moreEvents = True
            eventInd = 0
            while moreEvents:
                try:
                    event = c3dData.GetEvent(eventInd)
                    eventData['eventName'].append(event.GetContext() + ' ' + event.GetLabel())
                    eventData['eventFrame'].append(event.GetFrame())
                    eventData['eventTime'].append(event.GetTime())
                    eventInd += 1
                except:
                    moreEvents = False
            events = pd.DataFrame.from_dict(eventData)
            events.sort_values(by='eventFrame', inplace=True)
            events.reset_index(drop=True, inplace=True)

            # Identify limb contact for force plate 1
            v1_ind = np.argmax(vForce1 > forceThreshold)
            v1_t = forcesTable.getIndependentColumn()[v1_ind]

            # Get the strike foot based on time difference between foot strike and event
            strikeFoot = events['eventName'][np.argmin(np.array(np.abs(events['eventTime'] - v1_t)))].split(' ')[0]

            # Create external loads for plate 1
            grf1 = osim.ExternalForce()
            grf1.setName('grf1')
            if strikeFoot == 'Left':
                grf1.setAppliedToBodyName('calcn_l')
            elif strikeFoot == 'Right':
                grf1.setAppliedToBodyName('calcn_r')
            else:
                raise ValueError('Left or right foot not identified as strike foot...')
            grf1.setForceExpressedInBodyName('ground')
            grf1.setPointExpressedInBodyName('ground')
            grf1.setForceIdentifier('ground_force_1_v')
            grf1.setPointIdentifier('ground_force_1_p')
            grf1.setTorqueIdentifier('ground_force_1_m')
            forceXML.cloneAndAppend(grf1)

            # Identify limb contact for force plate 2
            v2_ind = np.argmax(vForce2 > forceThreshold)
            v2_t = forcesTable.getIndependentColumn()[v2_ind]

            # Get the strike foot based on time difference between foot strike and event
            strikeFoot = events['eventName'][np.argmin(np.array(np.abs(events['eventTime'] - v2_t)))].split(' ')[0]

            # Create external loads for plate 2
            grf2 = osim.ExternalForce()
            grf2.setName('grf2')
            if strikeFoot == 'Left':
                grf2.setAppliedToBodyName('calcn_l')
            elif strikeFoot == 'Right':
                grf2.setAppliedToBodyName('calcn_r')
            else:
                raise ValueError('Left or right foot not identified as strike foot...')
            grf2.setForceExpressedInBodyName('ground')
            grf2.setPointExpressedInBodyName('ground')
            grf2.setForceIdentifier('ground_force_2_v')
            grf2.setPointIdentifier('ground_force_2_p')
            grf2.setTorqueIdentifier('ground_force_2_m')
            forceXML.cloneAndAppend(grf2)

            # Set GRF datafile in external loads
            forceXML.setDataFileName(os.path.split(dynaFile)[-1].split('.')[0]+'_grf.mot')

            # Write to file
            forceXML.printToXML(os.path.join('..', 'data', participant, 'dynamic',
                                             os.path.split(dynaFile)[-1].split('.')[0]+'_grf.xml'))

            # Extract trial event timings
            # -------------------------------------------------------------------------

            # Get first and second foot strike timings
            firstFootStrike = np.min((v1_t, v2_t))
            secondFootStrike = np.max((v1_t, v2_t))

            # Find the toe off after the first foot strike as the start time for simulations
            # Find the subsequent toe off to indicate the corresponding end time for simulations
            offEvents = [events.iloc[ii]['eventTime'] for ii in range(len(events)) if
                         'Off' in events.iloc[ii]['eventName'] and events.iloc[ii]['eventTime'] > firstFootStrike]
            offEvents.sort()
            startTime = offEvents[0]
            endTime = offEvents[1]

            # Store timings in a dictionary and save
            eventData = {'startTime': startTime, 'endTime': endTime}
            with open(os.path.join('..', 'data', participant, 'events',
                                   os.path.split(dynaFile)[-1].split('.')[0]+'.pkl'), 'wb') as pklFile:
                pickle.dump(eventData, pklFile, protocol=pickle.HIGHEST_PROTOCOL)

            # Run IK on trial
            # -------------------------------------------------------------------------

            # Get trial name and condition
            condition = os.path.split(dynaFile)[-1].split('.')[0].split('_')[1]
            trialName = os.path.split(dynaFile)[-1].split('.')[0].replace(f'{participant}_','')

            # Create a folder for the condition and trial
            os.makedirs(os.path.join('..', 'data', participant, 'ik', trialName), exist_ok=True)

            # Set model
            ikTool[participant][condition].set_model_file(os.path.join('..', 'data', participant, 'scaling',
                                                                       f'{participant}_scaledModelAdjusted.osim'))

            # Set task set (consistent for all trials)
            for taskInd in range(ikTaskSet.getSize()):
                ikTool[participant][condition].getIKTaskSet().adoptAndAppend(ikTaskSet.get(taskInd))

            # Set to report marker locations
            ikTool[participant][condition].set_report_marker_locations(True)

            # Set the marker file (relative to setup file location)
            ikTool[participant][condition].setMarkerDataFileName(os.path.join('..', 'data', participant, 'dynamic',
                                                                              os.path.split(dynaFile)[-1].split('.')[0]+'.trc'))

            # Set times
            ikTool[participant][condition].setStartTime(eventData['startTime'])
            ikTool[participant][condition].setEndTime(eventData['endTime'])

            # Set output filename (relative to setup file location)
            ikTool[participant][condition].setOutputMotionFileName(os.path.join('..', 'data', participant, 'ik', trialName,
                                                                                f'{trialName}_ik.mot'))

            # Save IK tool to file
            ikTool[participant][condition].printToXML('ikSetup.xml')

            # Bring the tool back in and run it (this seems to avoid Python kernel crashing)
            ikRun = osim.InverseKinematicsTool('ikSetup.xml')
            ikRun.run()

            # Rename supplementary marker outputs
            shutil.move('_ik_marker_errors.sto',
                        os.path.join('..', 'data', participant, 'ik', trialName,
                                     f'{trialName}_ikMarkerErrors.sto'))
            shutil.move('_ik_model_marker_locations.sto',
                        os.path.join('..', 'data', participant, 'ik', trialName,
                                     f'{trialName}_ikModelMarkerLocations.sto'))
            shutil.move('ikSetup.xml',
                        os.path.join('..', 'data', participant, 'ik', trialName,
                                     f'{trialName}_ikSetup.xml'))

        # Print confirmatory output
        # -------------------------------------------------------------------------
        print(f'Finished processing data for participant {participant}...')

    else:

        # Print out confirmation that participant was skipped
        # -------------------------------------------------------------------------
        print(f'No data processed for participant {participant} due to invalid trials...')

# =========================================================================
# Finalise and exit
# =========================================================================

# Exit console to avoid exit code error
os._exit(00)

# %% ----- end of extractData.py ----- %% #