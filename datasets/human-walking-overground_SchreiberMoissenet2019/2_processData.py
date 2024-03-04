# -*- coding: utf-8 -*-
'''

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code processes the data from the SchreiberMoissenet2019 dataset through
    an OpenSim Moco marker tracking framework. A torque-driven model is used to
    produce kinematics that track markers consistent with the applied external
    loads.
    
'''

# %% Import packages

import opensim as osim
import os
import glob
import numpy as np
import pandas as pd
import random
import shutil
import btk

# %% Flags for running analyses

#Change these boolean flags to choose which analyses to run
runScaling = False #doesn't take long
runIK = False #takes a little time
runTracking = True #takes a long time

# %% Set-up
        
#Get home path
homeDir = os.getcwd()

#Try adding OpenSim geometry path
#This can avoid issues with body geometries constantly not showing up when loading models
#Note that this is defaulting to the standard C: drive path for Windows and Opensim 4.4
if os.path.isdir(os.path.join('C:', os.sep, 'OpenSim 4.4', 'Geometry')):
    osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.join('C:', os.sep, 'OpenSim 4.4', 'Geometry'))
else:
    print(f'OpenSim geometry folder not located at {os.path.join("C:", os.sep, "OpenSim 4.4", "Geometry")}')
    print('This code can be modified or you can just deal with the consequences of not adding a geometry path...')

#Set the participant list based on the codes in the raw folder
participantList = [ii for ii in os.listdir() if os.path.isdir(ii)]

#Read in participant anthropometrics
anthropometrics = pd.read_csv('anthropometrics.csv')

#Set the list of speed conditions to extract 
#Modify this if you want to extract different walking speeds
conditionList  = [
    'C1', #0-0.4 m/s
    # 'C2', #0.4-0.8 m/s
    'C3', #0.8-1.2 m/s
    # 'C4', #self-selected spontaneuous speed
    # 'C5', #self-selected fast speed
    ]

#Create dictionaries for tools to avoid over-writing
scaleTool = {}
ikTool = {}

#Create dictionary that sets optimal forces for torque model actuators
actForces = {'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz':1,
             'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1,
             'hip_flexion_r': 300, 'hip_adduction_r': 200, 'hip_rotation_r': 100,
             'knee_angle_r': 300, 'ankle_angle_r': 300,
             'hip_flexion_l': 300, 'hip_adduction_l': 200, 'hip_rotation_l': 100,
             'knee_angle_l': 300, 'ankle_angle_l': 300,
             'lumbar_extension': 300, 'lumbar_bending': 200, 'lumbar_rotation': 100,
             'arm_flex_r': 300, 'arm_add_r': 200, 'arm_rot_r': 100,
             'elbow_flex_r': 300, 'pro_sup_r': 100,
             'wrist_flex_r': 100, 'wrist_dev_r': 100,
             'arm_flex_l': 300, 'arm_add_l': 200, 'arm_rot_l': 100,
             'elbow_flex_l': 300, 'pro_sup_l': 100,
             'wrist_flex_l': 100, 'wrist_dev_l': 100,
             }

#Set weights for optimisations
globalMarkerTrackingWeight = 10

#Set vertical force threshold for identifying foot contact
forceThreshold = 20

# %% Create the measurement set to use in scaling

#Set the parameters for the measurement sets
measurementSetParams = {
    #Torso
    'torso_depth': {'markerPairs': [['CV7', 'SJN'],['TV10', 'SXS'],], 'bodyScale': ['torso'], 'axes': 'X'},
    'torso_height': {'markerPairs': [['CV7', 'TV10'],['SJN', 'SXS'],], 'bodyScale': ['torso'], 'axes': 'Y'},
    'torso_width': {'markerPairs': [['R_SAE', 'L_SAE'],['R_SAA', 'L_SAA'],], 'bodyScale': ['torso'], 'axes': 'Z'},
    #Pelvis
    'pelvis_depth': {'markerPairs': [['R_IAS', 'R_IPS'],['L_IAS', 'L_IPS'],], 'bodyScale': ['pelvis'], 'axes': 'X'},
    'pelvis_height': {'markerPairs': [['R_IAS', 'R_FTC'],['L_IAS', 'L_FTC'],], 'bodyScale': ['pelvis'], 'axes': 'Y'},
    'pelvis_width': {'markerPairs': [['R_IAS', 'L_IAS'],['R_IPS', 'L_IPS'],], 'bodyScale': ['pelvis'], 'axes': 'Z'},
    #Right thigh
    'r_thigh_length': {'markerPairs': [['R_FTC', 'R_FLE'],], 'bodyScale': ['femur_r'], 'axes': 'Y'},
    'r_thigh_breadth': {'markerPairs': [['R_FLE', 'R_FME'],], 'bodyScale': ['femur_r'], 'axes': 'XZ'},
    #Right patella
    'r_patella': {'markerPairs': [['R_FTC', 'R_FLE'],], 'bodyScale': ['patella_r'], 'axes': 'XYZ'},
    #Right shank
    'r_tibia_length': {'markerPairs': [['R_TTC', 'R_TAM'],], 'bodyScale': ['tibia_r', 'talus_r'], 'axes': 'Y'},
    'r_tibia_breadth': {'markerPairs': [['R_TAM', 'R_FAL'],], 'bodyScale': ['tibia_r', 'talus_r'], 'axes': 'XZ'},
    #Right foot
    'r_foot_length': {'markerPairs': [['R_FCC', 'R_FM2'],], 'bodyScale': ['calcn_r', 'toes_r'], 'axes': 'X'},
    'r_foot_width': {'markerPairs': [['R_FM1', 'R_FM5'],], 'bodyScale': ['calcn_r', 'toes_r'], 'axes': 'YZ'},
    #Left thigh
    'l_thigh_length': {'markerPairs': [['L_FTC', 'L_FLE'],], 'bodyScale': ['femur_l'], 'axes': 'Y'},
    'l_thigh_breadth': {'markerPairs': [['L_FLE', 'L_FME'],], 'bodyScale': ['femur_l'], 'axes': 'XZ'},
    #Left patella
    'l_patella': {'markerPairs': [['L_FTC', 'L_FLE'],], 'bodyScale': ['patella_l'], 'axes': 'XYZ'},
    #Left shank
    'l_tibia_length': {'markerPairs': [['L_TTC', 'L_TAM'],], 'bodyScale': ['tibia_l', 'talus_l'], 'axes': 'Y'},
    'l_tibia_breadth': {'markerPairs': [['L_TAM', 'L_FAL'],], 'bodyScale': ['tibia_l', 'talus_l'], 'axes': 'XZ'},
    #Left foot
    'l_foot_length': {'markerPairs': [['L_FCC', 'L_FM2'],], 'bodyScale': ['calcn_l', 'toes_l'], 'axes': 'X'},
    'l_foot_width': {'markerPairs': [['L_FM1', 'L_FM5'],], 'bodyScale': ['calcn_l', 'toes_l'], 'axes': 'YZ'},
    #Right upper arm
    'r_upperarm_length': {'markerPairs': [['R_SAA', 'R_UOA'],], 'bodyScale': ['humerus_r'], 'axes': 'Y'},
    'r_upperarm_breadth': {'markerPairs': [['R_HLE', 'R_HME'],], 'bodyScale': ['humerus_r'], 'axes': 'XZ'},
    #Right forearm
    'r_forearm_length': {'markerPairs': [['R_HLE', 'R_RSP'],['R_HME', 'R_UHE'],], 'bodyScale': ['radius_r', 'ulna_r'], 'axes': 'Y'},
    'r_forearm_breadth': {'markerPairs': [['R_RSP', 'R_UHE'],], 'bodyScale': ['radius_r', 'ulna_r'], 'axes': 'XZ'},
    #Right hand
    'r_hand_length': {'markerPairs': [['R_RSP', 'R_HM2'],['R_UHE', 'R_HM5'],], 'bodyScale': ['hand_r'], 'axes': 'Y'},
    'r_hand_breadth': {'markerPairs': [['R_HM2', 'R_HM5'],], 'bodyScale': ['hand_r'], 'axes': 'XZ'},
    #Left upper arm
    'l_upperarm_length': {'markerPairs': [['L_SAA', 'L_UOA'],], 'bodyScale': ['humerus_l'], 'axes': 'Y'},
    'l_upperarm_breadth': {'markerPairs': [['L_HLE', 'L_HME'],], 'bodyScale': ['humerus_l'], 'axes': 'XZ'},
    #Left forearm
    'l_forearm_length': {'markerPairs': [['L_HLE', 'L_RSP'],['L_HME', 'L_UHE'],], 'bodyScale': ['radius_l', 'ulna_l'], 'axes': 'Y'},
    'l_forearm_breadth': {'markerPairs': [['L_RSP', 'L_UHE'],], 'bodyScale': ['radius_l', 'ulna_l'], 'axes': 'XZ'},
    #Left hand
    'l_hand_length': {'markerPairs': [['L_RSP', 'L_HM2'],['L_UHE', 'L_HM5'],], 'bodyScale': ['hand_l'], 'axes': 'Y'},
    'l_hand_breadth': {'markerPairs': [['L_HM2', 'L_HM5'],], 'bodyScale': ['hand_l'], 'axes': 'XZ'},    
    }

#Create the measurement set
scaleMeasurementSet = osim.MeasurementSet()

#Append the measurements from parameters
for measureName in measurementSetParams.keys():
    #Create the measurement
    measurement = osim.Measurement()
    measurement.setName(measureName)
    #Append the marker pairs
    for ii in range(len(measurementSetParams[measureName]['markerPairs'])):
        measurement.getMarkerPairSet().cloneAndAppend(
            osim.MarkerPair(measurementSetParams[measureName]['markerPairs'][ii][0],
                            measurementSetParams[measureName]['markerPairs'][ii][1]))
    #Append the body scales
    for ii in range(len(measurementSetParams[measureName]['bodyScale'])):
        #Create body scale
        bodyScale = osim.BodyScale()
        bodyScale.setName(measurementSetParams[measureName]['bodyScale'][ii])
        #Create and set axis names
        axes = osim.ArrayStr()
        for jj in range(len(measurementSetParams[measureName]['axes'])):
            axes.append(measurementSetParams[measureName]['axes'][jj])
        bodyScale.setAxisNames(axes)
        #Apppend to body scale set
        measurement.getBodyScaleSet().cloneAndAppend(bodyScale)
    #Append the measurement to the set
    scaleMeasurementSet.cloneAndAppend(measurement)
    
# %% Create the IK task set for scaling and tracking

""" NOTE: Currently just defaulting to all 1's... """

#Set the parameters for the IK task sets
ikTaskSetParams = {
    #Torso
    'SJN': {'weight': 1}, 'SXS': {'weight': 1}, 'CV7': {'weight': 1}, 'TV10': {'weight': 1}, 'R_SAE': {'weight': 1}, 'L_SAE': {'weight': 1},
    #Pelvis 
    'L_IAS': {'weight': 1}, 'L_IPS': {'weight': 1}, 'R_IAS': {'weight': 1}, 'R_IPS': {'weight': 1},
    #Right thigh
    'R_FTC': {'weight': 1}, 'R_FLE': {'weight': 1}, 'R_FME': {'weight': 1}, 
    #Right shank
    'R_TTC': {'weight': 1}, 'R_FAX': {'weight': 1}, 'R_TAM': {'weight': 1}, 'R_FAL': {'weight': 1},
    #Right foot
    'R_FCC': {'weight': 1}, 'R_FM1': {'weight': 1}, 'R_FM2': {'weight': 1}, 'R_FM5': {'weight': 1},
    #Left thigh
    'L_FTC': {'weight': 1}, 'L_FLE': {'weight': 1}, 'L_FME': {'weight': 1}, 
    #Left shank
    'L_TTC': {'weight': 1}, 'L_FAX': {'weight': 1}, 'L_TAM': {'weight': 1}, 'L_FAL': {'weight': 1},
    #Left foot
    'L_FCC': {'weight': 1}, 'L_FM1': {'weight': 1}, 'L_FM2': {'weight': 1}, 'L_FM5': {'weight': 1},
    #Right upper arm
    'R_SAA': {'weight': 1}, 'R_HLE': {'weight': 1}, 'R_HME': {'weight': 1},
    #Right forearm
    'R_UOA': {'weight': 1}, 'R_RSP': {'weight': 1}, 'R_UHE': {'weight': 1},
    #Right hand
    'R_HM2': {'weight': 1}, 'R_HM5': {'weight': 1},
    #Left upper arm
    'L_SAA': {'weight': 1}, 'L_HLE': {'weight': 1}, 'L_HME': {'weight': 1},
    #Left forearm
    'L_UOA': {'weight': 1}, 'L_RSP': {'weight': 1}, 'L_UHE': {'weight': 1},
    #Left hand
    'L_HM2': {'weight': 1}, 'L_HM5': {'weight': 1},
    }

#Create the task set
ikTaskSet = osim.IKTaskSet()

#Append the tasks from the parameters
for taskName in ikTaskSetParams.keys():
    #Create the task and add details
    task = osim.IKMarkerTask()
    task.setName(taskName)
    task.setWeight(ikTaskSetParams[taskName]['weight'])
    #Append to task set
    ikTaskSet.cloneAndAppend(task)

# %% Loop through participants

for participant in participantList:
    
    #### TODO: check for OK paticipant and trials
    
    # %% Check and select participant trials
    
    """
    For a valid participant they need at least three trials from each listed 
    condition with two force plate contacts.
    
    """
    
    #Get the external loads files from each condition
    exLoadFiles = [glob.glob(os.path.join(participant, 'dynamic', f'*_{condition}*_grf.xml')) for condition in conditionList]
    exLoadFilesUse = []
    
    #Create a counter for the number of 'good' files in each condition
    conditionCounter = {condition: 0 for condition in conditionList}
    
    #Loop through conditions and files to get the count
    for ii in range(len(conditionList)):
        exLoadFilesUse.append([])
        for loadFile in exLoadFiles[ii]:
            #Load the file and count the external loads
            #If it is at least 2 then we're good
            if osim.ExternalLoads(loadFile, True).getSize() >= 2:
                #Add to counter
                conditionCounter[conditionList[ii]]  += 1
                #Add the file to the usable list
                exLoadFilesUse[ii].append(loadFile)
                
    #Check if all values are at least 3
    if (np.array([conditionCounter[condition] for condition in conditionList]) >= 3).all():
        #Set bool to process participant
        processParticipant = True        
        #Select a file list of 3 trials from each condition
        dynamicTrials = {condition: [] for condition in conditionList}
        for ii in range(len(conditionList)):
            #Set random seed for choice as participant id
            random.seed(int(participant))
            #Select from current condition
            conditionFiles = random.choices(exLoadFilesUse[ii], k = 3)
            #Partition the trial name out and drop the grf label
            #Put into dynamic files dictionary
            dynamicTrials[conditionList[ii]] = [os.path.split(file)[-1].replace('_grf.xml','') for file in conditionFiles]
    else:
        #Set bool to not process participant
        processParticipant = False
    
    #Get participants static trial
    staticTRC = glob.glob(os.path.join(participant, 'static', '*_ST.trc'))[0]
    
    # %% Process data
    
    if processParticipant:
    
        # %% Scale the generic model
        
        if runScaling:
        
            #Make scaling directory for participant
            os.makedirs(os.path.join(participant, 'scaling'), exist_ok = True)
            
            # =====================================================================
            #     Scale model
            # =====================================================================
            
            #Create scale tool for current participant
            #Perhaps check if tool exists to avoid kernel crahses
            if participant not in scaleTool.keys():
                scaleTool[participant] = osim.ScaleTool()
            
            #Set participant mass
            scaleTool[participant].setSubjectMass(
                anthropometrics[anthropometrics['subjectID'] == int(participant)]['mass'].values[0]
                )
            
            #Set generic model file
            scaleTool[participant].getGenericModelMaker().setModelFileName('Rajagopal2015_SchreiberMoissenet2019.osim')
            
            #Set measurement set in model scaler
            scaleTool[participant].getModelScaler().setMeasurementSet(scaleMeasurementSet)
            
            #Set scale tasks in tool
            for ii in range(ikTaskSet.getSize()):
                scaleTool[participant].getMarkerPlacer().getIKTaskSet().cloneAndAppend(ikTaskSet.get(ii))
                
            #Set marker file
            scaleTool[participant].getMarkerPlacer().setMarkerFileName(staticTRC)
            scaleTool[participant].getModelScaler().setMarkerFileName(staticTRC)
            
            #Set options
            scaleTool[participant].getModelScaler().setPreserveMassDist(True)
            scaleOrder = osim.ArrayStr()
            scaleOrder.set(0,'measurements')
            scaleTool[participant].getModelScaler().setScalingOrder(scaleOrder)
            
            #Set time ranges
            timeRange = osim.ArrayDouble()
            timeRange.set(0,0) #initial time
            timeRange.set(1,1) #final time
            scaleTool[participant].getMarkerPlacer().setTimeRange(timeRange)
            scaleTool[participant].getModelScaler().setTimeRange(timeRange)
            
            #Set output files
            scaleTool[participant].getModelScaler().setOutputModelFileName(os.path.join(participant, 'scaling', f'{participant}_scaledModel.osim'))
            scaleTool[participant].getModelScaler().setOutputScaleFileName(os.path.join(participant, 'scaling', f'{participant}_scaleSet.xml'))
            
            #Set marker adjustment parameters
            scaleTool[participant].getMarkerPlacer().setOutputMotionFileName(os.path.join(participant, 'scaling', f'{participant}_staticMotion.mot'))
            scaleTool[participant].getMarkerPlacer().setOutputModelFileName(os.path.join(participant, 'scaling', f'{participant}_scaledModelAdjusted.osim'))
            
            #Save and run scale tool
            scaleTool[participant].printToXML(os.path.join(participant, 'scaling', f'{participant}_scaleSetup.xml'))
            scaleTool[participant].run()
            
            # =====================================================================
            #     Adjust model
            # =====================================================================
            
            #Add marker locations underneath foot markers at floor level based on static motion
            #These may be useful later in determining foot-ground contact points
            
            #Set the list of markers to project to the floot
            floorMarkers = ['R_FM1', 'R_FM2', 'R_FM5', 'R_FCC',
                            'L_FM1', 'L_FM2', 'L_FM5', 'L_FCC']
            
            #Load the model
            scaledModel = osim.Model(os.path.join(participant, 'scaling', f'{participant}_scaledModelAdjusted.osim'))
            
            #Load the static motion
            staticMotion = osim.TimeSeriesTable(os.path.join(participant, 'scaling', f'{participant}_staticMotion.mot'))
            
            #Initialise the model state
            modelState = scaledModel.initSystem()
            
            #Set model states from static motion
            for stateName in staticMotion.getColumnLabels():
                scaledModel.setStateVariableValue(
                    modelState,
                    stateName,
                    staticMotion.getDependentColumn(stateName).to_numpy()[0])
                
            #Realize model to position stage
            scaledModel.realizePosition(modelState)
            
            #Loop through floor markers
            #Identify position in ground at floor level and translate back to associated body
            for marker in floorMarkers:
                
                #Get marker location in ground        
                markerInGround = scaledModel.updMarkerSet().get(marker).getLocationInGround(modelState)
            
                #Set the y-axis to be zero for the marker in ground position (i.e. floor level)
                markerInGround.set(1,0)
                
                #Translate this new position back to the markers associated body
                floorMarkerLoc = scaledModel.getGround().findStationLocationInAnotherFrame(
                    modelState, markerInGround, scaledModel.updMarkerSet().get(marker).getParentFrame()
                    )
                
                #Create new marker and append to model
                newMarker = osim.Marker()
                newMarker.setName(marker+'_ground')
                newMarker.setParentFrameName(scaledModel.updMarkerSet().get(marker).getParentFrameName())
                newMarker.set_location(floorMarkerLoc)
                scaledModel.getMarkerSet().cloneAndAppend(newMarker)
                
            #Finalise model connections
            scaledModel.finalizeConnections()
            
            #Print to file
            scaledModel.printToXML(os.path.join(participant, 'scaling', f'{participant}_scaledModelAdjusted_floorMarkers.osim'))
        
        # %% Run IK as a baseline for kinematics
        
        if runIK:
        
            #Make IK directory for participant
            os.makedirs(os.path.join(participant, 'ik'), exist_ok = True)
            
            # =====================================================================
            #     Run IK
            # =====================================================================
            
            #Create an IK tool for current participant
            #Perhaps check if tool exists to avoid kernel crahses
            if participant not in ikTool.keys():
                ikTool[participant] = osim.InverseKinematicsTool()
                
            #Set model (consistent for all trials)
            ikTool[participant].set_model_file(os.path.join(participant, 'scaling', f'{participant}_scaledModelAdjusted_floorMarkers.osim'))
            
            #Set task set (consistent for all trials)
            for taskInd in range(ikTaskSet.getSize()):
                ikTool[participant].getIKTaskSet().adoptAndAppend(ikTaskSet.get(taskInd))
                
            #Set to report marker locations
            ikTool[participant].set_report_marker_locations(True)
            
            #Loop through the conditions
            for condition in conditionList:
            
                #Loop through the dynamic files
                for trial in dynamicTrials[condition]:
                    
                    #Set the marker file (relative to setup file location)
                    # ikTool[participant].setMarkerDataFileName(os.path.join(participant, 'dynamic', f'{trial}.trc'))
                    ikTool[participant].setMarkerDataFileName(os.path.join('..', 'dynamic', f'{trial}.trc'))
                    
                    #Set times
                    ikTool[participant].setStartTime(osim.TimeSeriesTableVec3(os.path.join(participant, 'dynamic', f'{trial}.trc')).getIndependentColumn()[0])
                    ikTool[participant].setEndTime(osim.TimeSeriesTableVec3(os.path.join(participant, 'dynamic', f'{trial}.trc')).getIndependentColumn()[-1])
                    
                    #Set output filename (relative to setup file location)
                    # ikTool[participant].setOutputMotionFileName(os.path.join(participant, 'ik', f'{trial}_ik.mot'))
                    ikTool[participant].setOutputMotionFileName(f'{trial}_ik.mot')
            
                    #Save IK tool to file
                    ikTool[participant].printToXML(os.path.join(participant, 'ik', f'{trial}_ikSetup.xml'))
    
                    #Bring the tool back in and run it (this seems to avoid Python kernel crashing)
                    ikRun = osim.InverseKinematicsTool(os.path.join(participant, 'ik', f'{trial}_ikSetup.xml'))
                    ikRun.run()
                    
                    #Rename supplementary marker outputs
                    shutil.move(os.path.join(participant, 'ik', '_ik_marker_errors.sto'),
                                os.path.join(participant, 'ik', f'{trial}_ikMarkerErrors.sto'))
                    shutil.move(os.path.join(participant, 'ik', '_ik_model_marker_locations.sto'),
                                os.path.join(participant, 'ik', f'{trial}_ikModelMarkerLocations.sto'))
        
        # %% Run the marker tracking problem
        
        if runTracking:
    
            #Make tracking directory for participant
            os.makedirs(os.path.join(participant, 'tracking'), exist_ok = True)
            
            # =====================================================================
            #     Run marker tracking
            # =====================================================================
            
            #Loop through the conditions
            for condition in conditionList:
            
                #Loop through the dynamic files
                for trial in dynamicTrials[condition]:
                    
                    #Create directory for trial and files
                    os.makedirs(os.path.join(participant, 'tracking', trial), exist_ok = True)
                    
                    #Create and name a Moco track object
                    #These don't seem to have the same kernel crashing issues as other tools
                    track = osim.MocoTrack()
                    track.setName(f'{trial}_markerTracking')
                    
                    # =============================================================
                    #     Set-up model and files
                    # =============================================================
        
                    #Construct a model processor to use with the tool
                    modelProc = osim.ModelProcessor(os.path.join(participant, 'scaling', f'{participant}_scaledModelAdjusted_floorMarkers.osim'))
                    
                    #Append external loads
                    modelProc.append(osim.ModOpAddExternalLoads(os.path.join(participant, 'dynamic', f'{trial}_grf.xml')))
                    
                    #Weld locked joints
                    jointsToWeld = ['subtalar_r', 'subtalar_l','mtp_r', 'mtp_l']
                    #Create vector string object
                    weldVectorStr = osim.StdVectorString()
                    [weldVectorStr.append(joint) for joint in jointsToWeld]
                    #Append to model processor
                    modelProc.append(osim.ModOpReplaceJointsWithWelds(weldVectorStr))
                    
                    #Process the model to edit
                    trackingModel = modelProc.process()
                    
                    #Add coordinate actuators to model
                    for coordinate in actForces:
                        #Create actuator
                        actu = osim.CoordinateActuator()
                        #Set name
                        actu.setName(f'{coordinate}_actuator')
                        #Set coordinate
                        actu.setCoordinate(trackingModel.updCoordinateSet().get(coordinate))
                        #Set optimal force
                        actu.setOptimalForce(actForces[coordinate])
                        #Set min and max control
                        actu.setMinControl(np.inf*-1)
                        actu.setMaxControl(np.inf*1)
                        #Append to model force set
                        trackingModel.updForceSet().cloneAndAppend(actu)
        
                    #Finalise model connections
                    trackingModel.finalizeConnections()
        
                    #Print model to file in tracking directory
                    trackingModel.printToXML(os.path.join(participant, 'tracking', trial, f'{participant}_trackingModel.osim'))
                    
                    #Copy TRC and GRF files to use
                    #External loads (maybe not needed...)
                    shutil.copyfile(os.path.join(participant, 'dynamic', f'{trial}_grf.xml'),
                                    os.path.join(participant, 'tracking', trial, f'{trial}_grf.xml'))
                    #Ground reactions
                    shutil.copyfile(os.path.join(participant, 'dynamic', f'{trial}_grf.mot'),
                                    os.path.join(participant, 'tracking', trial, f'{trial}_grf.mot'))
                    #Marker data
                    shutil.copyfile(os.path.join(participant, 'dynamic', f'{trial}.trc'),
                                    os.path.join(participant, 'tracking', trial, f'{trial}.trc'))
                    
                    #Navigate to tracking folder
                    os.chdir(os.path.join(participant, 'tracking', trial))
                    
                    # =============================================================
                    #     Set-up tracking problem
                    # =============================================================
            
                    #Set model
                    track.setModel(osim.ModelProcessor(f'{participant}_trackingModel.osim'))
                    
                    #Set the markers reference file
                    track.setMarkersReferenceFromTRC(f'{trial}.trc')
                    track.set_allow_unused_references(True)
                    
                    #Set global marker tracking weight
                    track.set_markers_global_tracking_weight(globalMarkerTrackingWeight)
                    
                    #Set individual marker tracking weights
                    #If all are set at 1, then this probably does nothing
                    markerWeights = osim.MocoWeightSet()
                    for marker in ikTaskSetParams.keys():
                        markerWeights.cloneAndAppend(osim.MocoWeight(marker, ikTaskSetParams[marker]['weight']))
                    track.set_markers_weight_set(markerWeights)
                    
                    #Set times
                    #Timing is based on when the trail toe-off happens alongside the
                    #first force plate contact, up until the toe-on of the lead leg
                    #after the second force plate contact. This is to ensure that there
                    #are no periods where a foot is in contact with the ground when
                    #no force plate data is present.
                    
                    #Get the c3d data with events
                    reader = btk.btkAcquisitionFileReader() 
                    reader.SetFilename(os.path.join(homeDir, participant, 'dynamic', f'{trial}.c3d'))
                    reader.Update()
                    c3dData = reader.GetOutput()
                    
                    #Get events
                    eventData = {'eventName': [], 'eventFrame': [], 'eventTime': []}
                    moreEvents = True
                    eventInd = 0
                    while moreEvents:
                        try:
                            event = c3dData.GetEvent(eventInd)
                            eventData['eventName'].append(event.GetContext()+' '+event.GetLabel())
                            eventData['eventFrame'].append(event.GetFrame())
                            eventData['eventTime'].append(event.GetTime())
                            eventInd += 1
                        except:
                            moreEvents = False
                    events = pd.DataFrame.from_dict(eventData)
                    events.sort_values(by = 'eventFrame', inplace = True)
                    events.reset_index(drop = True, inplace = True)
                    
                    #Use external loads data to identify timing of first foot strike
                    #Final set times fill be from firstFootOff to lastFootStrike
                    
                    #Read in GRF data
                    grfData = osim.TimeSeriesTable(f'{trial}_grf.mot')
                    
                    #Examine vertical force data and find time of first foot strike
                    v1_t = grfData.getIndependentColumn()[np.argmax(grfData.getDependentColumn('ground_force_1_vy').to_numpy() > forceThreshold)]
                    v2_t = grfData.getIndependentColumn()[np.argmax(grfData.getDependentColumn('ground_force_2_vy').to_numpy() > forceThreshold)]
                    firstFootStrike = np.min((v1_t,v2_t))
                    
                    #Find the first foot off after the first foot strike
                    firstFootOff = np.min([events.iloc[ii]['eventTime'] for ii in range(len(events)) if 'Off' in events.iloc[ii]['eventName'] and events.iloc[ii]['eventTime'] > firstFootStrike])
                    
                    #Get index of second foot strike
                    secondFootStrike = np.max((v1_t,v2_t))
                    
                    #Find the next foot strike after this
                    lastFootStrike = np.min([events.iloc[ii]['eventTime'] for ii in range(len(events)) if 'Strike' in events.iloc[ii]['eventName'] and events.iloc[ii]['eventTime'] > secondFootStrike])
                    
                    #Set times
                    track.set_initial_time(firstFootOff)
                    track.set_final_time(lastFootStrike)
                    track.set_mesh_interval(0.05) #should result in a ~50 nodes
                    
                    #Initialise to a Moco study and problem
                    study = track.initialize()
                    problem = study.updProblem()
                    
                    # Get a reference to the MocoControlCost goal
                    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
                
                    #Put a large weight on the pelvis CoordinateActuators, which act as the
                    #residual, or 'hand-of-god', forces which we would like to keep as small
                    #as possible.
                    for pelvisCoord in ['pelvis_tx', 'pelvis_ty', 'pelvis_tz', 'pelvis_tilt', 'pelvis_list', 'pelvis_rotation']:
                        effort.setWeightForControl(f'/forceset/{pelvisCoord}_actuator', 10)
                    
                    #Solve the problem
                    solution = study.solve()

                    #Check if solution is sealed
                    if solution.isSealed():
                        solution.unseal()
                        
                    #WRite the solution to file
                    solution.write(f'{trial}_markerTracking_solution.sto')
                    
                    #Remove tracked markers file
                    os.remove(f'{trial}_markerTracking_tracked_markers.sto')
                    
                    #Return to home directory
                    os.chdir(homeDir)        
        
# %% ----- end of 2_processData.py ----- %% #