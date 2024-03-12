# -*- coding: utf-8 -*-
'''

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code processes the data from the HamnerDelp2013 dataset through both
    an OpenSim RRA and Moco states tracking framework. A torque-driven model is
    used to produce kinematics that track the kinematics while remaining consistent
    with external loads
    
'''

# %% Import packages

import opensim as osim
import os
import numpy as np
import pickle
import time
import re
import shutil
import osimFunctions as helper

# %% Flags for running analyses

#Change these boolean flags to choose which analyses to run
runRRA = True #doesn't take too long
runMoco = True #takes a little bit of time (~20 mins per sim)

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

#Set the participant list based on the codes in the folder
participantList = [ii for ii in os.listdir() if os.path.isdir(ii) and ii.startswith('subject')]

#Set run names list
#This is useful if you want to assess further running speeds that have been extracted
#Modify if changing the dataset to examine
runList  = [
    #'run2',
    'run3',
    #'run4',
    'run5',
    ]
runLabels  = [
    #'Run_2',
    'Run_3',
    #'Run_4',
    'Run_5',
    ]

#Set run cycle list
cycleList = ['cycle1',
             'cycle2',
             'cycle3']

#Create a dictionary of the coordinate tasks originally used in Hamner & Delp
rraTasks = {'pelvis_tx': 2.5e1, 'pelvis_ty': 1.0e2, 'pelvis_tz': 2.5e1,
            'pelvis_tilt': 7.5e2, 'pelvis_list': 2.5e2, 'pelvis_rotation': 5.0e1,
            'hip_flexion_r': 7.5e1, 'hip_adduction_r': 5.0e1, 'hip_rotation_r': 1.0e1,
            'knee_angle_r': 1.0e1, 'ankle_angle_r': 1.0e1,
            'hip_flexion_l': 7.5e1, 'hip_adduction_l': 5.0e1, 'hip_rotation_l': 1.0e1,
            'knee_angle_l': 1.0e1, 'ankle_angle_l': 1.0e1,
            'lumbar_extension': 7.5e1, 'lumbar_bending': 5.0e1, 'lumbar_rotation': 2.5e1,
            'arm_flex_r': 1.0e0, 'arm_add_r': 1.0e0, 'arm_rot_r': 1.0e0,
            'elbow_flex_r': 1.0e0, 'pro_sup_r': 1.0e0,
            'arm_flex_l': 1.0e0, 'arm_add_l': 1.0e0, 'arm_rot_l': 1.0e0,
            'elbow_flex_l': 1.0e0, 'pro_sup_l': 1.0e0
            }

#Create a dictionary of the optimal forces for actuators used in Hamner & Delp
rraActuators = {'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
                'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1,
                'hip_flexion_r': 1000, 'hip_adduction_r': 1000, 'hip_rotation_r': 1000,
                'knee_angle_r': 1000, 'ankle_angle_r': 1000,
                'hip_flexion_l': 1000, 'hip_adduction_l': 1000, 'hip_rotation_l': 1000,
                'knee_angle_l': 1000, 'ankle_angle_l': 1000,
                'lumbar_extension': 1000, 'lumbar_bending': 1000, 'lumbar_rotation': 1000,
                'arm_flex_r': 500, 'arm_add_r': 500, 'arm_rot_r': 500,
                'elbow_flex_r': 500, 'pro_sup_r': 500,
                'arm_flex_l': 500, 'arm_add_l': 500, 'arm_rot_l': 500,
                'elbow_flex_l': 500, 'pro_sup_l': 500
                }

#Create a dictionary for actuator limits used in Hamner & Delp
rraLimits = {'pelvis_tx': 10000, 'pelvis_ty': 10000, 'pelvis_tz': 10000,
             'pelvis_tilt': 10000, 'pelvis_list': 10000, 'pelvis_rotation': 10000,
             'hip_flexion_r': 1, 'hip_adduction_r': 1, 'hip_rotation_r': 1,
             'knee_angle_r': 1, 'ankle_angle_r': 1,
             'hip_flexion_l': 1, 'hip_adduction_l': 1, 'hip_rotation_l': 1,
             'knee_angle_l': 1, 'ankle_angle_l': 1,
             'lumbar_extension': 1, 'lumbar_bending': 1, 'lumbar_rotation': 1,
             'arm_flex_r': 1, 'arm_add_r': 1, 'arm_rot_r': 1,
             'elbow_flex_r': 1, 'pro_sup_r': 1,
             'arm_flex_l': 1, 'arm_add_l': 1, 'arm_rot_l': 1,
             'elbow_flex_l': 1, 'pro_sup_l': 1
             }

#Create a dictionary for kinematic boundary limits (+/- to max and min)
kinematicLimits = {'pelvis_tx': 0.2, 'pelvis_ty': 0.1, 'pelvis_tz': 0.2,
                   'pelvis_tilt': np.deg2rad(10), 'pelvis_list': np.deg2rad(10), 'pelvis_rotation': np.deg2rad(10),
                   'hip_flexion_r': np.deg2rad(10), 'hip_adduction_r': np.deg2rad(5), 'hip_rotation_r': np.deg2rad(5),
                   'knee_angle_r': np.deg2rad(15), 'ankle_angle_r': np.deg2rad(10),
                   'hip_flexion_l': np.deg2rad(10), 'hip_adduction_l': np.deg2rad(5), 'hip_rotation_l': np.deg2rad(5),
                   'knee_angle_l': np.deg2rad(15), 'ankle_angle_l': np.deg2rad(10),
                   'lumbar_extension': np.deg2rad(10), 'lumbar_bending': np.deg2rad(5), 'lumbar_rotation': np.deg2rad(5),
                   'arm_flex_r': np.deg2rad(5), 'arm_add_r': np.deg2rad(5), 'arm_rot_r': np.deg2rad(5),
                   'elbow_flex_r': np.deg2rad(10), 'pro_sup_r': np.deg2rad(5),
                   'arm_flex_l': np.deg2rad(5), 'arm_add_l': np.deg2rad(5), 'arm_rot_l': np.deg2rad(5),
                   'elbow_flex_l': np.deg2rad(10), 'pro_sup_l': np.deg2rad(5)
                   }

# %% Loop through participants

for participant in participantList:
    
    # %% Get and set info for current participant
    
    #Load in the subjects gait timing data
    with open(os.path.join(participant,'expData','gaitTimes.pkl'), 'rb') as openFile:
        gaitTimings = pickle.load(openFile)
        
    #Create an RRA directory in the subjects folder
    if runRRA:
        os.makedirs(os.path.join(participant,'rra'), exist_ok = True)        
        #Create run trial specific directories as well
        for runTrial in runList:
            os.makedirs(os.path.join(participant,'rra',runTrial), exist_ok = True)
            
    #Create a Moco directory in the subjects folder
    if runMoco:
        os.makedirs(os.path.join(participant,'moco'), exist_ok = True)        
        #Create run trial specific directories as well
        for runTrial in runList:
            os.makedirs(os.path.join(participant,'moco',runTrial), exist_ok = True)
            
    #Create dictionary to store timing data for RRA process
    if runRRA:
        rraRunTimeData = {run: {cyc: {'rraRunTime': []} for cyc in cycleList} for run in runList}

    #Create dictionary to store timing data
    if runMoco:
        mocoRunTimeData = {run: {cyc: {'mocoRunTime': [], 'nIters': [], 'solved': []} for cyc in cycleList} for run in runList}
    
    # %% Run RRA
    
    if runRRA:
        
        #Loop through run list
        for runTrial in runList:
            
            #Set the runName from list
            runName = runLabels[runList.index(runTrial)]
            
            # =================================================================
            #     Set-up for RRA
            # =================================================================
        
            #Change to rra directory for ease of use with tools
            os.chdir(os.path.join(participant,'rra',runTrial))
            
            #Add in opensim logger
            osim.Logger.removeFileSink()
            osim.Logger.addFileSink('rraLog.log')
            
            #Load the subject model to refer to body parameters
            osimModel = osim.Model(os.path.join('..','..','model',f'{participant}_adjusted_scaled.osim'))
        
            #Create dictionary to store mass adjustments
            bodyList = [osimModel.updBodySet().get(ii).getName() for ii in range(osimModel.updBodySet().getSize())]
            massAdjustmentData = {run: {cyc: {body: {'origMass': [], 'newMass': [], 'massChange': []} for body in bodyList} for cyc in cycleList} for run in runList}
            
            #Create the RRA actuators file
            
            #Create and set name
            rraForceSet = osim.ForceSet()
            rraForceSet.setName(f'{participant}_{runTrial}_RRA_Actuators')
            
            #Loop through coordinates and append to force set
            for actuator in rraActuators.keys():
                
                #Create the actuator. First we must check if point, torque or coordinate
                #actuators are required depending on the coordinate
                
                #Check for point actuator
                if actuator in ['pelvis_tx', 'pelvis_ty', 'pelvis_tz']:            
                    #Create the point actuator
                    pointActuator = osim.PointActuator()            
                    #Set the name to the residual coordinate
                    pointActuator.setName(f'F{actuator[-1].capitalize()}')            
                    #Set the max and min controls to those provided
                    pointActuator.set_min_control(rraLimits[actuator]*-1)
                    pointActuator.set_max_control(rraLimits[actuator])
                    #Set the force body as the pelvis
                    pointActuator.set_body('pelvis')
                    #Set the direction
                    pointActuator.set_direction(osim.Vec3(
                        (np.array([actuator[-1] == ii for ii in ['x','y','z']], dtype = int)[0],
                         np.array([actuator[-1] == ii for ii in ['x','y','z']], dtype = int)[1],
                         np.array([actuator[-1] == ii for ii in ['x','y','z']], dtype = int)[2])
                        ))            
                    #Set force to be global
                    pointActuator.set_point_is_global(True)
                    #Set the point from the model
                    pointActuator.set_point(osimModel.updBodySet().get('pelvis').get_mass_center())
                    #Set optimal force
                    pointActuator.set_optimal_force(rraActuators[actuator])
                    #Clone and append to force set
                    rraForceSet.cloneAndAppend(pointActuator)
                    
                #Check for torque actuator
                elif actuator in ['pelvis_list', 'pelvis_rotation', 'pelvis_tilt']:
                    #Create a torque actuator
                    torqueActuator = osim.TorqueActuator()
                    #Set the name to the residual coordinate            
                    torqueActuator.setName(f'M{[x for i, x in enumerate(["X","Y","Z"]) if [actuator == ii for ii in ["pelvis_list", "pelvis_rotation", "pelvis_tilt"]][i]][0]}')
                    #Set the max and min controls to those provided
                    torqueActuator.set_min_control(rraLimits[actuator]*-1)
                    torqueActuator.set_max_control(rraLimits[actuator])
                    #Set the torque to act on the pelvis relative to the ground
                    torqueActuator.set_bodyA('pelvis')
                    torqueActuator.set_bodyB('ground')
                    #Set the axis
                    torqueActuator.set_axis(osim.Vec3(
                        (np.array([actuator == ii for ii in ['pelvis_list', 'pelvis_rotation', 'pelvis_tilt']], dtype = int)[0],
                         np.array([actuator == ii for ii in ['pelvis_list', 'pelvis_rotation', 'pelvis_tilt']], dtype = int)[1],
                         np.array([actuator == ii for ii in ['pelvis_list', 'pelvis_rotation', 'pelvis_tilt']], dtype = int)[2])
                        ))
                    #Set torque to be global
                    torqueActuator.set_torque_is_global(True)
                    #Set optimal force
                    torqueActuator.set_optimal_force(rraActuators[actuator])
                    #Clone and append to force set
                    rraForceSet.cloneAndAppend(torqueActuator)
                    
                #Remaining should be coordinate actuators
                else:
                    #Create a coordinate actuator
                    coordActuator = osim.CoordinateActuator()
                    #Set name to coordinate
                    coordActuator.setName(actuator)
                    #Set coordinate
                    coordActuator.set_coordinate(actuator)
                    #Set min and max control to those provided
                    coordActuator.set_min_control(rraLimits[actuator]*-1)
                    coordActuator.set_max_control(rraLimits[actuator])
                    #Set optimal force
                    coordActuator.set_optimal_force(rraActuators[actuator])
                    #Clone and append to force set
                    rraForceSet.cloneAndAppend(coordActuator)
                    
            #Print the force set to file
            rraForceSet.printToXML(f'{participant}_{runTrial}_RRA_Actuators.xml')
            
            #Create the RRA tasks file
            
            #Create and set name
            rraTaskSet = osim.CMC_TaskSet()
            rraTaskSet.setName(f'{participant}_{runTrial}_RRA_Tasks')
            
            #Loop through coordinates and append to force set
            for task in rraTasks.keys():
                
                #Create the task
                cmcTask = osim.CMC_Joint()
                
                #Set the name to the coordinate
                cmcTask.setName(task)
                
                #Set task weight
                cmcTask.setWeight(rraTasks[task])
                
                #Set active parameters
                cmcTask.setActive(True, False, False)
                
                #Set kp and kv
                cmcTask.setKP(100)
                cmcTask.setKV(20)
                
                #Set coordinate
                cmcTask.setCoordinateName(task)
                
                #Clone and append to task set
                rraTaskSet.cloneAndAppend(cmcTask)
                    
            #Print the force set to file
            rraTaskSet.printToXML(f'{participant}_{runTrial}_RRA_Tasks.xml')
        
            # =================================================================
            #     Run the RRA
            # =================================================================
            
            #Create a generic RRA tool to manipulate for the 3 cycles
            rraTool = osim.RRATool()
            
            #Set the generic elements in the tool
            
            #Model file
            rraTool.setModelFilename(os.path.join('..','..','model',f'{participant}_adjusted_scaled.osim'))
            
            #Append the force set files
            forceSetFiles = osim.ArrayStr()
            forceSetFiles.append(f'{participant}_{runTrial}_RRA_Actuators.xml')
            rraTool.setForceSetFiles(forceSetFiles)
            rraTool.setReplaceForceSet(True)
            
            #External loads file
            rraTool.setExternalLoadsFileName(os.path.join('..','..','expData',f'{runName}_grf.xml'))
    
            #Kinematics file
            rraTool.setDesiredKinematicsFileName(os.path.join('..','..','ik',f'{runName}.mot'))
            
            #Cutoff frequency for kinematics
            rraTool.setLowpassCutoffFrequency(15.0)
            
            #Task set file
            rraTool.setTaskSetFileName(f'{participant}_{runTrial}_RRA_Tasks.xml')
            
            #Output precision
            rraTool.setOutputPrecision(20)
            
            #Loop through gait cycles
            for cycle in cycleList:
                
                #Create directory for cycle
                os.makedirs(cycle, exist_ok = True)
                
                #Add in opensim logger for cycle
                osim.Logger.removeFileSink()
                osim.Logger.addFileSink(os.path.join(cycle,f'{runTrial}_{cycle}_rraLog.log'))
                
                #Add in cycle specific details
                
                #Tool name
                rraTool.setName(f'{participant}_{runTrial}_{cycle}')
                
                #Start and end time
                rraTool.setInitialTime(gaitTimings[runTrial][cycle]['initialTime'])
                rraTool.setFinalTime(gaitTimings[runTrial][cycle]['finalTime'])
                
                #Results directory
                rraTool.setResultsDir(f'{cycle}/')
                
                #Output model file
                rraTool.setOutputModelFileName(f'{cycle}/{participant}_{runTrial}_{cycle}_rraAdjusted.osim')
                
                #Adjusted COM body
                rraTool.setAdjustCOMToReduceResiduals(True)
                rraTool.setAdjustedCOMBody('torso')
                
                #Print to file
                rraTool.printToXML(f'{participant}_{runTrial}_{cycle}_setupRRA.xml')
                
                #Load and run rra tool
                #For some reason rra works better when the tool is reloaded
                rraToolRun = osim.RRATool(f'{participant}_{runTrial}_{cycle}_setupRRA.xml')
                
                #Set-up start timer
                startRunTime = time.time()
                
                #Run tool        
                rraToolRun.run()
                
                #End timer and record
                rraRunTime = round(time.time() - startRunTime, 2)
                
                #Record run-time to dictionary
                rraRunTimeData[runTrial][cycle]['rraRunTime'] = rraRunTime
                
                # =============================================================
                #     Review RRA results
                # =============================================================
                
                #Mass adjustments
                #Stop the logger
                osim.Logger.removeFileSink()
                #Read in the log file
                fid = open(os.path.join(cycle,f'{runTrial}_{cycle}_rraLog.log'), 'r')
                fileText = fid.readlines()
                fid.close()
                #Loop through the bodies
                for body in bodyList:
                    #Search through log file lines for current body adjustment
                    for li in fileText:
                        if body in li and 'orig mass' in li and 'new mass' in li:
                            #Extract out the original mass and new mass
                            stringToGetOrig = re.search('orig mass = (.*),', li)
                            stringToGetNew = re.search('new mass = (.*)\n', li)
                            #Get the values and append to dictionary
                            massAdjustmentData[runTrial][cycle][body]['origMass'] = float(stringToGetOrig.group(1))
                            massAdjustmentData[runTrial][cycle][body]['newMass'] = float(stringToGetNew.group(1))
                            massAdjustmentData[runTrial][cycle][body]['massChange'] = float(stringToGetNew.group(1)) - float(stringToGetOrig.group(1))
                
                #Adjust mass in the newly created model
                #Load the model
                rraAdjustedModel = osim.Model(os.path.join(cycle,f'{participant}_{runTrial}_{cycle}_rraAdjusted.osim'))
                #Loop through the bodies and set the mass from the dictionary
                for body in bodyList:
                    #Get the new mass
                    newMass = massAdjustmentData[runTrial][cycle][body]['newMass']
                    #Update in the model
                    rraAdjustedModel.updBodySet().get(body).setMass(newMass)
                #Finalise the model connections
                rraAdjustedModel.finalizeConnections()
                #Re-save the model
                rraAdjustedModel.printToXML(os.path.join(cycle,f'{participant}_{runTrial}_{cycle}_rraAdjusted.osim'))
                                
            #Navigate back to home directory for next run trial or subject
            os.chdir(homeDir)
                
        #Save run time and mass adjustment data dictionaries
        with open(os.path.join(participant, 'rra', f'{participant}_rraRunTimeData.pkl'), 'wb') as writeFile:
            pickle.dump(rraRunTimeData, writeFile)
        with open(os.path.join(participant, 'rra', f'{participant}_massAdjustmentData.pkl'), 'wb') as writeFile:
            pickle.dump(massAdjustmentData, writeFile)
            
    # %% Run Moco
    
    if runMoco:
        
        #Loop through run list
        for runTrial in runList:
            
            #Set the runName from list
            runName = runLabels[runList.index(runTrial)]
            
            # =================================================================
            #     Set-up for Moco
            # =================================================================
        
            #Change to Moco directory for ease of use with tools
            os.chdir(os.path.join(participant,'moco',runTrial))
    
            #Add in opensim logger
            osim.Logger.removeFileSink()
            osim.Logger.addFileSink('mocoLog.log')
        
            #Copy external load files across as there are issues with using these out of
            #directory with Moco tools
            shutil.copyfile(os.path.join('..','..','expData',f'{runName}_grf.xml'),
                            f'{runName}_grf.xml')
            shutil.copyfile(os.path.join('..','..','expData',f'{runName}_grf.mot'),
                            f'{runName}_grf.mot')
            
            #Convert kinematics to states version for use with Moco
            helper.kinematicsToStates(kinematicsFileName = os.path.join('..','..','ik',f'{runName}.mot'),
                               osimModelFileName = os.path.join('..','..','model',f'{participant}_adjusted_scaled.osim'),
                               outputFileName = f'{runName}_coordinates.sto',
                               inDegrees = True, outDegrees = False,
                               filtFreq = 15.0)
        
            # =================================================================
            #     Run the Moco tracking
            # =================================================================
        
            #Create a generic tracking tool to manipulate for the 3 cycles
            mocoTrack = osim.MocoTrack()
            mocoTrack.setName('mocoResidualReduction')
            
            #Construct a ModelProcessor and set it on the tool.
            modelProcessor = osim.ModelProcessor(os.path.join('..','..','model',f'{participant}_adjusted_scaled.osim'))
            modelProcessor.append(osim.ModOpAddExternalLoads(f'{runName}_grf.xml'))
            modelProcessor.append(osim.ModOpRemoveMuscles())
            
            #Process model to edit
            mocoModel = modelProcessor.process()
            
            #Add in torque actuators that replicate the RRA actuators
            mocoModel = helper.addTorqueActuators(osimModel = mocoModel,
                                                  optForces = rraActuators,
                                                  controlLimits = rraLimits)
            
            #Set model in tracking tool
            mocoTrack.setModel(osim.ModelProcessor(mocoModel))
            
            #Construct a table processor to append to the tracking tool for kinematics
            #The kinematics can't be filtered here with the operator as it messes with
            #time stamps in a funky way. This however has already been done in the 
            #conversion to state coordinates
            tableProcessor = osim.TableProcessor(f'{runName}_coordinates.sto')
            mocoTrack.setStatesReference(tableProcessor)
            
            #Create a dictionary to set kinematic bounds
            #Create this based on maximum and minimum values in the kinematic data
            #plus/minus some generic values
            
            #Load the kinematics file as a table
            ikTable = osim.TimeSeriesTable(f'{runName}_coordinates.sto')
            
            #Create the bounds dictionary
            kinematicBounds = {}
            #Loop through the coordinates
            for coord in kinematicLimits.keys():
                #Get the coordinate path
                coordPath = mocoModel.updCoordinateSet().get(coord).getAbsolutePathString()+'/value'
                #Set bounds in dictionary
                kinematicBounds[coord] = [ikTable.getDependentColumn(coordPath).to_numpy().min() - kinematicLimits[coord],
                                          ikTable.getDependentColumn(coordPath).to_numpy().max() + kinematicLimits[coord]]
        
            #Set the global states tracking weight in the tracking problem
            mocoTrack.set_states_global_tracking_weight(1)
            
            #Set tracking tool to apply states to guess
            mocoTrack.set_apply_tracked_states_to_guess(True)
            
            #Provide the setting to ignore unused columns in kinematic data
            mocoTrack.set_allow_unused_references(True)
            
            #Set Moco to track the speed derivatives from kinematic data
            #This is probably a fair comparison to RRA given that it tracks accelerations
            #The case study test of this looked like it would help with smoothness as well
            mocoTrack.set_track_reference_position_derivatives(True)
            
            #Set tracking mesh interval time
            mocoTrack.set_mesh_interval(0.01) 
            
            #Set the coordinate reference task weights to match RRA
        
            #Create weight set for state tracking
            stateWeights = osim.MocoWeightSet()
            
            #Set a scaling factor for tracking the speeds
            speedsTrackingScale = 0.01
            
            #Loop through coordinates to apply weights
            for coordInd in range(mocoModel.updCoordinateSet().getSize()):
                
                #Get name and absolute path to coordinate
                coordName = mocoModel.updCoordinateSet().get(coordInd).getName()
                coordPath = mocoModel.updCoordinateSet().get(coordInd).getAbsolutePathString()
            
                #If a task weight is provided, add it in
                if coordName in list(rraTasks.keys()):
                    #Append state into weight set
                    #Track the coordinate value
                    stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/value',
                                                                rraTasks[coordName]))
                    #Don't track the Coordinate speed
                    stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/speed',
                                                                rraTasks[coordName] * speedsTrackingScale))
                    
            #Add state weights to the tracking tool
            mocoTrack.set_states_weight_set(stateWeights)
        
            #Loop through gait cycles
            for cycle in cycleList:
            
                #Create directory for cycle
                os.makedirs(cycle, exist_ok = True)
                
                #Add in opensim logger for cycle
                osim.Logger.removeFileSink()
                osim.Logger.addFileSink(os.path.join(cycle,f'{runTrial}_{cycle}_mocoLog.log'))
                
                #Add in cycle specific details
            
                #Set the gait timings in tracking tool
                mocoTrack.set_initial_time(gaitTimings[runTrial][cycle]['initialTime'])
                mocoTrack.set_final_time(gaitTimings[runTrial][cycle]['finalTime'])
                
                #Initialise the Moco study
                study = mocoTrack.initialize()
                problem = study.updProblem()
                
                #Set the parameters for the regularization term on MocoTrack problem
                #(minimize squared excitations)
                effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
                effort.setWeight(0.001)
                        
                #Lock time bounds to the IK data
                problem.setTimeBounds(gaitTimings[runTrial][cycle]['initialTime'],
                                      gaitTimings[runTrial][cycle]['finalTime'])
                
                #Set kinematic bounds using the dictionary values and experimental data
                for coordInd in range(mocoModel.updCoordinateSet().getSize()):
                    #First check if coordinate is in kinematic bounds dictionary
                    if mocoModel.updCoordinateSet().get(coordInd).getName() in list(kinematicBounds.keys()):                
                        #Get coordinate name and path
                        coordName = mocoModel.updCoordinateSet().get(coordInd).getName()
                        coordPath = mocoModel.updCoordinateSet().get(coordInd).getAbsolutePathString()+'/value'
                        #Set bounds in problem
                        problem.setStateInfo(coordPath,
                                             #Bounds set to model ranges
                                             [kinematicBounds[coordName][0], kinematicBounds[coordName][1]]
                                             )
                
                #Get the solver
                solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
                        
                #Set solver parameters
                solver.set_optim_max_iterations(1000)
                solver.set_optim_constraint_tolerance(1e-2)
                solver.set_optim_convergence_tolerance(1e-2)
                
                #Reset problem (required if changing to implicit mode)
                solver.resetProblem(problem)
                
                #Print to file
                study.printToXML(f'{participant}_{runTrial}_{cycle}_setupMoco.omoco')
                
                #Set-up start timer
                startRunTime = time.time()
                
                #Solve!       
                solution = study.solve()
                
                #End timer and record
                mocoRunTime = round(time.time() - startRunTime, 2)
                
                #Record run-time to dictionary
                mocoRunTimeData[runTrial][cycle]['mocoRunTime'] = mocoRunTime
                
                # =============================================================
                #     Review Moco results
                # =============================================================
                
                #Check need to unseal and store outcome
                if solution.isSealed():
                    solution.unseal()
                    mocoRunTimeData[runTrial][cycle]['solved'] = False
                else:
                    mocoRunTimeData[runTrial][cycle]['solved'] = True
                    
                #Store number of iterations
                mocoRunTimeData[runTrial][cycle]['nIters'] = solution.getNumIterations()
        
                #Write solution to file
                solution.write(os.path.join(cycle,f'{participant}_{runTrial}_{cycle}_mocoSolution.sto'))
                
                #Remove tracked states file
                os.remove('mocoResidualReduction_tracked_states.sto')
                
                #Stop the logger
                osim.Logger.removeFileSink()
                                
                #Convert the solution to a states table and back to standard kinematic coordinates for ease of use
                
                #Write states table to file
                osim.STOFileAdapter().write(solution.exportToStatesTable(),
                                            os.path.join(cycle,f'{participant}_{runTrial}_{cycle}_mocoStates.sto'))
                
                #Convert states back to kinematic coordinates with helper function
                helper.statesToKinematics(statesFileName = os.path.join(cycle,f'{participant}_{runTrial}_{cycle}_mocoStates.sto'),
                                          outputFileName = os.path.join(cycle,f'{participant}_{runTrial}_{cycle}_mocoKinematics.sto'))
                
            #Navigate back to home directory for next run trial or subject
            os.chdir(homeDir)
            
        #Save run time and mass adjustment data dictionaries
        with open(os.path.join(participant, 'moco', f'{participant}_mocoRunTimeData.pkl'), 'wb') as writeFile:
            pickle.dump(mocoRunTimeData, writeFile)
        
# %% ----- end of 2_processData.py ----- %% #