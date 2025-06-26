# -*- coding: utf-8 -*-
'''

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code processes the data from the Rankin2016 dataset through a Moco states
    tracking framework. A torque-driven model is used to produce kinematics that
    track the kinematics while remaining consistent with external loads.
    
    NOTES:
        > Moco solution seems to come back to ground level for some reason?
    
'''

# %% Import packages

import opensim as osim
import os
import shutil
import osimFunctions as helper

# %% Flags for running analyses

#Change these boolean flags to choose which analyses to run
runMoco = True #takes a little bit of time (~30 mins per trial)

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
    
#Also add the model specific geometry path
osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.join(os.getcwd(), 'models', 'Geometry'))

#Set trial names list
#This is useful if you want to assess the different trials in the dataset
#Modify if changing the dataset to examine
trialList = [
    #'Walk',
    'Run',
    ]

#Create a dictionary of the coordinate tasks to track
taskWeights = {
    'pelvis_tx': 2.5e1, 'pelvis_ty': 1.0e2, 'pelvis_tz': 2.5e1,
    'pelvis_tilt': 7.5e2, 'pelvis_list': 2.5e2, 'pelvis_rotation': 5.0e1,
    'R_hip_extension': 7.5e1, 'R_hip_abduction': 5.0e1, 'R_hip_rotation': 1.0e1,
    'R_knee_angle': 1.0e1, 'R_knee_varus': 5.0e0, 'R_knee_rotation': 5.0e0,
    'R_ankle_angle': 1.0e1,
    #'R_ankle_abduction': 1.0e1, 'R_ankle_rotation': 1.0e1,
    'R_mtp_angle': 5.0e0,
    #'R_mtp_abduction': 1.0e1, 'R_mtp_rotation': 1.0e1,
    'L_hip_extension': 7.5e1, 'L_hip_abduction': 5.0e1, 'L_hip_rotation': 1.0e1,
    'L_knee_angle': 1.0e1, 'L_knee_varus': 5.0e0, 'L_knee_rotation': 5.0e0,
    'L_ankle_angle': 1.0e1,
    'L_ankle_abduction': 5.0e0,
    'L_mtp_angle': 5.0e0,
    }
# taskWeights = {
#     'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
#     'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1,
#     'R_hip_extension': 1, 'R_hip_abduction': 1, 'R_hip_rotation': 1,
#     'R_knee_angle': 1, 'R_knee_varus': 1, 'R_knee_rotation': 1,
#     'R_ankle_angle': 1,
#     #'R_ankle_abduction': 1, 'R_ankle_rotation': 1,
#     'R_mtp_angle': 1,
#     #'R_mtp_abduction': 1, 'R_mtp_rotation': 1,
#     'L_hip_extension': 1, 'L_hip_abduction': 1, 'L_hip_rotation': 1,
#     'L_knee_angle': 1, 'L_knee_varus': 1, 'L_knee_rotation': 1,
#     'L_ankle_angle': 1,
#     'L_ankle_abduction': 1,
#     'L_mtp_angle': 1,
#     }

#Create a dictionary of the optimal forces for actuators
actuatorVals = {'pelvis_tx': 10, 'pelvis_ty': 10, 'pelvis_tz': 10,
                'pelvis_tilt': 25, 'pelvis_list': 25, 'pelvis_rotation': 25,
                'R_hip_extension': 500, 'R_hip_abduction': 500, 'R_hip_rotation': 500,
                'R_knee_angle': 500, 'R_knee_varus': 500, 'R_knee_rotation': 500,
                'R_ankle_angle': 500,
                #'R_ankle_abduction': 500, 'R_ankle_rotation': 500,
                'R_mtp_angle': 500,
                #'R_mtp_abduction': 500, 'R_mtp_rotation': 500,
                'L_hip_extension': 500, 'L_hip_abduction': 500, 'L_hip_rotation': 500,
                'L_knee_angle': 500, 'L_knee_varus': 500, 'L_knee_rotation': 500,
                'L_ankle_angle': 500,
                'L_ankle_abduction': 500,
                'L_mtp_angle': 500,
                }

#Set time-frames in dataset - manually identified based on force data
gaitTimings = {'Walk': [], 'Run': [0.137, 0.573]}
    
# %% Settings for running simulations

#Create a Moco directory in the subjects folder
if runMoco:
    os.makedirs('moco', exist_ok = True)
    #Create trial specific directories as well
    for trial in trialList:
        os.makedirs(os.path.join('moco', trial), exist_ok = True)
        
# %% Run Moco

if runMoco:
    
    #Loop through run list
    for trial in trialList:
        
        # =================================================================
        #     Set-up for Moco
        # =================================================================
    
        #Change to Moco directory for ease of use with tools
        os.chdir(os.path.join('moco', trial))

        #Add in opensim logger
        osim.Logger.removeFileSink()
        osim.Logger.addFileSink('mocoLog.log')
    
        #Copy external load files across as there are issues with using these out of
        #directory with Moco tools
        shutil.copyfile(os.path.join('..','..','data',f'Ostrich_{trial}_GRF_ground.xml'),
                        f'Ostrich_{trial}_GRF_ground.xml')
        shutil.copyfile(os.path.join('..','..','data',f'Ostrich_{trial}_GRF_ground.mot'),
                        f'Ostrich_{trial}_GRF_ground.mot')
        shutil.copyfile(os.path.join('..','..','data',f'Ostrich_{trial}_coordinates.mot'),
                        f'Ostrich_{trial}_coordinates.mot')
    
        # =================================================================
        #     Run the Moco tracking
        # =================================================================
    
        #Create a generic tracking tool to manipulate for the 3 cycles
        mocoTrack = osim.MocoTrack()
        mocoTrack.setName('mocoOstrich')
        
        #Construct a ModelProcessor and set it on the tool.
        modelProcessor = osim.ModelProcessor(os.path.join('..','..','models','Rankin2016_Ostrich_noMuscles.osim'))
        modelProcessor.append(osim.ModOpAddExternalLoads(f'Ostrich_{trial}_GRF_ground.xml'))
        
        #Process model to edit
        mocoModel = modelProcessor.process()
        
        #Add in torque actuators that replicate the RRA actuators
        mocoModel = helper.addTorqueActuators(osimModel = mocoModel,
                                              optForces = actuatorVals)
        
        #Set model in tracking tool
        mocoTrack.setModel(osim.ModelProcessor(mocoModel))
        
        #Construct a table processor to append to the tracking tool for kinematics
        #The kinematics can't be filtered here with the operator as it messes with
        #time stamps in a funky way. This however has already been done in the 
        #conversion to state coordinates
        tableProcessor = osim.TableProcessor(f'Ostrich_{trial}_coordinates.mot')
        mocoTrack.setStatesReference(tableProcessor)
        
        #Create a dictionary to set kinematic bounds
        #Create this based on maximum and minimum values in the kinematic data
        #plus/minus some generic values
        
        #### TODO: are kinematic bounds needed?
    
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
        mocoTrack.set_mesh_interval(0.01) ### TODO: is this appropriate for this motion?
        
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
            if coordName in list(taskWeights.keys()):
                #Append state into weight set
                #Track the coordinate value
                stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/value',
                                                            taskWeights[coordName]))
                #Don't track the Coordinate speed
                stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/speed',
                                                            taskWeights[coordName] * speedsTrackingScale))
                
        #Add state weights to the tracking tool
        mocoTrack.set_states_weight_set(stateWeights)
    
        #Set the gait timings in tracking tool
        mocoTrack.set_initial_time(gaitTimings[trial][0])
        mocoTrack.set_final_time(gaitTimings[trial][1])
        
        #Initialise the Moco study
        study = mocoTrack.initialize()
        problem = study.updProblem()
        
        #Set the parameters for the regularization term on MocoTrack problem
        #(minimize squared excitations)
        effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
        effort.setWeight(0.001)
                
        #Lock time bounds to the IK data
        problem.setTimeBounds(gaitTimings[trial][0], gaitTimings[trial][1])
        
        #Set kinematic bounds using model ranges
        
        #### TODO: are bounds needed? Model bounds don't cover experimental data...
        
        #Get the solver
        solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
                
        #Set solver parameters
        solver.set_optim_max_iterations(1000)
        solver.set_optim_constraint_tolerance(1e-2)
        solver.set_optim_convergence_tolerance(1e-2)
        
        #Reset problem (required if changing to implicit mode)
        solver.resetProblem(problem)
        
        #Print to file
        study.printToXML(f'{trial}_setupMoco.omoco')

        #Solve!       
        solution = study.solve()
        
        # =============================================================
        #     Review Moco results
        # =============================================================
                
        #Check need to unseal and store outcome
        if solution.isSealed():
            solution.unseal()

        #Write solution to file
        solution.write(f'{trial}_mocoSolution.sto')
        
        #Remove tracked states file
        os.remove('mocoOstrich_tracked_states.sto')
        
        #Stop the logger
        osim.Logger.removeFileSink()
                        
        #Convert the solution to a states table and back to standard kinematic coordinates for ease of use
        
        #Write states table to file
        osim.STOFileAdapter().write(solution.exportToStatesTable(), f'{trial}_mocoStates.sto')
        
        #Convert states back to kinematic coordinates with helper function
        helper.statesToKinematics(statesFileName = f'{trial}_mocoStates.sto',
                                  outputFileName = f'{trial}_mocoKinematics.sto')
        
        #Navigate back to home directory for next run trial or subject
        os.chdir(homeDir)
        
# %% ----- end of 2_processData.py ----- %% #