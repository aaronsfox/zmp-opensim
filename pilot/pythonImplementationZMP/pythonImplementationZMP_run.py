# -*- coding: utf-8 -*-
"""

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Test implementation of ZMP predicted GRF using OpenSim interface. See 
    README.MD for any extra details.
    
    Relevant papers:
        
        Dijkstra & Gutierrez-Farewik (2015). J Biomech, 48: 3776-3781.
        
        Xiang et al. (2009). Int J Numer Meth ENg, 79: 667-695.
            > Calculations for transforming and calculating COP


"""

# %% Import packages

import opensim as osim
import os
import matplotlib.pyplot as plt
import numpy as np

# %% Define helpful functions

# %% Add torque actuators

def addTorqueActuators(osimModel = None,
                       optForces = None,
					   controlLimits = None):
    
    """
    
    Convenience function for adding series of torque actuators to model
    
    Input:    osimModel - OpenSim model object for use
              optForces - dict of coordinates and their associated optimal forces to add
			  controlLimits - dict of coordinates and their associated max/min control limits to add
              
    Output:   osimModel - updated torque driven model
                  
    """
    
    #Check inputs
    if osimModel is None or optForces is None:
            raise ValueError('All inputs for this function are required!')

    #Intialise model system
    osimModel.initSystem()
    
    #Get coordinate list
    coordinatesList = list(optForces.keys())
    
    #Get coordinate set
    coordSet = osimModel.getCoordinateSet()
    
    #Loop through coordinates and add actuators
    for coordinate in coordinatesList:
        #Create actuator
        actu = osim.CoordinateActuator()
        #Set name
        actu.setName(f'{coordinate}_actuator')
        #Set coordinate
        actu.setCoordinate(coordSet.get(coordinate))
        #Set optimal force
        actu.setOptimalForce(optForces[coordinate])
        #Set min and max control
        actu.setMinControl(controlLimits[coordinate]*-1)
        actu.setMaxControl(controlLimits[coordinate]*1)
        #Append to model force set
        osimModel.updForceSet().cloneAndAppend(actu)
    
    #Finalise model connections
    osimModel.finalizeConnections()
    
    #Return model
    return osimModel

# %% Set-up

#Import model geometry
geomDir = os.path.join('C:', os.sep, 'OpenSim 4.4', 'Geometry')
osim.ModelVisualizer.addDirToGeometrySearchPaths(geomDir)

#Set matplotlib parameters
from matplotlib import rcParams
# rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['font.weight'] = 'bold'
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 16
rcParams['axes.linewidth'] = 1.5
rcParams['axes.labelweight'] = 'bold'
rcParams['legend.fontsize'] = 10
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.width'] = 1.5
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
rcParams['legend.framealpha'] = 0.0
rcParams['savefig.dpi'] = 300
rcParams['savefig.format'] = 'pdf'

#Create a dictionary of the optimal forces for actuators used in Hamner & Delp
actuatorForces = {'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
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
actuatorLimits = {'pelvis_tx': 10000, 'pelvis_ty': 10000, 'pelvis_tz': 10000,
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

#Set ground contact force calculation
forceThreshold = 50

#Set the right and left bodies to consider the applied forces to
rightBodyName = 'calcn_r'
leftBodyName = 'calcn_l'

# %% First, use the Inverse Dynamics tool to calculate body forces without GRFs

"""

TODO:
    
    > Only calculate when foot contact is detected
        >> This is done based on foot threshold at the moment, but could be done
           by something to do with foot position...

"""

##### ----- RUN INVERSE DYNAMICS TOOL ----- #####

#Create inverse dynamics tool
idTool = osim.InverseDynamicsTool(os.path.join('utilities','blank_id_setup.xml'))

#Set model
idTool.setModelFileName(os.path.join('data_run','subject01_adjusted_scaled.osim'))

#Set coordinates
idTool.setCoordinatesFileName(os.path.join('data_run','subject01_run5_cycle1_mocoKinematics.sto'))

#Set time ranges
idTool.setStartTime(osim.MocoTrajectory(os.path.join('data_run','subject01_run5_cycle1_mocoSolution.sto')).getInitialTime())
idTool.setEndTime(osim.MocoTrajectory(os.path.join('data_run','subject01_run5_cycle1_mocoSolution.sto')).getFinalTime())

#Set output file names
idTool.setResultsDir('outputs_run')
idTool.setOutputGenForceFileName('id_genForces.sto')

#Print and reload to avoid filepath issues
idTool.printToXML('setupID.xml')
osim.InverseDynamicsTool('setupID.xml').run()

#Rename body forces file
os.replace(os.path.join('outputs_run','body_forces_at_joints.sto'),
           os.path.join('outputs_run','id_bodyForces.sto'))

##### ----- PREPARE FILES FOR STORING ZMP GRFs ----- #####

#Read in and initialise model

#Load model for use
modelProcessor = osim.ModelProcessor(os.path.join('data_run','subject01_adjusted_scaled.osim'))

#Remove muscles to be consistent with torque driven solution
modelProcessor.append(osim.ModOpRemoveMuscles())

#Process model
osimModel = modelProcessor.process()

#Test adding actuators given that these are disabled later
osimModel = addTorqueActuators(osimModel = osimModel,
                               optForces = actuatorForces,
                               controlLimits = actuatorLimits)

#Unlock any coordinates to avoid issues with setting values
for cInd in range(osimModel.updCoordinateSet().getSize()):
    if osimModel.updCoordinateSet().get(cInd).get_locked():
        osimModel.updCoordinateSet().get(cInd).set_locked(False)
        
#Finalize connections
osimModel.finalizeConnections()

#Initialise model
osimModel.initSystem()

#Get model frame for transforming back to ground
ground = osimModel.getGround()
pelvisFrame = osimModel.getJointSet().get('ground_pelvis').getChildFrame()

#Get model states from motion
mocoTraj = osim.MocoTrajectory(os.path.join('data_run','subject01_run5_cycle1_mocoSolution.sto'))
statesTraj = mocoTraj.exportToStatesTrajectory(osimModel)

#Read in body forces data to transform to ground and calculate COP
idBodyForces = osim.TimeSeriesTable(os.path.join('outputs_run','id_bodyForces.sto'))

#Get number of times in ID results
nt = idBodyForces.getNumRows()

#Create storage file to append ZMP GRFs to
zmpResults = osim.Storage(nt)

#Create columns for ZMP force and COP calculations
#These are partitioned to left and right sides as per usual external load data
#(i.e. 6 x columns per side)
zmpLabels = osim.ArrayStr('time', 6*2+1)
zmpLabels.set(1, 'R_ground_force_vx')
zmpLabels.set(2, 'R_ground_force_vy')
zmpLabels.set(3, 'R_ground_force_vz')
zmpLabels.set(4, 'R_ground_force_px')
zmpLabels.set(5, 'R_ground_force_py')
zmpLabels.set(6, 'R_ground_force_pz')
zmpLabels.set(7, 'L_ground_force_vx')
zmpLabels.set(8, 'L_ground_force_vy')
zmpLabels.set(9, 'L_ground_force_vz')
zmpLabels.set(10, 'L_ground_force_px')
zmpLabels.set(11, 'L_ground_force_py')
zmpLabels.set(12, 'L_ground_force_pz')

##### ----- PREDICT GRFs AT EACH STATE USING ZMP ----- #####

#Loop through each time-point in states
for sInd in range(nt):
    
    #Set vector to store results for current state in (6 values per limb)
    #Set all values to zero
    zmpVec = osim.Vector(6*2, 0)
    
    #Initialise state
    s = statesTraj[sInd]
    osimModel.realizeDynamics(s)
    
    #Get the corresponding index in the ID body forces related to states time
    sTime = s.getTime()
    idInd = idBodyForces.getNearestRowIndexForTime(sTime)
    
    #Get the forces and torques from the pelvis body forces
    pelvisF = osim.Vec3(idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FX')[idInd],
                        idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FY')[idInd],
                        idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FZ')[idInd])
    pelvisM = osim.Vec3(idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_MX')[idInd],
                        idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_MY')[idInd],
                        idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_MZ')[idInd])
    
    #Check if vertical force greater than force threshold for calculations
    if pelvisF.get(1) > forceThreshold:
    
        #Get the transform between pelvis to ground at the current state
        groundTransform = pelvisFrame.getTransformInGround(s)
        
        #Get the translation component of this transform to get the pelvis position vector
        rp = groundTransform.T()
        # rr = groundTransform.R()
        
        #Take cross product of pelvis position and force vector to get moment at origin
        groundM = osim.Vec3((rp.get(1) * pelvisF.get(2)) - (rp.get(2) * pelvisF.get(1)),
                            -((rp.get(0) * pelvisF.get(2)) - (rp.get(2) * pelvisF.get(0))),
                            (rp.get(0) * pelvisF.get(1)) - (rp.get(1) * pelvisF.get(0)))
            
        # #Calculate force at ground origin (just equals force at pelvis given already in ground system)
        groundF = osim.Vec3(pelvisF.get(0), pelvisF.get(1), pelvisF.get(2)) #given these are in the same frame
            
        # #Review difference
        # for v in (range(3)):
        #     print(f'F{"XYZ"[v]} in {pelvisFrame.getName()}: {pelvisF.get(v)} vs. F{"XYZ"[v]} in Ground: {groundF.get(v)}')
        #     print(f'M{"XYZ"[v]} in {pelvisFrame.getName()}: {pelvisM.get(v)} vs. M{"XYZ"[v]} in Ground: {groundM.get(v)}')
        
        ##### TODO: calculations performed correctly?
        
        #Calculate X & Z cZMP, noting that yZMP is set as 0
        #Formulas come from Xiang et al.
        zmpCOP = osim.Vec3(groundM.get(2) / groundF.get(1),
                           0,
                           -groundM.get(0) / groundF.get(1))
        
        #Calculate the resultant active moment at ZMP along the y-axis
        #### TODO: this still seems wrong
        # myZMP = groundM.get(1) + (groundF.get(0) * zmpCOP.get(2)) - (groundF.get(2) * zmpCOP.get(0))
        
        #Get the position of the right and left body origins in the ground
        #This will determine which is closer to the predicted COP
        rightBodyPos = osimModel.updBodySet().get(rightBodyName).findStationLocationInGround(s, osim.Vec3(0,0,0))
        leftBodyPos = osimModel.updBodySet().get(leftBodyName).findStationLocationInGround(s, osim.Vec3(0,0,0))
        
        #Calculate distances from body to predicted ZMP
        rightBodyDist = (rightBodyPos.get(0) - zmpCOP.get(0))**2 + (rightBodyPos.get(1) - zmpCOP.get(1))**2 + (rightBodyPos.get(2) - zmpCOP.get(2))**2
        leftBodyDist = (leftBodyPos.get(0) - zmpCOP.get(0))**2 + (leftBodyPos.get(1) - zmpCOP.get(1))**2 + (leftBodyPos.get(2) - zmpCOP.get(2))**2
        
        #Find smallest distance and allocate to appropriate part of vector
        #### TODO: edge cases where identical?
        if rightBodyDist < leftBodyDist:
            #Right side forces
            zmpVec.set(0, groundF.get(0))
            zmpVec.set(1, groundF.get(1))
            zmpVec.set(2, groundF.get(2))
            #Right side COP
            zmpVec.set(3, zmpCOP.get(0))
            zmpVec.set(4, zmpCOP.get(1))
            zmpVec.set(5, zmpCOP.get(2))
        elif leftBodyDist < rightBodyDist:
            #Left side forces
            zmpVec.set(6, groundF.get(0))
            zmpVec.set(7, groundF.get(1))
            zmpVec.set(8, groundF.get(2))
            #Left side COP
            zmpVec.set(9, zmpCOP.get(0))
            zmpVec.set(10, zmpCOP.get(1))
            zmpVec.set(11, zmpCOP.get(2))
        
    #There is no need to change zmp vector as they start as zeros
    
    #Create state vector with time and values
    #Convert to state vector
    zmpStateVec = osim.StateVector(s.getTime(), zmpVec)
    
    #Append current vector to body forces results
    zmpResults.append(zmpStateVec)
        
#Set column labels in ZMP results storage
zmpResults.setColumnLabels(zmpLabels)

#Set name in ZMP results storage
zmpResults.setName('ZMP Predicted Ground Reaction Forces')

#Print ZMP results to file
osim.Storage().printResult(zmpResults, 'id_zmpForces', 'outputs_run', -1, '.sto')

"""

NOTES:
    
    > Forces are OK
    > COP is OK in the middle part of stance - but trails off badly at lower forces
        >> X is pretty good - Z is where large issues are, which relates to ground-X moments
        >> Possible issues with static ID approach - or something else?
        >> There's possibly a better approach to estimating COP?    

"""

# %% Estimate GRF using ZMP from kinematic states in Inverse Dynamics solver

"""

The idea from this approach is to use the coordinates from the states in the 
inverse dynamics solver to calculate body forces. Most of the code below comes
from the InverseDynamicsTool.cpp opensim-core code.

TODO:
    
    > Need a more creative way to recognise pelvis states in general forces
        >> Can probably come with tree order working properly
    
    > Do body forces need to be actually calculated, given that the ID solver
      gives the forces and moments of the free pelvis coordinates? When running
      ID through the tool - these gen and body forces are the same...
      
    > Flexible ordering of states vs. multibody tree in C++ implementation
    
    > Only calculate when foot is in contact with ground
        >> Currently does this based on a force threshold, which isn't applicable
           to later applications

"""

##### ----- SET-UP MODEL AND SOLVER ----- #####

#Load model for use
modelProcessor = osim.ModelProcessor(os.path.join('data_run','subject01_adjusted_scaled.osim'))

#Remove muscles to be consistent with torque driven solution
modelProcessor.append(osim.ModOpRemoveMuscles())

#Process model
osimModel = modelProcessor.process()

#Test adding actuators given that these are disabled later
osimModel = addTorqueActuators(osimModel = osimModel,
                               optForces = actuatorForces,
                               controlLimits = actuatorLimits)

#Unlock any coordinates to avoid issues with setting values
for cInd in range(osimModel.updCoordinateSet().getSize()):
    if osimModel.updCoordinateSet().get(cInd).get_locked():
        osimModel.updCoordinateSet().get(cInd).set_locked(False)
        
#Finalize connections
osimModel.finalizeConnections()

#Initialise model
osimModel.initSystem()

#Get model working copy
sWorkingCopy = osimModel.getWorkingState()

#Initialise an inverse dynamics solver with model
ivdSolver = osim.InverseDynamicsSolver(osimModel)

##### ----- SET-UP STATES AND STORAGE FILES ----- #####

#Get model states from motion
mocoTraj = osim.MocoTrajectory(os.path.join('data_run','subject01_run5_cycle1_mocoSolution.sto'))
statesTraj = mocoTraj.exportToStatesTrajectory(osimModel)

#In lieu of udot not working later, create an accelerations table
mocoTraj.generateAccelerationsFromSpeeds()
uDotTable = mocoTraj.exportToAccelerationsTable()

#Define number of times from states trajectory
nt = mocoTraj.getNumTimes()

#Define the number of coordinates in the model
nq = osimModel.updCoordinateSet().getSize()

#We can get the order of coordinates in tree order from the model
# coords = osimModel.getCoordinatesInMultibodyTreeOrder()
#### this doesn't seem to help in Python - potentially not wrapped?

# #Given the above doesn't work - do it manually
# #States Q and U order is in:
# orderStatesQU = ['pelvis_tilt', 'pelvis_list', 'pelvis_rotation',
#                  'pelvis_tx', 'pelvis_ty', 'pelvis_tz',
#                  'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r',
#                  'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l',
#                  'lumbar_extension', 'lumbar_bending', 'lumbar_rotation',
#                  'knee_angle_r', 'knee_angle_l',
#                  'arm_flex_r', 'arm_add_r', 'arm_rot_r',
#                  'arm_flex_l', 'arm_add_l', 'arm_rot_l',
#                  'ankle_angle_r', 'ankle_angle_l',
#                  'elbow_flex_r', 'elbow_flex_l',
#                  'subtalar_angle_r', 'subtalar_angle_l',
#                  'pro_sup_r', 'pro_sup_l',
#                  'mtp_angle_r', 'mtp_angle_l',
#                  'wrist_flex_r', 'wrist_dev_r',
#                  'wrist_flex_l', 'wrist_dev_l']

# #Best guess for multibody tree order is based on model jointset
# multiBodyOrder = ['pelvis_tilt', 'pelvis_list', 'pelvis_rotation',
#                   'pelvis_tx', 'pelvis_ty', 'pelvis_tz',
#                   'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r',
#                   'knee_angle_r',
#                   'ankle_angle_r', 'subtalar_angle_r', 'mtp_angle_r',
#                   'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l',
#                   'knee_angle_l',
#                   'ankle_angle_l', 'subtalar_angle_l', 'mtp_angle_l',
#                   'lumbar_extension', 'lumbar_bending', 'lumbar_rotation',
#                   'arm_flex_r', 'arm_add_r', 'arm_rot_r',
#                   'elbow_flex_r', 'pro_sup_r',
#                   'wrist_flex_r', 'wrist_dev_r',                  
#                   'arm_flex_l', 'arm_add_l', 'arm_rot_l',
#                   'elbow_flex_l', 'pro_sup_l',
#                   'wrist_flex_l', 'wrist_dev_l']

#### TODO: unsure exactly what needs to be mapped in a different order
#### NOTE: didn't reorder and somehow got the correct values somehow...?
    #### NOTE: this came from the coordinate functions being evaluated though
    #### Nontetheless, suggests that no re-ordering is required, it's just that udot is wrong
#Start by assuming udot comes out in multi-body order and needs to be re-arranged to states
# multiBodyReorder = []
# for statesBodyName in orderStatesQU:
#     #Find states index in proposed multibody order
#     multiBodyReorder.append(multiBodyOrder.index(statesBodyName))
    
# #Above did fix the general force trajectory actually creating non-zero vals --- but they looked wrong
# #Test out reversing the multibody order approach (i.e. multi-body to states)
# multiBodyReorder = []
# for multiBodyName in multiBodyOrder:
#     #Find states index in proposed multibody order
#     multiBodyReorder.append(orderStatesQU.index(multiBodyName))
    
# #Given the above ideas aren't working, test out coordinate functions approach to
# #see what inverse dynamics solver spits out with those
# coordinateValues = osim.Storage(os.path.join('data_run','subject01_run5_cycle1_mocoKinematics.sto'))
# coordFunctions = osim.FunctionSet()
# coordSplines = osim.GCVSplineSet(5, coordinateValues)
# for iq in range(nq):
#     coordFunctions.insert(iq, coordSplines.get(orderStatesQU[iq]))

#Create storage objects for gen and body forces to append data to
genForcesResults = osim.Storage(nt)
bodyForcesResults = osim.Storage(nt)

#Create coordinate labels for gen forces
genForceLabels = osim.ArrayStr('time', nq+1)
for iq in range(nq):
    genForceLabels.set(iq+1, osimModel.updCoordinateSet().get(iq).getName())

#Build the list of joints for computing and reporting equivalent body forces

#Set bodies to report forces for
jointsForReportingBodyForces = osim.ArrayStr()
jointsForReportingBodyForces.append('ground_pelvis')

#Create jointset for calculating equivalent body forces
jointsForEquivalentBodyForces = osim.JointSet()

#Get model joints
modelJoints = osimModel.getJointSet()

#Loop through reporting joint names, clone and append the corresponding model joint
for ij in range(jointsForReportingBodyForces.getSize()):
    ### NOTE: there is an option for "ALL" in C++ code
    k = modelJoints.getIndex(jointsForReportingBodyForces.get(ij))
    if k >= 0:
        jointsForEquivalentBodyForces.adoptAndAppend(modelJoints.get(k))
    else:
        print('Tool could not find joint named %s for reporting body forces...' % jointsForReportingBodyForces.get(ij))

#Get number of joints for reporting body forces
nj = jointsForEquivalentBodyForces.getSize()

#Create column labels for body forces storage file
#Starts with time plus the space for six components per joint to report on
bodyForceLabels = osim.ArrayStr('time', 6*nj+1)
stringXYZ = 'XYZ'
for ij in range(nj):
    #Start with the joint body label
    joint_body_label = jointsForEquivalentBodyForces.get(ij).getName()+'_'
    joint_body_label += jointsForEquivalentBodyForces.get(ij).getChildFrame().getName()
    #Set the XYZ components in array string for current joint/body
    for k in range(3):
        bodyForceLabels.set(6*ij+k+1, joint_body_label + '_F' + stringXYZ[k])
        bodyForceLabels.set(6*ij+k+3+1, joint_body_label + '_M' + stringXYZ[k])

#Create storage file to append ZMP GRFs to
zmpResults = osim.Storage(nt)

#Create columns for ZMP force and COP calculations
#These are partitioned to left and right sides as per usual external load data
#(i.e. 6 x columns per side)
zmpLabels = osim.ArrayStr('time', 6*2+1)
zmpLabels.set(1, 'R_ground_force_vx')
zmpLabels.set(2, 'R_ground_force_vy')
zmpLabels.set(3, 'R_ground_force_vz')
zmpLabels.set(4, 'R_ground_force_px')
zmpLabels.set(5, 'R_ground_force_py')
zmpLabels.set(6, 'R_ground_force_pz')
zmpLabels.set(7, 'L_ground_force_vx')
zmpLabels.set(8, 'L_ground_force_vy')
zmpLabels.set(9, 'L_ground_force_vz')
zmpLabels.set(10, 'L_ground_force_px')
zmpLabels.set(11, 'L_ground_force_py')
zmpLabels.set(12, 'L_ground_force_pz')
        
##### ----- LOOP THROUGH STATES AND CALCULATE BODY FORCES ----- #####

#Loop through times in states trajectory
for sInd in range(nt):

    #Get initial model state to test
    s = statesTraj[sInd]
    
    #Set modelling options for Actuators to be over-ridden
    for i in range(osimModel.updForceSet().getSize()):
        act = osim.ScalarActuator.safeDownCast(osimModel.updForceSet().get(i))
        act.overrideActuation(sWorkingCopy, True)
        
    #Set current states in working copy
    sWorkingCopy.setTime(s.getTime())
    sWorkingCopy.setQ(s.getQ())
    sWorkingCopy.setU(s.getU())
    # udot = sWorkingCopy.getUDot()

    #Realise model to acceleration stage
    # osimModel.realizeAcceleration(sWorkingCopy)
    osimModel.realizeDynamics(sWorkingCopy)
    
    # #Check if Q's are updated
    # for qi in range(sWorkingCopy.getQ().size()):
    #     print(sWorkingCopy.getQ().get(qi))
    
    #Get accelerations
    ### NOTE: only seems achievable via Moco trajectory data
    # udot = sWorkingCopy.getUDot()
    udot = osim.Vector(nq, 0)
    for iq in range(nq):
        udot.set(iq,
                 uDotTable.getRowAtIndex(uDotTable.getNearestRowIndexForTime(s.getTime())).getElt(0,iq)
                 )
    # for u in range(udot.size()):
    #     print(udot.get(u))
    
    # #Reorder udot to match states based on earlier created map
    # udotOrdered = osim.Vector(udot.size(), 0)
    # for u in range(udot.size()):
    #     udotOrdered.set(u, udot.get(multiBodyReorder[u]))
    
    #### TODO: important that udot order matches the states --- need checks in place for this...
    
    #Solve ID
    genForceTraj = ivdSolver.solve(sWorkingCopy, udot)
    # genForceTraj = ivdSolver.solve(sWorkingCopy, coordFunctions, s.getTime())
    # for u in range(genForceTraj.size()):
    #     print(genForceTraj.get(u))
    
    #The above results contain the generalised coordinate forces to generate the accelerations
    #based on the current state. Note that these aren't necessarily in the order of the coordinate
    #set, but rather the multibody tree order
    
    #Next we want to calculate the equivalent body forces related to these
    
    #Create a state vector to store the forces
    #Note that this allocates the time from the current state and 6 slots for the 
    #force and moment components for each body (with a zero allocated for each)
    # bodyForcesVec = osim.StateVector(s.getTime(), osim.Vector(6*nj, 0))
    genVec = osim.Vector(nq, 0)
    forcesVec = osim.Vector(6*nj, 0)
    
    #### TODO: mapping if q != u in index (i.e. tree vs. model)...
    
    #Calculate equivalent body force at joint for those listed
    for j in range(nj):
        equivalentBodyForceAtJoint = jointsForEquivalentBodyForces.get(j).calcEquivalentSpatialForce(s, genForceTraj)
        #Extract the Vec 3 components
        for k in range(3):
            #Body force components
            forcesVec.set(6*j+k, equivalentBodyForceAtJoint.get(1)[k])
            #Body torque components
            forcesVec.set(6*j+k+3, equivalentBodyForceAtJoint.get(0)[k])
            
    #Convert to state vector
    genForcesVec = osim.StateVector(s.getTime(), genForceTraj)
    bodyForcesVec = osim.StateVector(s.getTime(), forcesVec)
    
    #Append current vector to body forces results
    genForcesResults.append(genForcesVec)
    bodyForcesResults.append(bodyForcesVec)
    
    #Calculate the ZMP estimated GRFs and COP
    
    #Set vector to store results for current state in (6 values per limb)
    #Set all values to zero
    zmpVec = osim.Vector(6*2, 0)
    
    #Get the body forces acting at the pelvis
    #### NOTE: this is manually indexed --- which isn't great for use across different models...
    pelvisF = osim.Vec3(forcesVec.get(0), forcesVec.get(1), forcesVec.get(2))
    
    #Check if vertical force greater than force threshold for calculations
    if pelvisF.get(1) > forceThreshold:
        
        #Get the position of the pelvis from Q states
        #### NOTE: this is manually indexed --- which isn't great for use across different models...
        rp = osim.Vec3(s.getQ().get(3), s.getQ().get(4), s.getQ().get(5))
        
        #Take cross product of pelvis position and force vector to get moment at origin
        groundM = osim.Vec3((rp.get(1) * pelvisF.get(2)) - (rp.get(2) * pelvisF.get(1)),
                            -((rp.get(0) * pelvisF.get(2)) - (rp.get(2) * pelvisF.get(0))),
                            (rp.get(0) * pelvisF.get(1)) - (rp.get(1) * pelvisF.get(0)))
            
        #Calculate force at ground origin (just equals force at pelvis given already in ground system)
        groundF = osim.Vec3(pelvisF.get(0), pelvisF.get(1), pelvisF.get(2)) #given these are in the same frame
        
        #Calculate X & Z cZMP, noting that yZMP is set as 0
        #Formulas come from Xiang et al.
        zmpCOP = osim.Vec3(groundM.get(2) / groundF.get(1),
                           0,
                           -groundM.get(0) / groundF.get(1))
        
        #Calculate the resultant active moment at ZMP along the y-axis
        #### TODO: this still seems wrong
        # myZMP = groundM.get(1) + (groundF.get(0) * zmpCOP.get(2)) - (groundF.get(2) * zmpCOP.get(0))
        
        #Get the position of the right and left body origins in the ground
        #This will determine which is closer to the predicted COP
        rightBodyPos = osimModel.updBodySet().get(rightBodyName).findStationLocationInGround(s, osim.Vec3(0,0,0))
        leftBodyPos = osimModel.updBodySet().get(leftBodyName).findStationLocationInGround(s, osim.Vec3(0,0,0))
        
        #Calculate distances from body to predicted ZMP
        rightBodyDist = (rightBodyPos.get(0) - zmpCOP.get(0))**2 + (rightBodyPos.get(1) - zmpCOP.get(1))**2 + (rightBodyPos.get(2) - zmpCOP.get(2))**2
        leftBodyDist = (leftBodyPos.get(0) - zmpCOP.get(0))**2 + (leftBodyPos.get(1) - zmpCOP.get(1))**2 + (leftBodyPos.get(2) - zmpCOP.get(2))**2
        
        #Find smallest distance and allocate to appropriate part of vector
        #### TODO: edge cases where identical?
        if rightBodyDist < leftBodyDist:
            #Right side forces
            zmpVec.set(0, groundF.get(0))
            zmpVec.set(1, groundF.get(1))
            zmpVec.set(2, groundF.get(2))
            #Right side COP
            zmpVec.set(3, zmpCOP.get(0))
            zmpVec.set(4, zmpCOP.get(1))
            zmpVec.set(5, zmpCOP.get(2))
        elif leftBodyDist < rightBodyDist:
            #Left side forces
            zmpVec.set(6, groundF.get(0))
            zmpVec.set(7, groundF.get(1))
            zmpVec.set(8, groundF.get(2))
            #Left side COP
            zmpVec.set(9, zmpCOP.get(0))
            zmpVec.set(10, zmpCOP.get(1))
            zmpVec.set(11, zmpCOP.get(2))
        
    #There is no need to change zmp vector as they start as zeros
    
    #Create state vector with time and values
    #Convert to state vector
    zmpStateVec = osim.StateVector(s.getTime(), zmpVec)
    
    #Append current vector to ZMP results
    zmpResults.append(zmpStateVec)
        

##### ----- WRITE RESULTS FORCES TO FILE ----- #####

##### NOTE: currently all zeros in both files - probably issue with calculation...

#Set column labels in storage objects
genForcesResults.setColumnLabels(genForceLabels)
bodyForcesResults.setColumnLabels(bodyForceLabels)
zmpResults.setColumnLabels(zmpLabels)

#Set name in storage
genForcesResults.setName('Inverse Dynamics Generalized Forces')
bodyForcesResults.setName('Inverse Dynamics Body Forces at Specified Joints')
zmpResults.setName('ZMP Predicted Ground Reaction Forces')

#Write to file
osim.Storage().printResult(genForcesResults, 'manual_genForces', 'outputs_run', -1, '.sto')
osim.Storage().printResult(bodyForcesResults, 'manual_bodyForces', 'outputs_run', -1, '.sto')
osim.Storage().printResult(zmpResults, 'manual_zmpForces', 'outputs_run', -1, '.sto')

# %% Compare outputs from different approaches

"""

This just focuses on the manual ZMP calculations as they are near identical to
the ones generated via the ID tool

"""

##### ----- COMPARE GRFs ----- #####

#Create figure
fig,ax = plt.subplots(nrows = 1, ncols = 3, figsize = (10,3))
plt.subplots_adjust(left = 0.07, right = 0.975, top = 0.9, bottom = 0.15,
                    wspace = 0.3)

#Read in experimental and ZMP generated GRFs
expGRF = osim.TimeSeriesTable(os.path.join('data_run','Run_5_grf.mot'))
zmpGRF = osim.TimeSeriesTable(os.path.join('outputs_run','manual_zmpForces.sto'))

#Set start and end times for comparison
startTime = zmpGRF.getIndependentColumn()[0]
endTime = zmpGRF.getIndependentColumn()[-1]

#Get index in full experimental trial that matches the gait cycle
expStartInd = expGRF.getNearestRowIndexForTime(startTime)
expEndInd = expGRF.getNearestRowIndexForTime(endTime)

#Sum and plot XYZ ground reaction forces
for axesLabel in 'xyz':
    
    #Extract experimental GRFs to numpy array and sum left and right sides
    expData = expGRF.getDependentColumn(f'R_ground_force_v{axesLabel}').to_numpy()[expStartInd:expEndInd+1] + \
        expGRF.getDependentColumn(f'L_ground_force_v{axesLabel}').to_numpy()[expStartInd:expEndInd+1]
    expTime = expGRF.getIndependentColumn()[expStartInd:expEndInd+1]
    
    #Extract ZMP GRFs to numpy array and sum left and right sides
    zmpData = zmpGRF.getDependentColumn(f'R_ground_force_v{axesLabel}').to_numpy() + \
        zmpGRF.getDependentColumn(f'L_ground_force_v{axesLabel}').to_numpy()
    zmpTime = zmpGRF.getIndependentColumn()
    
    #Plot current axes data
    ax['xyz'.index(axesLabel)].plot(expTime, expData, ls = '-', lw = 1.5, c = 'black', label = 'Experimental')
    ax['xyz'.index(axesLabel)].plot(zmpTime, zmpData, ls = '--', lw = 1.5, c = 'red', label = 'ZMP Predicted')
    
    #Add labels and titles
    ax['xyz'.index(axesLabel)].set_xlabel('Time (s)', fontsize = 12, fontweight = 'bold')
    ax['xyz'.index(axesLabel)].set_ylabel('Force (N)', fontsize = 12, fontweight = 'bold')
    ax['xyz'.index(axesLabel)].set_title(f'{axesLabel.capitalize()}-Axis GRFs', fontsize = 14, fontweight = 'bold')
    
    #Add zero line
    ax['xyz'.index(axesLabel)].axhline(y = 0, ls = ':', lw = 1.0, c = 'dimgrey')
    
    #Add legend to last axis
    if axesLabel == 'z':
        ax['xyz'.index(axesLabel)].legend()
    
#Save figure
fig.savefig(os.path.join('outputs_run','grfComparison.png'), format = 'png', dpi = 300)

#Close figure
plt.close()    
    
##### ----- COMPARE COP ----- #####

#Set constant treadmill speed to update COP visualisation in X axis
treadmillSpeed = 5.0

#Create figure
fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (5,5))

#Extract experimental COP in X and Z axis to numpy array 
expCopX_r = expGRF.getDependentColumn('R_ground_force_px').to_numpy()[expStartInd:expEndInd+1]
expCopZ_r = expGRF.getDependentColumn('R_ground_force_pz').to_numpy()[expStartInd:expEndInd+1]
expCopX_l = expGRF.getDependentColumn('L_ground_force_px').to_numpy()[expStartInd:expEndInd+1]
expCopZ_l = expGRF.getDependentColumn('L_ground_force_pz').to_numpy()[expStartInd:expEndInd+1]
expTime = expGRF.getIndependentColumn()[expStartInd:expEndInd+1]

#Extract ZMP predicted COP in X and Z axis to numpy array
zmpCopX_r = expGRF.getDependentColumn('R_ground_force_px').to_numpy()
zmpCopZ_r = expGRF.getDependentColumn('R_ground_force_pz').to_numpy()
zmpCopX_l = expGRF.getDependentColumn('L_ground_force_px').to_numpy()
zmpCopZ_l = expGRF.getDependentColumn('L_ground_force_pz').to_numpy()
zmpTime = expGRF.getIndependentColumn()

#Replace zeros with NaN's as this means no contact is present
expCopX_r[expCopX_r == 0] = np.nan
expCopZ_r[expCopZ_r == 0] = np.nan
expCopX_l[expCopX_l == 0] = np.nan
expCopZ_l[expCopZ_l == 0] = np.nan
zmpCopX_r[zmpCopX_r == 0] = np.nan
zmpCopZ_r[zmpCopZ_r == 0] = np.nan
zmpCopX_l[zmpCopX_l == 0] = np.nan
zmpCopZ_l[zmpCopZ_l == 0] = np.nan

#Improve visualisation by adjusting X-COP for treadmill speed
#Experimental
for ii in range(1,len(expTime)):
    #Calculate how far treadmill will have moved in time-step based on speed and time
    timeDiff = expTime[ii] - expTime[ii-1]
    #Add treadmill displacement to COP
    expCopX_r[ii] = expCopX_r[ii] + (treadmillSpeed * timeDiff)
    expCopX_l[ii] = expCopX_l[ii] + (treadmillSpeed * timeDiff)
#ZMP predicted
for ii in range(1,len(zmpTime)):
    #Calculate how far treadmill will have moved in time-step based on speed and time
    timeDiff = zmpTime[ii] - zmpTime[ii-1]
    #Add treadmill displacement to COP
    zmpCopX_r[ii] = zmpCopX_r[ii] + (treadmillSpeed * timeDiff)
    zmpCopX_l[ii] = zmpCopX_l[ii] + (treadmillSpeed * timeDiff)

#Plot data
ax[0].plot(expCopZ_l, expCopX_l, ls = '-', lw = 1.0, c = 'black', label = '_Experimental')
ax[1].plot(expCopZ_r, expCopX_r, ls = '-', lw = 1.0, c = 'black', label = 'Experimental')
ax[0].plot(zmpCopZ_l, zmpCopX_l, ls = '--', lw = 1.0, c = 'red', label = '_ZMP Predicted')
ax[1].plot(zmpCopZ_r, zmpCopX_r, ls = '--', lw = 1.0, c = 'red', label = 'ZMP Predicted')

#Set equal axis
ax[0].axis('equal')
ax[1].axis('equal')

#Add labels and titles
ax[0].set_xlabel('X-Position (m)', fontsize = 12, fontweight = 'bold')
ax[1].set_xlabel('X-Position (m)', fontsize = 12, fontweight = 'bold')
ax[0].set_ylabel('Z-Position (m)', fontsize = 12, fontweight = 'bold')
fig.suptitle('Centre of Pressure', fontsize = 14, fontweight = 'bold')

#Add legend
ax[1].legend()

#Tight layout
plt.tight_layout()

#Save figure
fig.savefig(os.path.join('outputs_run','copComparison.png'), format = 'png', dpi = 300)

#Close figure
plt.close()


# %% ----- End of pythonImplementationZMP_run.py ----- %% #