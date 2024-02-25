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
        
        Xiang et al. (2009). Int J NUmer Meth ENg, 79: 667-695.
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

# #Create a dictionary of the optimal forces for actuators used in Hamner & Delp
# actuatorForces = {'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
#                   'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1,
#                   'hip_flexion_r': 1000, 'hip_adduction_r': 1000, 'hip_rotation_r': 1000,
#                   'knee_angle_r': 1000, 'ankle_angle_r': 1000,
#                   'hip_flexion_l': 1000, 'hip_adduction_l': 1000, 'hip_rotation_l': 1000,
#                   'knee_angle_l': 1000, 'ankle_angle_l': 1000,
#                   'lumbar_extension': 1000, 'lumbar_bending': 1000, 'lumbar_rotation': 1000,
#                   'arm_flex_r': 500, 'arm_add_r': 500, 'arm_rot_r': 500,
#                   'elbow_flex_r': 500, 'pro_sup_r': 500,
#                   'arm_flex_l': 500, 'arm_add_l': 500, 'arm_rot_l': 500,
#                   'elbow_flex_l': 500, 'pro_sup_l': 500
#                   }

# #Create a dictionary for actuator limits used in Hamner & Delp
# actuatorLimits = {'pelvis_tx': 10000, 'pelvis_ty': 10000, 'pelvis_tz': 10000,
#                   'pelvis_tilt': 10000, 'pelvis_list': 10000, 'pelvis_rotation': 10000,
#                   'hip_flexion_r': 1, 'hip_adduction_r': 1, 'hip_rotation_r': 1,
#                   'knee_angle_r': 1, 'ankle_angle_r': 1,
#                   'hip_flexion_l': 1, 'hip_adduction_l': 1, 'hip_rotation_l': 1,
#                   'knee_angle_l': 1, 'ankle_angle_l': 1,
#                   'lumbar_extension': 1, 'lumbar_bending': 1, 'lumbar_rotation': 1,
#                   'arm_flex_r': 1, 'arm_add_r': 1, 'arm_rot_r': 1,
#                   'elbow_flex_r': 1, 'pro_sup_r': 1,
#                   'arm_flex_l': 1, 'arm_add_l': 1, 'arm_rot_l': 1,
#                   'elbow_flex_l': 1, 'pro_sup_l': 1
#                   }

#Set ground contact force calculation
forceThreshold = 50

#Set the right and left bodies to consider the applied forces to
rightBodyName = 'calcn_r'
leftBodyName = 'calcn_l'

# %% Modify the model to be suitable in subsequent processes

#Read in model in processor
modelProc = osim.ModelProcessor(os.path.join('data_cut','scaledModelAdjusted.osim'))

#Remove muscles
modelProc.append(osim.ModOpRemoveMuscles())

#Process to access model
procModel = modelProc.process()

#Remove contact geometry and force set
procModel.updContactGeometrySet().clearAndDestroy()
procModel.updForceSet().clearAndDestroy()

#Print model to file
procModel.finalizeConnections()
procModel.printToXML(os.path.join('data_cut','processedModel.osim'))

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
idTool.setModelFileName(os.path.join('data_cut','processedModel.osim'))

#Set coordinates
idTool.setCoordinatesFileName(os.path.join('data_cut','rra1_Kinematics_q.sto'))

#Set time ranges
idTool.setStartTime(osim.Storage(os.path.join('data_cut','rra1_Kinematics_q.sto')).getFirstTime())
idTool.setEndTime(osim.Storage(os.path.join('data_cut','rra1_Kinematics_q.sto')).getLastTime())

#Set output file names
idTool.setResultsDir('outputs_cut')
idTool.setOutputGenForceFileName('id_genForces.sto')

#Print and reload to avoid filepath issues
idTool.printToXML('setupID.xml')
osim.InverseDynamicsTool('setupID.xml').run()

#Rename body forces file
os.replace(os.path.join('outputs_cut','body_forces_at_joints.sto'),
           os.path.join('outputs_cut','id_bodyForces.sto'))

##### ----- PREPARE FILES FOR STORING ZMP GRFs ----- #####

#Read in and initialise model

#Load model for use
osimModel = osim.Model(os.path.join('data_cut','processedModel.osim'))

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
statesTraj = osim.StatesTrajectory.createFromStatesStorage(osimModel,
                                                           os.path.join('data_cut','rra1_states.sto'))

#Read in body forces data to transform to ground and calculate COP
idBodyForces = osim.TimeSeriesTable(os.path.join('outputs_cut','id_bodyForces.sto'))

#Get number of times in ID results
nt = idBodyForces.getNumRows()

#Create storage file to append ZMP GRFs to
zmpResults = osim.Storage(nt)

#Create columns for ZMP force and COP calculations
#These are just allocated to 1 force plate rather than left & right sides
#(i.e. 6 x columns per force plate for this data)
zmpLabels = osim.ArrayStr('time', 6*1+1)
zmpLabels.set(1, 'ground_force_1_vx')
zmpLabels.set(2, 'ground_force_1_vy')
zmpLabels.set(3, 'ground_force_1_vz')
zmpLabels.set(4, 'ground_force_1_px')
zmpLabels.set(5, 'ground_force_1_py')
zmpLabels.set(6, 'ground_force_1_pz')

##### ----- PREDICT GRFs AT EACH STATE USING ZMP ----- #####

#Loop through each time-point in states
for sInd in range(nt):
    
    #Set vector to store results for current state in (6 values per limb)
    #Set all values to zero
    zmpVec = osim.Vector(6*1, 0)
    
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
        
        # #Get the position of the right and left body origins in the ground
        # #This will determine which is closer to the predicted COP
        # rightBodyPos = osimModel.updBodySet().get(rightBodyName).findStationLocationInGround(s, osim.Vec3(0,0,0))
        # leftBodyPos = osimModel.updBodySet().get(leftBodyName).findStationLocationInGround(s, osim.Vec3(0,0,0))
        
        # #Calculate distances from body to predicted ZMP
        # rightBodyDist = (rightBodyPos.get(0) - zmpCOP.get(0))**2 + (rightBodyPos.get(1) - zmpCOP.get(1))**2 + (rightBodyPos.get(2) - zmpCOP.get(2))**2
        # leftBodyDist = (leftBodyPos.get(0) - zmpCOP.get(0))**2 + (leftBodyPos.get(1) - zmpCOP.get(1))**2 + (leftBodyPos.get(2) - zmpCOP.get(2))**2
        
        #Allocate data to vector
        zmpVec.set(0, groundF.get(0))
        zmpVec.set(1, groundF.get(1))
        zmpVec.set(2, groundF.get(2))
        zmpVec.set(3, zmpCOP.get(0))
        zmpVec.set(4, zmpCOP.get(1))
        zmpVec.set(5, zmpCOP.get(2))
        
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
osim.Storage().printResult(zmpResults, 'id_zmpForces', 'outputs_cut', -1, '.sto')

"""

NOTES:
    
    > Forces are OK
    > COP is off again in the Z direction, OK in the X direction
    > A bit noisy at beginning and end of trial

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
    
NOTES:
    
    > Given the ID and this manual approach produce near identical results, have
      not repeated here for the cutting motion

"""

# %% Compare outputs from different approaches

"""

This just focuses on the ID ZMP calculations as the manual approach wasn't run for cutting

"""

##### ----- COMPARE GRFs ----- #####

#Create figure
fig,ax = plt.subplots(nrows = 1, ncols = 3, figsize = (11,5))
plt.subplots_adjust(left = 0.06, right = 0.975, top = 0.9, bottom = 0.125,
                    wspace = 0.3)

#Read in experimental and ZMP generated GRFs
expGRF = osim.TimeSeriesTable(os.path.join('data_cut','A11_grf.mot'))
zmpGRF = osim.TimeSeriesTable(os.path.join('outputs_cut','id_zmpForces.sto'))

#Set start and end times for comparison
startTime = zmpGRF.getIndependentColumn()[0]
endTime = zmpGRF.getIndependentColumn()[-1]

#Get index in full experimental trial that matches the gait cycle
expStartInd = expGRF.getNearestRowIndexForTime(startTime)
expEndInd = expGRF.getNearestRowIndexForTime(endTime)

#Sum and plot XYZ ground reaction forces
for axesLabel in 'xyz':
    
    #Extract experimental GRFs to numpy array
    expData = expGRF.getDependentColumn(f'ground_force_1_v{axesLabel}').to_numpy()[expStartInd:expEndInd+1]
    expTime = expGRF.getIndependentColumn()[expStartInd:expEndInd+1]
    
    #Extract ZMP GRFs to numpy array
    zmpData = zmpGRF.getDependentColumn(f'ground_force_1_v{axesLabel}').to_numpy() 
    zmpTime = zmpGRF.getIndependentColumn()
    
    #Plot current axes data
    ax['xyz'.index(axesLabel)].plot(expTime, expData, ls = '-', lw = 1.5, c = 'black', label = 'Experimental')
    ax['xyz'.index(axesLabel)].plot(zmpTime, zmpData, ls = '--', lw = 1.5, c = 'red', label = 'ZMP Predicted')
    
    #Add labels and titles
    ax['xyz'.index(axesLabel)].set_xlabel('Time (s)', fontsize = 12, fontweight = 'bold')
    ax['xyz'.index(axesLabel)].set_ylabel('Force (N)', fontsize = 12, fontweight = 'bold')
    ax['xyz'.index(axesLabel)].set_title(f'{axesLabel.capitalize()}-Axis GRFs', fontsize = 14, fontweight = 'bold')
    
    #Add legend to last axis
    if axesLabel == 'z':
        ax['xyz'.index(axesLabel)].legend()
    
#Save figure
fig.savefig(os.path.join('outputs_cut','grfComparison.png'), format = 'png', dpi = 300)

#Close figure
plt.close()    
    
##### ----- COMPARE COP ----- #####

#Create figure
fig,ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))

#Extract experimental COP in X and Z axis to numpy array 
expCopX = expGRF.getDependentColumn('ground_force_1_px').to_numpy()[expStartInd:expEndInd+1]
expCopZ = expGRF.getDependentColumn('ground_force_1_pz').to_numpy()[expStartInd:expEndInd+1]
# expTime = expGRF.getIndependentColumn()[expStartInd:expEndInd+1]

#Extract ZMP predicted COP in X and Z axis to numpy array
zmpCopX = zmpGRF.getDependentColumn('ground_force_1_px').to_numpy()
zmpCopZ = zmpGRF.getDependentColumn('ground_force_1_pz').to_numpy()
# zmpTime = zmpGRF.getIndependentColumn()

#Replace zeros with NaN's as this means no contact is present
expCopX[expCopX == 0] = np.nan
expCopZ[expCopZ == 0] = np.nan
zmpCopX[zmpCopX == 0] = np.nan
zmpCopZ[zmpCopZ == 0] = np.nan

#Plot data
ax.plot(expCopZ, expCopX, ls = '-', lw = 1.0, c = 'black', label = 'Experimental')
ax.plot(zmpCopZ, zmpCopX, ls = '--', lw = 1.0, c = 'red', label = 'ZMP Predicted')

#Set equal axis
ax.axis('equal')

#Add labels and titles
ax.set_xlabel('X-Position (m)', fontsize = 12, fontweight = 'bold')
ax.set_ylabel('Z-Position (m)', fontsize = 12, fontweight = 'bold')
ax.set_title('Centre of Pressure', fontsize = 14, fontweight = 'bold')

#Add legend
ax.legend()

#Tight layout
plt.tight_layout()

#Save figure
fig.savefig(os.path.join('outputs_cut','copComparison.png'), format = 'png', dpi = 300)

#Close figure
plt.close()

# %% ----- End of pythonImplementationZMP_cut.py ----- %% #