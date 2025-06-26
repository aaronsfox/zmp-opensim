# -*- coding: utf-8 -*-
"""

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Test implementation of ZMP predicted GRF using OpenSim interface. See 
    README.MD for any extra details.
    
    In particular this script has to deal with the double support phase of a
    walking movement.
    
    Relevant papers:
        
        Dijkstra & Gutierrez-Farewik (2015). J Biomech, 48: 3776-3781.
        
        Xiang et al. (2009). Int J NUmer Meth ENg, 79: 667-695.
            > Calculations for transforming and calculating COP
            
    TODO:
        
        > Consider Djikstra approach in subtracting accelerations from front limb etc.
            >> See section 2.2.3 & 2.2.4 in paper for details


"""

# %% Import packages

import opensim as osim
import os
import matplotlib.pyplot as plt
import numpy as np

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

#Set ground contact force calculation
forceThreshold = 20

#Set the right and left bodies to consider the applied forces to
rightBodyName = 'calcn_r'
leftBodyName = 'calcn_l'

#Set the markers on each foot related to ground contact
rightContactPts = ['R_FM1_ground', 'R_FM2_ground', 'R_FM5_ground', 'R_FCC_ground']
leftContactPts = ['L_FM1_ground', 'L_FM2_ground', 'L_FM5_ground', 'L_FCC_ground']

#Define the points of contact when left and right foots are in contact with ground
#Right ground contact = 0.150 - 0.732
#Left ground contact = 0.611 - 1.120
rightOn = 0.150
rightOff = 0.723
leftOn = 0.611
leftOff = 1.120

# %% Convert states file to kinematics

#Read in moco solution
mocoSolution = osim.MocoTrajectory(os.path.join('data_walk','2014001_C3_05_markerTracking_solution.sto'))

#Export to states table
statesTable = mocoSolution.exportToStatesTable()

#Get column labels
colLabels = list(statesTable.getColumnLabels())

#Loop through and rename column labels to joint angles where appropriate
newColLabels = []
keepCol = []
for col in colLabels:
    if col.endswith('/value'):
        newColLabels.append(col.split('/')[3])
        keepCol.append(True)
    else:
        newColLabels.append(col)
        keepCol.append(False)
        
#Remove undesired columns
for col in newColLabels:
    if col.endswith('/speed'):
        statesTable.removeColumn(col)

#Rename column labels
statesTable.setColumnLabels([xx for ii, xx in enumerate(newColLabels) if keepCol[ii]])

#Print to file
osim.STOFileAdapter.write(statesTable, os.path.join('data_walk','kinematics.sto'))


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
idTool.setModelFileName(os.path.join('data_walk','2014001_trackingModel.osim'))

#Set coordinates
idTool.setCoordinatesFileName(os.path.join('data_walk','kinematics.sto'))

#Set time ranges
idTool.setStartTime(osim.Storage(os.path.join('data_walk','kinematics.sto')).getFirstTime())
idTool.setEndTime(osim.Storage(os.path.join('data_walk','kinematics.sto')).getLastTime())

#Set output file names
idTool.setResultsDir('outputs_walk')
idTool.setOutputGenForceFileName('id_genForces.sto')

#Print and reload to avoid filepath issues
idTool.printToXML('setupID.xml')
osim.InverseDynamicsTool('setupID.xml').run()

#Rename body forces file
os.replace(os.path.join('outputs_walk','body_forces_at_joints.sto'),
           os.path.join('outputs_walk','id_bodyForces.sto'))

##### ----- PREPARE FILES FOR STORING ZMP GRFs ----- #####

#Read in and initialise model

#Load model for use
osimModel = osim.Model(os.path.join('data_walk','2014001_trackingModel.osim'))

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
mocoTraj = osim.MocoTrajectory(os.path.join('data_walk','2014001_C3_05_markerTracking_solution.sto'))
statesTraj = mocoTraj.exportToStatesTrajectory(osimModel)

#Read in body forces data to transform to ground and calculate COP
idBodyForces = osim.TimeSeriesTable(os.path.join('outputs_walk','id_bodyForces.sto'))

#Get number of times in ID results
nt = idBodyForces.getNumRows()

#Create storage file to append ZMP GRFs to
zmpResults = osim.Storage(nt)

#Create columns for ZMP force and COP calculations
#These are just allocated to 2 forces for left and right sides
#(i.e. 6 x columns per force plate for this data)
zmpLabels = osim.ArrayStr('time', 6*2+1)
zmpLabels.set(1, 'ground_force_left_vx')
zmpLabels.set(2, 'ground_force_left_vy')
zmpLabels.set(3, 'ground_force_left_vz')
zmpLabels.set(4, 'ground_force_left_px')
zmpLabels.set(5, 'ground_force_left_py')
zmpLabels.set(6, 'ground_force_left_pz')
zmpLabels.set(7, 'ground_force_right_vx')
zmpLabels.set(8, 'ground_force_right_vy')
zmpLabels.set(9, 'ground_force_right_vz')
zmpLabels.set(10, 'ground_force_right_px')
zmpLabels.set(11, 'ground_force_right_py')
zmpLabels.set(12, 'ground_force_right_pz')

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
    # idInd = idBodyForces.getNearestRowIndexForTime(sTime)
    
    #Get the forces and torques from the pelvis body forces
    pelvisF = osim.Vec3(idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FX')[sInd],
                        idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FY')[sInd],
                        idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_FZ')[sInd])
    pelvisM = osim.Vec3(idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_MX')[sInd],
                        idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_MY')[sInd],
                        idBodyForces.getDependentColumn('ground_pelvis_pelvis_offset_MZ')[sInd])
    
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
        
        #Note that there are some difficulties in this simulation with markers dragging
        #through the ground --- need to fix...
        
        #Single vs. double support was defined manually earlier...
                                                                                             
        #Identify whether right, left or both are in contact
        if rightOn <= sTime <= rightOff:
            rightFootOn = True
        else:
            rightFootOn = False
        if leftOn <= sTime <= leftOff:
            leftFootOn = True
        else:
            leftFootOn = False
        
        #Allocate data to vector
        #If both feet are on, the force and CoP needs to be split
        #### TODO: case when force is recognised and feet aren't on ground...?
        if rightFootOn and leftFootOn:
            
            #Identify the leading foot by which contact is in front
            #This would usually require identifying pelvis orientation for direction
            #Cheating here as we know the model is moving in +z direction
            rightContactLoc = osimModel.updBodySet().get(rightBodyName).findStationLocationInGround(s, osim.Vec3(0,0,0))
            leftContactLoc = osimModel.updBodySet().get(leftBodyName).findStationLocationInGround(s, osim.Vec3(0,0,0))
            if rightContactLoc.get(2) > leftContactLoc.get(2):
                
                #Get the position of the right-leading heel in the ground frame
                heelContactLoc = osimModel.updBodySet().get('calcn_r').findStationLocationInGround(s, osim.Vec3(0,0,0))
                
                #Get the position of the left-trailing toes in the ground frame
                toeContactLoc = osimModel.updBodySet().get('toes_l').findStationLocationInGround(s, osim.Vec3(0,0,0))
                
                #Calculate the 2D distance between the above points and the ZMP CoP (ignoring vertical y-axis)
                #This will be used to distribute forces based on distance
                leadDist = np.sqrt(((heelContactLoc.get(0) - zmpCOP.get(0))**2) + ((heelContactLoc.get(2) - zmpCOP.get(2))**2))
                trailDist = np.sqrt(((toeContactLoc.get(0) - zmpCOP.get(0))**2) + ((toeContactLoc.get(2) - zmpCOP.get(2))**2))
                
                #Partition the GRF based on distance between COP and body locations
                #See Xiang et al. 2009
                leadGRF = (trailDist / (trailDist + leadDist)) * groundF.to_numpy()
                trailGRF = (leadDist / (trailDist + leadDist)) * groundF.to_numpy()
                
            elif leftContactLoc.get(2) > rightContactLoc.get(2):
                
                #### NOTE: These might not to be all individually split...
                
                #Get the position of the left-leading heel in the ground frame
                heelContactLoc = osimModel.updBodySet().get('calcn_l').findStationLocationInGround(s, osim.Vec3(0,0,0))
                
                #Get the position of the right-trailing toes in the ground frame
                toeContactLoc = osimModel.updBodySet().get('toes_r').findStationLocationInGround(s, osim.Vec3(0,0,0))
                
                #Calculate the 2D distance between the above points and the ZMP CoP (ignoring vertical y-axis)
                #This will be used to distribute forces based on distance
                leadDist = np.sqrt(((heelContactLoc.get(0) - zmpCOP.get(0))**2) + ((heelContactLoc.get(2) - zmpCOP.get(2))**2))
                trailDist = np.sqrt(((toeContactLoc.get(0) - zmpCOP.get(0))**2) + ((toeContactLoc.get(2) - zmpCOP.get(2))**2))
                
                #Partition the GRF based on distance between COP and body locations
                #See Xiang et al. 2009
                leadGRF = (trailDist / (trailDist + leadDist)) * groundF.to_numpy()
                trailGRF = (leadDist / (trailDist + leadDist)) * groundF.to_numpy()
                leadGRM = (trailDist / (trailDist + leadDist)) * groundM.to_numpy()
                trailGRM = (leadDist / (trailDist + leadDist)) * groundM.to_numpy()
                
                #Partition the COP based on position vectors between points
                #See Xiang et al. 2009
                leadF = leadGRF
                leadM = leadGRM + (heelContactLoc.to_numpy() - zmpCOP.to_numpy()) * leadGRF                
                trailF = trailGRF                
                trailM = trailGRM + (toeContactLoc.to_numpy() - zmpCOP.to_numpy()) * trailGRF
            
            
            # ### TODO: calculations --- are these correct? Probably not...
            # leadCOP = osim.Vec3(leadM[2] / leadF[1],
            #                    0,
            #                    -leadM[0] / leadF[1])
            
            # trailCOP = osim.Vec3(trailM[2] / trailF[1],
            #                    0,
            #                    -trailM[0] / trailF[1])
            
            
        elif rightFootOn and not leftFootOn:
            
            #Right side forces
            zmpVec.set(0, groundF.get(0))
            zmpVec.set(1, groundF.get(1))
            zmpVec.set(2, groundF.get(2))
            
            #Right side COP
            zmpVec.set(3, zmpCOP.get(0))
            zmpVec.set(4, zmpCOP.get(1))
            zmpVec.set(5, zmpCOP.get(2))
            
        elif leftFootOn and not rightFootOn:
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
        
# #Set column labels in ZMP results storage
# zmpResults.setColumnLabels(zmpLabels)

# #Set name in ZMP results storage
# zmpResults.setName('ZMP Predicted Ground Reaction Forces')

# #Print ZMP results to file
# osim.Storage().printResult(zmpResults, 'id_zmpForces', 'outputs_cut', -1, '.sto')

"""

NOTES:
    
    > Overall forces look OK
    > Partitioning not quite complete
    > Probably need accurate moments for this --- not sure these are accurate
      though, they seem off by a scale of 10x. Could be a moments issue problem?
    
    
"""


# %% TODO: states-based implementation?

# %% ----- End of pythonImplementationZMP_walk.py ----- %% #