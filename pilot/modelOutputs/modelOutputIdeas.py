# -*- coding: utf-8 -*-
"""

@author:

    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Script that tinkers with the idea of implementing Model outputs functions in
    Model.cpp and Model.h

"""


# %% Import packages

import opensim as osim
import os
import numpy as np

# %% Set-up

#Import model geometry for cleanliness
geomDir = os.path.join('C:', os.sep, 'OpenSim 4.4', 'Geometry')
osim.ModelVisualizer.addDirToGeometrySearchPaths(geomDir)

#Load model for use
osimModel = osim.Model('osimModel.osim')
osimModel.initSystem()

#Read in moco trajectory
mocoTraj = osim.MocoTrajectory('subject01_run5_cycle1_mocoSolution.sto')
statesTraj = mocoTraj.exportToStatesTrajectory(osimModel)
mocoTraj.generateAccelerationsFromSpeeds() #only for testing
uDotTable = mocoTraj.exportToAccelerationsTable() #only for testing

# %% Function with state, free body joint, single contact body, force threshold

""" TODO: how do determine when body is in contact...? """
""" TODO: multi-body force order based on correct q ordering...? """

#Set some sample inputs to test with
s = statesTraj[27] #t = 0.683
freeBodyJointName = 'ground_pelvis' ### is there an easy way to find this without input?
contactBodyname = 'calcn_r'
forceThreshold = 20

# def model_calcZeroMomentPointGroundReactions(s, freeBodyJoint, contactBodyName, forceThreshold):

#Initialise an inverse dynamics solver with model
ivdSolver = osim.InverseDynamicsSolver(osimModel)

#Realise model to acceleration stage with state
osimModel.realizeAcceleration(s)

# #Check if Q's are updated
# for qi in range(s.getQ().size()):
#     print(s.getQ().get(qi))

#Get accelerations
### NOTE: this seems to be where results differ from getting MocoTraj accelerations?
udot = s.getUDot() ### TODO: is this right way to get accelerations?
udot2 = osim.Vector(osimModel.updCoordinateSet().getSize(), 0)
for iq in range(osimModel.updCoordinateSet().getSize()):
    udot2.set(iq,
              uDotTable.getRowAtIndex(uDotTable.getNearestRowIndexForTime(s.getTime())).getElt(0,iq)
              )
for u in range(udot.size()):
    print('getUDot: '+str(udot.get(u))+' / mocoTraj: '+str(udot2.get(u)))
    
# osimModel.realizeAcceleration(s)
    
#Solve ID
""" Gives all zeros when using udot for some reason??? """
genForceTraj = ivdSolver.solve(s, udot2)
# for u in range(genForceTraj.size()):
#     print(genForceTraj.get(u))

#Calculate equivalent body force at free joint

#Get the free body joint based on name
freeJoint = osimModel.getJointSet().get(freeBodyJointName)

#Calculate equivalent force
equivalentBodyForceAtJoint = freeJoint.calcEquivalentSpatialForce(s, genForceTraj)

#Get the free torque components
freeBodyTorque = equivalentBodyForceAtJoint.get(0)
freeBodyForce = equivalentBodyForceAtJoint.get(1)

#Check if vertical force greater than force threshold for ZMP calculations
#Assumes vertical axis is Y
#Can probably just allocate zeros otherwise
if freeBodyForce.get(1) > forceThreshold:
    
    #Get the translational positions of the free joint in the ground
    #Assumes that the origin of the child frame is the best point
    rp = freeJoint.getChildFrame().findStationLocationInGround(s, osim.Vec3(0,0,0))
    
    #Take cross product of pelvis position and force vector to get moment at origin
    groundM = osim.Vec3((rp.get(1) * freeBodyForce.get(2)) - (rp.get(2) * freeBodyForce.get(1)),
                        -((rp.get(0) * freeBodyForce.get(2)) - (rp.get(2) * freeBodyForce.get(0))),
                        (rp.get(0) * freeBodyForce.get(1)) - (rp.get(1) * freeBodyForce.get(0)))
        
    #Calculate X & Z cZMP, noting that yZMP is set as 0
    #Again, assumes y-axis is vertical
    #Could identify this based on model gravity?
    #Formulas come from Xiang et al.
    zmpCOP = osim.Vec3(groundM.get(2) / freeBodyForce.get(1),
                       0,
                       -groundM.get(0) / freeBodyForce.get(1))
    
    #Transform the moment at the origin to the ZMP CoP
    
    
    #Calculate the resultant active moment at ZMP along the y-axis
    #Scaling by 1000 makes sense to get in Nmm?
    myZMP = (groundM.get(1) + (freeBodyForce.get(0) * zmpCOP.get(2)) - (freeBodyForce.get(2) * zmpCOP.get(0))) / 1000
   
    """ TODO: body contact with ground:
        > How to determine?
        > Needed in this single contact body iteration? """


# %% ----- End of modelOutputIdeas.py ----- %% #