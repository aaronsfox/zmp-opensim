# -*- coding: utf-8 -*-
'''

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code is a preliminary step in cleaning the relevant data needed from the
    Rankin2016 dataset. Not a lot is done here, just some editing of the model to
    make it torque-driven, adusting the experimental GRFs to be in the ground
    frame, and creating XML files for external loads. Perhaps the biggest change
    here is converting the plane of the model ground to equate to y = 0, as the
    original data seemingly has the ground plane higher than this. This isn't a
    problem, except it may cause some issues with later ZMP calculations that will
    likely assume the ground is at y = 0.
    
'''

# %% Import packages

import opensim as osim
import os
import numpy as np
import osimFunctions as helper

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

# %% Edit model

#Load model into processor
modelProc = osim.ModelProcessor(os.path.join('models', 'Rankin2016_Ostrich.osim'))

#Remove muscles
modelProc.append(osim.ModOpRemoveMuscles())

#Process model to access
osimModel = modelProc.process()
        
#Finalise model connections
osimModel.finalizeConnections()

#Save model to file
osimModel.printToXML(os.path.join('models', 'Rankin2016_Ostrich_noMuscles.osim'))

# %% Clean up experimental files

# =============================================================================
# Convert kinematics to states version for use with Moco
# =============================================================================

#Walk
helper.kinematicsToStates(kinematicsFileName = os.path.join('data', 'Ostrich_Walk.mot'),
                          osimModelFileName = os.path.join('models', 'Rankin2016_Ostrich_noMuscles.osim'),
                          outputFileName = os.path.join('data', 'Ostrich_Walk_coordinates.mot'),
                          inDegrees = True, outDegrees = False)

#Run
helper.kinematicsToStates(kinematicsFileName = os.path.join('data', 'Ostrich_Run.mot'),
                          osimModelFileName = os.path.join('models', 'Rankin2016_Ostrich_noMuscles.osim'),
                          outputFileName = os.path.join('data', 'Ostrich_Run_coordinates.mot'),
                          inDegrees = True, outDegrees = False)

# =============================================================================
# Convert GRF .mot file to be in ground frame
# =============================================================================

#1_ground_* is for the L_pes segment
#ground_ is for the R_pes segment

#Initialise model state
s = osimModel.initSystem()

#Loop through trial types
for trial in ['Walk', 'Run']:

    #Read in the GRF file as a table
    grfTable = osim.TimeSeriesTable(os.path.join('data', f'Ostrich_{trial}_GRF.mot'))
    
    #Read in the kinematic data to pose model
    kinematicTable = osim.TimeSeriesTable(os.path.join('data', f'Ostrich_{trial}_coordinates.mot'))
    
    #Initialise a numpy array to store new COP positions
    L_pes_cop = np.zeros((grfTable.getNumRows(), 3))
    R_pes_cop = np.zeros((grfTable.getNumRows(), 3))
    
    #Loop through GRF time points
    for ii in range(grfTable.getNumRows()):
        
        #Get GRF time for index
        t = grfTable.getIndependentColumn()[ii]
        
        #Get closest index in kinematic data
        jj = kinematicTable.getNearestRowIndexForTime(t)
        
        #Loop through model state coordinates and set
        for stateLabel in kinematicTable.getColumnLabels():
            osimModel.setStateVariableValue(s, stateLabel,
                                            kinematicTable.getDependentColumn(stateLabel).to_numpy()[jj])
            
        #Realize to position
        osimModel.realizePosition(s)
        
        #Get the vertical force to determine whether COP needs to be transformed
        
        #Left pes
        if grfTable.getDependentColumn('1_ground_force_vy').to_numpy()[ii] > 0:
            
            #Get the current COP in the body frame
            localCOP = osim.Vec3(grfTable.getDependentColumn('1_ground_force_px').to_numpy()[ii],
                                 grfTable.getDependentColumn('1_ground_force_py').to_numpy()[ii],
                                 grfTable.getDependentColumn('1_ground_force_pz').to_numpy()[ii])
            
            #Get new COP position for this segment
            newCOP = osimModel.updBodySet().get('L_pes').findStationLocationInGround(s, localCOP)
            
            #Set in appropriate array
            L_pes_cop[ii,0] = newCOP.get(0)
            L_pes_cop[ii,1] = newCOP.get(1)
            L_pes_cop[ii,2] = newCOP.get(2)
            
        #Right pes
        if grfTable.getDependentColumn('ground_force_vy').to_numpy()[ii] > 0:
            
            #Get the current COP in the body frame
            localCOP = osim.Vec3(grfTable.getDependentColumn('ground_force_px').to_numpy()[ii],
                                 grfTable.getDependentColumn('ground_force_py').to_numpy()[ii],
                                 grfTable.getDependentColumn('ground_force_pz').to_numpy()[ii])
            
            #Get new COP position for this segment
            newCOP = osimModel.updBodySet().get('R_pes').findStationLocationInGround(s,localCOP)
            
            #Set in appropriate array
            R_pes_cop[ii,0] = newCOP.get(0)
            R_pes_cop[ii,1] = newCOP.get(1)
            R_pes_cop[ii,2] = newCOP.get(2)
            
    #Re-loop through GRF time points to set rows with new COP data
    for ii in range(grfTable.getNumRows()):
            
        #Create numpy array in appropriate ordering for GRF file row
        #Here we enforce all vertical COP data to be 0 for correcting ground level to zero
        newDataArr = np.array([
            grfTable.getDependentColumn('ground_force_vx').to_numpy()[ii],
            grfTable.getDependentColumn('ground_force_vy').to_numpy()[ii],
            grfTable.getDependentColumn('ground_force_vz').to_numpy()[ii],
            # R_pes_cop[ii,0], R_pes_cop[ii,1], R_pes_cop[ii,2],
            R_pes_cop[ii,0], 0, R_pes_cop[ii,2],
            grfTable.getDependentColumn('1_ground_force_vx').to_numpy()[ii],
            grfTable.getDependentColumn('1_ground_force_vy').to_numpy()[ii],
            grfTable.getDependentColumn('1_ground_force_vz').to_numpy()[ii],
            # L_pes_cop[ii,0], L_pes_cop[ii,1], L_pes_cop[ii,2],
            L_pes_cop[ii,0], 0, L_pes_cop[ii,2],
            grfTable.getDependentColumn('ground_torque_x').to_numpy()[ii],
            grfTable.getDependentColumn('ground_torque_y').to_numpy()[ii],
            grfTable.getDependentColumn('ground_torque_z').to_numpy()[ii],
            grfTable.getDependentColumn('1_ground_torque_x').to_numpy()[ii],
            grfTable.getDependentColumn('1_ground_torque_y').to_numpy()[ii],
            grfTable.getDependentColumn('1_ground_torque_z').to_numpy()[ii],
            ])
        
        #Create row from array
        newDataRow = osim.RowVector_createFromMat(newDataArr)
        
        #Set in table
        grfTable.setRowAtIndex(ii, newDataRow)
        
    #Write to file
    osim.STOFileAdapter().write(grfTable, os.path.join('data', f'Ostrich_{trial}_GRF_ground.mot'))
    
    #Use vertical COP data to correct vertical pelvis positioning to ensure that
    #model kinematics reflect ground being equal to zero
    
    #Create a variable for new pelvis ty values
    new_pelvis_ty = np.zeros((kinematicTable.getNumRows()))
    
    #Loop through kinematic table rows
    for jj in range(kinematicTable.getNumRows()):
        
        #Get kinematic time for index
        t = kinematicTable.getIndependentColumn()[jj]
        
        #Get closest index in grf data
        #Check for data outside of range with exception
        try:
            ii = grfTable.getNearestRowIndexForTime(t)
        except:
            if t > grfTable.getIndependentColumn()[-1]:
                #Use last value as there are some slight differences in timing
                ii = grfTable.getNumRows() - 1
            elif t < grfTable.getIndependentColumn()[0]:
                #Use first index as there are some slight differences in timing
                ii = 0
        
        #Get associated vertical COP values from GRF data index
        L_vert_cop = L_pes_cop[ii,1]
        R_vert_cop = R_pes_cop[ii,1]
        
        #Identify vertical correction based on average of COP if both are in contact,
        #or just the one value if only one foot is in contact
        if L_vert_cop > 0 and R_vert_cop > 0:
            vertCorrection = (L_vert_cop + R_vert_cop) / 2
        elif L_vert_cop > 0 and R_vert_cop == 0:
            vertCorrection = L_vert_cop
        elif  L_vert_cop == 0 and R_vert_cop > 0:
            vertCorrection = R_vert_cop
        
        #Calculate new pelvis_ty value
        new_pelvis_ty[jj] = kinematicTable.getDependentColumn('/jointset/pelvis/pelvis_ty/value').to_numpy()[jj] - vertCorrection
        
    #Get corrected pelvis ty value in the kinematics table
    #It's easier to set a data column in a storage object, so reload here
    kinematicStorage = osim.Storage(os.path.join('data', f'Ostrich_{trial}_coordinates.mot'))
    
    #Create the array double from the new pelvis position data
    pelvisArrayDouble = osim.ArrayDouble()
    for jj in range(kinematicTable.getNumRows()):
        pelvisArrayDouble.set(jj, new_pelvis_ty[jj])
        
    #Get index of desired column in storage object
    #Subtract 1 as this includes the time in column labels
    pelvis_ty_ind = [kinematicStorage.getColumnLabels().get(ii) for ii in range(kinematicStorage.getColumnLabels().getSize())].index('/jointset/pelvis/pelvis_ty/value') - 1
    
    #Set data column for new pelvis ty
    kinematicStorage.setDataColumn(pelvis_ty_ind, pelvisArrayDouble)
    
    #Write to file
    kinematicStorage.printToXML(os.path.join('data', f'Ostrich_{trial}_coordinates.mot'))

# =============================================================================
# Create .xml external loads files
# =============================================================================

#Loop through trial types
for trial in ['Walk', 'Run']:
    
    #Create external loads
    extLoads = osim.ExternalLoads()
    
    #Create the external forces
    
    #Left pes
    extForceLeft = osim.ExternalForce()
    extForceLeft.set_applied_to_body('L_pes')
    extForceLeft.set_force_expressed_in_body('ground')
    extForceLeft.set_point_expressed_in_body('ground')
    extForceLeft.set_force_identifier('1_ground_force_v')
    extForceLeft.set_point_identifier('1_ground_force_p')
    extForceLeft.set_torque_identifier('1_ground_torque_')
    extLoads.cloneAndAppend(extForceLeft)
    
    #Right pes
    extForceRight = osim.ExternalForce()
    extForceRight.set_applied_to_body('R_pes')
    extForceRight.set_force_expressed_in_body('ground')
    extForceRight.set_point_expressed_in_body('ground')
    extForceRight.set_force_identifier('ground_force_v')
    extForceRight.set_point_identifier('ground_force_p')
    extForceRight.set_torque_identifier('ground_torque_')
    extLoads.cloneAndAppend(extForceRight)
    
    #Set the data file name
    extLoads.setDataFileName(f'Ostrich_{trial}_GRF_ground.mot')
    
    #Print to file
    extLoads.printToXML(os.path.join('data', f'Ostrich_{trial}_GRF_ground.xml'))

# %% ----- end of 1_cleanData.py ----- %% #