# -*- coding: utf-8 -*-
'''

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code is a preliminary step in extracting just the relevant data needed
    from the SchreiberMoissenet2019 dataset. There are a number of C3D files in
    each participants folder, however we only extract two speeds here. The c3d
    files are converted to OpenSim friendly formats, those being TRC and MOT
    file types. Each participant needs to have their data within a 'raw' folder
    and then in a folder named with the participant code (e.g. '2014001') for
    this script to work.
    
'''

# %% Import packages

import opensim as osim
import os
import glob
import shutil
import numpy as np

# %% Set-up
        
#Get home path
homeDir = os.getcwd()

#Set the participant list based on the codes in the raw folder
participantList = [ii for ii in os.listdir('raw') if os.path.isdir(ii)]

#Set the list of speed conditions to extract 
#Modify this if you want to extract different walking speeds
conditionList  = [
    'C1', #0-0.4 m/s
    # 'C2', #0.4-0.8 m/s
    'C3', #0.8-1.2 m/s
    # 'C4', #self-selected spontaneuous speed
    # 'C5', #self-selected fast speed
    ]

# %% Loop through participants to extract data

#Start loop
for participant in participantList:
        
    # =========================================================================
    #     Create folders for the participant
    # =========================================================================
    
    #Starting directory
    os.makedirs(participant, exist_ok = True)
    #Static files directory
    os.makedirs(os.path.join(participant, 'static'), exist_ok = True)
    #Dynamic trials file directory
    os.makedirs(os.path.join(participant, 'dynamic'), exist_ok = True)

    # =========================================================================
    #     Get the static file
    # =========================================================================

    #Get the path to the static file
    staticFile = glob.glob(os.path.join('raw', participant, '*_ST.c3d*'))[0]
    
    #Copy c3d across to static directory
    shutil.copyfile(staticFile, os.path.join(participant, 'static', os.path.split(staticFile)[-1]))
    
    #Construct opensim 3d object
    c3dFile = osim.C3DFileAdapter()
    c3dFile.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)
    
    #Read in the static trial
    staticC3D = c3dFile.read(staticFile)
    
    #Get markers table
    staticMarkers = c3dFile.getMarkersTable(staticC3D)
    
    #Rotate marker data
    
    #Create the two rotations needed
    markerRot1 = osim.Rotation(np.deg2rad(-90), osim.Vec3(0,0,1))
    markerRot2 = osim.Rotation(np.deg2rad(-90), osim.Vec3(1,0,0))
    
    #Rotate the data
    for iRow in range(staticMarkers.getNumRows()):
        #Apply the two rotations
        staticMarkers.setRowAtIndex(iRow, markerRot1.multiply(staticMarkers.getRowAtIndex(iRow)))
        staticMarkers.setRowAtIndex(iRow, markerRot2.multiply(staticMarkers.getRowAtIndex(iRow)))
        
    #There's a potentially annoying bug to deal with in scaling when only one row of marker data exists
    #Add a pseudo row of identical marker data to avoid this
    staticMarkers.appendRow(1, staticMarkers.getRow(0))
        
    #Write static markers to TRC file
    osim.TRCFileAdapter().write(staticMarkers, os.path.join(participant, 'static', os.path.split(staticFile)[-1].split('.')[0]+'.trc'))
    
    # =========================================================================
    #     Get the dynamic files
    # =========================================================================
    
    #Identify dynamic files
    dynamicFiles = glob.glob(os.path.join('raw', participant, '*_C*.c3d*'))
    
    #Extract those with the desired condition label
    extractFiles = [dynaFile for dynaFile in dynamicFiles if os.path.split(dynaFile)[-1].split('_')[1] in conditionList]
    
    #Loop through files to copy and convert
    for dynaFile in extractFiles:
        
        #Copy c3d across to dynamic directory
        shutil.copyfile(dynaFile, os.path.join(participant, 'dynamic', os.path.split(dynaFile)[-1]))
        
        #Construct opensim 3d object
        c3dFile = osim.C3DFileAdapter()
        c3dFile.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)
        
        #Read in the dynamic trial
        dynaC3D = c3dFile.read(dynaFile)
        
        # =====================================================================
        #     Marker data
        # =====================================================================
        
        #Get markers table
        dynaMarkers = c3dFile.getMarkersTable(dynaC3D)
        
        #Rotate marker data (se same rotations as earlier)
        for iRow in range(dynaMarkers.getNumRows()):
            #Apply the two rotations
            dynaMarkers.setRowAtIndex(iRow, markerRot1.multiply(dynaMarkers.getRowAtIndex(iRow)))
            dynaMarkers.setRowAtIndex(iRow, markerRot2.multiply(dynaMarkers.getRowAtIndex(iRow)))
            
        #Write dynamic markers to TRC file
        osim.TRCFileAdapter().write(dynaMarkers, os.path.join(participant, 'dynamic', os.path.split(dynaFile)[-1].split('.')[0]+'.trc'))
        
        # =====================================================================
        #     Forces data
        # =====================================================================
        
        #Get forces table
        dynaForces = c3dFile.getForcesTable(dynaC3D)
        
        #Rotate forces data
        #Use the same rotations as earlier
        for iRow in range(dynaForces.getNumRows()):
            dynaForces.setRowAtIndex(iRow, markerRot1.multiply(dynaForces.getRowAtIndex(iRow)))
            dynaForces.setRowAtIndex(iRow, markerRot2.multiply(dynaForces.getRowAtIndex(iRow)))
            
        #Flatten forces data
        forcesFlat = dynaForces.flatten()
        
        #Convert to numpy array
        #Pre-allocate numpy array based on data size
        dataArray = np.zeros((forcesFlat.getNumRows(),
                              forcesFlat.getNumColumns()))
        #Extract data
        for forceInd in range(forcesFlat.getNumColumns()):
            dataArray[:,forceInd] = forcesFlat.getDependentColumn(forcesFlat.getColumnLabels()[forceInd]).to_numpy()
            
        #Replace nan's for COP and moment data with zeros
        np.nan_to_num(dataArray, copy = False, nan = 0.0)
        
        #Convert force point data from mm to m
        for forceName in list(forcesFlat.getColumnLabels()):
            if forceName.startswith('p') or forceName.startswith('m'):
                #Get force index
                forceInd = list(forcesFlat.getColumnLabels()).index(forceName)
                #Convert to m units in data array
                dataArray[:,forceInd] = dataArray[:,forceInd] / 1000

        #Build the new time series table
        forcesStorage = osim.Storage()
        
        #Get the time data
        time = forcesFlat.getIndependentColumn()
        
        #Create maps to replace text from force labels with
        #Force plate and type identifiers
        forceType = {}
        for ii in range(1,int(len(forcesFlat.getColumnLabels()) / 9)+1):
            forceType[f'f{ii}'] = f'ground_force_{ii}_v'
            forceType[f'p{ii}'] = f'ground_force_{ii}_p'
            forceType[f'm{ii}'] = f'ground_force_{ii}_m'
        #Axis identifiers
        forceAxis = {'1': 'x',
                     '2': 'y',
                     '3': 'z'}
        
        #Set labels in table
        newLabels = osim.ArrayStr()
        newLabels.append('time')
        for forceLabel in forcesFlat.getColumnLabels():
            #Split the label to get parts
            labelSplit = forceLabel.split('_')
            #Create new label
            forceLabel = f'{forceType[labelSplit[0]]}{forceAxis[labelSplit[1]]}'
            #Append to labels vector
            newLabels.append(forceLabel)
        forcesStorage.setColumnLabels(newLabels)
        
        #Add data
        for iRow in range(dataArray.shape[0]):
            row = osim.ArrayDouble()
            for iCol in range(dataArray.shape[1]):
                row.append(dataArray[iRow,iCol])
            #Add data to storage
            forcesStorage.append(time[iRow], row)
            
        #Set name for storage object
        forcesStorage.setName(os.path.split(dynaFile)[-1].split('.')[0]+'_grf')
        
        #Write to file
        forcesStorage.printResult(forcesStorage,
                                  os.path.split(dynaFile)[-1].split('.')[0]+'_grf',
                                  os.path.join(participant, 'dynamic'),
                                  0.001, '.mot')
        
        # =====================================================================
        #     External loads data
        # =====================================================================
        
        #Note that not all trials have force plate contacts on both plates
        #This is accounted for in applying external loads
        
        #Create the external loads file
        forceXML = osim.ExternalLoads()
        
        #Convert forces to time-series table for easier use
        forcesTable = osim.TimeSeriesTable(os.path.join(participant, 'dynamic', os.path.split(dynaFile)[-1].split('.')[0]+'_grf.mot'))
        
        #Extract the vertical force data from the two plates
        vForce1 = forcesTable.getDependentColumn('ground_force_1_vy').to_numpy()
        vForce2 = forcesTable.getDependentColumn('ground_force_2_vy').to_numpy()
        
        #Identify contact indices and times on plates based on force threshold > 20
        forceThreshold = 20
        
        #Check for no force data
        if np.sum(vForce1 > forceThreshold) > 0:
            v1_ind = np.argmax(vForce1 > forceThreshold)
            #Check if index is at least past the first frame
            #If it isn't then the foot is already on the plate
            if v1_ind == 0:
                v1_ind = None            
        else:
            v1_ind = None
        if np.sum(vForce2 > forceThreshold) > 0:
            v2_ind = np.argmax(vForce2 > forceThreshold)
            #Check if index is at least past the first frame
            #If it isn't then the foot is already on the plate
            if v2_ind == 0:
                v2_ind = None  
        else:
            v2_ind = None
        
        #Grab the heel marker at each force plate contact
        #Determine which is closer to the force plate CoP to determine which
        #foot is in contact with the plate
        
        #Force plate 1
        if v1_ind is not None:
        
            #Get the time on the plate
            v1_t = forcesTable.getIndependentColumn()[v1_ind]
            
            #Get the COP
            cop1 = np.array((forcesTable.getDependentColumn('ground_force_1_px').to_numpy()[v1_ind],
                             forcesTable.getDependentColumn('ground_force_1_py').to_numpy()[v1_ind],
                             forcesTable.getDependentColumn('ground_force_1_pz').to_numpy()[v1_ind]))
            
            #Get the two markers at the closest timestamp
            #Divide by 1000 to convert from mm to m
            m1_ind = dynaMarkers.getNearestRowIndexForTime(v1_t)
            rightMkr = np.array((dynaMarkers.flatten().getDependentColumn('R_FCC_1').to_numpy()[m1_ind],
                                 dynaMarkers.flatten().getDependentColumn('R_FCC_2').to_numpy()[m1_ind],
                                 dynaMarkers.flatten().getDependentColumn('R_FCC_3').to_numpy()[m1_ind])) / 1000
            leftMkr = np.array((dynaMarkers.flatten().getDependentColumn('L_FCC_1').to_numpy()[m1_ind],
                                dynaMarkers.flatten().getDependentColumn('L_FCC_2').to_numpy()[m1_ind],
                                dynaMarkers.flatten().getDependentColumn('L_FCC_3').to_numpy()[m1_ind])) / 1000
            
            #Calculate the distances
            rightDist = np.sqrt(np.sum((cop1 - rightMkr)**2, axis = 0))
            leftDist = np.sqrt(np.sum((cop1 - leftMkr)**2, axis = 0))
            
            #Create external loads for plate 1
            grf1 = osim.ExternalForce()
            grf1.setName('grf1')
            if rightDist < leftDist:
                grf1.setAppliedToBodyName('calcn_r')
            elif leftDist < rightDist:
                grf1.setAppliedToBodyName('calcn_l')
            grf1.setForceExpressedInBodyName('ground')
            grf1.setPointExpressedInBodyName('ground')
            grf1.setForceIdentifier('ground_force_1_v')
            grf1.setPointIdentifier('ground_force_1_p')
            grf1.setTorqueIdentifier('ground_force_1_m')
            forceXML.cloneAndAppend(grf1)
        
        #Force plate 2
        if v2_ind is not None:
            
            #Get the time on the plate
            v2_t = forcesTable.getIndependentColumn()[v2_ind]
        
            #Get the COP
            cop2 = np.array((forcesTable.getDependentColumn('ground_force_2_px').to_numpy()[v2_ind],
                             forcesTable.getDependentColumn('ground_force_2_py').to_numpy()[v2_ind],
                             forcesTable.getDependentColumn('ground_force_2_pz').to_numpy()[v2_ind]))
            
            #Get the two markers at the closest timestamp
            #Divide by 1000 to convert from mm to m
            m2_ind = dynaMarkers.getNearestRowIndexForTime(v2_t)
            rightMkr = np.array((dynaMarkers.flatten().getDependentColumn('R_FCC_1').to_numpy()[m2_ind],
                                 dynaMarkers.flatten().getDependentColumn('R_FCC_2').to_numpy()[m2_ind],
                                 dynaMarkers.flatten().getDependentColumn('R_FCC_3').to_numpy()[m2_ind])) / 1000
            leftMkr = np.array((dynaMarkers.flatten().getDependentColumn('L_FCC_1').to_numpy()[m2_ind],
                                dynaMarkers.flatten().getDependentColumn('L_FCC_2').to_numpy()[m2_ind],
                                dynaMarkers.flatten().getDependentColumn('L_FCC_3').to_numpy()[m2_ind])) / 1000
            
            #Calculate the distances
            rightDist = np.sqrt(np.sum((cop2 - rightMkr)**2, axis = 0))
            leftDist = np.sqrt(np.sum((cop2 - leftMkr)**2, axis = 0))
            
            #Create external loads for plate 2
            grf2 = osim.ExternalForce()
            grf2.setName('grf2')
            if rightDist < leftDist:
                grf2.setAppliedToBodyName('calcn_r')
            elif leftDist < rightDist:
                grf2.setAppliedToBodyName('calcn_l')
            grf2.setForceExpressedInBodyName('ground')
            grf2.setPointExpressedInBodyName('ground')
            grf2.setForceIdentifier('ground_force_2_v')
            grf2.setPointIdentifier('ground_force_2_p')
            grf2.setTorqueIdentifier('ground_force_2_m')
            forceXML.cloneAndAppend(grf2)
        
        #Set GRF datafile in external loads
        forceXML.setDataFileName(os.path.split(dynaFile)[-1].split('.')[0]+'_grf.mot')
        
        #Write to file
        forceXML.printToXML(os.path.join(participant, 'dynamic', os.path.split(dynaFile)[-1].split('.')[0]+'_grf.xml'))
        
# %% ----- end of 1_extractData.py ----- %% #