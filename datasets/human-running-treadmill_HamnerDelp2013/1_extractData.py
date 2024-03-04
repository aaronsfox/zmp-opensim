# -*- coding: utf-8 -*-
'''

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This code is a preliminary step in extracting just the relevant data needed
    from the HamnerDelp2013 dataset. There are numerous result and set-up files
    that we don't need, so we just keep the necessary things that we need (e.g.
    scaled model, IK results etc.). Each subject needs to have their data within
    a 'raw' folder and then in a folder named with the subject code (e.g. 
    'subject01') for this script to work.
    
'''

# %% Import packages

import opensim as osim
import os
import glob
import shutil
import re
import pickle

# %% Set-up

#Try adding OpenSim geometry path
#This can avoid issues with body geometries constantly not showing up when loading models
#Note that this is defaulting to the standard C: drive path for Windows and Opensim 4.4
if os.path.isdir(os.path.join('C:', os.sep, 'OpenSim 4.4', 'Geometry')):
    osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.join('C:', os.sep, 'OpenSim 4.4', 'Geometry'))
else:
    print(f'OpenSim geometry folder not located at {os.path.join("C:", os.sep, "OpenSim 4.4", "Geometry")}')
    print('This code can be modified or you can just deal with the consequences of not adding a geometry path...')
        
#Get home path
homeDir = os.getcwd()

#Set subject list
#Modify this list if you only want to use specific subject folders
subList = [
    'subject01',
    # 'subject02',
    # 'subject03',
    # 'subject04',
    # 'subject08',
    # 'subject10',
    # 'subject11',
    # 'subject17',
    # 'subject19',
    # 'subject20'
    ]

#Set run names list
#Modify this if you want to extract different run speeds
runList  = [
    #'run2',
    'run3',
    #'run4',
    'run5'
    ]

#Set run cycle list
#This labels the three gait cycles appropriately
cycleList = ['cycle1',
             'cycle2',
             'cycle3']

# %% Loop through subjects to extract data

#Start loop
for subject in subList:
        
    # =========================================================================
    #     Create folders for the subject
    # =========================================================================
    
    #Starting directory
    os.makedirs(subject, exist_ok = True)
    #Model directory
    os.makedirs(os.path.join(subject, 'model'), exist_ok = True)
    #Kinematics directory
    os.makedirs(os.path.join(subject, 'ik'), exist_ok = True)
    #Experimental data directory
    os.makedirs(os.path.join(subject, 'expData'), exist_ok = True)
      
    #Create dictionary to store timing data
    gaitTimings = {run: {cyc: {'initialTime': [], 'finalTime': []} for cyc in cycleList} for run in runList}
        
    #Set the scale and IK directories to extract data from
    scaleDir = os.path.join('raw', subject, 'scale')
    ikDir = os.path.join('raw', subject, 'ik')
    
    #Identify the name of the generic model, scale setup file and static marker data
    modelFile = glob.glob(os.path.join('raw',subject,f'*{subject}.osim'))[0]
    scaleSetupFile = glob.glob(os.path.join(f'{scaleDir}','*_setup_scale*.xml'))[0]
    staticFile = glob.glob(os.path.join('raw', subject, 'ExportedData', 'Static*.trc'))[0]   
    
    #Copy the generic model file, scale tool and static file to the model directory
    shutil.copyfile(modelFile, os.path.join(subject, 'model', 'genericModel.osim'))
    shutil.copyfile(scaleSetupFile, os.path.join(subject, 'model', 'setupScale.xml'))
    shutil.copyfile(staticFile, os.path.join(subject, 'model', 'static.trc'))
    
    # =========================================================================
    #     Scale the generic model
    # =========================================================================
    
    #Navigate to model directory
    os.chdir(os.path.join(subject, 'model'))
    
    #Load the scale tool
    scaleTool = osim.ScaleTool('setupScale.xml')
    
    #Alter parameters in scale tool
    scaleTool.getGenericModelMaker().setModelFileName('genericModel.osim')
    scaleTool.getModelScaler().setMarkerFileName('static.trc')
    scaleTool.getMarkerPlacer().setMarkerFileName('static.trc')
    
    #Ensure all the file labels are appropriate for subject
    scaleTool.getMarkerPlacer().setOutputMotionFileName(f'{subject}_static_output.mot')
    scaleTool.getMarkerPlacer().setOutputModelFileName(f'{subject}_adjusted_scaled.osim')
    scaleTool.getModelScaler().setOutputScaleFileName(f'{subject}_scaleSet_applied.xml')
    
    #Check that scaling times are appropriate
    
    #Model scaler
    initialTime = scaleTool.getModelScaler().getTimeRange().get(0)
    finalTime = scaleTool.getModelScaler().getTimeRange().get(1)
    #Check and fix if necessary
    if initialTime >= finalTime:
        newTimes = osim.ArrayDouble()
        newTimes.append(initialTime - 0.1)
        newTimes.append(finalTime)
        scaleTool.getModelScaler().setTimeRange(newTimes)
        
    #Marker placer
    initialTime = scaleTool.getMarkerPlacer().getTimeRange().get(0)
    finalTime = scaleTool.getMarkerPlacer().getTimeRange().get(1)
    #Check and fix if necessary
    if initialTime >= finalTime:
        newTimes = osim.ArrayDouble()
        newTimes.append(initialTime - 0.1)
        newTimes.append(finalTime)
        scaleTool.getMarkerPlacer().setTimeRange(newTimes)
    
    #Run the scale tool
    scaleTool.run()
    
    #Navigate back to home directory
    os.chdir(homeDir)
    
    # =========================================================================
    #     Identify the IK files
    # =========================================================================
    
    #Find the folder that starts with results
    ikResultsDir = glob.glob(os.path.join(ikDir, 'Results*'))[0]
    
    #Identify .mot files in this directory
    ikFilesList = glob.glob(os.path.join(ikResultsDir,'*.mot'))
    
    #Loop through files to copy and rename
    for ikFile in ikFilesList:
        
        #Split up fileparts
        fileName = os.path.split(ikFile)[-1]
        
        #Identify what speed the file is
        runSpeed = fileName.split('Run_')[-1][0]
        
        #Copy and rename the IK file if it is in the run list
        if 'run'+runSpeed in runList:
            
            #Copy the file
            shutil.copyfile(ikFile, os.path.join(subject, 'ik', f'Run_{runSpeed}.mot'))
        
            #Copy and rename the trc file
            trcFile = glob.glob(os.path.join('raw', subject, 'ExportedData', f'Run_{runSpeed}*.trc'))[0]
            shutil.copyfile(trcFile, os.path.join(subject, 'expData', f'Run_{runSpeed}.trc'))
        
            #Find the relevant RRA directory
            rraDir = glob.glob(os.path.join('raw', subject, 'rra_multipleSteps', f'*Run_{runSpeed}*'))[0]
        
            #Get the trial specific cycle list
            rraSetupList = glob.glob(os.path.join(rraDir, f'*Setup_RRA_Run_{runSpeed}*.xml'))
            subCycleList = ['cycle'+os.path.split(rraSetupList[ii])[-1].split('cycle')[-1][0] for ii in range(len(rraSetupList))]
        
            #Find the relevant CMC directory
            cmcDir = glob.glob(os.path.join('raw', subject, 'cmc_multipleSteps*', f'CMC_Results_*Run_{runSpeed}*'))[0]
        
            #While we're doing this, extract the cycle timings from the RRA files
            #Here we can also grab the external loads file info
            #Loop through cycles
            for cycle in subCycleList:
                
                #Find the relevant RRA setup file for the current cycle
                rraSetupFile = glob.glob(os.path.join(rraDir, f'*Setup_RRA_Run_{runSpeed}*_{cycle}*.xml'))[0]
                
                #Find the relevant RRA setup file for the current cycle
                cmcSetupFile = glob.glob(os.path.join(cmcDir, f'*Setup_CMC_Run_{runSpeed}*_{cycle}*.xml'))[0]
                
                #Read in RRA setup file as text
                fid = open(rraSetupFile, 'r')
                lines = fid.readlines()
                fid.close()
                
                #Get the relevant data from the text lines
                for li in lines:
                    #Initial time
                    if '<initial_time>' in li:
                        #Get info in between tags
                        stringToGet = re.search('<initial_time>(.*)</initial_time>', li)
                        initialTime = float(stringToGet.group(1).strip(' '))
                    #Final time
                    if '<final_time>' in li:
                        #Get info in between tags
                        stringToGet = re.search('<final_time>(.*)</final_time>', li)
                        finalTime = float(stringToGet.group(1).strip(' '))
                        
                #Store in gait timings dictionary
                gaitTimings[f'run{runSpeed}'][cycleList[subCycleList.index(cycle)]]['initialTime'] = initialTime
                gaitTimings[f'run{runSpeed}'][cycleList[subCycleList.index(cycle)]]['finalTime'] = finalTime
                
                #Read in CMC setup file as text
                fid = open(cmcSetupFile, 'r')
                lines = fid.readlines()
                fid.close()
                        
                #Get the relevant data from the text lines
                #Note this only needs to be done on the first cycle as the same .mot
                #file is used for all cycles
                if subCycleList.index(cycle) == 0:
                    
                    #Identify external loads file
                    for li in lines:
                        #External loads file
                        if '<external_loads_file>' in li:
                            #Get info in between tags
                            stringToGet = re.search('<external_loads_file>(.*)</external_loads_file>', li)
                            externalLoadsFile = stringToGet.group(1).strip(' ')
                        
                    #Search for the GRF filename in the external loads file
                    fid = open(os.path.join(cmcDir, externalLoadsFile), 'r')
                    for li in fid.readlines():
                        if '<datafile>' in li:
                            #Extract out the datafile
                            stringToGet = re.search('<datafile>(.*)</datafile>', li)
                            grfDataFile = os.path.split(stringToGet.group(1))[-1]
                    fid.close()
                    
                    #Read in GRF mot file and write to new file
                    grfTable = osim.TimeSeriesTable(os.path.join('raw', subject, 'ExportedData', grfDataFile))
                    osim.STOFileAdapter().write(grfTable, os.path.join(subject, 'expData', f'Run_{runSpeed}_grf.mot'))
                    
                    #Create associated external loads file
                    extLoads = osim.ExternalLoads()
                    
                    #Create the external forces
                    #Left foot
                    extForceLeft = osim.ExternalForce()
                    extForceLeft.set_applied_to_body('calcn_l')
                    extForceLeft.set_force_expressed_in_body('ground')
                    extForceLeft.set_point_expressed_in_body('ground')
                    extForceLeft.set_force_identifier('L_ground_force_v')
                    extForceLeft.set_point_identifier('L_ground_force_p')
                    extForceLeft.set_torque_identifier('L_ground_torque_')
                    extLoads.cloneAndAppend(extForceLeft)
                    #Right foot
                    extForceRight = osim.ExternalForce()
                    extForceRight.set_applied_to_body('calcn_r')
                    extForceRight.set_force_expressed_in_body('ground')
                    extForceRight.set_point_expressed_in_body('ground')
                    extForceRight.set_force_identifier('R_ground_force_v')
                    extForceRight.set_point_identifier('R_ground_force_p')
                    extForceRight.set_torque_identifier('R_ground_torque_')
                    extLoads.cloneAndAppend(extForceRight)
                    
                    #Set the data file name
                    extLoads.setDataFileName(f'Run_{runSpeed}_grf.mot')
                    
                    #Print to file
                    extLoads.printToXML(os.path.join(subject, 'expData', f'Run_{runSpeed}_grf.xml'))
                
    # =========================================================================
    #     Save gait timings dictionary
    # =========================================================================
    
    #Save pickle file
    with open(os.path.join(subject, 'expData', 'gaitTimes.pkl'), 'wb') as writeFile:
        pickle.dump(gaitTimings, writeFile)
    # with open(os.path.join(subject, 'expData', 'gaitTimes.pkl'), 'rb') as openFile:
    #     loadedDict = pickle.load(openFile)
                
# %% ----- end of 1_extractData.py ----- %% #