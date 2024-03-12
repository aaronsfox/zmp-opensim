# -*- coding: utf-8 -*-
"""

@author: 
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    Build of functions to assist with running simulations
    
"""


import opensim as osim

# %% Function to add set of torque actuators to model

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

# %% Function to convert IK coordinates to states

def kinematicsToStates(kinematicsFileName = None, osimModelFileName = None,
                       outputFileName = 'coordinates.sto',
                       inDegrees = True, outDegrees = False,
                       filtFreq = None):
    
    # Convenience function for converting IK results to a states storage.
    #
    # Input:    kinematicsFileName - file containing kinematic data. Header should only be coordinates name, rather than path to state
    #           osimModelFileName - opensim model filename that corresponds to kinematic data
    #           outputFileName - optional filename to output to (defaults to coordinates.sto)
    #           inDegrees - set to true if kinematics file is in degrees (defaults to True)
    #           outDegrees - set to true if desired output is in degrees (defaults to False)
    #           filtFreq - lowpass filter frequency for kinematics

    if kinematicsFileName is None:
        raise ValueError('Filename for kinematics is required')
    if osimModelFileName is None:
        raise ValueError('OpenSim model filename is required')

    #Load in the kinematic data
    kinematicsStorage = osim.Storage(kinematicsFileName)
    
    #Create a copy of the kinematics data to alter the column labels in
    statesStorage = osim.Storage(kinematicsFileName)
    
    #Check to filter data
    if filtFreq is not None:
        #Filter both storage objects
        #Note that this resamples time stamps so eliminates the need to do so
        kinematicsStorage.lowpassFIR(4, filtFreq)
        kinematicsStorage.lowpassFIR(4, filtFreq)
    else:
        #Resample the data points linearly to avoid any later issues with matching
        #time points. Use a time stamp for 250 Hz
        kinematicsStorage.resampleLinear(1/250)
        statesStorage.resampleLinear(1/250)
    
    #Get the column headers for the storage file
    angleNames = kinematicsStorage.getColumnLabels()
    
    #Get the corresponding full paths from the model to rename the
    #angles in the kinematics file
    kinematicModel = osim.Model(osimModelFileName)
    for ii in range(angleNames.getSize()):
        currAngle = angleNames.get(ii)
        if currAngle != 'time':
            #Try getting the full path to coordinate
            #This may fail due to their being marker data included in these files
            try:
                #Loo for full coordinate path
                fullPath = kinematicModel.updCoordinateSet().get(currAngle).getAbsolutePathString()+'/value'
                #Set angle name appropriately using full path
                angleNames.set(ii,fullPath)
            except:
                #Print out that current column isn't a coordinate
                print(f'{currAngle} not a coordinate...skipping name conversion...')
                #Set to the same as originaly
                angleNames.set(ii,currAngle)
    
    #Set the states storage object to have the updated column labels
    statesStorage.setColumnLabels(angleNames)
    
    #Appropriately set output in degrees or radians
    if inDegrees and not outDegrees:
        #Convert degrees values to radians for consistency with the current
        #file label (defaults back to inDegrees=no). Radians seem to work
        #better with the Moco process as well.
        kinematicModel.initSystem()
        kinematicModel.getSimbodyEngine().convertDegreesToRadians(statesStorage)
    elif inDegrees and outDegrees:
        #Change the storage label back to specifying indegrees=yes
        statesStorage.setInDegrees(True)
    elif not inDegrees and outDegrees:
        #Convert radians to degrees
        kinematicModel.initSystem()
        kinematicModel.getSimbodyEngine().convertRadiansToDegrees(statesStorage)
        #Reset labeling for degrees
        statesStorage.setInDegrees(True)
    
    #Write the states storage object to file
    statesStorage.printToXML(outputFileName)
    
# %% Function to convert states kinematics back to standard kinematics
    
def statesToKinematics(statesFileName = None,
                       outputFileName = 'coordinates.sto',
                       inDegrees = False, outDegrees = True):
    
    # Convenience function for converting IK results to a states storage.
    #
    # Input:    statesFileName - file containing kinematic data. Header should only be coordinates name, rather than path to state
    #           osimModelFileName - opensim model filename that corresponds to kinematic data
    #           outputFileName - optional filename to output to (defaults to coordinates.sto)
    #           inDegrees - set to true if kinematics file is in degrees (defaults to False)
    #           outDegrees - set to true if desired output is in degrees (defaults to True)

    if statesFileName is None:
        raise ValueError('Filename for states is required')

    #Load in the states data as a table
    statesTable = osim.TimeSeriesTable(statesFileName)
    
    # #Create a copy of the states data to alter the column labels in
    # kinematicsStorage = osim.Storage(statesFileName)
    
    #Get the column headers for the states file
    stateNames = statesTable.getColumnLabels()
    
    #Loop through the column names and identify which columns need to be removed
    #given they aren't kinematic values
    for state in stateNames:
        if '/value' not in state:
            statesTable.removeColumn(state)
    
    #Create new set of column labels that removes all but the coordinate
    #This is based off the fact that kinematic states are presented as:
    # /jointset/joint_name/coordinate_name/value
    
    #Get the selected state names and create a string array to store names in
    selectStateNames = statesTable.getColumnLabels()
    newStateNames = osim.StdVectorString()
    
    #Loop through state names to alter
    for state in selectStateNames:
        #Split the string by the / and get the 4th output
        #Append this to the vector string
        newStateNames.append(state.split('/')[3])
        
    #Set new names in table
    statesTable.setColumnLabels(newStateNames)
    
    # #Appropriately set output in degrees or radians
    # if inDegrees and not outDegrees:
    #     #Convert degrees values to radians for consistency with the current
    #     #file label (defaults back to inDegrees=no). Radians seem to work
    #     #better with the Moco process as well.
    #     kinematicModel.initSystem()
    #     kinematicModel.getSimbodyEngine().convertDegreesToRadians(statesStorage)
    # elif inDegrees and outDegrees:
    #     #Change the storage label back to specifying indegrees=yes
    #     statesStorage.setInDegrees(True)
    # elif not inDegrees and outDegrees:
    #     #Convert radians to degrees
    #     kinematicModel.initSystem()
    #     kinematicModel.getSimbodyEngine().convertRadiansToDegrees(statesStorage)
    #     #Reset labeling for degrees
    #     statesStorage.setInDegrees(True)
    
    #Write to file
    osim.STOFileAdapter().write(statesTable, outputFileName)

# %% ----- End of osimFunctions.py -----