# -*- coding: utf-8 -*-
"""

@author:

    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script implements a ground contact identification method. See
    README in this folder for details

"""

# =========================================================================
# Import packages
# =========================================================================

import opensim as osim
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Rectangle
import numpy as np
import os

# =========================================================================
# Set-up
# =========================================================================

# Set matplotlib parameters
# rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['font.weight'] = 'bold'
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 16
rcParams['axes.linewidth'] = 1.5
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False
rcParams['legend.fontsize'] = 10
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.width'] = 1.5
rcParams['legend.framealpha'] = 0.0
rcParams['savefig.dpi'] = 300
rcParams['savefig.format'] = 'pdf'

# Add the utility geometry path for model visualisation
osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.join('utilities', 'Geometry'))

# Read in the model file
osimModel = osim.Model('2014001_C4_05_trackingModel.osim')
osimModel.initSystem()

# Read in the tracking solution and convert to a states trajectory
trackingSolution = osim.MocoTrajectory('2014001_C4_05_trackingSolution.sto')
states = trackingSolution.exportToStatesTrajectory(osimModel)

# Set the contact times for each limb
rightContact = (trackingSolution.getInitialTime(), trackingSolution.getFinalTime())
leftContact = (0.546, trackingSolution.getFinalTime())

# Set the list of markers designating floor level on the models feet
floorMarkers = ['R_FM1_ground', 'R_FM2_ground', 'R_FM5_ground',
                'R_FM1_mid_ground', 'R_FM2_mid_ground', 'R_FM5_mid_ground', 'R_FCC_ground',
                'L_FM1_ground', 'L_FM2_ground', 'L_FM5_ground',
                'L_FM1_mid_ground', 'L_FM2_mid_ground', 'L_FM5_mid_ground', 'L_FCC_ground']

# Create array to get vertical position and velocity of markers from states
markerPos = {marker: np.zeros(trackingSolution.getNumTimes()) for marker in floorMarkers}
markerVel = {marker: np.zeros(trackingSolution.getNumTimes()) for marker in floorMarkers}

# =========================================================================
# Review contact times with respect to marker position
# =========================================================================

# Loop through states and extract marker data
for ii in range(trackingSolution.getNumTimes()):

    # Get the current state
    s = states[ii]

    # Set state in model
    # Probably only need to realize to velocity stage here
    osimModel.realizeVelocity(s)

    # Get the position and velocity for each marker in the ground frame
    for marker in floorMarkers:

        # Get y-axis marker position in ground frame
        markerPos[marker][ii] = osimModel.updBodySet().get(
            osimModel.getMarkerSet().get(marker).getParentFrameName().split('/')[-1]).findStationLocationInGround(
            s, osimModel.getMarkerSet().get(marker).get_location()).get(1)

        # Get marker velocity in ground frame
        # Get the Vec3 velocity component
        v = osimModel.updBodySet().get(
            osimModel.getMarkerSet().get(marker).getParentFrameName().split('/')[-1]).findStationVelocityInGround(
            s, osimModel.getMarkerSet().get(marker).get_location())
        # Calculate the absolute marker velocity
        markerVel[marker][ii] = np.sqrt(v.get(0)**2 + v.get(1)**2 + v.get(2)**2)

# Create a plot that shows the relationship between when ground contact occurs vs. the marker positions and velocity
# -------------------------------------------------------------------------
fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (6,6))

# Plot marker positions for the left and right on the top plots
for marker in floorMarkers:
    if marker.startswith('L_'):
        ax[0, 0].plot(trackingSolution.getTime(), markerPos[marker],
                     c = 'black', lw = 1)
    elif marker.startswith('R_'):
        ax[0, 1].plot(trackingSolution.getTime(), markerPos[marker],
                      c='black', lw=1)

# Plot marker velocities for the left and right on the top plots
for marker in floorMarkers:
    if marker.startswith('L_'):
        ax[1, 0].plot(trackingSolution.getTime(), markerVel[marker],
                      c = 'black', lw = 1, zorder = 2)
    elif marker.startswith('R_'):
        ax[1, 1].plot(trackingSolution.getTime(), markerVel[marker],
                      c = 'black', lw = 1, zorder = 2)

# Add axes titles
ax[0,0].set_title('Left Foot Marker Positions', fontsize = 10, fontweight = 'bold')
ax[0,1].set_title('Right Foot Marker Positions', fontsize = 10, fontweight = 'bold')
ax[1,0].set_title('Left Foot Marker Velocities', fontsize = 10, fontweight = 'bold')
ax[1,1].set_title('Right Foot Marker Velocities', fontsize = 10, fontweight = 'bold')

# Add patches representing when ground contact is occurring
# Left contact on
ax[0,0].add_patch(Rectangle((leftContact[0], 0), leftContact[1]-leftContact[0], ax[0,0].get_ylim()[1],
                  facecolor = 'green', alpha = 0.3, zorder = 1))
ax[1,0].add_patch(Rectangle((leftContact[0], 0), leftContact[1]-leftContact[0], ax[1,0].get_ylim()[1],
                  facecolor = 'green', alpha = 0.3, zorder = 1))
# Left contact off
ax[0,0].add_patch(Rectangle((ax[0,0].get_xlim()[0], 0), leftContact[0]-ax[0,0].get_xlim()[0], ax[0,0].get_ylim()[1],
                  facecolor = 'red', alpha = 0.3, zorder = 1))
ax[1,0].add_patch(Rectangle((ax[1,0].get_xlim()[0], 0), leftContact[0]-ax[0,0].get_xlim()[0], ax[1,0].get_ylim()[1],
                  facecolor = 'red', alpha = 0.3, zorder = 1))
# Right contact on
ax[0,1].add_patch(Rectangle((rightContact[0], 0), rightContact[1]-rightContact[0], ax[0,1].get_ylim()[1],
                  facecolor = 'green', alpha = 0.3, zorder = 1))
ax[1,1].add_patch(Rectangle((rightContact[0], 0), rightContact[1]-rightContact[0], ax[1,1].get_ylim()[1],
                  facecolor = 'green', alpha = 0.3, zorder = 1))
# Right contact off
ax[0,1].add_patch(Rectangle((ax[0,1].get_xlim()[0], 0), rightContact[0]-ax[0,1].get_xlim()[0], ax[0,1].get_ylim()[1],
                  facecolor = 'red', alpha = 0.3, zorder = 1))
ax[1,1].add_patch(Rectangle((ax[1,1].get_xlim()[0], 0), rightContact[0]-ax[1,1].get_xlim()[0], ax[1,1].get_ylim()[1],
                  facecolor = 'red', alpha = 0.3, zorder = 1))

# Set axes limits
ax[0,0].set_xlim([trackingSolution.getTime().to_numpy()[0],trackingSolution.getTime().to_numpy()[-1]])
ax[0,1].set_xlim([trackingSolution.getTime().to_numpy()[0],trackingSolution.getTime().to_numpy()[-1]])
ax[1,0].set_xlim([trackingSolution.getTime().to_numpy()[0],trackingSolution.getTime().to_numpy()[-1]])
ax[1,1].set_xlim([trackingSolution.getTime().to_numpy()[0],trackingSolution.getTime().to_numpy()[-1]])

# Add axes labels
ax[0,0].set_ylabel('Vertical Position (m)', fontsize = 10, fontweight = 'bold')
ax[0,1].set_ylabel('Velocity (m.s)', fontsize = 10, fontweight = 'bold')
ax[1,0].set_ylabel('Vertical Position (m)', fontsize = 10, fontweight = 'bold')
ax[1,1].set_ylabel('Velocity (m.s)', fontsize = 10, fontweight = 'bold')

# Tight layout
plt.tight_layout()

# =========================================================================
# Estimate ground body contact based on point position and velocity
# =========================================================================

# Set the points to consider for each body
bodyPoints = {
    'right': ['R_FM1_ground', 'R_FM2_ground', 'R_FM5_ground',
              'R_FM1_mid_ground', 'R_FM2_mid_ground', 'R_FM5_mid_ground', 'R_FCC_ground'],
    'left': ['L_FM1_ground', 'L_FM2_ground', 'L_FM5_ground',
             'L_FM1_mid_ground', 'L_FM2_mid_ground', 'L_FM5_mid_ground', 'L_FCC_ground']
}

# Set the distance and velocity thresholds to consider
# Karcnik suggests 8cm and 1.5cm/sec --- but unsure if the units are correct (feel like it's 1.5m/sec?)
distThreshold = 0.08  # 8 cm
velThreshold = 1.5  # 1.5 m/sec

# Loop through each time index and determine whether any markers for body are under threshold
# A point will be added at the bottom of the relevant plot if contact is identified
for ii in range(trackingSolution.getNumTimes()):

    # Left contact body
    # -------------------------------------------------------------------------

    # Loop through markers and check if they meet thresholds
    inContact = [markerPos[marker][ii] < distThreshold and markerVel[marker][ii] < velThreshold for marker in bodyPoints['left']]

    # Check if any are in contact
    if any(inContact):

        # Add a point to plots
        ax[0,0].scatter(trackingSolution.getTime().to_numpy()[ii], 0, marker = 's', fc = 'blue', s = 10, zorder = 3)
        ax[1,0].scatter(trackingSolution.getTime().to_numpy()[ii], 0, marker='s', fc='blue', s=10, zorder=3)

    # Right contact body
    # -------------------------------------------------------------------------

    # Loop through markers and check if they meet thresholds
    inContact = [markerPos[marker][ii] < distThreshold and markerVel[marker][ii] < velThreshold for marker in bodyPoints['right']]

    # Check if any are in contact
    if any(inContact):

        # Add a point to plots
        ax[0,1].scatter(trackingSolution.getTime().to_numpy()[ii], 0, marker='s', fc='blue', s=10, zorder=3)
        ax[1,1].scatter(trackingSolution.getTime().to_numpy()[ii], 0, marker='s', fc='blue', s=10, zorder=3)

"""

The above works OK, and the thresholds probably change with different
scenarios (e.g. gait speeds, tasks etc.). There are some slight errors
around foot contact, but it isn't too bad. THese could be addressed by
changing the thresholds.

"""

# =========================================================================
# Combine ground contact and ZMP methods to estimate GRFs and CoP
# =========================================================================

# Read in model - note that the external loads were removed manually
osimModel = osim.Model('2014001_C4_05_trackingModel.osim')
osimModel.initSystem()

# Initialise an inverse dynamics solver with model
ivdSolver = osim.InverseDynamicsSolver(osimModel)

# Extract accelerations from Moco solution
trackingSolution.generateAccelerationsFromSpeeds()
udotTable = trackingSolution.exportToAccelerationsTable()
osim.STOFileAdapter().write(udotTable, 'accelerations.sto')

# Export states trajectory from Moco solution
states = trackingSolution.exportToStatesTrajectory(osimModel)
# states = trackingSolution.exportToStatesTable()

# Define number of times
# nt = trackingSolution.getNumTimes()
nt = states.getSize()
# nt = states.getNumRows()

# Set contact body names
contactBodyNames = osim.ArrayStr()
contactBodyNames.append('calcn_r')
contactBodyNames.append('calcn_l')

# Create a storage object for the calculated forces
zmpResults = osim.Storage(nt)
# zmpResults = osim.TimeSeriesTable()

# Create the columns for storing force data based on the body names
# Necessary order is Fx, Fy, Fz, Px, Py, Pz, Mx, My, Mz
zmpLabels = osim.ArrayStr('time', 9 * contactBodyNames.size() + 1)
for nb in range(contactBodyNames.size()):
    zmpLabels.set(nb * 9 + 1, contactBodyNames.get(nb) + '_force_vx')
    zmpLabels.set(nb * 9 + 2, contactBodyNames.get(nb) + '_force_vx')
    zmpLabels.set(nb * 9 + 3, contactBodyNames.get(nb) + '_force_vx')
    zmpLabels.set(nb * 9 + 4, contactBodyNames.get(nb) + '_force_px')
    zmpLabels.set(nb * 9 + 5, contactBodyNames.get(nb) + '_force_py')
    zmpLabels.set(nb * 9 + 6, contactBodyNames.get(nb) + '_force_pz')
    zmpLabels.set(nb * 9 + 7, contactBodyNames.get(nb) + '_torque_x')
    zmpLabels.set(nb * 9 + 8, contactBodyNames.get(nb) + '_torque_y')
    zmpLabels.set(nb * 9 + 9, contactBodyNames.get(nb) + '_torque_z')

# Set column labels in ZMP storage
zmpResults.setColumnLabels(zmpLabels)

# Set the free joint and force threshold parameters
freeJoint = 'ground_pelvis'
forceThreshold = 20

# Set the points to consider for each body
# Right side
rightPointNames = osim.ArrayStr()
rightPointNames.append('R_FM1_ground')
rightPointNames.append('R_FM2_ground')
rightPointNames.append('R_FM5_ground')
rightPointNames.append('R_FM1_mid_ground')
rightPointNames.append('R_FM2_mid_ground')
rightPointNames.append('R_FM5_mid_ground')
rightPointNames.append('R_FCC_ground')
# Left side
leftPointNames = osim.ArrayStr()
leftPointNames.append('L_FM1_ground')
leftPointNames.append('L_FM2_ground')
leftPointNames.append('L_FM5_ground')
leftPointNames.append('L_FM1_mid_ground')
leftPointNames.append('L_FM2_mid_ground')
leftPointNames.append('L_FM5_mid_ground')
leftPointNames.append('L_FCC_ground')

# Build the list of joints for computing and reporting equivalent body forces

# Set bodies to report forces for
jointsForReportingBodyForces = osim.ArrayStr()
jointsForReportingBodyForces.append(freeJoint)

# Create jointset for calculating equivalent body forces
jointsForEquivalentBodyForces = osim.JointSet()

# Get model joints
modelJoints = osimModel.getJointSet()

# Loop through reporting joint names, clone and append the corresponding model joint
for ij in range(jointsForReportingBodyForces.getSize()):
    # Note that there is an option for "ALL" in C++ code
    k = modelJoints.getIndex(jointsForReportingBodyForces.get(ij))
    if k >= 0:
        jointsForEquivalentBodyForces.adoptAndAppend(modelJoints.get(k))
    else:
        print('Tool could not find joint named %s for reporting body forces...' % jointsForReportingBodyForces.get(ij))

# Get number of joints for reporting body forces
nj = jointsForEquivalentBodyForces.getSize()

# Create column labels for body forces storage file
# Starts with time plus the space for six components per joint to report on
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

# Define the number of coordinates in the model
nq = osimModel.updCoordinateSet().getSize()

# Set the distance and velocity thresholds to consider for contact
distThreshold = 0.05  # Original 8 cm seemed too high. Even 5cm seems a little high
velThreshold = 1.5

# Loop through states and calculate forces
for ii in range(nt):

    # Get the current state
    s = states[ii]

    # Realize to accelerations stage
    osimModel.realizeAcceleration(s)

    # Get accelerations from table
    # Note that using the Moco trajectory still seems to resolve things better than calculating udot manually
    # Note that it needs to be a vector instead of row vector
    udot = osim.Vector().createFromMat(udotTable.getRowAtIndex(ii).to_numpy())

    # Use the ID solver with the current state
    genForceTraj = ivdSolver.solve(s, udot)

    # The above results contain the generalised coordinate forces to generate the accelerations
    # based on the current state. Note that these aren't necessarily in the order of the coordinate
    # set, but rather the multibody tree order
    # TODO: mapping if q != u in index (i.e. tree vs. model)...

    # Next we want to calculate the equivalent body forces related to these

    # Create a state vector to store the forces
    # Note that this allocates the time from the current state and 6 slots for the
    # force and moment components for each body (with a zero allocated for each)
    forcesVec = osim.Vector(6 * nj, 0)

    # Calculate equivalent body force at joint for those listed (i.e. the free joint)
    for j in range(nj):
        equivalentBodyForceAtJoint = jointsForEquivalentBodyForces.get(j).calcEquivalentSpatialForce(s, genForceTraj)
        # Extract the Vec 3 components
        for k in range(3):
            # Body force components
            forcesVec.set(6 * j + k, equivalentBodyForceAtJoint.get(1)[k])
            # Body torque components
            forcesVec.set(6 * j + k + 3, equivalentBodyForceAtJoint.get(0)[k])

    # Get the body forces acting at the listed free joint
    # TODO: this is manually indexed --- which isn't great for use across different models...
    freeJointF = osim.Vec3(forcesVec.get(0), forcesVec.get(1), forcesVec.get(2))

    # Check ground contact points
    # -------------------------------------------------------------------------

    # Create a vector to store if body is in contact
    bodyInContact = osim.Vector(contactBodyNames.size(),0)

    # Note that we could traditionally loop through the bodies here based on model
    # Here a manual approach is used to do the right and them left body
    # Similarly the points wouldn't be markers, but rather contact points in model component

    # Right body
    # -------------------------------------------------------------------------

    # Create a variable to check contact of points
    inContact1 = osim.Vector(rightPointNames.size(), 0)

    # Loop through contact points
    for cp in range(rightPointNames.size()):

        # Get point name
        pointName = rightPointNames.get(cp)

        # Get body marker is attached to
        bodyName = osimModel.getMarkerSet().get(pointName).getParentFrameName().split('/')[-1]

        # Get point location
        pointLocation = osimModel.getMarkerSet().get(pointName).get_location()

        # Get vertical point position in ground frame
        vertPointPos = osimModel.updBodySet().get(bodyName).findStationLocationInGround(s, pointLocation).get(1)

        # Get marker velocity in ground frame
        v = osimModel.updBodySet().get(bodyName).findStationVelocityInGround(s, pointLocation)
        pointVel = np.sqrt(v.get(0) ** 2 + v.get(1) ** 2 + v.get(2) ** 2)

        # Check if it meets distance and velocity criteria
        if vertPointPos < distThreshold and pointVel < velThreshold:
            inContact1.set(cp, 1)

    # Check if sum of contact vector is greater than zero
    # If it is then the body is deemed in contact
    if inContact1.sum() > 0:
        bodyInContact.set(0,1)  # TODO: manual allocation here

    # Left body
    # -------------------------------------------------------------------------

    # Create a variable to check contact of points
    inContact2 = osim.Vector(leftPointNames.size(), 0)

    # Loop through contact points
    for cp in range(leftPointNames.size()):

        # Get point name
        pointName = leftPointNames.get(cp)

        # Get body marker is attached to
        bodyName = osimModel.getMarkerSet().get(pointName).getParentFrameName().split('/')[-1]

        # Get point location
        pointLocation = osimModel.getMarkerSet().get(pointName).get_location()

        # Get vertical point position in ground frame
        vertPointPos = osimModel.updBodySet().get(bodyName).findStationLocationInGround(s, pointLocation).get(1)

        # Get marker velocity in ground frame
        v = osimModel.updBodySet().get(bodyName).findStationVelocityInGround(s, pointLocation)
        pointVel = np.sqrt(v.get(0) ** 2 + v.get(1) ** 2 + v.get(2) ** 2)

        # Check if it meets distance and velocity criteria
        if vertPointPos < distThreshold and pointVel < velThreshold:
            inContact2.set(cp, 1)

    # Check if sum of contact vector is greater than zero
    # If it is then the body is deemed in contact
    if inContact2.sum() > 0:
        bodyInContact.set(1, 1)  # TODO: manual allocation here

    # Calculate ZMP results for current state
    # -------------------------------------------------------------------------

    # Set a vector to store the state results in
    # Default values to zero
    zmpVec = osim.Vector(9 * contactBodyNames.size(), 0)

    # TODO: note that these don't need to actually be variables --- the vectors could just be repeatedly referenced

    # Check if any contact has been identified
    if bodyInContact.sum() > 0:
        contactFound = True
    else:
        contactFound = False

    # Check if unilateral contact has been found
    if contactFound:
        if bodyInContact.sum() == 1:
            unilateralContact = True
        else:
            unilateralContact = False

    # Do a check to see if contact is found and the force threshold has been met
    # Assumes vertical force is in the y-axis
    # -------------------------------------------------------------------------
    if contactFound and freeJointF.get(1) > forceThreshold:

        # Get the position of the free joint in the ground frame at the current state
        rp = osimModel.getJointSet().get(freeJoint).getChildFrame().findStationLocationInGround(s, osim.Vec3(0, 0, 0))

        # Take cross product of free joint body position and force vector to get moment at origin
        groundM = osim.Vec3((rp.get(1) * freeJointF.get(2)) - (rp.get(2) * freeJointF.get(1)),
                            -((rp.get(0) * freeJointF.get(2)) - (rp.get(2) * freeJointF.get(0))),
                            (rp.get(0) * freeJointF.get(1)) - (rp.get(1) * freeJointF.get(0)))

        # Calculate X & Z cZMP, noting that yZMP is set as 0
        # Formulas come from Xiang et al.
        zmpCOP = osim.Vec3(groundM.get(2) / freeJointF.get(1),
                           0,
                           -groundM.get(0) / freeJointF.get(1))

        # Calculate the active moment at the ZMP along the y-axis
        # See equation 20 in Xiang et al.
        # TODO: still unsure how correct this is?
        zmpMy = groundM.get(1) + (freeJointF.get(0) * zmpCOP.get(2)) - (freeJointF.get(2) * zmpCOP.get(1))

        # Check for unilateral contact option
        # If this is the case then all force gets allocated to the contacting body
        if unilateralContact:

            # Figure out which body is in contact to allocate the forces and moments to
            for bb in range(bodyInContact.size()):
                if bodyInContact.get(bb) == 1:
                    contactInd = bb
                    break

            # Set the values at the appropriate indices in the ZMP vector
            # Note standardised order of Fx, Fy, Fz, Px, Py, Pz, Mx, My, Mz
            # ZMP assumes zero for Mx and Mz, so only My is allocated
            # TODO: this uses the ground moment at the origin and not the centre of pressure --- how to convert to CoP location?
            zmpVec.set(contactInd * 9 + 0, freeJointF.get(0))
            zmpVec.set(contactInd * 9 + 1, freeJointF.get(1))
            zmpVec.set(contactInd * 9 + 2, freeJointF.get(2))
            zmpVec.set(contactInd * 9 + 3, zmpCOP.get(0))
            zmpVec.set(contactInd * 9 + 4, zmpCOP.get(1))
            zmpVec.set(contactInd * 9 + 5, zmpCOP.get(2))
            # zmpVec.set(contactInd * 9 + 6, groundM.get(0))
            zmpVec.set(contactInd * 9 + 7, zmpMy)
            # zmpVec.set(contactInd * 9 + 8, groundM.get(2))

        else:

            # Split the contact between the bodies
            # -------------------------------------------------------------------------
            # TODO: this currently assumes 2 bodies and perhaps needs to be specified as the upper limit of tool
            # TODO: this is probably also only appropriate for straight line gait?

            # Right body
            # -------------------------------------------------------------------------

            # Get the positions of the points in contact with the ground
            contactLoc1 = osim.ArrayVec3()

            # Loop through points
            for ptInd in range(inContact1.size()):

                # Check if point is in contact
                if inContact1.get(ptInd) == 1:

                    # Get point name
                    pointName = rightPointNames.get(ptInd)

                    # Get body marker is attached to
                    bodyName = osimModel.getMarkerSet().get(pointName).getParentFrameName().split('/')[-1]

                    # Get point location in ground frame
                    pointLocation = osimModel.updBodySet().get(bodyName).findStationLocationInGround(
                        s, osimModel.getMarkerSet().get(pointName).get_location())

                    # Append to vector
                    contactLoc1.append(pointLocation)

            # Get the average position of the extracted points
            # TODO: this is a Python way to get it. Will need a C++ alternative

            # Average the X, Y and Z elements
            xMu = np.array([contactLoc1.get(cInd).get(0) for cInd in range(contactLoc1.size())]).mean()
            # yMu = np.array([contactLoc1.get(cInd).get(1) for cInd in range(contactLoc1.size())]).mean()
            zMu = np.array([contactLoc1.get(cInd).get(2) for cInd in range(contactLoc1.size())]).mean()

            # Put average location into Vec3 format
            # Assumes a zero y location
            avgContact1 = osim.Vec3(xMu, 0, zMu)

            # Left body
            # -------------------------------------------------------------------------

            # Get the positions of the points in contact with the ground
            contactLoc2 = osim.ArrayVec3()

            # Loop through points
            for ptInd in range(inContact2.size()):

                # Check if point is in contact
                if inContact2.get(ptInd) == 1:

                    # Get point name
                    pointName = leftPointNames.get(ptInd)

                    # Get body marker is attached to
                    bodyName = osimModel.getMarkerSet().get(pointName).getParentFrameName().split('/')[-1]

                    # Get point location in ground frame
                    pointLocation = osimModel.updBodySet().get(bodyName).findStationLocationInGround(
                        s, osimModel.getMarkerSet().get(pointName).get_location())

                    # Append to vector
                    contactLoc2.append(pointLocation)

            # Get the average position of the extracted points
            # TODO: this is a Python way to get it. Will need a C++ alternative

            # Average the X, Y and Z elements
            xMu = np.array([contactLoc2.get(cInd).get(0) for cInd in range(contactLoc2.size())]).mean()
            # yMu = np.array([contactLoc1.get(cInd).get(1) for cInd in range(contactLoc1.size())]).mean()
            zMu = np.array([contactLoc2.get(cInd).get(2) for cInd in range(contactLoc2.size())]).mean()

            # Put average location into Vec3 format
            # Assumes a zero y location
            avgContact2 = osim.Vec3(xMu, 0, zMu)

            # Calculate distances between average contact points to ZMP
            # Note this is calcualated as 2D distance as y is always equal to 0
            # -------------------------------------------------------------------------

            # Contact body one to ZMP
            d1 = np.sqrt(((zmpCOP.get(0) - avgContact1.get(0))**2) + ((zmpCOP.get(2) - avgContact1.get(2))**2))

            # Contact body two to ZMP
            d2 = np.sqrt(((zmpCOP.get(0) - avgContact2.get(0)) ** 2) + ((zmpCOP.get(2) - avgContact2.get(2)) ** 2))

            # TODO: up to here --- how to partition GRF between the two points
            # (d2 / (d1+d2)) * freeJointF.to_numpy()
            # (d1 / (d1 + d2)) * freeJointF.to_numpy()

    # Create a state vector with the time and ZMP values
    # Note this will remain as zeros if no force allocated
    zmpStateVec = osim.StateVector(s.getTime(), zmpVec)

    # Append the current state vector to the ZMP results
    zmpResults.append(zmpStateVec)

# Write the ZMP results to file
# -------------------------------------------------------------------------

# Set name in ZMP storage
zmpResults.setName('ZMP Predicted Ground Reactions')

# # Set column labels in ZMP storage
# zmpResults.setColumnLabels(zmpLabels)

# Write to file
osim.Storage().printResult(zmpResults, 'zmpPredictedGroundReactions', '.', -1, '.sto')



# %% ---------- end of groundContact.py ---------- %% #
