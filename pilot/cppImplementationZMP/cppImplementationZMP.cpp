/* -------------------------------------------------------------------------- *
 * OpenSim: cppImplementationZMP.cpp                                          *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2017-19 Stanford University and the Authors                  *
 *                                                                            *
 * Author(s): Aaron Fox		                                                  *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0          *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

 /// This example estimates the GRFs using the Zero Moment Point method for a
 /// running gait cycle. It is an example implementation of the required
 /// calculations in C++ as a means to transfer into a probe implementation.
 
 /// NOTE: zeros come out when only realizing to acceleration stage, values when going to dynamics
 /// Despite getting values, they are different to Python implementation --- why???
 /// Perhaps use of s instead of sWorkingCopy in calculations?
 /// Seems like it stems from the ivdSolver given that the body forces are different in the first place?

#include <OpenSim/OpenSim.h>
#include <OpenSim/Simulation/InverseDynamicsSolver.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimulationUtilities.h>
#include <OpenSim/Moco/MocoTrajectory.h>
#include <OpenSim/Common/IO.h>

using namespace OpenSim;
using namespace std;
using namespace SimTK;

int main() {

    // ----- SETTINGS ----- //
    //* TODO: these could be probe inputs at some point? *//

    // Set vertical force threshold
    int forceThreshold = 50;

    // Set right and left body ground contact names
    string rightBodyName = "calcn_r";
    string leftBodyName = "calcn_l";

    // ----- DEFINE MODEL ----- //
    cout << "Loading model...\n";

    // Import model geometries to avoid warning print-outs --- full path needed?
    //ModelVisualizer.addDirToGeometrySearchPaths("C:\\OpenSim 4.4\\Geometry");

	// Define the model --- need to specify full path for some reason...?
	//Model osimModel("osimModel.osim");
    Model osimModel("C:\\+GitRepos+\\zmp-opensim\\pilot\\cppImplementationZMP\\osimModel.osim");

 	// Initialize the model's underlying computational system and get its default state.
    State& sWorkingCopy = osimModel.initSystem();

    // Get number of forces in model
    int nf = osimModel.updForceSet().getSize();

    // OpenSim::Coordinates represent degrees of freedom for a model.
        // Each Coordinate's value and speed maps to an index
        // in the model's underlying SimTK::State (value to a slot in the
        // State's q, and speed to a slot in the State's u).
        // So we need to map each OpenSim::Coordinate value and speed to the
        // corresponding SimTK::State's q and u indices, respectively.
    auto coords = osimModel.getCoordinatesInMultibodyTreeOrder();
    int nq = sWorkingCopy.getNQ();
    int nu = sWorkingCopy.getNU();
    int nCoords = (int)coords.size(); // TODO: needed?
    // int intUnusedSlot = -1; // TODO: needed?

    // Initialise an inverse dynamics solver with model
    InverseDynamicsSolver ivdSolver(osimModel);

    // Finish section
    cout << "Model loaded...\n";

    // ----- SET-UP STATES AND STORAGE FILES ----- //
    cout << "Loading files and setting storage...\n";

    // Get model states from Moco solution
    MocoTrajectory mocoTraj = MocoTrajectory("C:\\+GitRepos+\\zmp-opensim\\pilot\\cppImplementationZMP\\subject01_run5_cycle1_mocoSolution.sto");
    StatesTrajectory statesTraj = mocoTraj.exportToStatesTrajectory(osimModel);

    //// In lieu of udot seemingly not working later, create an accelerations table
    //mocoTraj.generateAccelerationsFromSpeeds();
    //TimeSeriesTable uDotTable = mocoTraj.exportToAccelerationsTable();

    // Define number of times from states trajectory
    int nt = mocoTraj.getNumTimes();

    // Create storage objects for gen and body forces to append data to
    Storage genForcesResults = Storage(nt);
    Storage bodyForcesResults = Storage(nt);

    // Create coordinate labels for gen forces
    Array<string> genForceLabels("time", nq + 1);
    for (int iq = 0; iq < nq; ++iq) {
        genForceLabels.set(iq + 1, osimModel.updCoordinateSet().get(iq).getName());
    }
    
    // Finish section
    cout << "Files and storage set...\n";

    // ----- SET REPORTING OF BODY & JOINT FORCES --- //
    cout << "Setting body and joint forces to report...\n";

    //// Create array string for reporting body forces
    //Array<string> jointsForReportingBodyForces("ground_pelvis", 1);

    // Set the free joint to report body forces for
    string freeJoint = "ground_pelvis";

    //// Create jointset for calculating equivalent body forces
    //JointSet jointsForEquivalentBodyForces; // TODO: needed...

    //// Get the model joints
    //JointSet modelJoints = osimModel.getJointSet();

    // Loop through desired body forces and append for reporting

    /* NOTE: there is a problem here, whereby the jointsForEquivalentBodyForces variable is not appending a joint...*/

    //for (int i = 0; i < jointsForReportingBodyForces.getSize(); ++i) {

        /*Note that an 'ALL' option could be included here*/

        //// Get index for joint in model
        //int k = modelJoints.getIndex(jointsForReportingBodyForces[i]);

        //// Adopt and append if joint identified
        //if (k >= 0) {
        //    jointsForEquivalentBodyForces.adoptAndAppend(&modelJoints[k]);
        //}
        //else {
        //    log_warn("Could not find Joint named '{}' to report body forces.", jointsForReportingBodyForces[i]);
        //}
    //}

    // Get number of joints for reporting (should be just 1 given above implementation)
    //int nj = jointsForEquivalentBodyForces.getSize();
    int nj = 1; // manually set

    // Create column labels for body force storage file
    // Starts with time plus the space for six components per joint to report on
    Array<string> bodyForceLabels("time", 6 * nj + 1);
    string XYZ = "XYZ";
    string joint_body_label = osimModel.getJointSet().get(freeJoint).getName() + "_";
    joint_body_label += osimModel.getJointSet().get(freeJoint).getChildFrame().getName(); // does this still cause error?
    for (int k = 0; k < 3; ++k) {
        bodyForceLabels[k + 1] = joint_body_label + "_F" + XYZ[k]; //first label is time
        bodyForceLabels[k + 3 + 1] = joint_body_label + "_M" + XYZ[k]; //first label is time
    }

    // Finish section
    cout << "Body and joint forces for reporting set...\n";

    // ----- CREATE STORAGE FOR ZMP CALCULATIONS ----- //
    cout << "Setting storage for ZMP calculations...\n";

    // Create storage file to append ZMP GRFs to
    Storage zmpResults = Storage(nt);

    // Create columns for ZMP force and COP calculations
    // These are currently partitioned to left and right sides as per
    // usual external loads data
    Array<string> zmpLabels("time", 6 * 2 + 1);
    zmpLabels.set(1, "R_ground_force_vx");
    zmpLabels.set(2, "R_ground_force_vy");
    zmpLabels.set(3, "R_ground_force_vz");
    zmpLabels.set(4, "R_ground_force_px");
    zmpLabels.set(5, "R_ground_force_py");
    zmpLabels.set(6, "R_ground_force_pz");
    zmpLabels.set(7, "L_ground_force_vx");
    zmpLabels.set(8, "L_ground_force_vy");
    zmpLabels.set(9, "L_ground_force_vz");
    zmpLabels.set(10, "L_ground_force_px");
    zmpLabels.set(11, "L_ground_force_py");
    zmpLabels.set(12, "L_ground_force_pz");

    // Finish section
    cout << "ZMP storage set...\n";

    // ----- LOOP THROUGH STATES AND CALCULATE BODY FORCES --- //
    cout << "Calculating forces at individual states...\n";

    // Loop through times
    for (int i = 0; i < nt; i++) {

        // Get the current state
        State s = statesTraj[i];

        // Set the current states in the working copy
        sWorkingCopy.setTime(s.getTime());
        sWorkingCopy.setQ(s.getQ());
        sWorkingCopy.setU(s.getU());

        // Overide actuation of forces
        for (int j = 0; j < nf; j++) {
            ScalarActuator* act = dynamic_cast<ScalarActuator*>(&osimModel.updForceSet().get(j));
            if (act) {
                act->setOverrideActuation(sWorkingCopy, 1);
            }
        }

        // Realize to accelerations stage
        //osimModel.getMultibodySystem().realize(sWorkingCopy, Stage::Acceleration);
        osimModel.getMultibodySystem().realize(sWorkingCopy, Stage::Dynamics);

        // Compute accelerations of current state
        Vector udot = osimModel.getMatterSubsystem().getUDot(sWorkingCopy);

        // Solve inverse dynamics given current states and udot
        // The output vector contains the generalised coordinate forces
        // to generate the accelerations based on the current state.
        // Note that these aren't necessarily in the order of the
        // coordinate set, but rather the multibody tree order.
        Vector genForceTraj = ivdSolver.solve(sWorkingCopy, udot);

        // Create a state vector to store the forces
        // Note that this allocates the time from the current state and 6 slots for the
        // force and moment components for each body (with a zero allocated for each)
        Vector genVec = Vector(nq); // TODO: used?
        Vector forcesVec = Vector(6 * nj);

        ///* TODO: mapping if q != u in index(i.e.tree vs.model)... */

        // Calculate the equivalent body force at joint for those listed
        for (int j = 0; j < nj; j++) {
            
            // Calculate the body force at the current joint
            SpatialVec equivalentBodyForceAtJoint = osimModel.getJointSet().get(freeJoint).calcEquivalentSpatialForce(s, genForceTraj); // problems...

            // Extract the Vec3 components and set in forces vector
            for (int k = 0; k < 3; k++) {
                // Body force component
                forcesVec.set(6 * j + k, equivalentBodyForceAtJoint.get(1)[k]);
                //Body torque component
                forcesVec.set(6 * j + k + 3, equivalentBodyForceAtJoint.get(0)[k]);
            }

        }

        // Convert to state vectors
        StateVector genForcesVec = StateVector(s.getTime(), genForceTraj);
        StateVector bodyForcesVec = StateVector(s.getTime(), forcesVec);

        // Append current state vectors to results storage
        genForcesResults.append(genForcesVec);
        bodyForcesResults.append(bodyForcesVec);

        // ----- CALCULATE ZMP ESTIMATED GRFs & COPs FOR CURRENT STATE ----- //

        // Set vector to store results for current state (6 values per limb)
        Vector zmpVec = Vector(6 * 2);

        // Get the body forces acting at the free body (generally pelvis)
        //* TODO: this is manually indexed to get pelvis --- needs to be fixed for different models... *//
        Vec3 freeBodyF = Vec3(forcesVec.get(0), forcesVec.get(1), forcesVec.get(2));

        // Check if vertical force is greater than threshold to perform calculations
        // There shouldn't be a need to change the non-used side as they should be zeros
        //* TODO: this should account for different vertical axes... *//
        if ( freeBodyF.get(1) > forceThreshold ) {

            // Get th position of the free body (generally pelvis) from Q states
            //* TODO: this is manually indexed to get pelvis --- needs to be fixed for different models... *//
            Vec3 rp = Vec3(s.getQ().get(3), s.getQ().get(4), s.getQ().get(5));
            
            // Take cross product of free body position and force vector to get moment at origin
            Vec3 groundM = Vec3((rp.get(1) * freeBodyF.get(2)) - (rp.get(2) * freeBodyF.get(1)),
                -((rp.get(0) * freeBodyF.get(2)) - (rp.get(2) * freeBodyF.get(0))),
                (rp.get(0) * freeBodyF.get(1)) - (rp.get(1) * freeBodyF.get(0)));

            // Calculate force at ground origin (should just equal force at free body given already in ground system)
            //* TODO: should just be free body force if expressed in ground --- perhaps need to check this...? *//

            // Calculate X & Z cZMP, noting that yZMP is set as 0
            // Formulas come from Xiang et al. (2009), Int J Numer Meth ENg, 79: 667-695.
            Vec3 zmpCOP = Vec3(groundM.get(2) / freeBodyF.get(1), 0, -groundM.get(0) / freeBodyF.get(1));

            // Get the position of the rightand left body origins in the ground
            // This will determine which is closer to the predicted COP
            Vec3 rightBodyPos = osimModel.updBodySet().get(rightBodyName).findStationLocationInGround(s, Vec3(0, 0, 0));
            Vec3 leftBodyPos = osimModel.updBodySet().get(leftBodyName).findStationLocationInGround(s, Vec3(0, 0, 0));

            // Calculate distances from body to predicted ZMP
            double rightBodyDist = pow((rightBodyPos.get(0) - zmpCOP.get(0)), 2) + pow((rightBodyPos.get(1) - zmpCOP.get(1)), 2) + pow((rightBodyPos.get(2) - zmpCOP.get(2)), 2);
            double leftBodyDist = pow((leftBodyPos.get(0) - zmpCOP.get(0)), 2) + pow((leftBodyPos.get(1) - zmpCOP.get(1)), 2) + pow((leftBodyPos.get(2) - zmpCOP.get(2)), 2);

            // Find smaller distance and allocate to appropriate part of ZMP vector
            // There shouldn't be a need to change the non-used side as they should be zeros
            //* TODO: edge cases where identical? It won't allocate anything if that happens... *//
            if (rightBodyDist < leftBodyDist) {
                // Right side forces
                zmpVec.set(0, freeBodyF.get(0));
                zmpVec.set(1, freeBodyF.get(1));
                zmpVec.set(2, freeBodyF.get(2));
                // Right side COP
                zmpVec.set(3, zmpCOP.get(0));
                zmpVec.set(4, zmpCOP.get(1));
                zmpVec.set(5, zmpCOP.get(2));
            }
            else if (leftBodyDist < rightBodyDist) {
                // Left side forces
                zmpVec.set(6, freeBodyF.get(0));
                zmpVec.set(7, freeBodyF.get(1));
                zmpVec.set(8, freeBodyF.get(2));
                // Left side COP
                zmpVec.set(9, zmpCOP.get(0));
                zmpVec.set(10, zmpCOP.get(1));
                zmpVec.set(11, zmpCOP.get(2));
            }

        }

        // Create state vector with time and values
        // If the force threshold isn't met, then the ZMP vector will simply be zeros
        StateVector zmpStateVec = StateVector(s.getTime(), zmpVec);

        // Append current vector to ZMP results
        zmpResults.append(zmpStateVec);

    }

    // Finish section
    cout << "Calculations made for all states...\n";

    // ----- WRITE RESULTS TO FILE ----- //
    cout << "Writing results to file...\n";

    // Set column labels in storage objects
    genForcesResults.setColumnLabels(genForceLabels);
    bodyForcesResults.setColumnLabels(bodyForceLabels);
    zmpResults.setColumnLabels(zmpLabels);

    //Set name in storage
    genForcesResults.setName("Inverse Dynamics Generalized Forces");
    bodyForcesResults.setName("Inverse Dynamics Body Forces at Specified Joints");
    zmpResults.setName("ZMP Predicted Ground Reaction Forces");
	
	// Create directory to write file (where this goes, who knows...)
	IO::makeDir("outputs");

    // Write to file
    Storage::printResult(&genForcesResults, "manual_genForces", "outputs", -1, ".sto");
    Storage::printResult(&bodyForcesResults, "manual_bodyForces", "outputs", -1, ".sto");
    Storage::printResult(&zmpResults, "manual_zmpForces", "outputs", -1, ".sto");

    // Finish section
    cout << "Results written to file...\n";

    // Print some stuff to know that we got this far
	std::cout << "Code completed successfully.\n";
    std::cin.get();
    return 0;

}