/* -------------------------------------------------------------------------- *
 * OpenSim: test.cpp												          *
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

// This code includes some tests for how ZMP predictions of GRF could be included
// in an analyze tool plugin. The idea is to take a prescribed motion from coordinates
// and run the ZMP calculations on this by using the Inverse Dynamics solver.

#include <OpenSim/OpenSim.h>
#include <OpenSim/Simulation/InverseDynamicsSolver.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimulationUtilities.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/FunctionSet.h>
#include <OpenSim/Common/GCVSplineSet.h>
#include <OpenSim/Common/IO.h>

using namespace OpenSim;
using namespace std;
using namespace SimTK;

int main() {

    // ----- SETTINGS ----- //
    //* TODO: these could be analysis inputs at some point? *//

    // Set vertical force threshold
    int forceThreshold = 50;

    // Set right and left body ground contact names
    string rightBodyName = "calcn_r";
    string leftBodyName = "calcn_l";

    // Set time ranges
    double start_time = 0.549;
    double final_time = 1.143;

    // Set the free joint
    string freeJointName = "ground_pelvis";

    // ----- DEFINE MODEL ----- //
    cout << "Loading model...\n";

    // Import model geometries to avoid warning print-outs --- full path needed?
    //ModelVisualizer.addDirToGeometrySearchPaths("C:\\OpenSim 4.4\\Geometry");

    // Define the model --- need to specify full path for some reason...?
    //Model osimModel("osimModel.osim");
    Model model("C:\\+GitRepos+\\zmp-opensim\\pilot\\analyzeImplementation\\osimModel.osim");

    // Initialize the model's underlying computational system and get its default state.
    State& s = model.initSystem();

    // OpenSim::Coordinates represent degrees of freedom for a model.
        // Each Coordinate's value and speed maps to an index
        // in the model's underlying SimTK::State (value to a slot in the
        // State's q, and speed to a slot in the State's u).
        // So we need to map each OpenSim::Coordinate value and speed to the
        // corresponding SimTK::State's q and u indices, respectively.
    auto coords = model.getCoordinatesInMultibodyTreeOrder();
    int nq = s.getNQ();
    int nu = s.getNU();
    int nCoords = (int)coords.size();
    int intUnusedSlot = -1;

    // Create a vector mapCoordinateToQ whose i'th element provides
        // the index in vector 'coords' that corresponds with the i'th 'q' value.
        // Do the same for mapCoordinateToU, which tracks the order for
        // each 'u' value.
    auto coordMap = createSystemYIndexMap(model);
    vector<int> mapCoordinateToQ(nq, intUnusedSlot);
    vector<int> mapCoordinateToU(nu, intUnusedSlot);
    for (const auto& c : coordMap) {
        // SimTK::State layout is [q u z] 
        // (i.e., all "q"s first then "u"s second).
        if (c.second < nq + nu) {
            string svName = c.first;

            // The state names corresponding to q's and u's will be of 
            // the form:
            // /jointset/(joint)/(coordinate)/(value or speed).
            // So the corresponding coordinate name is second from the end.
            ComponentPath svPath(svName);
            string lastPathElement = svPath.getComponentName();
            string coordName = svPath.getSubcomponentNameAtLevel(
                svPath.getNumPathLevels() - 2);

            for (int i = 0; i < nCoords; ++i) {
                if (coordName == coords[i]->getName()) {
                    if (lastPathElement == "value") {
                        mapCoordinateToQ[c.second] = i;
                        break;
                    }

                    // Shift State/System indices by nq since u's follow q's
                    else if (lastPathElement == "speed") {
                        mapCoordinateToU[c.second - nq] = i;
                        break;
                    }

                    else {
                        throw OpenSim::Exception("Last element in state variable "
                            " name " + svName + " is neither "
                            "'value' nor 'speed'");
                    }
                }
                if (i == nCoords - 1) {
                    throw OpenSim::Exception("Coordinate " + coordName +
                        " not found in model.");
                }
            }

        }
    }

    // Make sure that order of coordFunctions (which define splines for 
        // State's q's) is in the same order as the State's q order.
        // Also make a new vector (coordinatesToSpeedsIndexMap) that, for each
        // u in the State, gives the corresponding index in q (which is same
        // order as  coordFunctions). This accounts for cases where qdot != u.
    FunctionSet coordFunctions;
    coordFunctions.ensureCapacity(nq);
    vector<int> coordinatesToSpeedsIndexMap(nu, intUnusedSlot);

    // Finish section
    cout << "Model loaded...\n";

    // ----- LOAD COORDINATES ----- //
    cout << "Loading coordinates...\n";

    // Define the coordinates --- need to specify full path for some reason...?
    Storage coordinateValues("C:\\+GitRepos+\\zmp-opensim\\pilot\\analyzeImplementation\\coordinates.sto");

    // Convert degrees to radians if indicated
    if (coordinateValues.isInDegrees()) {
        model.getSimbodyEngine().convertDegreesToRadians(coordinateValues);
    }

    // Create differentiable splines of the coordinate data
    GCVSplineSet coordSplines(5, &coordinateValues);

    // Functions must correspond to model coordinates.
            // Solver needs the order of Function's to be the same as order
            // in State's q's.
    for (int i = 0; i < nq; i++) {
        int coordInd = mapCoordinateToQ[i];

        // unused q slot
        if (coordInd == intUnusedSlot) {
            coordFunctions.insert(i, new Constant(0));
            continue;
        }

        const Coordinate& coord = *coords[coordInd];
        if (coordSplines.contains(coord.getName())) {
            coordFunctions.insert(i, coordSplines.get(coord.getName()));
        }
        else {
            coordFunctions.insert(i, new Constant(coord.getDefaultValue()));
            log_info("InverseDynamicsTool: coordinate file does not "
                "contain coordinate '{}'. Assuming default value.",
                coord.getName());
        }

        // Fill in coordinatesToSpeedsIndexMap as we go along to make
        // sure we know which function corresponds to State's u's.
        for (int j = 0; j < nu; ++j) {
            if (mapCoordinateToU[j] == coordInd) {
                coordinatesToSpeedsIndexMap[j] = i;
            }
        }
    }

    if (coordFunctions.getSize() > nq) {
        coordFunctions.setSize(nq);
    }

    // TODO: disable model forces?

    //// Get times from coordinate values
    //double first_time = _coordinateValues->getFirstTime();
    //double last_time = _coordinateValues->getLastTime();

    // Determine time indices
    int start_index = coordinateValues.findIndex(start_time);
    int final_index = coordinateValues.findIndex(final_time);

    // Finish section
    cout << "Coordinates data loaded...\n";

    // ----- RUN INVERSE DYNAMICS ----- //
    cout << "Running inverse dynamics...\n";

    // create the solver given the input data
    InverseDynamicsSolver ivdSolver(model);

    // Set number of times in data for solver
    int nt = final_index - start_index + 1;

    // Get time data for solver
    Array_<double> times(nt, 0.0);
    for (int i = 0; i < nt; i++) {
        times[i] = coordinateValues.getStateVector(start_index + i)->getTime();
    }

    // Preallocate results
    Array_<Vector> genForceTraj(nt, Vector(nCoords, 0.0));

    // solve for the trajectory of generalized forces that correspond to the 
        // coordinate trajectories provided
    ivdSolver.solve(s, coordFunctions, coordinatesToSpeedsIndexMap, times,
        genForceTraj);

    // Finish section
    cout << "Inverse dynamics successful...\n";

    // ----- EXTRACT RELEVANT FORCES ----- //
    cout << "Extracting forces to trajectories...\n";

    // Set the free joint name in array
    Array<string> jointNames(freeJointName, 1);

    // Create the jointset for extracting body forces
    const JointSet& modelJoints = model.getJointSet();

    //  Create the joints for equivalent body forces to fill
    JointSet jointsForEquivalentBodyForces;

    // Identify joints in model for reporting body forces
    // There is an 'ALL' option here, but it should be unsued
    /* The search for individual group or force names IS case-sensitive BUT keywords are not*/
    for (int i = 0; i < jointNames.getSize(); ++i) {
        //Check for keywords first starting with ALL
        if (IO::Uppercase(jointNames[i]) == "ALL") {
            for (int j = 0; j < modelJoints.getSize(); ++j) {
                jointsForEquivalentBodyForces.adoptAndAppend(&modelJoints[j]);
            }
            break;
        }

        int k = modelJoints.getIndex(jointNames[i]);
        if (k >= 0) {
            jointsForEquivalentBodyForces.adoptAndAppend(&modelJoints[k]);
        }
        else {
            log_warn("InverseDynamicsTool could not find Joint named '{}' to "
                "report body forces.", jointNames[i]);
        }
    }
    jointsForEquivalentBodyForces.setMemoryOwner(false);

    // Get number of joints for reporting
    int nj = jointsForEquivalentBodyForces.getSize();

    // Generalized forces from ID Solver are in MultibodyTree order and not
        // necessarily in the order of the Coordinates in the Model.
        // We can get the Coordinates in Tree order from the Model.
    Array<string> labels("time", nCoords + 1);
    for (int i = 0; i < nCoords; i++) {
        labels[i + 1] = coords[i]->getName();
        labels[i + 1] += (coords[i]->getMotionType() == Coordinate::Rotational) ?
            "_moment" : "_force";
    }

    Array<string> body_force_labels("time", 6 * nj + 1);
    string XYZ = "XYZ";
    for (int i = 0; i < nj; i++) {
        string joint_body_label = jointsForEquivalentBodyForces[i].getName() + "_";
        joint_body_label += jointsForEquivalentBodyForces[i].getChildFrame().getName();
        for (int k = 0; k < 3; ++k) {
            body_force_labels[6 * i + k + 1] = joint_body_label + "_F" + XYZ[k]; //first label is time
            body_force_labels[6 * i + k + 3 + 1] = joint_body_label + "_M" + XYZ[k];
        }
    }

    // Set storage locations
    Storage genForceResults(nt);
    Storage bodyForcesResults(nt);
    Storage zmpResults(nt);
    SpatialVec equivalentBodyForceAtJoint;

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

    // Loop through times to extract forces
    for (int i = 0; i < nt; i++) {

        StateVector
            genForceVec(times[i], genForceTraj[i]);
        genForceResults.append(genForceVec);

        // if there are joints requested for equivalent body forces then calculate them
        if (nj > 0) {
            Vector forces(6 * nj, 0.0);
            StateVector bodyForcesVec(times[i],
                SimTK::Vector_<double>(6 * nj,
                    &forces[0]));

            s.updTime() = times[i];
            Vector& q = s.updQ();
            Vector& u = s.updU();

            // Account for cases where qdot != u with coordinatesToSpeedsIndexMap
            for (int j = 0; j < nq; ++j) {
                q[j] = coordFunctions.evaluate(j, 0, times[i]);
            }
            for (int j = 0; j < nu; ++j) {
                u[j] = coordFunctions.evaluate(
                    coordinatesToSpeedsIndexMap[j], 1, times[i]);
            }


            for (int j = 0; j < nj; ++j) {
                equivalentBodyForceAtJoint = jointsForEquivalentBodyForces[j].calcEquivalentSpatialForce(s, genForceTraj[i]);
                for (int k = 0; k < 3; ++k) {
                    // body force components
                    bodyForcesVec.setDataValue(6 * j + k, equivalentBodyForceAtJoint[1][k]);
                    // body torque components
                    bodyForcesVec.setDataValue(6 * j + k + 3, equivalentBodyForceAtJoint[0][k]);
                }
            }
            bodyForcesResults.append(bodyForcesVec);

        }

        // Calculate ZMP forces for current time state

        // Set vector to store results for current state (6 values per limb)
        Vector zmpVec(6 * 2);

        // Get the body forces acting at the free body (generally pelvis)
        //* TODO: this is manually indexed to get pelvis --- this would need to be fixed if more joints are in there... *//
        Vec3 freeBodyF = Vec3(equivalentBodyForceAtJoint[1][0],
            equivalentBodyForceAtJoint[1][1],
            equivalentBodyForceAtJoint[1][2]);

        // Check if vertical force is greater than threshold to perform calculations
        // There shouldn't be a need to change the non-used side as they should be zeros
        //* TODO: this should account for different vertical axes... *//
        if (freeBodyF.get(1) > forceThreshold) {

            // Get the position of the free body (generally pelvis) from Q states
            //* TODO: this is manually indexed to get pelvis --- needs to be fixed for different models... *//
            Vec3 rp = Vec3(coordinateValues.getStateVector(i)->getData().get(3),
                coordinateValues.getStateVector(i)->getData().get(4),
                coordinateValues.getStateVector(i)->getData().get(5));

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
            /*TODO: this is using the default state to identify body location...*/
            /* NOTE: this may be more easily accountable when working within the analyze tool framework...*/
            Vec3 rightBodyPos = model.updBodySet().get(rightBodyName).findStationLocationInGround(s, Vec3(0, 0, 0));
            Vec3 leftBodyPos = model.updBodySet().get(leftBodyName).findStationLocationInGround(s, Vec3(0, 0, 0));

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
        StateVector zmpStateVec = StateVector(times[i], zmpVec);

        // Append current vector to ZMP results
        zmpResults.append(zmpStateVec);

    }

    // Set names in gen force results
    genForceResults.setColumnLabels(labels);
    genForceResults.setName("Inverse Dynamics Generalized Forces");

    // Set names in body force results
    bodyForcesResults.setColumnLabels(body_force_labels);
    bodyForcesResults.setName("Inverse Dynamics Body Forces at Specified Joints");

    // Set names in ZMP results
    zmpResults.setColumnLabels(zmpLabels);
    zmpResults.setName("ZMP Predicted Ground Reaction Forces");

    // Finish section
    cout << "Forces extracted...\n";

    // ----- PRINT RESULTS TO FILE ----- //
    cout << "Printing results to file...\n";

    // Create directory to write file (where this goes, who knows...)
    IO::makeDir("outputs");

    // Write to file
    Storage::printResult(&genForceResults, "manual_genForces", "outputs", -1, ".sto");
    Storage::printResult(&bodyForcesResults, "manual_bodyForces", "outputs", -1, ".sto");
    Storage::printResult(&zmpResults, "manual_zmpForces", "outputs", -1, ".sto");

    // Finish section
    cout << "Results printed to file...\n";

    // Print some stuff to know that we got this far
    std::cout << "Code completed successfully.\n";
    std::cin.get();
    return 0;

}