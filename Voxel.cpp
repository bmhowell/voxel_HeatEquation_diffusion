// Class ==> Function Definitions
//
// Created by Brian Howell on 11/9/21.
//

#include "Voxel.h"

// Default Constructor - sets member variables to their default (null) state
//Voxel::Voxel() {
//    /* double colon :: is referred to as the scope resolution
//     *
//     * */
//}

// Overload Constructor
//   - Sets the input variables to whatever we pass through the Class.
Voxel::Voxel(int node_, float tFinal_, double dt_){

    // MEMBER VARIABLES

    // Particle parameters
    thetaWall = 300.;                                   // |    K    |  temperature at the wall
    theta0 = 300.;                                      // |    K    |  initial temperature
    intensity = 2 * pow(10, 7);           // |  W/m^2  |  initial laser intensity
    absorb = 25.;                                       // |   1/m   |  absorption constant
    lenBlock = 0.05;                                    // |    m    |  sample length
    intThick = 1.;                                      // |   ___   |  number of nodes

    // Material parameters
    kP = 903;                                           // |  W/m-K  |  thermal conductivity of particle
    kM = 1;                                             // |  W/m-K  |  thermal conductivity of material
    cP = 903;                                           // | J/kg-K  |  heat capacity of particle
    cM = 156;                                           // | J/kg-K  |  heat capacity of material
    rho0P = 2700;                                       // | kg/m^3  |  initial particle density
    rho0M = 1000;                                       // | kg/m^3  |  initial material density
    vP = 0.3;                                           // |   ---   |  volume fraction of particles
    rParticle = lenBlock / 10;                          // |    m    |  radius of the particles

    // simulation parameters
    tFinal = tFinal_;                                   // |    s    |  final simulation time
    dt = dt_;                                           // |    s    |  initial time discretization
    node = node_;                                       // |   ---   |  number of nodes
    h = (float) lenBlock / (node - 1);                  // |    m    |  physical discretization length
    sizeTime = (int) tFinal / dt;                       // |   ---   |  total number of time steps

    // misc. parameters
    SIZEA3 = pow(node, 3);                      // |    ---   |  total number of nodes
    SIZEA2 = pow(node, 2);                      // |    ---   |  number of nodes on a single face

    // initialize material parameters
    for (int i=0; i<SIZEA3; i++){
        density.push_back(rho0M);                      //  | kg/m^3  |  initialize as density of material
        heatCap.push_back(cM);                         //  | J/kg-K  |  initialize as heat capacity of material
        thermCond.push_back(kM);                       //  | J/kg-K  |  initialize as thermal conductivity of material
        laserValues.push_back(0.0);                    //  |  W/m^2  |  initialize as 0 at each node ( update in laserProfile() )
        theta.push_back(theta0);                       //  |    K    |  initial temperature of material
    }



//    // vectors and arrays for particle generation
//    std::vector<int> particlesInd;                               // vector holding total indices for each particle
//    std::vector<int> particlesInterInd;                          // interfacial distance indices
//
//    // initialize material properties for each node

//    ** change to vectors **
//    double density[SIZEA3];
//    double heatCap[SIZEA3];
//    double thermCond[SIZEA3];
//    std::fill_n(density, SIZEA3, rho0M);
//    std::fill_n(heatCap, SIZEA3, cM);
//    std::fill_n(thermCond, SIZEA3, kM);
//
//    // initialize temperature
//    std::vector< std::vector<double> > temperature;               // store temperature all time steps
//    std::vector<double> theta(SIZEA3, theta0);                 // store temperature for individual time steps

}

// Destructor
Voxel::~Voxel() {
}

// Accessor functions

// ** currently throws error: reference to non-static member function must be called; did you
//    mean to call it with no arguments?
//double Voxel::testVector(std::vector<double> &myVector) {
//    return myVector.size;
//}

int Voxel::get_sizeA3() const{
    return SIZEA3;
}

double Voxel::get_thetaWall() const{
    return thetaWall;
}

double Voxel::get_theta0() const{
    return theta0;
}

double Voxel::get_intensity() const{
    return intensity;
}

double Voxel::get_absorb() const{
    return absorb;
}

double Voxel::get_lenBlock() const{
    return lenBlock;
}

double Voxel::get_intThick() const{
    return intThick;
}

int Voxel::get_kP() const{
    return kP;
}

int Voxel::get_kM() const{
    return kM;
}

int Voxel::get_cM() const{
    return cM;
}

int Voxel::get_cP() const{
    return cP;
}

double Voxel::get_vP() const{
    return vP;
}

int Voxel::get_node() const{
    return node;
}

double Voxel::get_tFinal() const{
    return tFinal;
}

void Voxel::get_densityVec(){
    std::cout << "density: " << std::endl;
    for (int i=0; i<density.size(); i++){
        std::cout << density[i] << std::endl;
    }
    std::cout << std::endl;
}

void Voxel::get_particleIndVec(){
    std::cout << "particleIndVec: " << std::endl;
    for (int i=0; i<particlesInd.size(); i++){
        std::cout << particlesInd[i] << std::endl;
    }
}

void Voxel::get_bNodes(){
    std::cout << "bNodes: " << std::endl;
    for (int i=0; i<bNodes.size(); i++){
        std::cout << "bNodes[i] size: " << bNodes[i].size() << std::endl;
        std::cout << "\n   bNodes[i]: " << i << std::endl;
        for (int j=0; j<bNodes[i].size(); j++){
            std::cout << "    " << bNodes[i][j] << std::endl;
        }
    }
}

void Voxel::get_nonbNodes(){
    std::cout << "nonbNodes.size() = " << nonbNodes.size() << std::endl;
    std::cout << "nonbNodes: " << std::endl;
    for (int i=0; i<nonbNodes.size(); i++){
        std::cout << nonbNodes[i] << std::endl;
    }
}

/* Mutator functions
 *
 * Allows us to edit/modify each of the member variables one at a time.
 * Does not return anything. Acts like overload constructor, but they only
 * modify one member variable at a time.
 *
 * */

// HELPER MUTATOR FUNCTIONS
void Voxel::set_thetaWall(double thetaWall_){
    thetaWall = thetaWall_;
}

void Voxel::set_node(int node_){
    node = node_;
}

void Voxel::set_dt(double dt_){
    dt = dt_;
}

/* Simulation mutator functions
 *
 * Simulation functions required for solving temperature distribution
 * problem.
 *
 * */

// helper functions
void Voxel::uniqueVec(std::vector<int> &vec){
    // return only unique nodes in particleInd
    // https://www.geeksforgeeks.org/stdunique-in-cpp/

//    std::cout << "before: " << std::endl;
//    for (int i = 0; i < vec.size(); i++){
//        std::cout << vec[i] << " ";
//    }
//    std::cout << std::endl;

    // begin removing duplicate indices
    std::vector<int>::iterator ip;
    ip = std::unique(vec.begin(), vec.begin() + vec.size());
    vec.resize(std::distance(vec.begin(), ip));

//    std::cout << "after: " << std::endl;
//    for (ip = vec.begin(); ip != vec.end(); ++ip) {
//        std::cout << *ip << " ";
//    }
//    std::cout << std::endl;
}

void Voxel::density2file(){

    // write to file
    printDensity.open("/Users/brianhowell/Desktop/Berkeley/MSOL/voxelTemperature/density.dat");
    printDensity << "X Y Z D C K" << std::endl;

    for (int i=0; i<SIZEA3; i++){
        printDensity << cubeCoord[i][1] << " " << cubeCoord[i][2] << " " << cubeCoord[i][3]
                     << " " << density[i] << " " << heatCap[i] << " " << thermCond[i] << std::endl;
    }
    printDensity.close();
}

void Voxel::temp2file(){
    // TO BE COMPLETED!!!
}

void Voxel::lastTemp2file(){
    // write to file
    printTemp.open("/Users/brianhowell/Desktop/Berkeley/MSOL/voxelTemperature/lastTemp.dat");
    printTemp << "X Y Z T" << std::endl;

    for (int i = 0; i < SIZEA3; i++){
        printTemp << cubeCoord[i][1] << " " << cubeCoord[i][2] << " " << cubeCoord[i][3]
                  << " " << theta[i] << std::endl;
    }
    printTemp.close();
}

// functions
void Voxel::computeCoord(){
    /* stores x-y-z coordinates for each node:
     *
     * @updateVec cubeCoord = [ [node #, x, y, z],
     *                          [      ...      ],
     *                          [node n, x, y, z] ]
     */

    // required variables
    double coordDist = 0;               // distance incrementer
    std::vector<double> x;              // node numbering in x direction
    std::vector<double> y;              // node numbering in y direction
    std::vector<double> z;              // node numbering in z direction
    std::vector<double> xyz;            // node numbering in xyz direction

    std::cout << "--- Constructing cubeCoord array ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    // populate vectors x-y-z
    for (int i = 0; i < node; i +=1 ){
//        xyz.push_back(coordDist);
        x.push_back(coordDist);         // add spatial increment to x vector
        y.push_back(coordDist);         // add spatial increment to y vector
        z.push_back(coordDist);         // add spatial increment to z vector
        coordDist += h;                 // increment by spatial discretization
    }

    // populate array coordCube with x-y-z directions
    int counter = 0;
    for (int i=0; i<node; i++){
        for (int j=0; j<node; j++){
            std::vector<double> temp;
            for (int k=0; k<node; k++){
                temp.push_back(counter);
                temp.push_back(x[k]);
                temp.push_back(y[j]);
                temp.push_back(z[i]);
//                std::cout << xyz[k] << " " << xyz[j]  << " " << xyz[i] << std::endl;
                counter++;
                cubeCoord.push_back(temp);
                temp.clear();
            }

        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "    computational time: " << duration << "s" << std::endl;

}

void Voxel::computeBoundary(){
    /* Stores the node numbers of all boundary nodes.
     *      - Dirichlet on the sides of the wall
     *      - Neumann on the tops and bottoms (requires adjacent nodes).
     *
     * @updateVec bNodes - [ [node #: top boundary node],
     *                       [node #: top adjacent bNode],
     *                       [node #: bottom boundary],
     *                       [node #: bottom adjacent bNode],
     *                       [node #: wall boundary] ]
     */

    std::cout << "--- Constructing bNodes array ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> topBoundary;                                   // boundary node along top
    std::vector<int> topBoundary_h;                                 // adjacent bNode along top
    std::vector<int> bottomBoundary;                                // boundary node along bottom
    std::vector<int> bottomBoundary_h;                              // adjacent bNode along bottom
    std::vector<int> wallBoundary;                                  // boundary node along wall

    int counter = 0;
    for (int i = 0; i < node; i++){
        for (int j = 0; j < node; j++){
            for (int k = 0; k < node; k++){
                if (i == 0){
                    // log the two boundary nodes
                    bottomBoundary.push_back(counter);              // boundary node
                    bottomBoundary_h.push_back(counter + SIZEA2);   // adjacent boundary node
                }
                else if (i == node - 1){
                    topBoundary.push_back(counter);                 // boundary node
                    topBoundary_h.push_back(counter - SIZEA2);      // adjacent boundary node
                }
                else if (j == 0 or j == node - 1 or k == 0 or k == node - 1){
                    wallBoundary.push_back(counter);
                }else{
                    nonbNodes.push_back(counter);
                }
                counter++;                                          // increments node counter
            }
        }
    }
    bNodes.push_back(topBoundary);
    bNodes.push_back(topBoundary_h);
    bNodes.push_back(bottomBoundary);
    bNodes.push_back(bottomBoundary_h);
    bNodes.push_back(wallBoundary);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "    compute time boundaries: " << duration << "s\n" << std::endl;
}

void Voxel::computeAsparse(){
    /* computeAsparse - Builds a sparse A matrix
     *
     * structure:
     *      A = [ [node 1, i, j, value],
     *            [node ..., i, j, value],
     *            [node SIZEA3, i, j, value] ]
     *
     * @updateVec Asparse - |  2D Vector  |
     */

    std::cout << "\n--- Constructing A sparse array ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    // Loop through every single node
    for (int i = 0; i < SIZEA3; i++){
        // Loop through only nodes that
        // are not on the bottom or top faces
        if (SIZEA2 < i and i < SIZEA3) {
            std::vector<int> tempVec;
            for (int j = 0; j < SIZEA3; j++) {
                if ((j == i - SIZEA2) or (j == i - node) or (j == i - 1)) {
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(1);
                } else if (j == i) {
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(-6);
                } else if ((j == i + 1) or (j == i + node) or (j == i + SIZEA2)) {
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(1);
                } else {
                    continue;
                }
                ASparse.push_back(tempVec);
                tempVec.clear();
            }
        }else{
            continue;
        }
    }

    /*
     * enforce dirichlet boundary conditions along the walls
     * loop through boundary wall nodes bNodes[4][i] and set the value of
     * ASparse[node][3] equal to zero
     */
    for (int i = 0; i < bNodes[4].size(); i++){
        for(int j = 0; j < ASparse.size(); j++) {
            if (ASparse[j][0] == bNodes[4][i]) {
                ASparse[j][2] = 0;
            }
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "    compute time A matrix: " << duration << "s" << std::endl;
}

void Voxel::computeParticles(){
    /* computeParticles - Function computes random location within Voxel for a particle.
     *                    Adjacent nodes that fall within radius of a particle are then
     *                    marked to be a part of the particle.
     *
     * @paramVector cubeCoord - coordinates for each node
     *
     * @updateVector particlesInd - vector holding total indices for each particle
     * @updateVector particlesInterInd - interfacial distance indices
     *
     */

    int nParticleNode = std::round(SIZEA3 * vP);     // total number of host nodes for particles
    double partDist, nodeParticle, randLoc, testVar;

    int counter1 = 0;
    std::cout << "\n---- GENERATING PARTICLES ----" << std::endl;
    while ((particlesInd.size() < nParticleNode) and (counter1 < 10000)){

        // choose a random node to generate particle
        // https://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::uniform_real_distribution<> dist{0, 1};
        randLoc = dist(gen);          // ensure diameter is positive

        // Generate random seed location for particle
//        std::cout << "\n-- Generating particle " << counter1 + 1 << std::endl;
        nodeParticle = ceil(SIZEA3 * randLoc);

        // find nodes that are within the distance of the seed location
        for (int i = 0; i < SIZEA3; i++){

            // compute distance between generated particle and all other particles
            partDist = sqrt(pow(cubeCoord[i][1] - cubeCoord[nodeParticle][1], 2) +
                            pow(cubeCoord[i][2] - cubeCoord[nodeParticle][2], 2) +
                            pow(cubeCoord[i][3] - cubeCoord[nodeParticle][3], 2));

            // determine if distance is within the range of another particle
            if (partDist <= rParticle){
                particlesInd.push_back(i);
            }else if (intThick != 0 and partDist <= rParticle + intThick * h){
                particlesInterInd.push_back(i);

            }else{
                continue;
            }
        }

        // ensure vectors contain no duplicates of nodes
        uniqueVec(particlesInd);
        uniqueVec(particlesInterInd);

        // print total number of nodes in both vectors:
//        std::cout << "particleInd_.size(): " << particlesInd.size() << std::endl;
//        std::cout << "particleInterInd_.size(): " << particlesInterInd.size() << std::endl;

        // assign interfacial material properties
        for (int i = 0; i < particlesInterInd.size(); i++){
            density[particlesInterInd[i]] = (rho0M + rho0P) / 2.;
            heatCap[particlesInterInd[i]] = (cP + cM) / 2.;
            thermCond[particlesInterInd[i]] = (kP + kM) / 2.;
        }

        // assign particle material properties
        for (int i = 0; i < particlesInd.size(); i++){
            density[particlesInd[i]] = rho0P;
            heatCap[particlesInd[i]] = cP;
            thermCond[particlesInd[i]] = kP;
        }
        counter1++;
    }
    std::cout << "--> total particles generated: " << counter1 << std::endl;
}

void Voxel::laserProfile(){
    /* laserProfile - updates "laserValues" vector with corresponding intensity values
     *                as computed by Beer-Lambert
     *
     * @updateVec laserValues - intensity at each node from laser beam
     */

    for (int i=0; i < SIZEA3; i++){
        if (0.0199 <= cubeCoord[i][1] and cubeCoord[i][1] <= 0.0301 and 0.0199 <= cubeCoord[i][2] and cubeCoord[i][2] <= 0.0301){
            laserValues[i] = absorb * intensity * exp(-absorb * (lenBlock - cubeCoord[i][3]));
        }
    }
}

void Voxel::solutionSchemeMatrix(){
    /* solutionScheme - using A matrix for FDM computation
     *
     * @paramVec - Asparse          | 2D vector |  A matrix for FDM stencil
     * @paramVec - bNodes           | 2D vector |  boundary nodes
     *
     * @paramVec - density          | 1D vector |  density at each node
     * @paramVec - heatCapacity     | 1D vector |  heat capacity at each node
     *
     * @updateVec - temperature     | 2D vector |  temperature for all time steps
     * @updateVec - theta           | 1D vector |  temperature for individual time steps
     */


    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    // initialize required variables for solution scheme
    double prefix1, prefix2, averageTemp;
    int row, col, val;
    int nnz = ASparse.size();

//    // initialize temperature (room temperature theta0)
//    std::vector<double> theta(SIZEA3, theta0);

    // begin time stepping
    std::cout << "---- SIMULATING: A MATRIX METHOD ----" << std::endl;
    for (int t = 0; t < sizeTime; t++){
        // print output information
        averageTemp = std::accumulate(theta.begin(), theta.end(), 0.0) / theta.size();
//        std::cout << "time: " << t + 1 << " / " << sizeTime;
//        std::cout << " -> average temperature: " << averageTemp << std::endl;

        // coordinate-wise sparse-matrix-vector multiplication
        // http://www.mathcs.emory.edu/~cheung/Courses/561/Syllabus/3-C/sparse.html
        double AtimesTheta[SIZEA3];
        std::fill_n(AtimesTheta, SIZEA3, 0);
        for (int i = 0; i < nnz; i++){
            row = ASparse[i][0];
            col = ASparse[i][1];
            val = ASparse[i][2];
            AtimesTheta[row] = AtimesTheta[row] + val * theta[col];
        }

        // using solve for temperatures at the next time step
        for (int j = 0; j < SIZEA3; j++){
            theta[j] = theta[j] + dt / density[j] / heatCap[j] *
                                  (thermCond[j] / pow(h, 2.0) * AtimesTheta[j] + laserValues[j]);
        }

        // enforce Neumann boundary conditions on the top boundary
        for (int k = 0; k < bNodes[0].size(); k++){
            theta[bNodes[0][k]] = theta[bNodes[1][k]];
        }

        // enforce Neumann boundary conditions on the bottom boundary
        for (int l = 0; l < bNodes[2].size(); l++){
            theta[bNodes[2][l]] = theta[bNodes[3][l]];
        }

        // store temperature results (every 100 steps)
        if (t % 100 == 0){
            temperature.push_back(theta);
        }else{
            continue;
        }
    }

//    std::cout << "time: " << t + 1 << " / " << sizeTime;
    std::cout << " -> average temperature: " << averageTemp << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "simulation time: " << duration << "s" << std::endl;

}

void Voxel::solutionSchemeFAST(){
    /* solutionScheme - using A matrix for FDM computation
     *
     * @paramVec - bNodes           | 2D vector |  boundary nodes
     *                              [ [node #: top boundary node],
     *                                [node #: top adjacent bNode],
     *                                [node #: bottom boundary],
     *                                [node #: bottom adjacent bNode],
     *                                [node #: wall boundary] ]
     *
     * @paramVec - density          | 1D vector |  density at each node
     * @paramVec - heatCapacity     | 1D vector |  heat capacity at each node
     *
     * @updateVec - temperature     | 2D vector |  temperature for all time steps
     * @updateVec - theta           | 1D vector |  temperature for individual time steps
     */

    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    // initialize required variables for solution scheme
    double prefix1, prefix2, Avalue;
    double averageTemp;
//    int row, col, val;
//    int nnz = ASparse.size();

    // begin time stepping
    std::cout << "\n---- SIMULATING: FAST METHOD ----" << std::endl;
    for (int t=0; t<sizeTime; t++){
        // print output information
        averageTemp = std::accumulate(theta.begin(), theta.end(), 0.0) / theta.size();
//        std::cout << "time: " << t + 1 << " / " << sizeTime;
//        std::cout << " -> average temperature: " << averageTemp << std::endl;

        /* HANDLING BOUNDARY CONDITIONS:
         *
         * - Dirichlet BCs on the walls: Enforce by not updating the temperature.
         * - Neumann BCs on the top and bottom: Enforce be setting adjacent nodes equal to the top and bottom
         *
         * Computational Approach:
         * - First, loop through all internal nodes (automatically skips boundary nodes)
         * - Second, loop through all adjacent nodes under the top boundary -> bNodes[1]
         * - Third, loop through all adjacent nodes under the bottom boundary -> bNodes[3]
         */

        // loop through internal nodes (non boundary nodes)
        for (int i=0; i<nonbNodes.size(); i++){

            // HEAT EQUATION:
            prefix1 = dt / (density[nonbNodes[i]] * heatCap[nonbNodes[i]]);
            prefix2 = thermCond[nonbNodes[i]] / pow(h, 2);
            Avalue = (theta[nonbNodes[i] - SIZEA2] + theta[nonbNodes[i] - node] + theta[nonbNodes[i] - 1]
                      - 6 * theta[nonbNodes[i]]
                      + theta[nonbNodes[i] + 1] + theta[nonbNodes[i] + node] + theta[nonbNodes[i] + SIZEA2]);

            theta[nonbNodes[i]] = theta[nonbNodes[i]] + prefix1 * (prefix2 * (Avalue) + laserValues[nonbNodes[i]]);
        }

        // loop through adjacent nodes to top boundary -> bNodes[1]
        for (int i=0; i<bNodes[1].size(); i++){
            // set node equals
//            theta[bNodes[1][i]] = theta[i + SIZEA2];
            // alternatively:
            theta[bNodes[1][i]] = theta[bNodes[0][i]];
        }

        // loop through adjacent nodes to top boundary -> bNodes[3]
        for (int i=0; i<bNodes[3].size(); i++){
//            theta[bNodes[3][i]] = theta[i - SIZEA2];
            // alternatively: ** (might be in issue -- check if one-to-one mapping between bNodes[2] and bNodes[3]
            theta[bNodes[3][i]] = theta[bNodes[2][i]];
        }

        // store temperature results (every 100 steps)
        if (t % 100 == 0){
            temperature.push_back(theta);
        }else{
            continue;
        }
    }

    std::cout << " -> average temperature: " << averageTemp << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "simulation time: " << duration << "s" << std::endl;
}

void Voxel::laserSimulation(){
    // laserSimulation - runs required functions to compute temperature evolution of
    //                   a material exposed to a laser

    computeCoord();                 // compute the coordinates for each node
    computeBoundary();              // compute which nodes are on the sides, top and bottom of cube
    computeAsparse();               // compute FDM mesh A matrix
    computeParticles();             // compute random particles embedded in voxel
    laserProfile();                 // compute laser intensity values at each node
    solutionSchemeMatrix();         // compute solution via A matrix
}

void Voxel::laserSimulationFAST(){
    // laserSimulation - runs required functions to compute temperature evolution of
    //                   a material exposed to a laser

    computeCoord();                 // compute the coordinates for each node
    computeBoundary();              // compute which nodes are on the sides, top and bottom of cube
    computeParticles();             // compute random particles embedded in voxel
    laserProfile();                 // compute laser intensity values at each node
    solutionSchemeFAST();         // compute solution via A matrix
}






