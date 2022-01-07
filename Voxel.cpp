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

    // misc. parameters
    SIZEA3 = pow(node, 3);                      // |    ---   |  total number of nodes
    SIZEA2 = pow(node, 2);                      // |    ---   |  number of nodes on a single face

    // initialize material parameters
    for (int i=0; i<SIZEA3; i++){
        density.push_back(rho0M);                      //  | kg/m^3  |  initialize as density of material
        heatCap.push_back(cM);                         //  | J/kg-K  |  initialize as heat capacity of material
        thermCond.push_back(kM);                       //  | J/kg-K  |  initialize as thermal conductivity of material
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

// functions
void Voxel::computeCoord(){
    /* structure cubeCoord:
     * cubeCoord = [node, x, y, z]
     */

    // required variables
    double coordDist = 0;               // distance incrementer
    std::vector<double> x;              // node numbering in x direction
    std::vector<double> y;              // node numbering in y direction
    std::vector<double> z;              // node numbering in z direction
    std::vector<double> xyz;            // node numbering in xyz direction

    std::cout << "--- Initializing computation --- " << std::endl;
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
    std::cout << "computational time: " << duration << "s\n" << std::endl;

}

void Voxel::computeBoundary(){
    std::cout << "--- Initializing computation --- " << std::endl;
    std::cout << "--- Constructing bNodes array ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> topBoundary;
    std::vector<int> topBoundary_h;
    std::vector<int> bottomBoundary;
    std::vector<int> bottomBoundary_h;
    std::vector<int> wallBoundary;

    int counter = 0;
    for (int i = 0; i < node; i++){
        for (int j = 0; j < node; j++){
            for (int k = 0; k < node; k++){
                if (i == 0){
                    // log the two boundary nodes
                    bottomBoundary.push_back(counter);              // boundary node
                    bottomBoundary_h.push_back(counter + SIZEA2);   // adjacent boundary node
                }
                if (i == node - 1){
                    topBoundary.push_back(counter);                 // boundary node
                    topBoundary_h.push_back(counter - SIZEA2);      // adjacent boundary node
                }
                if (j == 0 or j == node - 1 or k == 0 or k == node - 1){
                    wallBoundary.push_back(counter);
                }
                counter++;
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
    std::cout << "compute time boundaries: " << duration << "s\n" << std::endl;
}

void Voxel::computeAsparse(){
    /* Build a sparse A matrix
     *  A = [N^3, 4]
     *  [node, i, j, value]
     */
    std::cout << "--- Initializing computation --- " << std::endl;
    std::cout << "--- Constructing A sparse array ---" << std::endl;
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
    std::cout << "compute time A matrix: " << duration << "s" << std::endl;
}

void Voxel::computeParticles(std::vector< std::vector<double> > &cubeCoord,
                      std::vector<int> &particlesInd,
                      std::vector<int> &particlesInterInd){

    int nParticleNode = std::round(SIZEA3 * vP);     // total number of host nodes for particles
    double partDist, nodeParticle, randLoc, testVar;

    int counter1 = 0;
    while ((particlesInd.size() < nParticleNode) and (counter1 < 10000)){

        // choose a random node to generate particle
        // https://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::uniform_real_distribution<> dist{0, 1};
        randLoc = dist(gen);          // ensure diameter is positive

        // Generate random seed location for particle
        std::cout << "\n-- Generating particle " << counter1 + 1 << std::endl;
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

        uniqueVec(particlesInd);
        uniqueVec(particlesInterInd);

//        std::cout << "particleInd_.size(): " << particlesInd.size() << std::endl;
//        std::cout << "particleInterInd_.size(): " << particlesInterInd.size() << std::endl;

        // assign interfacial material properties
//        for (int i = 0; i < particlesInterInd.size(); i++){
//            density[particlesInterInd[i]] = (rho0M + rho0P) / 2.;
//            heatCap[particlesInterInd[i]] = (cP + cM) / 2.;
//            thermCond[particlesInterInd[i]] = (kP + kM) / 2.;
//        }
//
//        // assign particle material properties
//        for (int i = 0; i < particlesInd.size(); i++){
//            density[particlesInd[i]] = rho0P;
//            heatCap[particlesInd[i]] = cP;
//            thermCond[particlesInd[i]] = kP;
//        }
        counter1++;
    }
}








