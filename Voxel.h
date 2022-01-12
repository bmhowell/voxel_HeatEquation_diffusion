//
// Created by Brian Howell on 11/9/21.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <random>
#include <vector>
#include <algorithm>

#ifndef VOXELTEMPERATURE_VOXEL_H
#define VOXELTEMPERATURE_VOXEL_H


class Voxel {

private:
    /* Encapsulation allows us to privatize certain information inside of our program
     * so it can't be accessed from other files (preventing tampering).
     * Sometimes, a function declaration goes here, but typically it is just the member variables.
     *
     * */

    // MEMBER VARIABLES

    // particle parameters
    float thetaWall;                                  // |    K    |  temperature at the wall
    float theta0;                                     // |    K    |  initial temperature
    double intensity;                                 // |  W/m^2  |  initial laser intensity
    float absorb;                                     // |   1/m   |  absorption constant
    float lenBlock;                                   // |    m    |  sample length
    double intThick;                                  // |   ---   |  interfacial thickness parameter

    // material parameters
    int kP;                                           // |  W/m-K  |  thermal conductivity of particle
    int kM;                                           // |  W/m-K  |  thermal conductivity of material
    int cP;                                           // | J/kg-K  |  heat capacity of particle
    int cM;                                           // | J/kg-K  |  heat capacity
    int rho0P;                                        // | kg/m^3  |  initial particle density
    int rho0M;                                        // | kg/m^3  |  initial material density
    double vP;                                        // |   ---   |  volume fraction of particles
    double rParticle;                                 // |    m    |  radius of the particles

    // simulation parameters
    double tFinal;                                    // |    s    |  final simulation time
    double dt;                                        // |    s    |  initial time discretization
    int node;                                         // |   ---   |  number of nodes
    double h;                                         // |    m    |  spatial discretization
    int sizeTime;                                     // |   ---   |  total number of time steps

    // misc. parameters
    int SIZEA3;                                       // |   ---   |  total number of nodes
    int SIZEA2;                                       // |   ---   |  number of nodes on a single face

    // initialize cubeCoord, bNode, and A matrix vectors
    std::vector< std::vector<double> > cubeCoord;
    std::vector< std::vector<int> > bNodes;
    std::vector<int> nonbNodes;
    std::vector< std::vector<int> > ASparse;

    // vectors and arrays for particle generation
    std::vector<int> particlesInd;                     // vector holding total indices for each particle
    std::vector<int> particlesInterInd;                // interfacial distance indices

    // initialize material properties for each node
    std::vector<int> testVec;
    std::vector<float> density;                        // density at each node in voxel
    std::vector<float> heatCap;                        // heat capacity at each node in voxel
    std::vector<float> thermCond;                      // thermal conductivity at each node in voxel
    std::vector<double> laserValues;                   // laser intensity at each node in voxel

    // initialize temperature vectors
    std::vector< std::vector<double> > temperature;    // store temperature for all time steps
    std::vector<double> theta;                         // store temperature at a single times

    // data outputs
    std::ofstream printDensity;
    std::ofstream printTemp;


public:
    // Typically function declarations go under "Public"

    /* Default Constructor */
    Voxel()
        : thetaWall(300.), theta0(300.), intensity(2 * pow(10, 7)), absorb(25.),
          lenBlock(0.05), node(31), intThick(1.), kP(903), kM(1), cM(156), rho0P(2700), rho0M(1000),
          vP(0.3), tFinal(2.), dt(pow(10, -4)),
          SIZEA3(pow(node, 3)),  SIZEA2(pow(node, 2))
    {

        // initialize cubeCoord, bNode, and A matrix vectors
        std::vector< std::vector<double> > cubeCoord;
        std::vector< std::vector<int> > bNodes;
        std::vector< std::vector<int> > ASparse;

        // vectors and arrays for particle generation
        std::vector<int> particlesInd;                               // vector holding total indices for each particle
        std::vector<int> particlesInterInd;                          // interfacial distance indices

        // initialize material properties for each node
        double density[SIZEA3];
        double heatCap[SIZEA3];
        double thermCond[SIZEA3];
        std::fill_n(density, SIZEA3, rho0M);
        std::fill_n(heatCap, SIZEA3, cM);
        std::fill_n(thermCond, SIZEA3, kM);

        // initialize temperature
        std::vector< std::vector<double> > temperature;               // store temperature all time steps
        std::vector<double> theta(SIZEA3, theta0);                 // store temperature for individual time steps

        // compute necessary matrices
//        computeCoord(cubeCoord);                                     // compute the x-y-z coordinates
//        computeBoundary(bNodes);                                     // find the boundary nodes
//        computeAsparse(ASparse, bNodes);                             // compute FDM mesh A matrix

    };

    /* Overload Constructor */
    Voxel(int node_, float tFinal_, double dt_);

    /* Destructor */
    ~Voxel();

    // Accessor Functions
    //   - Returns optimization parameters

//    double testVector(std::vector<double> &);
//    // testVector: returns size of vector

    int get_sizeA3() const;
    // get_sizeA3: Returns total number of nodes in 3D voxel

    double get_thetaWall() const;
    // get_thetaWall: Returns temperature of the wall

    double get_theta0() const;
    // get_theta0: Returns initial temperature

    double get_intensity() const;
    // get_intensity(): return intensity of the laser

    double get_absorb() const;
    // get_absorb(): return absorption constant

    double get_lenBlock() const;
    // get_lenBlock(): returns the length sample material

    double get_intThick() const;
    // get_intThick(): returns interstitial thickness

    int get_kP() const;
    // get_kP(): returns thermal conductivity of particle

    int get_kM() const;
    // get_kM(): returns thermal conductivity of material

    int get_cM() const;
    // get_cM(): returns heat capacity of material

    int get_cP() const;
    // get_cP(): returns heat capacity of particle

    double get_vP() const;
    // get_vP(): returns particle volume fraction

    int get_node() const;
    // get_node: Returns number of nodes in one direction

    double get_tFinal() const;
    // get_tFinal: Returns total simulation time

    void get_densityVec();
    // get_floatVec: Returns all values within myVec

    void get_particleIndVec();
    // get_particleIndVec(): Returns all indices within Voxel for particles properties

    void get_bNodes();
    // get_bNodes(): Returns nodes that are on the boundary

    void get_nonbNodes();
    // get_nonbNodes(): Returns nodes that are not on the boundary

    ///////////////////////////////////////////////////////////////////////////////////
    /* Mutator functions
     *
     * Allows us to edit/modify each of the member variables one at a time.
     * Does not return anything. Acts like overload constructor, but they only
     * modify one member variable at a time.
     * */
    /////////////////////////////////////////
    // HELPER MUTATOR FUNCTIONS
    void set_thetaWall(double thetaWall_);
    // set_thetaWall(): sets thetaWall parameter
    // @param double - new thetaWall parameter

    void set_node(int node_);
    // get_node(): sets node parameter
    // @param int - new node parameter

    void set_dt(double dt_);
    // set_dt(): sets dt parameter
    // @param double - new dt parameter

    /////////////////////////////////////////
    /* Simulation mutator functions
     *
     * Simulation functions required for solving temperature distribution
     * problem.
     *
     * */

    // helper function declarations
    void uniqueVec(std::vector<int> &);
    // uniqueVec - modifies input vector to only contain unique indices
    // @param - vec: vector to be passed in

    void density2file();
    // density2file - prints material parameters of each node to density.dat

    void temp2file();
    // temp2file - prints temperatures at all timesteps to folder
    void lastTemp2file();
    // lastTemp2file - prints the final temperature at each node to lastTemp.dat

    // function declarations
    void computeCoord();
    // computeCoord - compute the coords in the x-y-z directions
    // @param vector< vector<double> > - the 2D vector to be filled
    // @param double N -  Number of nodes
    // @param double L -  Sample length

    void computeBoundary();
    // boundaryNodes - computes which nodes are on the sides, top and bottom of the cube
    // @params vector< vector<double> > - input 2D vector to be filled
    // @params double - number of nodes
    // @return vector< vector<double> > - vector of vectors containing each side

    void computeAsparse();
    // computeAsparse - computes the A matrix, in the sparse variety (see comment below for structure)
    // @param - modifies 2D vector into a [n x 4] vector

    void computeParticles();
    // computeParticles - compute random particles within the medium of material
    // @paramVector cubeCoord - coordinates for each node
    // @updateVector particlesInd - vector holding total indices for each particle
    // @updateVector particlesInterInd - interfacial distance indices

    void laserProfile();
    // laserProfile - computes energy supplied to each node based on Beer-Lambert equation
    // @paramVector - cubeCoord
    // @paramVector - laserValues

    void solutionSchemeMatrix();
    // solutionScheme - using A matrix for FDM computation
    // @paramVec - ASparse: FDM mesh matrix
    // @paramVec - bNodes: nodes of the boundaries
    // @paramVec - density: array for density of each node
    // @paramVec - heatCap: array for heat capacity of each node
    // @paramVec - thermCond: array for thermal conductivity
    // @paramVec - laserValues: energy input from laser at each node
    // @updateVec - temperature: solution vector for every time step
    // @updateVec - theta: solution vector for a single time step

    void solutionSchemeFAST();
    // solutionScheme - using FAST FDM computation
    // @paramVec - bNodes: nodes of the boundaries
    // @paramVec - density: array for density of each node
    // @paramVec - heatCap: array for heat capacity of each node
    // @paramVec - thermCond: array for thermal conductivity
    // @paramVec - laserValues: energy input from laser at each node
    // @updateVec - temperature: solution vector for every time step
    // @updateVec - theta: solution vector for a single time step

    void laserSimulation();
    // laserSimulation - computes the temperature evolution of a material exposed to a laser

    void laserSimulationFAST();
    // laserSimulationFAST -- compute a fast version without matrix multiplication


};


#endif //VOXELTEMPERATURE_VOXEL_H
