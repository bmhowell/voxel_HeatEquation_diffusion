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

    // misc. parameters
    int SIZEA3;                                 // |   ---   |  total number of nodes
    int SIZEA2;                                 // |   ---   |  number of nodes on a single face

    // initialize cubeCoord, bNode, and A matrix vectors
    std::vector< std::vector<double> > cubeCoord;
    std::vector< std::vector<int> > bNodes;
    std::vector< std::vector<int> > ASparse;

    // vectors and arrays for particle generation
    std::vector<int> particlesInd;                               // vector holding total indices for each particle
    std::vector<int> particlesInterInd;                          // interfacial distance indices

    // initialize material properties for each node
    std::vector<int> testVec;
    std::vector<float> density;
    std::vector<float> heatCap;
    std::vector<float> thermCond;

//    double density[SIZEA3];
//    double heatCap[SIZEA3];
//    double thermCond[SIZEA3];
//    std::fill_n(density.begin(), SIZEA3, rho0M);
//    std::fill_n(heatCap.begin(), SIZEA3, CM);
//    std::fill_n(thermCond, SIZEA3, KM);

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
    // get_intThick(): returns interstial thickness

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

    /* Mutator functions
     *
     * Allows us to edit/modify each of the member variables one at a time.
     * Does not return anything. Acts like overload constructor, but they only
     * modify one member variable at a time.
     * */

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

    void computeParticles(std::vector< std::vector<double> >&,
                          std::vector<int>&,
                          std::vector<int>&);
    // computeParticles - compute random particles within the medium of material
    // @param - modifies 2D vector and fills with particles









};


#endif //VOXELTEMPERATURE_VOXEL_H
