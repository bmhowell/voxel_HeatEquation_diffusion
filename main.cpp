//
// Voxel code CPP FAST
// Created by Brian Howell on 11/09/21.
//

#include <iostream>
#include "Voxel.h"

/* Object-oriented simulation of voxel computed heat equation
 * Brian Howell
 * MSOL UC Berkeley
 *
 *
 * https://www.youtube.com/watch?v=vz1O9nRyZaY&list=PL318A5EB91569E29A&index=20
 */

int main() {

    // data outputs
    std::ofstream printDensity;
    std::ofstream printTemp;


    // Declare variables
    int NODE = 21;
    double TFINAL = 2.;
    double DT = pow(10, -4);

    // Create VoxelSystem
//    Voxel VoxelSystem;                        // Using the default constructor
    Voxel VoxelSystem1(NODE, TFINAL, DT);      // Using overload constructor
    Voxel VoxelSystem2(NODE, TFINAL, DT);

    // Testing the instance of the voxel class

    std::cout << "\n----- TESTING SLOW METHOD -----" << std::endl;
    VoxelSystem1.laserSimulation();              // run laser simulation
    std::cout << "\n----- TESTING FAST METHOD -----" << std::endl;
    VoxelSystem2.laserSimulationFAST();

//    VoxelSystem2.lastTemp2file();                // print final temperature to lastTemp.dat
//    VoxelSystem2.density2file();                 // print material parameters for each node to density.dat


    /*
     * next steps:
     *
     * -code return functions for multidimensional vectors: DONE
     * -confirm computeCoord, computeBoundary and computeAsparse are working correctly: DONE
     * -modify compute particles to handle vectors rather than arrays: DONE
     * -finish coding member functions for simulation: DONE
     * -write one single member function to compute all the steps above: DONE
     *
     * -add function to print temperature at every time step: DONE
     * -replace A matrix method with direct computation: DEBUG -> set random seed for particle gen.
     * -replace forward euler with adaptive trapezoidal method, possibly AB methods
     * -add solution for diffusion equation
     */



    return 0;
}
