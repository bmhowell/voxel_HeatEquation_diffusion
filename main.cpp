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
 * The following simulation is my first real attempt in writing simulations in C++
 * and contains overly detailed comments directed for my own learning of this language.
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
    Voxel VoxelSystem(NODE, TFINAL, DT);      // Using overload constructor

    // Testing the instance of the voxel class
    std::cout << "Hello, World!" << std::endl;
    std::cout << "Wall temperature: " << VoxelSystem.get_thetaWall() << std::endl;
    std::cout << "Simulation time: " << VoxelSystem.get_tFinal() << std::endl;
    std::cout << "Total Nodes: " << VoxelSystem.get_node() << std::endl;
    std::cout << "Total nodes in 3D domain: " << VoxelSystem.get_sizeA3() << std::endl;

    VoxelSystem.laserSimulation();              // run laser simulation
    VoxelSystem.lastTemp2file();                // print final temperature to lastTemp.dat
    VoxelSystem.density2file();                 // print material parameters for each node to density.dat

//    VoxelSystem.set_thetaWall(500.);
//    std::cout << "Wall temperature: " << VoxelSystem.get_thetaWall() << std::endl;
//    std::cout << "SIZEA3: " << VoxelSystem.get_sizeA3() << std::endl;

    /*
     * next steps:
     *
     * -code return functions for multidimensional vectors: DONE
     * -confirm computeCoord, computeBoundary and computeAsparse are working correctly: DONE
     * -modify compute particles to handle vectors rather than arrays: DONE
     * -finish coding member functions for simulation: DONE
     * -write one single member function to compute all the steps above: DONE
     *
     * -add function to print temperature at every time step
     * -replace A matrix method with direct computation
     * -replace forward euler with adaptive trapezoidal method
     * -add solution for diffusion equation
     */



    return 0;
}
