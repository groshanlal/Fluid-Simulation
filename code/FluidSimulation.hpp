#ifndef FLUID_SIMULATION_H
#define FLUID_SIMULATION_H

#include "Grid.hpp"

using namespace std;

class FluidSimulation
{
  public:
    Grid pressure, density;
    Grid vx, vy;
    float timestep;

    FluidSimulation(float origin_x, float origin_y, int xn, float dx, float t);
    void computeNextStep();
    void applyAcceleration(float ax, float ay);
    void advect();
    void pressureSolve();
};

#endif
