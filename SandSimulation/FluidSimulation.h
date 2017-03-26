#pragma once

#include "Grid3f.h"
#include "pcg_solver.h"
#include "BoundarySet.h"

#define CFL_C 3

class FluidSimulation
{
public:
	/// Where the fluid simulation is in space
	Vector3f		systemOrigin_;
	/// How big is it - this will decide dx, dy, dz for each grid
	Vector3f		extent_;
	/// Grid dimensions - same as the respective quantites in centrally sampled grids
	int				m_;
	int				n_;
	int				o_;

	/// Pressure
	Grid3f			pressure_;
	/// Density
	Grid3f			density_, backupDensity_;
	/// Velocity
	Grid3f			velocity_[3], backupVelocity_[3];
	/// Convienient references to velocity components - I've often found them useful and cleaner
	Grid3f&			u_;
	Grid3f&			v_;
	Grid3f&			w_;
	/// Fluid Mass
	Grid3f			mass_[3];
	/// Current time
	float			currentTime_;
	/// Time Step
	float			timeStep_;

	/// Variables for pressure solve
	SparseMatrixf	A;
	vector<float>	rhs;

	BoundarySet		simulationLimits;

public:
	FluidSimulation(void);
	~FluidSimulation(void);
	FluidSimulation(int m, int n, int o, Vector3f& origin, Vector3f& extent);

	int&		m();
	int&		n();
	int&		o();
	Vector3f&	origin();
	Vector3f&	extent();
	Grid3f&		pressure();
	Grid3f&		density();
	Grid3f&		velocity(int i);
	Grid3f&		u();
	Grid3f&		v();
	Grid3f&		w();
	Grid3f&		mass(int i);
	float&		currentTime();
	float&		currentTimeStep();
	float		getVelocityDivergence();
	float		pressureSolveCheck();

	/// Single use function, useAdaptiveTimeStep chooses the next timestep based on the CFL condition
	void		computeNextStep(Vector3f& acceleration, bool useAdaptiveTimeStep = true);
	void		doPressureSolve(bool useAdaptiveTimeStep = true);
};

void advect(Grid3f& src, Grid3f& dest, Grid3f velocity[3], Grid3f mass[3], float timeStep, float d[3], bool advectDensity = false, float direction = -1, int nSteps=10);

void pressureSolve(Grid3f velocity[3], Grid3f& pressure, Grid3f mass[3], SparseMatrixf& A, vector<float>& rhs, float timeStep, Grid3f& density, float d[3]);

void applyAcceleration(Grid3f velocity[3], Grid3f mass[3], Vector3f& acceleration, float timeStep);