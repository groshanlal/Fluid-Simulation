#include "FluidSimulation.h"
#include "common.h"

FluidSimulation::FluidSimulation(void): u_(velocity_[0]), v_(velocity_[1]), w_(velocity_[2])
{
	m_ = n_ = o_ = 0;
	timeStep_ = 0.001f;
	currentTime_ = 0;
}

FluidSimulation::~FluidSimulation(void)
{
}

FluidSimulation::FluidSimulation(int m, int n, int o, Vector3f &origin, Vector3f &extent) : systemOrigin_(origin), extent_(extent), u_(velocity_[0]), v_(velocity_[1]), w_(velocity_[2])
{
	m_ = m;
	n_ = n;
	o_ = o;
	float d[3] = {extent.x() / m, extent.y() / n, extent.z() / o};
	timeStep_ = 0.001f;
	currentTime_ = 0;

	pressure_ = Grid3f(m, n, o, origin + Vector3f(d[0]/2, d[1]/2, d[2]/2));
	density_ = Grid3f(m, n, o, origin + Vector3f(d[0]/2, d[1]/2, d[2]/2));
	backupDensity_ = Grid3f(m, n, o, origin + Vector3f(d[0]/2, d[1]/2, d[2]/2));
	
	velocity_[0] = Grid3f(m+1, n, o, origin + Vector3f(0, d[1]/2, d[2]/2));
	backupVelocity_[0] = Grid3f(m+1, n, o, origin + Vector3f(0, d[1]/2, d[2]/2));
	mass_[0] = Grid3f(m+1, n, o, origin + Vector3f(0, d[1]/2, d[2]/2));
	
	velocity_[1] = Grid3f(m, n+1, o, origin + Vector3f(d[0]/2, 0, d[2]/2));
	backupVelocity_[1] = Grid3f(m, n+1, o, origin + Vector3f(d[0]/2, 0, d[2]/2));
	mass_[1] = Grid3f(m, n+1, o, origin + Vector3f(d[0]/2, 0, d[2]/2));
	
	velocity_[2] = Grid3f(m, n, o+1, origin + Vector3f(d[0]/2, d[1]/2, 0));
	backupVelocity_[2] = Grid3f(m, n, o+1, origin + Vector3f(d[0]/2, d[1]/2, 0));
	mass_[2] = Grid3f(m, n, o+1, origin + Vector3f(d[0]/2, d[1]/2, 0));

	for ( int i = 0; i < 3; i++ ) {
		velocity_[i].dx_ = backupVelocity_[i].dx_ = mass_[i].dx_ = d[0];
		velocity_[i].dy_ = backupVelocity_[i].dy_ = mass_[i].dy_ = d[1];
		velocity_[i].dz_ = backupVelocity_[i].dz_ = mass_[i].dz_ = d[2];
	}
	pressure_.dx_ = density_.dx_ = backupDensity_.dx_ = d[0];
	pressure_.dy_ = density_.dy_ = backupDensity_.dy_ = d[1];
	pressure_.dz_ = density_.dz_ = backupDensity_.dz_ = d[2];

	int mno = m * n * o;
	A.resize(mno);
	rhs.resize(mno);

	/// Default initializations
	/// We want a default walled cubical region, so mass values at boundary = 0, next to boundary = 0.5, rest = 1
	/// Scratch that, we have a proper boundary system now
	mass_[0].setAllOne(); mass_[1].setAllOne(); mass_[2].setAllOne(); density_.setAllOne();
	float offset = 0.5f;
	simulationLimits.clearSet();
	simulationLimits.addBoundary(Vector3f(offset * d[0], 0, 0),Vector3f(1, 0, 0));
	simulationLimits.addBoundary(Vector3f( (m - offset) * d[0],0,0),Vector3f(-1, 0, 0));
	simulationLimits.addBoundary(Vector3f(0, offset * d[1], 0), Vector3f(0, 1, 0));
	simulationLimits.addBoundary(Vector3f(0, (n - offset) * d[1], 0), Vector3f(0, -1, 0));
	simulationLimits.addBoundary(Vector3f(0, 0, offset * d[2]), Vector3f(0, 0, 1));
	simulationLimits.addBoundary(Vector3f(0, 0, (o - offset) * d[2]), Vector3f(0, 0, -1));
	
	/// Set masses x, y, z
	for ( int i = 0; i < 3; i++) {
		int im = mass_[i].m_, in = mass_[i].n_, io = mass_[i].o_;
		Vector3f cellStep(mass_[i].dx_ / 2, mass_[i].dy_ / 2, mass_[i].dz_ / 2);
		for ( int x = 0; x < im; x++ ) {
			for ( int y = 0; y < in; y++ ) {
				for ( int z = 0; z < io; z++ ) {
					Vector3f cellCenter = mass_[i].convert(x, y, z);
					mass_[i](x, y, z) = simulationLimits.getPartialMass(cellCenter - cellStep, cellCenter + cellStep);
				}
			}
		}
	}
	/*
	for ( int y = 0; y < n; y++ ) {
		for ( int z = 0; z < o; z++ ) {
			mass_[0](0, y, z) = 0; mass_[0](m, y, z) = 0;
			mass_[0](1, y, z) = 0.5f; mass_[0](m-1, y, z) = 0.5f;
		}
	}

	for ( int x = 0; x < m; x++ ) {
		for ( int z = 0; z < o; z++ ) {
			mass_[1](x, 0, z) = 0; mass_[1](x, n, z) = 0;
			mass_[1](x, 1, z) = 0.5f; mass_[1](x, n-1, z) = 0.5f;
		}
	}

	for ( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			mass_[2](x, y, 0) = 0; mass_[2](x, y, o) = 0;
			mass_[2](x, y, 1) = 0.5f; mass_[2](x, y, o-1) = 0.5f;
		}
	}

	/// One last pass to clean things up, Slightly redundant loops - values per cell may be set more than once
	/// But leads to simpler loops
	for ( int y = 0; y < n; y++ ) {
		for ( int z = 0; z < o; z++ ) {
			mass_[1](0, y, z) = 0; mass_[1](m-1, y, z) = 0; mass_[2](0, y, z) = 0; mass_[2](m-1, y, z) = 0;
		}
	}
	for ( int x = 0; x < m; x++ ) {
		for ( int z = 0; z < o; z++ ) {
			mass_[0](x, 0, z) = 0; mass_[0](x, n-1, z) = 0; mass_[2](x, 0, z) = 0; mass_[2](x, n-1, z) = 0;
		}
	}
	for ( int x = 0; x < m; x++ ) {
		for ( int y =0; y < n; y++ ) {
			mass_[0](x, y, 0) = 0; mass_[0](x, y, o-1) = 0; mass_[1](x, y, 0) = 0; mass_[1](x, y, o-1) = 0;
		}
	}
	*/
	mass_[0].writeMatlab("debugMass1.m","mass1");
	mass_[1].writeMatlab("debugMass2.m","mass2");
	mass_[2].writeMatlab("debugMass3.m","mass3");
}

int& FluidSimulation::m() {
	return m_;
}

int& FluidSimulation::n() {
	return n_;
}

int& FluidSimulation::o() {
	return o_;
}

Vector3f& FluidSimulation::origin() {
	return systemOrigin_;
}

Vector3f& FluidSimulation::extent() {
	return extent_;
}

Grid3f& FluidSimulation::pressure() {
	return pressure_;
}

Grid3f& FluidSimulation::density() {
	return density_;
}

Grid3f& FluidSimulation::velocity(int i) {
	assert ( (i==0) || (i==1) || (i == 2) );
	return velocity_[i];
}

Grid3f& FluidSimulation::u() {
	return u_;
}

Grid3f& FluidSimulation::v() {
	return v_;
}

Grid3f& FluidSimulation::w() {
	return w_;
}

float& FluidSimulation::currentTime() {
	return currentTime_;
}

float& FluidSimulation::currentTimeStep() {
	return timeStep_;
}

Grid3f& FluidSimulation::mass(int i) {
	assert ( (i==0) || (i==1) || (i == 2) );
	return mass_[i];
}

float FluidSimulation::getVelocityDivergence() {
	return pressureSolveCheck();
}

float FluidSimulation::pressureSolveCheck() {
	int m = pressure_.m_, n = pressure_.n_, o = pressure_.o_;
	float minDivergence = 1 / EPSILON, maxDivergence = 0;
	float cumulativeDivergence = 0;
	int minPos[3]={-1,-1,-1}, maxPos[3]={-1,-1,-1};
	for ( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			for ( int z = 0; z < o; z++ ) {
				/// Shouldn't need this bypass anymore
				//if(mass_[0](x, y, z) + mass_[0](x+1, y, z) + mass_[1](x, y, z) + mass_[1](x, y+1, z) + mass_[2](x, y, z) + mass_[2](x, y, z+1) < EPSILON)
					//continue;
				float cellDivergence = fabsf ( ( mass_[0](x+1, y, z) * velocity_[0](x+1, y, z) - mass_[0](x, y, z) * velocity_[0](x, y, z) ) + ( mass_[1](x, y+1, z) * velocity_[1](x, y+1, z) - mass_[1](x, y, z) * velocity_[1](x, y, z) ) + ( mass_[2](x, y, z+1) * velocity_[2](x, y, z+1) - mass_[2](x, y, z) * velocity_[2](x, y, z) ) );
				if(cellDivergence<minDivergence) {
					minPos[0] = x; minPos[1] = y; minPos[2] = z;
				}
				else if(cellDivergence > maxDivergence) {
					maxPos[0] = x; maxPos[1] = y; maxPos[2] = z;
				}
				minDivergence = min<float>( minDivergence, cellDivergence );
				maxDivergence = max<float>( maxDivergence, cellDivergence );
				cumulativeDivergence += cellDivergence;
			}
		}
	}
#ifdef SHOW_DEBUG_OUTPUT
	fprintf(stderr,"Velocity Divergence: Cumulative= %f, Range= (%f,%f) at (%d,%d,%d),(%d,%d,%d)\n",cumulativeDivergence, minDivergence, maxDivergence,minPos[0],minPos[1],minPos[2], maxPos[0],maxPos[1],maxPos[2]);
#endif
	return cumulativeDivergence;
}

void advect(Grid3f& src, Grid3f& dest, Grid3f velocity[3], Grid3f mass[3], float timeStep, float d[3], bool advectDensity, float direction, int nSteps) {
	int m = src.m(), n = src.n(), o = src.o();
	//Vector3f delta(d);
	float velocityMultiplier = direction * timeStep / nSteps;
	for ( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			for ( int z = 0; z < o; z++ ) {
				if(advectDensity && (mass[0](x , y, z) + mass[0](x+1, y, z) + mass[1](x, y, z) + mass[1](x, y+1, z) + mass[2](x, y, z) + mass[2](x, y, z+1) < EPSILON)) {
					/// Solid cell, don't need to advect density
					dest(x, y, z) = src(x, y, z);
					continue;
				}
				Vector3f position = src.convert(x, y, z), nextPosition;
				for ( int step = 0; step < nSteps; step++ ) {
					Vector3f currentVelocity = Vector3f(velocity[0](position), velocity[1](position), velocity[2](position));
					/// Only move forward if the next position is still in fluid!
					saxpby(nextPosition, position, 1, currentVelocity, velocityMultiplier);
					if ( mass[0](nextPosition) < 0.5f || mass[1](nextPosition) < 0.5f || mass[2](nextPosition) < 0.5f) {
						/// Your are going to go to into an obstacle, don't
						break;
					}
					position = nextPosition;
				}
				/// We've reached the position, now set the new value
				dest(x, y, z) = src(position);
			}
		}
	}
}

void pressureSolve(Grid3f velocity[3], Grid3f& pressure, Grid3f mass[3], SparseMatrixf& A, vector<float>& rhs, float timeStep, Grid3f& density, float d[3]) {
	A.zero();
	int m = pressure.m(), n = pressure.n(), o = pressure.o();

	/// Currently restricted to cubical cells
	assert ( (d[0] == d[1]) && (d[1] == d[2]) );

	/// Static multipliers for LHS and RHS
	float staticLeftMultiplier = timeStep / d[0];
	float staticRightMultiplier = 1;

	for( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			for ( int z = 0; z < o; z++ ) {
				/// Each of this loop fills in line z + o*( y + n*x )
				/// For each line, we look at 6 cells {x+-1,y,z},{x,y+-1,z},{x,y,z+-1} apart from {x,y,z}
				
				/// Get fluid masses, they will always exist (no boundary problems), and (x,y,z) on each of these grids means a different location in space. Sorry for the naming convention, it's simple :)
				/// Each mass variable stores fraction of fluid in cell, 0 implies no fluid, 1 implies full of fluid
				float fluidMass[6] = {mass[0](x,y,z), mass[0](x+1,y,z), mass[1](x,y,z), mass[1](x,y+1,z), mass[2](x,y,z), mass[2](x,y,z+1)};

				/// If summation of all these mass values is close to 0, we're in a solid cell
				if( fluidMass[0] + fluidMass[1] + fluidMass[2] + fluidMass[3] + fluidMass[4] + fluidMass[5] < EPSILON) {
					continue;
				}

				/// Thankfully density shares the same numbering
				float cellDensity = density(x,y,z);

				float xyzCoefficient = 0;
				float rhsCoefficient = 0;
				int currentRow = pressure.iconvert(x, y, z);

				/// Now iterate over all the neighboring cells
				for ( int i = 0; i < 6; i++ ) {
					/// Coordinate of this neighbor cell
					int neighborPos[3] = { x + (i == 1) - (i == 0), y + (i == 3) - (i == 2), z + (i == 5) - (i == 4)};
					/// Check whether this coordinate is valid
					if ( fluidMass[i] < EPSILON || pressure.outOfLimits(neighborPos[0], neighborPos[1], neighborPos[2]) ) {
					//if ( fluidMass[i] < EPSILON || pressure.outOfLimits(neighborPos[0] - 1, neighborPos[1] - 1, neighborPos[2] - 1) || pressure.outOfLimits(neighborPos[0] + 1, neighborPos[1] + 1, neighborPos[2] + 1) ) {
						/// This cell boundary is with so little fluid, it won't make a difference, or it's outside the grid
						continue;
					}

					/// Need separate coordinates for velocity, because instead of %-1 and %+1, velocity requires %, %+1
					int neighborVelocityPos[3] = { x + (i == 1), y + (i == 3), z + (i == 5)};
					int relevantVelocity = i / 2;
					float neighborDensityInverse = 2 / ( cellDensity + density(neighborPos[0], neighborPos[1], neighborPos[2]) );

					/// Add component to xyzCoefficient
					xyzCoefficient += fluidMass[i] * neighborDensityInverse;
					/// Old code, further derivation showed that 1/density can be cancelled on both sides
					//xyzCoefficient += fluidMass[i] * SQR(neighborDensityInverse);
					/// Add corresponding term to RHS
					int dynamicRHSMultiplier = (x - neighborPos[0]) + (y - neighborPos[1]) + (z - neighborPos[2]);
					assert ( dynamicRHSMultiplier == 1 || dynamicRHSMultiplier == -1 );
					rhsCoefficient += dynamicRHSMultiplier * fluidMass[i] * velocity[relevantVelocity](neighborVelocityPos[0], neighborVelocityPos[1], neighborVelocityPos[2]);
					/// Old code, further derivation showed that 1/density can be cancelled on both sides
					//rhsCoefficient += dynamicRHSMultiplier * fluidMass[i] * neighborDensityInverse * velocity[relevantVelocity](neighborVelocityPos[0], neighborVelocityPos[1], neighborVelocityPos[2]);

					/// Add non-diagonal component
					A.set_element( currentRow, pressure.iconvert(neighborPos[0], neighborPos[1], neighborPos[2]), -fluidMass[i] * staticLeftMultiplier * (neighborDensityInverse) );
					//A.set_element( currentRow, pressure.iconvert(neighborPos[0], neighborPos[1], neighborPos[2]), -fluidMass[i] * staticLeftMultiplier * SQR(neighborDensityInverse) );
				}

				/// All neighborhood cells are done, now we can set the xyz cell, and rhs
				A.set_element( currentRow, currentRow, xyzCoefficient * staticLeftMultiplier );
				rhs[currentRow] = staticRightMultiplier * rhsCoefficient;
			}
		}
	}

	/// Matrix and RHS are set, now we can solve the pressure system
#ifdef DO_DEBUG_TESTS
	fprintf(stderr,"Pressure solve debug tests:\nDumping A to file\n");
	A.write_matlab(fstream("debugA.m",ios::out),"A");
	FILE *temporaryFile = fopen("debugRHS.m","w");
	fprintf(temporaryFile,"rhs=[");
	for ( int i = 0; i < (int)rhs.size(); i++ ) fprintf(temporaryFile,"%f ",rhs[i]);
	fprintf(temporaryFile,"];\n");
	fclose(temporaryFile);
#endif
	PCGSolver<float> equationSolver;
	float residual = 0;
	int iterations = 0;
	equationSolver.solve(A, rhs, pressure.data_, residual, iterations);
#ifdef SHOW_DEBUG_OUTPUT
	fprintf(stderr,"Solver complete in %d iterations, with residual = %f\n", iterations, residual);
#endif

	/// Solver complete, update velocities
	float velocityUpdateMultiplier = timeStep;
	for ( int x = 0; x < m + 1; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			for ( int z = 0; z < o; z++ ) {
				velocity[0](x, y, z) -= velocityUpdateMultiplier / d[0] * ( pressure(x, y, z) - pressure(x-1, y, z) ) * 2 / ( density(x, y, z) + density(x-1, y, z) );
			}
		}
	}

	for ( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n + 1; y++ ) {
			for ( int z = 0; z < o; z++ ) {
				velocity[1](x, y, z) -= velocityUpdateMultiplier / d[1] * ( pressure(x, y, z) - pressure(x, y-1, z) ) * 2 / ( density(x, y, z) + density(x, y-1, z) );
			}
		}
	}

	for ( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			for ( int z = 0; z < o + 1; z++ ) {
				velocity[2](x, y, z) -= velocityUpdateMultiplier / d[2] * ( pressure(x, y, z) - pressure(x, y, z-1) ) * 2 / ( density(x, y, z) + density(x, y, z-1) );
			}
		}
	}
}

void applyAcceleration(Grid3f velocity[3], Grid3f mass[3], Vector3f& acceleration, float timeStep) {
	/// Still assume that the 1st and last column planes are static velocity
	int m = velocity[0].m(), n = velocity[0].n(), o = velocity[0].o();
	float acc = acceleration.x();
	for ( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			for ( int z = 0; z < o; z++ ) {
				if(mass[0](x, y, z) > EPSILON)
					velocity[0](x, y, z) += acc * timeStep;
				//else velocity[0](x, y, z) = 0;
			}
		}
	}
	m = velocity[1].m(); n = velocity[1].n();
	acc = acceleration.y();
	for ( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			for ( int z = 0; z < o; z++ ) {
				if(mass[1](x, y, z) > EPSILON)
					velocity[1](x, y, z) += acc * timeStep;
				//else velocity[1](x, y, z) = 0;
			}
		}
	}
	n = velocity[2].n(); o = velocity[2].o();
	acc = acceleration.z();
	for ( int x = 0; x < m; x++ ) {
		for ( int y = 0; y < n; y++ ) {
			for ( int z = 0; z < o; z++ ) {
				if(mass[2](x, y, z) > EPSILON)
					velocity[2](x, y, z) += acc * timeStep;
				//else velocity[2](x, y, z) = 0;
			}
		}
	}
}

void FluidSimulation::computeNextStep(Vector3f &acceleration, bool useAdaptiveTimeStep) {
#ifdef SHOW_DEBUG_OUTPUT
	float valmin, valmax;
	fprintf(stderr,"Stage 1: Starting advection\n");
#ifdef DO_DEBUG_TESTS
	fprintf(stderr,"Stage 1: Pre-advection Ranges: ");
	for ( int i = 0; i < 3; i++) {
		velocity_[i].getMinMax(valmin,valmax);
		fprintf(stderr,"%d(%f, %f) ",i,valmin,valmax);
	}
	fprintf(stderr,"\n");
	velocity_[0].writeMatlab("debug_ubase.m","baseu");
	velocity_[1].writeMatlab("debug_vbase.m","basev");
	velocity_[2].writeMatlab("debug_wbase.m","basew");
#endif
#endif
	float d[3] = {pressure_.dx_, pressure_.dy_, pressure_.dz_}; 
	pressureSolve(velocity_, pressure_, mass_, A, rhs, timeStep_, density_, d);
	for ( int i = 0; i < 3; i++ )
		advect(velocity_[i], backupVelocity_[i], velocity_, mass_, timeStep_, d);
	for ( int i = 0; i < 3; i++ )
		velocity_[i].data_.swap(backupVelocity_[i].data_);

#ifdef SHOW_DEBUG_OUTPUT
#ifdef DO_DEBUG_TESTS
	fprintf(stderr,"Stage 1: Ranges: ");
	for ( int i = 0; i < 3; i++) {
		velocity_[i].getMinMax(valmin,valmax);
		fprintf(stderr,"%d(%f, %f) ",i,valmin,valmax);
	}
	fprintf(stderr,"\n");
#endif
	fprintf(stderr,"Stage 1: Advection of velocity complete, doing density advection\n");
#ifdef DO_DEBUG_TESTS
	density_.getMinMax(valmin,valmax);
	fprintf(stderr,"Pre advection density range: (%f, %f)\n",valmin,valmax);
#endif
#endif
	advect(density_, backupDensity_, velocity_, mass_, timeStep_, d, true);
	density_.data_.swap(backupDensity_.data_);
#ifdef SHOW_DEBUG_OUTPUT
#ifdef DO_DEBUG_TESTS
	density_.getMinMax(valmin,valmax);
	fprintf(stderr,"Post advection density range: (%f, %f)\n",valmin,valmax);
#endif
	fprintf(stderr,"Stage 1: Complete\nStage 2: Starting body forces\n");
#endif

#ifdef DO_DEBUG_TESTS
	velocity_[0].writeMatlab("debug_uba.m","bua");
	velocity_[1].writeMatlab("debug_vba.m","bva");
	velocity_[2].writeMatlab("debug_wba.m","bwa");
#endif
	applyAcceleration(velocity_, mass_, acceleration, timeStep_);
#ifdef SHOW_DEBUG_OUTPUT
	fprintf(stderr,"Stage 2: Divergence test:\n");
	getVelocityDivergence();
	fprintf(stderr,"Stage 2: Complete\nStage 3: Starting pressure solve\n");
#endif

#ifdef DO_DEBUG_TESTS
	velocity_[0].writeMatlab("debug_ub.m","bu");
	velocity_[1].writeMatlab("debug_vb.m","bv");
	velocity_[2].writeMatlab("debug_wb.m","bw");
	mass_[0].writeMatlab("debug_mu.m","mu");
	mass_[1].writeMatlab("debug_mv.m","mv");
	mass_[2].writeMatlab("debug_mw.m","mw");
#endif
	pressureSolve(velocity_, pressure_, mass_, A, rhs, timeStep_, density_, d);
#ifdef DO_DEBUG_TESTS
	velocity_[0].writeMatlab("debug_ua.m","au");
	velocity_[1].writeMatlab("debug_va.m","av");
	velocity_[2].writeMatlab("debug_wa.m","aw");
#endif

#ifdef SHOW_DEBUG_OUTPUT
#ifdef DO_DEBUG_TESTS
	fprintf(stderr,"Post Pressure solve Ranges: ");
	for ( int i = 0; i < 3; i++) {
		velocity_[i].getMinMax(valmin,valmax);
		fprintf(stderr,"%d(%f, %f) ",i,valmin,valmax);
	}
	fprintf(stderr,"\n");
#endif
	fprintf(stderr,"Stage 3: Complete\n");
#endif
#ifdef DO_DEBUG_TESTS
	getVelocityDivergence();
#endif
	currentTime_ += timeStep_;
	if(useAdaptiveTimeStep) {
		/// Find norm(vel) - approximate by Vector3f(max(u),max(v),max(w))
		float vmax[3], vmin[3];
		for ( int i = 0; i < 3; i++ ) {
			velocity_[i].getMinMax(vmin[i],vmax[i]);
		}

		timeStep_ = min<float>(CFL_C * pressure_.dx_ / ( Vector3f(vmax).length() + EPSILON ), 0.1f);
	}
}