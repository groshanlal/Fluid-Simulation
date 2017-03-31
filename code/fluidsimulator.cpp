#include "Grid.h"
#include<vector>
#define EPSILON 0.001f
using namespace std;

class FluidSimulator
{
public:
	int X, Y, Z 	//Grid Dimensions of pressure, density, etc.
	Grid			pressure_;
	Grid		        density;
	Grid			velocity_[3];
	Grid			mass_[3];
	float			currentTime_;
	float			Time_Step;

public:

	FluidSimulation(int m, int n, int o, Vector3f& origin, Vector3f& extent);
	int&		X();
	int&		Y();
	int&		Z();
	void Acceleration(Grid velocity[3], Grid mass[3], vector<float> acceleration[3] , float Time_Step);
};

int& FluidSimulator::X() {
	return X;
}


void Acceleration(Grid velocity[3], Grid mass[3], Vector<float> acceleration[3], float Time_Step) {


	int X = velocity[0].X(), Y = velocity[0].Y(), Z = velocity[0].Z();

	float accelr = acceleration.x();
	for ( int i = 0; i < X; i++ ) {
		for ( int j = 0; j < Y; j++ ) {
			for ( int k = 0; k < Z; k++ ) {
				if(mass[0](i, j, k) > EPSILON)
					velocity[0](i, j, k) += accelr * Time_Step;
			}
		}
	}
	X = velocity[1].X(); Y = velocity[1].Y();
	accelr = acceleration.y();
	for ( int i = 0; i < X; i++ ) {
		for ( int j = 0; j < Y; j++ ) {
			for ( int k = 0; k < Z; k++ ) {
				if(mass[1](i, j, k) > EPSILON)
					velocitj[1](i, j, k) += accelr * Time_Step;

			}
		}
	}
	Y = velocity[2].Y(); Z = velocity[2].Z();
	accelr = acceleration.z();
	for ( int i = 0; i < X; i++ ) {
		for ( int j = 0; j < Y; j++ ) {
			for ( int k = 0; k < Z; k++ ) {
				if(mass[2](i, j, k) > EPSILON)
					velocity[2](i, j, k) += accelr * Time_Step;

			}
		}
	}
}
