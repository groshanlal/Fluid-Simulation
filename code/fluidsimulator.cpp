#include "Grid.h"

using namespace std;

class FluidSimulator
{public:
	Grid			pressure_;
	Grid		density;
	Grid			velocity_[3];
	Grid3f			mass_[3];
	float			currentTime_;
	float			Time_Step;

public:
	FluidSimulation(void);
	~FluidSimulation(void);
	FluidSimulation(int m, int n, int o, Vector3f& origin, Vector3f& extent);

	int&		X();
	int&		Y();
	int&		Z();
	Grid&		density();
	Grid&		velocity(int i);
	Grid&		mass(int i);
  void Acceleration(Grid velocity[3], Grid mass[3], vector<float> acceleration[3] , float Time_Step);
};
int& FluidSimulator::X() {
	return X;
}

int& FluidSimulator::Y() {
	return Y;
}

int& FluidSimulator::Z() {
	return Z;
}
void Acceleration(Grid velocity[3], Grid mass[3], Vector<float> acceleration[3], float Time_Step) {
	
	int X = velocity[0].X(), Y = velocity[0].Y(), Z = velocity[0].Z();
	float accelr = acceleration.x();
	for ( int i = 0; i < X; i++ ) {
		for ( int y = 0; j < Y; j++ ) {
			for ( int k = 0; k < Z; k++ ) {
				if(mass[0](i, j, k) > EPSILON)
					velocity[0](i, j, k) += accelr * Time_Step;
			}
		}
	}
	X = velocity[1].m(); Y = velocity[1].n();
	accelr = acceleration.y();
	for ( int i = 0; i < X; i++ ) {
		for ( int j = 0; j < Y; j++ ) {
			for ( int k = 0; k < Z; k++ ) {
				if(mass[1](i, j, k) > EPSILON)
					velocitj[1](i, j, k) += accelr * Time_Step;
				
			}
		}
	}
	Y = velocity[2].n(); Z = velocity[2].o();
	accelr = acceleration.z();
	for ( int i = 0; i < X; i++ ) {
		for ( int j = 0; j < Y; j++ ) {
			for ( int k = 0; k < Z; k++ ) {
				if(mass[2](i, j, z) > EPSILON)
					velocity[2](i, j, k) += accelr * Time_Step;
			
			}
		}
	}
}

