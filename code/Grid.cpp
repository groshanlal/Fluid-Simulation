#include<vector>
#include "Grid.h"

using namespace std;

Grid::Grid(int m, int n, int o, float (&orig)[3])
{
  X = m;
  Y = n;
  Z = o;
  dx = 1;
  dy = 1;
  dz = 1;
  origin[0] = orig[0];
  origin[1] = orig[1];
  origin[2] = orig[2];
}

float Grid::getVolume()
{
  return X*Y*Z;
}

float Grid::getX()
{
  return X;
}

float Grid::getY()
{
  return Y;
}

float Grid::getZ()
{
  return Z;
}

float* Grid::getOrigin()
{
  return origin;
}
