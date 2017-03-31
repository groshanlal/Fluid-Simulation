#include<vector>
#include "Grid.h"

using namespace std;

Grid::Grid(int m, int n, int o)
{
  X = m;
  Y = n;
  Z = o;
  dx = 1;
  dy = 1;
  dz = 1;
}

int Grid::getVolume()
{
  return X*Y*Z;
}

int Grid::getX()
{
  return X;
}

int Grid::getY()
{
  return Y;
}

int Grid::getZ()
{
  return Z;
}
