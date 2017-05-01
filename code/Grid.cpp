#include "Grid.hpp"

Grid::Grid(float origin_x, float origin_y, int xn, int yn, float dx, float dy)
{
  this->origin_x = origin_x;
  this->origin_y = origin_y;
  this->xn = xn;
  this->yn = yn;
  this->dx = dx;
  this->dy = dy;

  for(int i=0;i<xn;i++)
  {
    for(int j=0;j<yn;j++)
    {
      data.push_back(0);
    }
  }
}

Grid::Grid()
{
  this->origin_x = 10;
  this->origin_y = 10;
  this->xn = 3;
  this->yn = 2;
  this->dx = 1;
  this->dy = 1;

  for(int i=0;i<xn;i++)
  {
    for(int j=0;j<yn;j++)
    {
      data.push_back(0);
    }
  }
}

float Grid::getdata(int i, int j)
{
  return this->data[i*this->yn + j];
}

void Grid::setdata(float val, int i, int j)
{
  this->data[i*this->yn + j]  = val;
}

void Grid::setdata_block(float val, int a, int b, int len)
{
  for(int i=a;i<a+len;i++)
  {
    for(int j=b;j<b+len;j++)
    {
      this->data[i*this->yn + j] = val;
    }
  }
}

Grid Grid::clone()
{
  Grid c = Grid(origin_x, origin_y, xn, yn, dx, dy);
  return c;
}
