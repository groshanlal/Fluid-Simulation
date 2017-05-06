#include "Particle.hpp"
#include <stdlib.h>
#include <time.h>

Particle::Particle(float x, float y)
{
  this->x = x;
  this->y = y;

  this->vx = 0;
  this->vy = 0;

}

Particle Particle_around(float x,float y)
{
  float r, rx, ry;

  r = (rand() % 998);
  r = r + 1;
  rx = r/1000;

  r = (rand() % 998);
  r = r + 1;
  ry = r/1000;

  return Particle(x + rx - 0.5, y + ry - 0.5);


}

void Particle::set_velocity(float vx, float vy)
{
  this->vx = vx;
  this->vy = vy;
}


void Particle::set_position(float x, float y)
{
  this->x = x;
  this->y = y;
}

float Particle::kernel(float x, float y)
{
  float valx, valy;
  if((x>=0)&&(x<=1))
    valx = 1 - x;
  else if((x>=-1)&&(x<=0))
    valx = 1 + x;
  else
    valx = 0;

  if((y>=0)&&(y<=1))
      valy = 1 - y;
  else if((y>=-1)&&(y<=0))
      valy = 1 + y;
  else
      valy = 0;

  return valx*valy;
}

float Particle::get_vx(float x, float y)
{
  return this->vx * this->kernel(x,y);
}

float Particle::get_vy(float x, float y)
{
  return this->vy * this->kernel(x,y);
}

Particle Particle::clone()
{
  Particle  p = Particle(this->x, this->y);
  p.set_velocity(this->vx, this->vy);
  return p;
}
