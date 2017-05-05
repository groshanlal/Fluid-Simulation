#ifndef PARTICLE_H
#define PARTICLE_H

class Particle
{
  public:
    float x,y;
    float vx,vy;

    Particle(float x,float y);
    void set_velocity(float vx, float vy);
    float kernel(float x, float y);
    float get_vx(float x=0, float y=0);
    float get_vy(float x=0, float y=0);


};

Particle Particle_around(float x,float y);

#endif
