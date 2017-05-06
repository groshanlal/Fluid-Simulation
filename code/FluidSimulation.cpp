#include <iostream>
#include <cmath>
#include "FluidSimulation.hpp"
#include "Grid.hpp"
#include "Particle.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stdlib.h>
#include <time.h>

using namespace Eigen;
using namespace std;

FluidSimulation::FluidSimulation(float origin_x, float origin_y, int xn, float dx, float t)
{
  int yn = xn;
  float dy = dx;

  pressure = Grid(origin_x, origin_y, xn, yn, dx, dy);
  density  = Grid(origin_x, origin_y, xn, yn, dx, dy);
  vx       = Grid(origin_x - dx/2, origin_y, xn+1, yn, dx, dy);
  vy       = Grid(origin_x, origin_y - dy/2, xn, yn+1, dx, dy);
  old_vx   = Grid(origin_x - dx/2, origin_y, xn+1, yn, dx, dy);
  old_vy   = Grid(origin_x, origin_y - dy/2, xn, yn+1, dx, dy);

  density.setdata_block(1,(int)xn/3,20,(int)xn/3,(int)xn/3);

  timestep = t;

  srand (time(NULL));

  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if(density.getdata(i,j)!=1)
        continue;

      for(int k=0;k<40;k++)
        ParticleSystem.push_back(Particle_around(i,j));
    }
  }
}

void FluidSimulation::computeNextStep()
{

  float cfl = this->CFL();
  float threshold = 2.0;

  if(cfl>threshold)
    cout<<"CFL greater than "<<threshold<<": "<<cfl<<endl;

  this->advect();
  this->applyAcceleration(0,-6);
  this->pressureSolve();


}

float FluidSimulation::CFL()
{
  float u_max=0;
  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if(density.getdata(i,j)!=1)
          continue;

      float ux = (vx.getdata(i,j) + vx.getdata(i+1,j))/2;
      float uy = (vy.getdata(i,j) + vy.getdata(i,j+1))/2;

      float u = ux*ux + uy*uy;
      if(u>u_max)
        u_max=u;
    }
  }

  u_max = sqrt(u_max);
  return u_max*timestep/density.dx;
}

void FluidSimulation::applyAcceleration(float ax, float ay)
{
  //cout<<"acc\n";

  for(int i=0;i<vx.xn;i++)
  {
    for(int j=0;j<vx.yn;j++)
    {
      if( ((i-1>=0) && (density.getdata(i-1,j) == 1)) || ((i<density.xn))&&(density.getdata(i,j)==1) )
      {
        vx.setdata(vx.getdata(i,j)  + ax*timestep, i, j);
      }
    }
  }

  for(int i=0;i<vy.xn;i++)
  {
    for(int j=0;j<vy.yn;j++)
    {
      if( ((j-1>=0) && (density.getdata(i,j-1) == 1)) || ((j<density.yn))&&(density.getdata(i,j)==1) )
      {
        vy.setdata(vy.getdata(i,j)  + ay*timestep, i, j);
      }
    }
  }
}

void FluidSimulation::advect()
{
  for(vector<Particle>::iterator it = ParticleSystem.begin(); it != ParticleSystem.end() ; it++)
  {
    float px = it->x;
    float py = it->y;

    float gx = floor(px + 0.5);
    float gy = floor(py + 0.5);

    //PIC
    float pvx_pic = (gx + 0.5 - px)*vx.getdata(gx,gy) + (px - gx + 0.5)*vx.getdata(gx+1,gy);
    float pvy_pic = (gy + 0.5 - py)*vy.getdata(gx,gy) + (py - gy + 0.5)*vy.getdata(gx,gy+1);

    //FLIP
    float pvx_flip = it->vx + (gx + 0.5 - px)*(vx.getdata(gx,gy) - old_vx.getdata(gx,gy)) + (px - gx + 0.5)*(vx.getdata(gx+1,gy) - old_vx.getdata(gx+1,gy));
    float pvy_flip = it->vy + (gy + 0.5 - py)*(vy.getdata(gx,gy) - old_vy.getdata(gx,gy)) + (py - gy + 0.5)*(vy.getdata(gx,gy+1) - old_vy.getdata(gx,gy+1));

    float alpha = 0.95;
    //float alpha = 0.00;
    float pvx = alpha*pvx_flip + (1-alpha)*pvx_pic;
    float pvy = alpha*pvy_flip + (1-alpha)*pvy_pic;

    it->set_velocity(pvx, pvy);
  }

  vector<Particle> ps;
  Grid backup_density = density.clone();
  Grid backup_vx = vx.clone();
  Grid backup_vy = vy.clone();
  Grid backup_wx = vx.clone();
  Grid backup_wy = vy.clone();

  for(vector<Particle>::iterator it = ParticleSystem.begin(); it != ParticleSystem.end(); it++)
  {
    Particle p = it->clone();
    float pvx = it->vx;
    float pvy = it->vy;
    float px  = it->x;
    float py  = it->y;
    px = px + pvx*timestep/density.dx;
    py = py + pvy*timestep/density.dy;

    while(px<0)
      px++;
    while(py<0)
      py++;
    while(px>density.xn-1)
      px--;
    while(py>density.yn-1)
      py--;

    p.set_position(px,py);
    ps.push_back(p);

    float gx = floor(px + 0.5);
    float gy = floor(py + 0.5);

    backup_density.setdata(1,gx,gy);
  }

  for(vector<Particle>::iterator it = ps.begin(); it != ps.end(); it++)
  {
    float px = it->x;
    float py = it->y;

    float gx = floor(px + 0.5);
    float gy = floor(py + 0.5);

    backup_vx.setdata(backup_vx.getdata(gx,gy) + it->get_vx(gx-0.5-px, gy-py),gx,gy);
    backup_vx.setdata(backup_vx.getdata(gx+1,gy) + it->get_vx(gx+0.5-px, gy-py),gx+1,gy);

    backup_vy.setdata(backup_vy.getdata(gx,gy) + it->get_vy(gx-px, gy-0.5-py),gx,gy);
    backup_vy.setdata(backup_vy.getdata(gx,gy+1) + it->get_vy(gx-px, gy+0.5-py),gx,gy+1);

    backup_wx.setdata(backup_wx.getdata(gx,gy) + it->kernel(gx-0.5-px, gy-py),gx,gy);
    backup_wx.setdata(backup_wx.getdata(gx+1,gy) + it->kernel(gx+0.5-px, gy-py),gx+1,gy);

    backup_wy.setdata(backup_wy.getdata(gx,gy) + it->kernel(gx-px, gy-0.5-py),gx,gy);
    backup_wy.setdata(backup_wy.getdata(gx,gy+1) + it->kernel(gx-px, gy+0.5-py),gx,gy+1);
  }

  for(int i=0;i<backup_vx.xn;i++)
  {
    for(int j=0;j<backup_vx.yn;j++)
    {
      if(backup_wx.getdata(i,j)==0)
        continue;

      backup_vx.setdata(backup_vx.getdata(i,j)/backup_wx.getdata(i,j),i,j);

    }
  }

  for(int i=0;i<backup_vy.xn;i++)
  {
    for(int j=0;j<backup_vy.yn;j++)
    {
      if(backup_wy.getdata(i,j)==0)
        continue;

      backup_vy.setdata(backup_vy.getdata(i,j)/backup_wy.getdata(i,j),i,j);

    }
  }

  ParticleSystem = ps;
  density = backup_density;
  vx = backup_vx;
  vy = backup_vy;
  old_vx = backup_vx;
  old_vy = backup_vy;
}

void FluidSimulation::pressureSolve()
{
  //cout<<"pressure-solve\n";
  int n=0;

  typedef Eigen::Triplet<double> T;
  vector<T> tripletList;
  vector<T> mapList;

  SparseMatrix<int> map(density.xn,density.yn);
  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if(density.getdata(i,j)==1)
      {
        mapList.push_back(T(i,j,n));
        n++;
      }
    }
  }
  map.setFromTriplets(mapList.begin(), mapList.end());

  float rho = 1.00;
  //cout<<"No. of grid cells, density: "<<n<<", "<<rho<<endl;

  VectorXd x(n), b(n);
  SparseMatrix<double> A(n,n);


  float scale = 1 / density.dx;
  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if(density.getdata(i,j)!=1)
        continue;

      b[map.coeff(i,j)] = -scale * (vx.getdata(i+1,j)-vx.getdata(i,j)+vy.getdata(i,j+1)-vy.getdata(i,j));

      if(i==0)
        b[map.coeff(i,j)] = b[map.coeff(i,j)] - scale * (vx.getdata(i,j) - 0);

      if(i+1==density.xn)
        b[map.coeff(i,j)] = b[map.coeff(i,j)] + scale * (vx.getdata(i+1,j) - 0);

      if(j==0)
        b[map.coeff(i,j)] = b[map.coeff(i,j)] - scale * (vy.getdata(i,j) - 0);

      if(j+1==density.yn)
        b[map.coeff(i,j)] = b[map.coeff(i,j)] + scale * (vy.getdata(i,j+1) - 0);
    }
  }

  scale = timestep / (rho * density.dx * density.dx);
  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if((density.getdata(i,j)==1)&&(i+1<density.xn)&&(density.getdata(i+1,j)==1))
      {
        tripletList.push_back(T(map.coeff(i,j),map.coeff(i,j),scale));
        tripletList.push_back(T(map.coeff(i+1,j),map.coeff(i+1,j),scale));
        tripletList.push_back(T(map.coeff(i,j),map.coeff(i+1,j), -scale));
        tripletList.push_back(T(map.coeff(i+1,j),map.coeff(i,j), -scale));
      }
      else if((density.getdata(i,j)==1)&&(i+1<density.xn)&&(density.getdata(i+1,j)==0))
        tripletList.push_back(T(map.coeff(i,j),map.coeff(i,j),scale));
      else if((density.getdata(i,j)==0)&&(i+1<density.xn)&&(density.getdata(i+1,j)==1))
        tripletList.push_back(T(map.coeff(i+1,j),map.coeff(i+1,j),scale));

      if((density.getdata(i,j)==1)&&(j+1<density.yn)&&(density.getdata(i,j+1)==1))
      {
        tripletList.push_back(T(map.coeff(i,j),map.coeff(i,j),scale));
        tripletList.push_back(T(map.coeff(i,j+1),map.coeff(i,j+1),scale));
        tripletList.push_back(T(map.coeff(i,j),map.coeff(i,j+1), -scale));
        tripletList.push_back(T(map.coeff(i,j+1),map.coeff(i,j), -scale));
      }
      else if((density.getdata(i,j)==1)&&(j+1<density.yn)&&(density.getdata(i,j+1)==0))
        tripletList.push_back(T(map.coeff(i,j),map.coeff(i,j),scale));
      else if((density.getdata(i,j)==0)&&(j+1<density.yn)&&(density.getdata(i,j+1)==1))
        tripletList.push_back(T(map.coeff(i,j+1),map.coeff(i,j+1),scale));

    }
  }

  A.setFromTriplets(tripletList.begin(), tripletList.end());

  ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
  //cg.setMaxIterations(10);
  cg.compute(A);
  x = cg.solve(b);
  //cout << "#iterations:     " << cg.iterations() << endl;
  //cout << "estimated error: " << cg.error()      << endl;
  //cout << "A: \n" << A << endl;
  //cout << "pressure: \n" << x << endl;
  //cout << "b: \n" << b << endl;

  for(int i=0;i<pressure.xn;i++)
  {
    for(int j=0;j<pressure.yn;j++)
    {
      if(density.getdata(i,j)!=1)
        continue;

      pressure.setdata(x(map.coeff(i,j)),i,j);
    }
  }


  scale = timestep / (rho*density.dx);
  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if(density.getdata(i,j)!=1)
        continue;

      //cout<<vy.getdata(i,j)<<", "<<scale<<", "<<pressure.getdata(i,j)<<endl;

      vx.setdata(vx.getdata(i,j) -  scale*pressure.getdata(i,j) ,i,j);
      vx.setdata(vx.getdata(i+1,j) +  scale*pressure.getdata(i,j) ,i+1,j);
      vy.setdata(vy.getdata(i,j) -  scale*pressure.getdata(i,j) ,i,j);
      vy.setdata(vy.getdata(i,j+1) +  scale*pressure.getdata(i,j) ,i,j+1);

      if(i==0)
        vx.setdata(0,i,j);

      if(i==density.xn-1)
        vx.setdata(0,i+1,j);

      if(j==0)
        vy.setdata(0,i,0);

      if(j==density.yn-1)
        vy.setdata(0,i,j+1);



    }
  }

  //cout<<"A: \n"<<A<<endl;
  //cout<<"x: \n"<<x<<endl;
  //cout<<"b: \n"<<b<<endl;

}
