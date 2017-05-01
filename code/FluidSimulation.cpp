#include <iostream>
#include <cmath>
#include "FluidSimulation.hpp"
#include "Grid.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

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

  density.setdata_block(1,(int)xn/3,(int)xn/3,(int)xn/3);

  timestep = t;
}

void FluidSimulation::computeNextStep()
{
  this->advect();
  this->applyAcceleration(0,-1.333);
  this->pressureSolve();

}
void FluidSimulation::applyAcceleration(float ax, float ay)
{
  for(int i=0;i<vx.xn;i++)
  {
    for(int j=0;j<vx.yn;j++)
    {
      if( ((i-1>0) && (density.getdata(i-1,j) == 1)) || ((i<density.xn))&&(density.getdata(i,j)==1) )
      {
        vx.setdata(vx.getdata(i,j)  + ax*timestep, i, j);
      }
    }
  }

  for(int i=0;i<vy.xn;i++)
  {
    for(int j=0;j<vy.yn;j++)
    {
      if( ((j-1>0) && (density.getdata(i,j-1) == 1)) || ((j<density.yn))&&(density.getdata(i,j)==1) )
      {
        vy.setdata(vy.getdata(i,j)  + ay*timestep, i, j);
      }
    }
  }
}

void FluidSimulation::advect()
{
  Grid backup_density = density.clone();
  //Grid cent_vx = density.clone();
  //Grid cent_vy = density.clone();
  Grid backup_vx = vx.clone();
  Grid backup_vy = vy.clone();

  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if(density.getdata(i,j) != 1)
        continue;

      float grid_vel_x = vx.getdata(i,j) + vx.getdata(i+1,j);
      grid_vel_x = grid_vel_x/2;

      float grid_vel_y = vy.getdata(i,j) + vy.getdata(i,j+1);
      grid_vel_y = grid_vel_y/2;

      float pos_x = i + grid_vel_x*timestep/vx.dx;
      float pos_y = j + grid_vel_y*timestep/vy.dy;

      int a = floor(pos_x+0.5);
      int b = floor(pos_y+0.5);

      //cout<<i<<", "<<j<<" goes to "<<pos_y<<", "<<b<<endl;

      backup_density.setdata(1, a, b);
    }
  }

  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if(density.getdata(i,j) != 1)
        continue;
      if(backup_density.getdata(i,j) != 1)
        continue;

      float grid_vel_x = vx.getdata(i,j) + vx.getdata(i+1,j);
      grid_vel_x = grid_vel_x/2;

      float grid_vel_y = vy.getdata(i,j) + vy.getdata(i,j+1);
      grid_vel_y = grid_vel_y/2;

      float pos_x = i - grid_vel_x*timestep/vx.dx;
      float pos_y = j - grid_vel_y*timestep/vy.dy;


      int a = floor(pos_x+0.5);
      int b = floor(pos_y+0.5);

      float vvxx = ((a+0.5)-pos_x)*vx.getdata(a,b) + (pos_x-(a-0.5))*vx.getdata(a+1,b);
      float vvyy = ((b+0.5)-pos_y)*vy.getdata(a,b) + (pos_y-(b-0.5))*vy.getdata(a,b+1);

      //cout<<i<<", "<<j<<" gets vel from "<<pos_y<<", "<<b<<": "<<vvyy<<endl;

      backup_vx.setdata(backup_vx.getdata(i,j) + vvxx,i,j);
      backup_vx.setdata(backup_vx.getdata(i+1,j) + vvxx,i+1,j);
      backup_vy.setdata(backup_vy.getdata(i,j) + vvyy,i,j);
      backup_vy.setdata(backup_vy.getdata(i,j+1) + vvyy,i,j+1);
    }
  }

  for(int i=0;i<backup_vx.xn;i++)
  {
    for(int j=0;j<backup_vx.yn;j++)
    {
      int count=0;
      if((i-1>0) && (backup_density.getdata(i-1,j) == 1))
        count++;
      if((i<backup_density.xn)&&(backup_density.getdata(i,j)==1))
        count++;
      if(count>0)
        backup_vx.setdata(backup_vx.getdata(i,j)/count, i, j);
    }
  }

  for(int i=0;i<backup_vy.xn;i++)
  {
    for(int j=0;j<backup_vy.yn;j++)
    {
      int count=0;
      if((j-1>0) && (backup_density.getdata(i,j-1) == 1))
        count++;
      if((j<backup_density.yn)&&(backup_density.getdata(i,j)==1))
        count++;
      if(count>0)
        backup_vy.setdata(backup_vy.getdata(i,j)/count, i, j);
    }
  }



  density = backup_density;
  vx = backup_vx;
  vy = backup_vy;


}

void FluidSimulation::pressureSolve()
{
  float rho = 1;
  int n = density.xn * density.xn;
  VectorXd x(n), b(n);
  SparseMatrix<double> A(n,n);

  typedef Eigen::Triplet<double> T;
  vector<T> tripletList;


  float rhs[density.xn][density.yn];
  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      rhs[i][j]=0;
    }
  }


  float scale = 1 / density.dx;
  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if(density.getdata(i,j)!=1)
        continue;
      rhs[i][j] = -scale * (vx.getdata(i+1,j)-vx.getdata(i,j)+vy.getdata(i,j+1)-vy.getdata(i,j));

      if(i==0)
        rhs[i][j] = rhs[i][j] - scale * (vx.getdata(i,j) - 0);

      if(i+1==density.xn)
        rhs[i][j] = rhs[i][j] + scale * (vx.getdata(i+1,j) - 0);

      if(j==0)
        rhs[i][j] = rhs[i][j] - scale * (vy.getdata(i,j) - 0);

      if(j+1==density.yn)
        rhs[i][j] = rhs[i][j] - scale * (vy.getdata(i,j+1) - 0);
    }
  }

  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      b[i*density.xn  + j] = rhs[i][j];
    }
  }

  scale = timestep / (rho * density.dx * density.dx);
  for(int i=0;i<density.xn;i++)
  {
    for(int j=0;j<density.yn;j++)
    {
      if((density.getdata(i,j)==1)&&(i+1<density.xn)&&(density.getdata(i+1,j)==1))
      {
        tripletList.push_back(T(i*density.xn+j,i*density.xn+j,scale));
        tripletList.push_back(T((i+1)*density.xn+j,(i+1)*density.xn+j,scale));
        tripletList.push_back(T(i*density.xn+j,(i+1)*density.xn+j, -scale));
        tripletList.push_back(T((i+1)*density.xn+j,i*density.xn+j, -scale));
      }
      else if((density.getdata(i,j)==1)&&(i+1<density.xn)&&(density.getdata(i+1,j)==0))
        tripletList.push_back(T(i*density.xn+j,i*density.xn+j,scale));

      if((density.getdata(i,j)==1)&&(j+1<density.yn)&&(density.getdata(i,j+1)==1))
      {
        tripletList.push_back(T(i*density.xn+j,i*density.xn+j,scale));
        tripletList.push_back(T(i*density.xn+(j+1),i*density.xn+(j+1),scale));
        tripletList.push_back(T(i*density.xn+j,i*density.xn+(j+1), -scale));
        tripletList.push_back(T(i*density.xn+(j+1),i*density.xn+j, -scale));
      }
      else if((density.getdata(i,j)==1)&&(j+1<density.yn)&&(density.getdata(i,j+1)==0))
        tripletList.push_back(T(i*density.xn+j,i*density.xn+j,scale));
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

      pressure.setdata(x[i*density.xn+j],i,j);
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

      if(i+1==density.xn)
        vx.setdata(0,i,j);

      if(j==0)
        vy.setdata(0,i,j);

      if(j+1==density.yn)
        vy.setdata(0,i,j);

    }
  }

}
