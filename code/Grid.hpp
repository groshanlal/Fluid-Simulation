#ifndef GRID_H
#define GRID_H

#include<vector>

using namespace std;

class Grid
{
  public:
    float origin_x, origin_y;
    float dx, dy;
    int xn, yn;
    vector<float> data;

    Grid();
    Grid(float origin_x, float origin_y, int xn, int yn, float dx, float dy);
    float getdata(int i, int j);
    void setdata(float val, int i, int j);
    void setdata_block(float val, int a, int b, int len);
    Grid clone();
};

#endif
