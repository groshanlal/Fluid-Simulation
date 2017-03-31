#include<vector>

using namespace std;

class Grid
{
    public:
        float origin[3];
        vector<float> data;
        float dx, dy, dz;
        float X, Y, Z;

    public:
      Grid(int m, int n, int o, float (&orig)[3]);
      float getVolume();
      float getX();
      float getY();
      float getZ();
      float* getOrigin();

};
