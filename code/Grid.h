#include<vector>

using namespace std;

class Grid
{
    public:
        vector<float> origin[3];
        vector<float> data;
        float dx, dy, dz;
        float X, Y, Z;

    public:
      Grid(int m, int n, int o);
      int getVolume();
      int getX();
      int getY();
      int getZ();
};
