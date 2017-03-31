#include<iostream>
#include "Grid.h"

using namespace std;

int main()
{
    cout<<"Fluid Simulation\n";
    Grid gr(5,4,1);
    cout<<gr.getVolume()<<endl;
    cout<<gr.getX()<<", "<<gr.getY()<<", "<<gr.getZ()<<endl;
    return 0;
}
