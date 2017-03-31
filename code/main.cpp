#include<iostream>
#include "Grid.h"

using namespace std;

int main()
{
    cout<<"Fluid Simulation\n";
    float o[3]={0,0,0};
    Grid gr(5,4,1,o);
    cout<<gr.getVolume()<<endl;
    cout<<gr.getX()<<", "<<gr.getY()<<", "<<gr.getZ()<<endl;
    cout<<gr.getOrigin()[0]<<", "<<gr.getOrigin()[0]<<", "<<gr.getOrigin()[0]<<endl;
    return 0;
}
