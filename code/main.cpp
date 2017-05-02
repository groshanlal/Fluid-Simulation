#include <iostream>
#include "GL/glut.h"
#include "Grid.hpp"
#include "time.h"
#include "FluidSimulation.hpp"

using namespace std;

//Grid t = Grid(4, 4, 4, 3, 10, 10);
//Grid t;
float FPS = 10.0;
//FluidSimulation f = FluidSimulation(8, 8, 9, 10, 1/FPS);
FluidSimulation f = FluidSimulation(8, 8, 80, 1, 1/FPS);


void draw(Grid t)
{
  glPushMatrix();
  glScalef(0.01,0.01,0.01);

  for(int i=0;i<t.xn;i++)
  {
    for(int j=0;j<t.yn;j++)
    {
      if(t.getdata(i,j)==1)
      {
        glColor3f (0.2, 0.2, 1.0);
        glBegin(GL_POLYGON);
          glVertex2f(t.origin_x + i* t.dx, t.origin_y + j* t.dy);
          glVertex2f(t.origin_x + (i+1)* t.dx, t.origin_y + j* t.dy);
          glVertex2f(t.origin_x + (i+1)* t.dx, t.origin_y + (j+1)* t.dy);
          glVertex2f(t.origin_x + i* t.dx, t.origin_y + (j+1)* t.dy);
        glEnd();
      }
    }
  }


  for(int i=0;i<=t.xn;i++)
  {
    glColor3f (0.0, 0.0, 0.0);
    glBegin(GL_LINES);
      glVertex2f(t.origin_x + i*t.dx, t.origin_y);
      glVertex2f(t.origin_x + i*t.dx, t.origin_y + t.yn*t.dy);
    glEnd();
  }

  for(int i=0;i<=t.yn;i++)
  {
    glColor3f (0.0, 0.0, 0.0);
    glBegin(GL_LINES);
      glVertex2f(t.origin_x , t.origin_y + i*t.dy);
      glVertex2f(t.origin_x + t.xn*t.dx , t.origin_y + i*t.dy);
    glEnd();
  }

  glPopMatrix();


}
void display()
{
    glClear (GL_COLOR_BUFFER_BIT);

    draw(f.density);
    glFlush ();
}

void init (void)
{
/*  select clearing (background) color       */
    glClearColor (1.0, 1.0, 1.0, 1.0);

/*  initialize viewing values  */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
}


void idle()
{
    clock_t current_time = clock();
    clock_t end_time = current_time + CLOCKS_PER_SEC/FPS;


    while(current_time<end_time)
        current_time=clock();

    f.computeNextStep();
    glutPostRedisplay();
}

void print_grid(Grid t)
{
  for(int i=t.yn-1;i>=0;i--)
  {
    for(int j=0;j<t.xn;j++)
    {
      cout<<t.getdata(j,i)<<", ";
    }
    cout<<endl;
  }
  cout<<endl;
}


/*
 *  Declare initial window size, position, and display mode
 *  (single buffer and RGBA).  Open window with "hello"
 *  in its title bar.  Call initialization routines.
 *  Register callback function to display graphics.
 *  Enter main loop and process events.
 */
int main(int argc, char** argv)
{
/*
    for(int i=0;i<35;i++)
    {
      print_grid(f.vy);
      print_grid(f.pressure);
      f.advect();
      f.applyAcceleration(0,-1);
      //f.pressureSolve();
      cout<<endl;
    }
    return 0;
*/

    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (600, 600);
    glutInitWindowPosition (10, 10);
    glutCreateWindow ("Fluid Simulation");
    init ();
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMainLoop();
    return 0;
}
