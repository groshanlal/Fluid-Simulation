#include <iostream>
#include <vector>
#include <cmath>
#include "GL/glut.h"
#include "Grid.hpp"
#include "time.h"
#include "FluidSimulation.hpp"
#include "Particle.hpp"

using namespace std;

//Grid t = Grid(4, 4, 4, 3, 10, 10);
//Grid t;
float FPS = 10;
//FluidSimulation f = FluidSimulation(8, 8, 9, 10, 1/FPS);
FluidSimulation f = FluidSimulation(8, 8, 81, 1, 1/FPS);

int moving, startx, starty;
float zoom = 1.0;
int mouseButton = 0;
float panx = 0.0, pany = 0.0;
bool show_grid = false;
bool show_grid_lines = false;
bool show_particles = true;

void mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
	mouseButton = button;
    moving = 1;
    startx = x;
    starty = y;
  }
  if (state == GLUT_UP) {
	mouseButton = button;
    moving = 0;
  }
}

void motion(int x, int y)
{
  if (moving)
  {
  	if(mouseButton==GLUT_LEFT_BUTTON)
  		zoom += ((y-starty)*0.001);

    startx = x;
    starty = y;
  	glutPostRedisplay();
  }
}

void handleSpecialKeypress(int key, int x, int y)
{
 switch (key)
 {
   case GLUT_KEY_LEFT:
   panx = panx - 0.001;
   break;

   case GLUT_KEY_RIGHT:
   panx = panx + 0.001;
   break;

   case GLUT_KEY_UP:
   pany = pany + 0.001;
   break;

   case GLUT_KEY_DOWN:
   pany = pany - 0.001;
   break;
 }
}

void keyboard (unsigned char key, int x, int y)
{
	switch (key) {
    case 'Q':
		case 'q':
			exit(0);
		break;

    case 'P':
    case 'p':
      show_particles = !show_particles;
    break;

    case 'G':
    case 'g':
      show_grid = !show_grid;
    break;

    case 'L':
    case 'l':
      show_grid_lines = !show_grid_lines;
    break;

		default:
		break;
	}

	glutPostRedisplay();
}


void draw(Grid t, vector<Particle> p)
{
  glPushMatrix();
  //glTranslatef(-5,-5,0);
  glScalef(0.01,0.01,0.01);

  glLineWidth(5.0);
  glBegin(GL_LINE_LOOP);
  glColor3f (0.85, 0.65, 0.12);
    glVertex2f(t.origin_x, t.origin_y);
    glVertex2f(t.origin_x + t.dx * t.xn, t.origin_y);
    glVertex2f(t.origin_x + t.dx * t.xn, t.origin_y + t.dy * t.yn);
    glVertex2f(t.origin_x , t.origin_y + t.dy * t.yn);
  glEnd();
  glLineWidth(1.0);



  if(show_grid)
  {
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
  }

  if(show_grid_lines)
  {
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
  }

  if(show_particles)
  {
    for(vector<Particle>::iterator it = p.begin(); it != p.end(); it++)
    {
      glColor3f (1.0, 0.2, 0.2);
      glBegin(GL_POLYGON);
        glVertex2f(t.origin_x + (it->x+0.5-0.1)* t.dx, t.origin_y + (it->y+0.5-0.1)* t.dy);
        glVertex2f(t.origin_x + (it->x+0.5+0.1)* t.dx, t.origin_y + (it->y+0.5-0.1)* t.dy);
        glVertex2f(t.origin_x + (it->x+0.5+0.1)* t.dx, t.origin_y + (it->y+0.5+0.1)* t.dy);
        glVertex2f(t.origin_x + (it->x+0.5-0.1)* t.dx, t.origin_y + (it->y+0.5+0.1)* t.dy);
      glEnd();

    }
  }


  glPopMatrix();

}

void display()
{
    glClear (GL_COLOR_BUFFER_BIT);
    glClearColor(1.0, 1.0, 1.0, 1.0);

    glShadeModel(GL_SMOOTH);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslatef(panx,pany,0);
    glScalef(zoom,zoom,zoom);


    draw(f.density, f.ParticleSystem);
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
    for(int i=0;i<40;i++)
    {
      f.advect();
      f.applyAcceleration(0,-1);
      cout<<"Iteration: "<<i<<endl;
      print_grid(f.vy);
      f.pressureSolve();
      print_grid(f.pressure);
      print_grid(f.vy);
      cout<<endl;
    }
    return 0;
*/
/*
    Particle p = Particle_around(10,20);
    cout<<p.x<<", "<<p.y<<endl;

    p.set_velocity(4,5);
    cout<<p.get_vx()<<", "<<p.get_vy()<<endl;
    cout<<p.get_vx(0.1,0.1)<<", "<<p.get_vy(0.1,0.1)<<endl;
    cout<<p.get_vx(0.9,0.9)<<", "<<p.get_vy(0.9,0.9)<<endl;
    cout<<p.kernel(0,0)<<endl;
    cout<<p.kernel(0.5,0.25)<<endl;
    cout<<p.kernel(1.5,0.25)<<endl;
    cout<<p.kernel(0.5,1.25)<<endl;
    return 1;
*/
    cout<<"Fluid Simulator\n";
    cout<<"Press arrow keys to pan\n";
    cout<<"Drag left mouse-button to zoom\n";
    cout<<"Press P to show/hide fluid-particles\n";
    cout<<"Press G to show/hide fluid-grids\n";
    cout<<"Press L to show/hide grid-lines\n";
    cout<<"Press Q to exit\n";

    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (600, 600);
    glutInitWindowPosition (10, 10);
    glutCreateWindow ("Fluid Simulation");
    init ();
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouse);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(handleSpecialKeypress);
    glutMotionFunc(motion);
    glutMainLoop();
    return 0;
}
