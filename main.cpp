#include <stdlib.h>
#include <stdio.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <time.h>
#include <stdlib.h>
#include <cmath>

int frame_iter=0;
time_t start, stop;
float seconds;
char framerate[8];

/* macros */

#define IX(i,j) ((i)+(N+2)*(j))
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

/* external definitions (from solver.c) */

extern void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt );
extern void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt );
extern bool *obstacle;

/* global variables */

static int N;
static float dt, diff, visc;
static float force, source, gravity;
static int dvel;

static float * u, * v, * u_prev, * v_prev;
static float * dens, * dens_prev;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;


/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/


static void free_data ( void )
{
    if ( u ) free ( u );
    if ( v ) free ( v );
    if ( u_prev ) free ( u_prev );
    if ( v_prev ) free ( v_prev );
    if ( dens ) free ( dens );
    if ( dens_prev ) free ( dens_prev );
}

static void clear_data ( void )
{
    int i, size=(N+2)*(N+2);

    for ( i=0 ; i<size ; i++ ) {
        u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = obstacle[i] = 0.0f;
    }
}

static int allocate_data ( void )
{
    int size = (N+2)*(N+2);

    u			= (float *) malloc ( size*sizeof(float) );
    v			= (float *) malloc ( size*sizeof(float) );
    u_prev		= (float *) malloc ( size*sizeof(float) );
    v_prev		= (float *) malloc ( size*sizeof(float) );
    dens		= (float *) malloc ( size*sizeof(float) );
    dens_prev	= (float *) malloc ( size*sizeof(float) );
    obstacle    = (bool *)  malloc ( size*sizeof(bool)  );

    if ( !u || !v || !u_prev || !v_prev || !dens || !dens_prev || !obstacle ) {
        fprintf ( stderr, "cannot allocate data\n" );
        return ( 0 );
    }

    return ( 1 );
}


/*
  ----------------------------------------------------------------------
   OpenGL specific drawing routines
  ----------------------------------------------------------------------
*/

static void pre_display ( void )
{
    glViewport ( 0, 0, win_x, win_y );
    glMatrixMode ( GL_PROJECTION );
    glLoadIdentity ();
    gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT );
}

static void post_display ( void )
{
    glutSwapBuffers ();
}

static void draw_velocity ( void )
{
    int i, j;
    float x, y, h;

    h = 1.0f/N;

    glColor3f ( 1.0f, 1.0f, 1.0f );
    glLineWidth ( 1.0f );

    glBegin ( GL_LINES );

    for ( i=1 ; i<=N ; i++ ) {
        x = ((float)i-0.5f)*h;
        for ( j=1 ; j<=N ; j++ ) {
            y = (j-0.5f)*h;
            glColor3f(1,0,0);
            glVertex2f ( x, y );
            glColor3f(0,3.f*fabs(u[IX(i,j)] + v[IX(i,j)]), 0);
            glVertex2f ( x+u[IX(i,j)], y+v[IX(i,j)] );
        }
    }

    glEnd ();
}

static void draw_density ( void )
{
    int i, j;
    float x, y, h, d00, d01;

    h = 1.0f/N;

    glBegin ( GL_POINTS );

    for ( i=0 ; i<=N ; i++ ) {
        x = (i-0.5f)*h;
        for ( j=0 ; j<=N ; j++ ) {
            y = (j-0.5f)*h;
            d00 = dens[IX(i,j)];
            d01 = fabs(u[IX(i,j)] + v[IX(i,j)])/2.f;
            glColor4f (2.f-d00, d01*d00, d01/5.f, d00); glVertex2f ( x, y );
        }
    }

    glEnd ();
}

static void draw_obstacle ( void )
{
    int i, j;
    float x, y, h, o00;
    h = 1.0f/N;

    glBegin ( GL_POINTS );

    for ( i=0 ; i<=N ; i++ ) {
        x = (i-0.5f)*h;
        for ( j=0 ; j<=N ; j++ ) {
            y = (j-0.5f)*h;
            o00 = obstacle[IX(i,j)];
            glColor4f ( 0, 0, o00-0.5f, o00); glVertex2f ( x, y );
        }
    }

    glEnd ();
}

/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/

static void get_from_UI ( float * d, float * u, float * v )
{
    int i, j, size = (N+2)*(N+2);

    for ( i=0 ; i<size ; i++ ) {
        u[i] = v[i] = d[i] = 0.0f;
    }

    if ( !mouse_down[0] && !mouse_down[2] && !mouse_down[1]) return;

    i = (int)((       mx /(float)win_x)*N+1);
    j = (int)(((win_y-my)/(float)win_y)*N+1);

    if ( i<1 || i>N || j<1 || j>N ) return;

    if ( mouse_down[0] ) {
        u[IX(i,j)] = force * (mx-omx);
        v[IX(i,j)] = force * (omy-my);
    }

    if ( mouse_down[1] ) {
        obstacle[IX(i,j)] = true;
    }

    if ( mouse_down[2] ) {
        if(i+N/6<N)
        {
            for(int k=i; k<i+N/6; k++)
                d[IX(k,j)] = source/5;
        }
        //d[IX(i,j)]=source;
    }

    omx = mx;
    omy = my;

    return;
}

static void apply_gravity(float *d, float *u, float *v)
{
    int size = (N)*(N);
    int i,j;
    for(int i=N/8; i<N-N/8; i++)
        v[IX(i,N/8)]+=100*gravity;

}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/

static void key_func ( unsigned char key, int x, int y )
{
    switch ( key )
    {
    case 'c':
    case 'C':
        clear_data ();
        break;

    case 'q':
    case 'Q':
        free_data ();
        exit ( 0 );
        break;

    case 'v':
    case 'V':
        dvel = !dvel;
        break;
    }
}

static void mouse_func ( int button, int state, int x, int y )
{
    omx = mx = x;
    omx = my = y;

    mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y )
{
    mx = x;
    my = y;
}

static void reshape_func ( int width, int height )
{
    glutSetWindow ( win_id );
    glutReshapeWindow ( width, height );

    win_x = width;
    win_y = height;
}

static void idle_func ( void )
{
    glPointSize((glutGet(GLUT_WINDOW_WIDTH)+glutGet(GLUT_WINDOW_WIDTH))/(1.2f*N));
    get_from_UI ( dens_prev, u_prev, v_prev );
    apply_gravity(dens_prev, u_prev, v_prev);

    vel_step ( N, u, v, u_prev, v_prev, visc, dt );
    dens_step ( N, dens, dens_prev, u, v, diff, dt );

    glutSetWindow ( win_id );
    glutPostRedisplay ();
}

static void show_framerate()
{
    if(frame_iter==0)
        start=clock();
    if(frame_iter<15)
        frame_iter++;
    else
    {
        frame_iter=0;
        stop=clock();
        seconds = (float)(stop - start) / CLOCKS_PER_SEC / 15.f;
        seconds = 1.f/seconds;
        //dt=1.f/seconds; //optional adjustment of physics frames to time rate to keep it realistic
        snprintf(framerate, 8, "%f", seconds);
    }


    glLoadIdentity();
    glColor3f(1, 1, 0.f);     //<-- this line controls the color (now text is yellow)
    glRasterPos2f(-1, 0.95);
    for(int z=0; z<8; z++)
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, (int)framerate[z]);
    glutSwapBuffers();
    glRasterPos2f(0, 0);
}

static void display_func ( void )
{

    show_framerate();
    pre_display ();

    if ( dvel ) draw_velocity ();
    else
    {
        draw_density ();
        draw_obstacle ();
    }
    post_display ();

}


/*
  ----------------------------------------------------------------------
   open_glut_window --- open a glut compatible window and set callbacks
  ----------------------------------------------------------------------
*/

static void open_glut_window ( void )
{
    glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

    glutInitWindowPosition ( 0, 0 );
    glutInitWindowSize ( win_x, win_y );
    win_id = glutCreateWindow ( "Fluid Sim" );

    glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();

    pre_display ();

    glutKeyboardFunc ( key_func );
    glutMouseFunc ( mouse_func );
    glutMotionFunc ( motion_func );
    glutReshapeFunc ( reshape_func );
    glutIdleFunc ( idle_func );
    glutDisplayFunc ( display_func );
}


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{
    glutInit ( &argc, argv );

    if ( argc != 1 && argc != 7 ) {
        fprintf ( stderr, "usage : %s N dt diff visc force source\n", argv[0] );
        fprintf ( stderr, "where:\n" );\
        fprintf ( stderr, "\t N       : grid resolution\n" );
        fprintf ( stderr, "\t dt      : time step\n" );
        fprintf ( stderr, "\t diff    : diffusion rate of the density\n" );
        fprintf ( stderr, "\t visc    : viscosity of the fluid\n" );
        fprintf ( stderr, "\t force   : scales the mouse movement that generate a force\n" );
        fprintf ( stderr, "\t source  : amount of density that will be deposited\n" );
        fprintf ( stderr, "\t gravity : magnitude of gravitational force\n" );
        exit ( 1 );
    }

    if ( argc == 1 ) {
        N = 500;
        dt = 0.005f;
        diff = 0.00001f;
        visc = 0.000001f;
        force = 800.0f;
        source = 1000.0f;
        gravity = 15.f;
        fprintf ( stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g gravity=%g\n",
                  N, dt, diff, visc, force, source, gravity );
    } else {
        N = atoi(argv[1]);
        dt = atof(argv[2]);
        diff = atof(argv[3]);
        visc = atof(argv[4]);
        force = atof(argv[5]);
        source = atof(argv[6]);
        gravity = atof(argv[7]);
    }

    printf ( "\n\nHow to use this demo:\n\n" );
    printf ( "\t Add densities with the right mouse button\n" );
    printf ( "\t Add velocities with the left mouse button and dragging the mouse\n" );
    printf ( "\t Toggle density/velocity display with the 'v' key\n" );
    printf ( "\t Clear the simulation by pressing the 'c' key\n" );
    printf ( "\t Quit by pressing the 'q' key\n" );

    dvel = 0;

    if ( !allocate_data () ) exit ( 1 );
    clear_data ();

    win_x = 512;
    win_y = 512;
    open_glut_window ();

    glutMainLoop ();

    exit ( 0 );
}
