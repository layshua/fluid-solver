#include "solver.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define IX(i,j) ((i)+(N+2)*(j))

#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}



/** Adds new source of either velocity or density from UI */
void add_source ( int N, float * x, float * s, float dt )
{
    int i, size=(N+2)*(N+2);
    for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}

/** Sets boundaries for whatever is passed as argument. This is used to
 * apply obstacles such as lattice boundaries or any other objects */
void set_bnd ( int N, int b, float * x)
{
    int i;
    int j;

    FOR_EACH_CELL
            //bounce off obstacles
            if(obstacle[IX(i-1, j)])
            x[IX(i ,j)] = b==1 ? -x[IX(i+1 ,j)] : x[IX(i+1 ,j)];
    if(obstacle[IX(i+1, j)])
        x[IX(i ,j)] = b==1 ? -x[IX(i-1 ,j)] : x[IX(i-1 ,j)];
    if(obstacle[IX(i, j-1)])
        x[IX(i ,j)] = b==2 ? -x[IX(i ,j+1)] : x[IX(i ,j+1)];
    if(obstacle[IX(i, j+1)])
        x[IX(i ,j)] = b==2 ? -x[IX(i ,j-1)] : x[IX(i ,j-1)];

    //bounce off domain boundaries
    x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
    x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
    x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
    x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
    END_FOR
}

/** Finds inverse matrix using jacobi method
* TODO: Replace with successive overrelaxation (O(N^(3/2)) instead of O(N^2). in paralell sqrt(N)
* instead of N as for jacobi. It can be paralellized directly)
* Also try to improve it to FFT based method (O(N * log(N)))
* http://www.cs.berkeley.edu/~demmel/cs267/lecture24/lecture24.html
*/
void lin_solve ( int N, int b, float * x, float * x0, float a, float c )
{
    int i, j, k;
    float w = 2.f/(1.f + sinf(3.1415926535/((N+2)+1.f)));
    float* nx = (float*)calloc(sizeof(float), (N+2)*(N+2));
    //SOR solver, doesn't work :(
    for ( k=0 ; k<5 ; k++ )
    {
        //imagine it's a checkerboard and here we do it for blacks
        for ( i=1 ; i<N ; i+=2 )
        {
            for ( j=1 ; j<N ; j+=1 )
            {
                //nx[IX(i,j)] = (a*(nx[IX(i-1,j)]+x[IX(i+1,j)]+nx[IX(i,j-1)]+x[IX(i,j+1)]) + x0[IX(i,j)] )/c; //worked with this so at least paralelization is achieved
                nx[IX(i,j)]=x[IX(i,j)] + w*((x[IX(i-1, j)] + x[IX(i+1, j)] + x[IX(i, j-1)] + x[IX(i, j+1)])*a + x0[IX(i,j)] - 4.f*x[IX(i,j)])/c;
            }
        }
        //and here for whites
        for ( i=2 ; i<N ; i+=2 )
        {
            for ( j=2 ; j<N ; j+=1 )
            {
                //nx[IX(i,j)] = (a*(nx[IX(i-1,j)]+x[IX(i+1,j)]+nx[IX(i,j-1)]+x[IX(i,j+1)]) + x0[IX(i,j)] )/c;
                nx[IX(i,j)]=x[IX(i,j)] + w*((x[IX(i-1, j)] + x[IX(i+1, j)] + x[IX(i, j-1)] + x[IX(i, j+1)])*a + x0[IX(i,j)] - 4.f*x[IX(i,j)])/c;
            }
        }
        memcpy(x, nx, (N+2)*(N+2)*sizeof(float));
        set_bnd ( N, b, nx );

    }


    //this works
    //k can be increased for higher precission, however this is very time consuming!
    //    for ( k=0 ; k<5 ; k++ ) {
    //        FOR_EACH_CELL
    //                nx[IX(i,j)] = (a*(nx[IX(i-1,j)]+x[IX(i+1,j)]+nx[IX(i,j-1)]+x[IX(i,j+1)]) + x0[IX(i,j)] )/c;
    //        END_FOR

    //        memcpy(x, nx, (N+2)*(N+2)*sizeof(float));
    //        set_bnd ( N, b, x );
    //    }
}

/** Simply diffuses whatever is passed as the argument */
void diffuse ( int N, int b, float * x, float * x0, float diff, float dt )
{
    float a=dt*diff*N*N;
    lin_solve ( N, b, x, x0, a, 1+4*a );
}

/** performs advection (movement of particles) by means of linear interpolation.
 * This function interpolates how much of current value of current particle should
 * propagate to another particle whose coordinates are determined by the linear
 * interpolation */
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{\
    //TODO: Replace this shit with MacCormack method that performs two
    //intermediate semi-Lagrangian advection steps.
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt*N;

    FOR_EACH_CELL
            //determine where to propagate current value
            x = i-dt0*u[IX(i,j)];
    y = j-dt0*v[IX(i,j)];

    //Move no further than lattice boundaries
    if (x<0.5f)
        x=0.5f;
    if (x>N+0.5f)
        x=N+0.5f;
    if (y<0.5f)
        y=0.5f;
    if (y>N+0.5f)
        y=N+0.5f;

    //linear interpolation
    i0=(int)x;
    i1=i0+1;
    j0=(int)y;
    j1=j0+1;

    //how much of values should move (s0 start, s1 end for x and same for y with t)
    s1 = x-i0;
    s0 = 1-s1;
    t1 = y-j0;
    t0 = 1-t1;

    //set new values
    d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)])+
            s1*(t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);
    END_FOR
            set_bnd ( N, b, d );
}

/**
 * To obtain velocity field for incompressible flow, we subtract
 * velocity gradient field from current velocities.
*/
void project ( int N, float * u, float * v, float * p, float * div )
{
    int i, j;

    FOR_EACH_CELL
            div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
    p[IX(i,j)] = 0;
    END_FOR

            //neccessary to prevent simulation from blowing up
            set_bnd ( N, 0, div );
    set_bnd ( N, 0, p );

    //find inverse matrix to obtain velocity gradient field
    lin_solve ( N, 0, p, div, 1, 4 );

    //subtract gradient field from current velocities to
    FOR_EACH_CELL
            u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
    v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
    END_FOR

            //neccessary to prevent simulation from blowing up
            set_bnd ( N, 1, u );
    set_bnd ( N, 2, v );
}

/** Performs velocity step for the simulation (this is performed first) */
void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt )
{
    add_source ( N, u, u0, dt );
    add_source ( N, v, v0, dt );

    SWAP ( u0, u );
    SWAP ( v0, v );

    diffuse ( N, 1, u, u0, visc, dt );
    diffuse ( N, 2, v, v0, visc, dt );
    project ( N, u, v, u0, v0 );

    SWAP ( u0, u );
    SWAP ( v0, v );

    advect ( N, 1, u, u0, u0, v0, dt );
    advect ( N, 2, v, v0, u0, v0, dt );
    project ( N, u, v, u0, v0 );
}

/** Performs density step for the simulation (performed after velocity step) */
void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt )
{
    add_source ( N, x, x0, dt );
    SWAP ( x0, x );

    diffuse ( N, 0, x, x0, diff, dt );
    SWAP ( x0, x );

    advect ( N, 0, x, x0, u, v, dt );
}


