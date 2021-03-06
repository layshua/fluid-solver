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
void set_bnd ( int N, int b, float * x, bool velocity)
{
    int i;
    int j;

    FOR_EACH_CELL
    //bounce off obstacles

    if(b!=0)
    {
        if(obstacle[IX(i-1, j)])
            x[IX(i ,j)] = b==1 ? -x[IX(i+1 ,j)] : x[IX(i+1 ,j)];
        if(obstacle[IX(i+1, j)])
            x[IX(i ,j)] = b==1 ? -x[IX(i-1 ,j)] : x[IX(i-1 ,j)];
        if(obstacle[IX(i, j-1)])
            x[IX(i ,j)] = b==2 ? -x[IX(i ,j+1)] : x[IX(i ,j+1)];
        if(obstacle[IX(i, j+1)])
            x[IX(i ,j)] = b==2 ? -x[IX(i ,j-1)] : x[IX(i ,j-1)];

        //use also Moore neigbourhood (decreases framerate by 8%!)
        //...however improves flows past obstacles greatly
        if(obstacle[IX(i+1, j+1)])
            x[IX(i ,j)] =  -x[IX(i-1 ,j-1)];
        if(obstacle[IX(i+1, j-1)])
            x[IX(i ,j)] =  -x[IX(i-1 ,j+1)];
        if(obstacle[IX(i-1, j-1)])
            x[IX(i ,j)] =  -x[IX(i+1 ,j+1)];
        if(obstacle[IX(i-1, j+1)])
            x[IX(i ,j)] =  -x[IX(i+1 ,j-1)];
    }
    //bounce off domain boundaries
//    x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];   //left
//    x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];   //right
//    x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];  //bottom
//    x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];  //top

    END_FOR
}

/** Finds inverse matrix using successive overrelaxation method
*/
void lin_solve ( int N, int b, float * x, float * x0, float a, float c )
{
    int i, j, k;
    //magintude. Fixed value is set but the formula can be used to find more optimal one.
    //it must be between 1 and 2
    float w = 1.333f;

    //A = x
    //b = x0
    for ( k=0 ; k<8 ; k++ ) {
        FOR_EACH_CELL
                //x[IX(i,j)]=x[IX(i,j)] + w*((x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c -x[IX(i,j)] );
                x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]) )/c; //old jacobi for comparison
        END_FOR

        set_bnd ( N, b, x, true );
    }
}

/** Simply diffuses whatever is passed as the argument */
void diffuse ( int N, int b, float * x, float * x0, float diff, float dt )
{
    float a=dt*diff*N*N;
    float c=1.f+4.f*a;
    int i,j;

    float w = 1.333f;

    for(int k = 0; k<8; k++)
    FOR_EACH_CELL
           x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]) )/c; //old jacobi for comparison
    END_FOR
}

/** performs BFECC advection basing on semi-lagrangian method (see function advect) */
void advect_BFECC ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    bool vel = true;

    if(b==0)
        vel=false;

    dt0 = dt*N;

    float *fi_ = (float*)calloc((N+2)*(N+2),sizeof(float));
    float *fi_t = (float*)calloc((N+2)*(N+2),sizeof(float));
    float *fi_n = d0;
    float *fi_n1 = d;

    advect(N, b, fi_t, fi_n, u, v, dt);

    //Step 2: determine fi_
    FOR_EACH_CELL
    //determine from which cell we should propagate density to current cell..
    x = i-dt0*(-u[IX(i,j)]);
    y = j-dt0*(-v[IX(i,j)]);

    //...but no further than lattice boundaries
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

    //s0*t0 is the area of IX(i0,j0) that is taken into account and so on
    s1 = x-i0;
    s0 = 1-s1;
    t1 = y-j0;
    t0 = 1-t1;

    //set new values
    fi_[IX(i,j)] = s0*(t0*fi_t[IX(i0,j0)] + t1*fi_t[IX(i0,j1)])+
                   s1*(t0*fi_t[IX(i1,j0)] + t1*fi_t[IX(i1,j1)]);
    END_FOR
    set_bnd ( N, b, fi_, vel);

    //subtract error comparing backward and forward differentiation
    FOR_EACH_CELL
            fi_t[IX(i,j)] = fi_n[IX(i,j)] + (fi_n[IX(i,j)] - fi_n[IX(i,j)])/2.f;
    END_FOR
    //set_bnd ( N, b, fi_t, vel);

    //Final step - integrate
    advect(N, b, d, fi_t, u, v, dt);
    set_bnd ( N, b, d, vel);

    free(fi_);
    free(fi_t);
}

//OBSOLETE//
/** performs advection (movement of particles) by means of semi lagrangian method.
 * This function traces velocity vectors backwards in order to see how much of each
 * value should propagate to current cell */
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{

    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt*N;

    FOR_EACH_CELL
    //determine from which cell we should propagate density to current cell..
    x = i-dt0*u[IX(i,j)];
    y = j-dt0*v[IX(i,j)];

    //...but no further than lattice boundaries
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

    //s0*t0 is the area of IX(i0,j0) that is taken into account and so on
    s1 = x-i0;
    s0 = 1-s1;
    t1 = y-j0;
    t0 = 1-t1;

    //set new values
    d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)])+
                 s1*(t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);
    END_FOR

    set_bnd ( N, b, d, true );
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
    set_bnd ( N, 0, div, true );
    set_bnd ( N, 0, p, true );

    //find inverse matrix to obtain velocity gradient field
    lin_solve ( N, 0, p, div, 1, 4 );

    //subtract gradient field from current velocities to
    FOR_EACH_CELL
            u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
            v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
    END_FOR

    //neccessary to prevent simulation from blowing up
    set_bnd ( N, 1, u, true);
    set_bnd ( N, 2, v, true );
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

    advect_BFECC ( N, 1, u, u0, u0, v0, dt );
    advect_BFECC ( N, 2, v, v0, u0, v0, dt );
    project ( N, u, v, u0, v0 );
}

/** Performs density step for the simulation (performed after velocity step) */
void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt )
{
    add_source ( N, x, x0, dt );
    SWAP ( x0, x );

    diffuse ( N, 0, x, x0, diff, dt );
    SWAP ( x0, x );
    advect_BFECC ( N, 0, x, x0, u, v, dt );
}
