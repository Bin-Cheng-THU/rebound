/**
 * Self-gravitating disc
 *
 * A self-gravitating disc is integrated using
 * the leap frog integrator. Collisions are not resolved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"


void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_create_simulation();
    // Setup constants
    r->integrator    = REB_INTEGRATOR_LEAPFROG;
    //r->collision        = REB_COLLISION_TREE;
    //r->collision_resolve = reb_collision_resolve_hardsphere;
    r->gravity    = REB_GRAVITY_TREE;
    r->boundary    = REB_BOUNDARY_OPEN;
    r->opening_angle2    = 1.5;        // This constant determines the accuracy of the tree code gravity estimate.
    r->G         = 1;        
    r->softening     = 0.02;        // Gravitational softening length
    r->dt         = 3e-2;        // Timestep
    const double boxsize = 10.2;
    reb_configure_box(r,boxsize,1,1,1);

    // Setup particles
    double disc_mass = 2e-1;    // Total disc mass
    int N = 10000;            // Number of particles
    // Initial conditions
    struct reb_particle star = {0};
    star.m         = 1;
    reb_add(r, star);
    for (int i=0;i<N;i++){
        struct reb_particle pt = {0};
        double a    = reb_random_powerlaw(r, boxsize/10.,boxsize/2./1.2,-1.5);
        double phi     = reb_random_uniform(r, 0,2.*M_PI);
        pt.x         = a*cos(phi);
        pt.y         = a*sin(phi);
        pt.z         = a*reb_random_normal(r, 0.001);
        double mu     = star.m + disc_mass * (pow(a,-3./2.)-pow(boxsize/10.,-3./2.))/(pow(boxsize/2./1.2,-3./2.)-pow(boxsize/10.,-3./2.));
        double vkep     = sqrt(r->G*mu/a);
        pt.vx         =  vkep * sin(phi);
        pt.vy         = -vkep * cos(phi);
        pt.vz         = 0;
        pt.m         = disc_mass/(double)N;
        reb_add(r, pt);
    }

    r->heartbeat = heartbeat;
    reb_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* const r){
    if (reb_output_check(r,10.0*r->dt)){
        reb_output_timing(r,0);
    }
}
