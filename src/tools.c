/**
 * @file 	tools.c
 * @brief 	Tools for creating distributions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdarg.h>
#include "particle.h"
#include "rebound.h"
#include "tools.h"


void reb_tools_init_srand(struct reb_simulation* r){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	r->rand_seed = tim.tv_usec + getpid();
}

double reb_random_uniform(struct reb_simulation* r, double min, double max){
	return ((double)rand_r(&(r->rand_seed)))/((double)(RAND_MAX))*(max-min)+min;
}


double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope){
	double y = reb_random_uniform(r, 0., 1.);
	if(slope == -1) return exp(y*log(max/min) + log(min));
    else return pow( (pow(max,slope+1.)-pow(min,slope+1.))*y+pow(min,slope+1.), 1./(slope+1.));
}

double reb_random_normal(struct reb_simulation* r, double variance){
	double v1,v2,rsq=1.;
	while(rsq>=1. || rsq<1.0e-12){
		v1=2.*((double)rand_r(&(r->rand_seed)))/((double)(RAND_MAX))-1.0;
		v2=2.*((double)rand_r(&(r->rand_seed)))/((double)(RAND_MAX))-1.0;
		rsq=v1*v1+v2*v2;
	}
	// Note: This gives another random variable for free, but we'll throw it away for simplicity and for thread-safety.
	return 	v1*sqrt(-2.*log(rsq)/rsq*variance);
}

double reb_random_rayleigh(struct reb_simulation* r, double sigma){
	double y = reb_random_uniform(r, 0.,1.);
	return sigma*sqrt(-2*log(y));
}

/// Other helper routines
double reb_tools_energy(const struct reb_simulation* const r){
    const int N = r->N;
    const int N_var = r->N_var;
    const int _N_active = (r->N_active==-1)?(N-N_var):r->N_active;
    const struct reb_particle* restrict const particles = r->particles;
    double e_kin = 0.;
    double e_pot = 0.;
    int N_interact = (r->testparticle_type==0)?_N_active:(N-N_var);
    for (int i=0;i<N_interact;i++){
        struct reb_particle pi = particles[i];
        e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
    }
    for (int i=0;i<_N_active;i++){
        struct reb_particle pi = particles[i];
        for (int j=i+1;j<N_interact;j++){
            struct reb_particle pj = particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double dz = pi.z - pj.z;
            e_pot -= r->G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    
    return e_kin + e_pot + r->energy_offset;
}

struct reb_vec3d reb_tools_angular_momentum(const struct reb_simulation* const r){
	const int N = r->N;
	const struct reb_particle* restrict const particles = r->particles;
	const int N_var = r->N_var;
    struct reb_vec3d L = {0};
    for (int i=0;i<N-N_var;i++){
		struct reb_particle pi = particles[i];
        L.x += pi.m*(pi.y*pi.vz - pi.z*pi.vy);
        L.y += pi.m*(pi.z*pi.vx - pi.x*pi.vz);
        L.z += pi.m*(pi.x*pi.vy - pi.y*pi.vx);
	}
	return L;
}

void reb_serialize_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
    const int N_real = r->N - r->N_var;
    struct reb_particle* restrict const particles = r->particles;
    for (int i=0;i<N_real;i++){
        if (hash){
            hash[i] = particles[i].hash;
        }
        if (m){
            m[i] = particles[i].m;
        }
        if (radius){
            radius[i] = particles[i].r;
        }
        if (xyz){
            xyz[i][0] = particles[i].x;
            xyz[i][1] = particles[i].y;
            xyz[i][2] = particles[i].z;
        }
        if (vxvyvz){
            vxvyvz[i][0] = particles[i].vx;
            vxvyvz[i][1] = particles[i].vy;
            vxvyvz[i][2] = particles[i].vz;
        }
        if (xyzvxvyvz){
            xyzvxvyvz[i][0] = particles[i].x;
            xyzvxvyvz[i][1] = particles[i].y;
            xyzvxvyvz[i][2] = particles[i].z;
            xyzvxvyvz[i][3] = particles[i].vx;
            xyzvxvyvz[i][4] = particles[i].vy;
            xyzvxvyvz[i][5] = particles[i].vz;
        }
    }
}

void reb_set_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
    const int N_real = r->N - r->N_var;
    struct reb_particle* restrict const particles = r->particles;
    for (int i=0;i<N_real;i++){
        if (hash){
           particles[i].hash = hash[i];
        }
        if (m){
            particles[i].m = m[i];
        }
        if (radius){
            particles[i].r = radius[i] ;
        }
        if (xyz){
            particles[i].x = xyz[i][0];
            particles[i].y = xyz[i][1];
            particles[i].z = xyz[i][2];
        }
        if (vxvyvz){
            particles[i].vx = vxvyvz[i][0];
            particles[i].vy = vxvyvz[i][1];
            particles[i].vz = vxvyvz[i][2];
        }
        if (xyzvxvyvz){
            particles[i].x = xyzvxvyvz[i][0];
            particles[i].y = xyzvxvyvz[i][1];
            particles[i].z = xyzvxvyvz[i][2];
            particles[i].vx = xyzvxvyvz[i][3];
            particles[i].vy = xyzvxvyvz[i][4];
            particles[i].vz = xyzvxvyvz[i][5];
        }
    }
}
	
void reb_tools_init_plummer(struct reb_simulation* r, int _N, double M, double R) {
	// Algorithm from:	
	// http://adsabs.harvard.edu/abs/1974A%26A....37..183A
	
	double E = 3./64.*M_PI*M*M/R;
	for (int i=0;i<_N;i++){
		struct reb_particle star = {0};
		double _r = pow(pow(reb_random_uniform(r, 0,1),-2./3.)-1.,-1./2.);
		double x2 = reb_random_uniform(r, 0,1);
		double x3 = reb_random_uniform(r, 0,2.*M_PI);
		star.z = (1.-2.*x2)*_r;
		star.x = sqrt(_r*_r-star.z*star.z)*cos(x3);
		star.y = sqrt(_r*_r-star.z*star.z)*sin(x3);
		double x5,g,q;
		do{
			x5 = reb_random_uniform(r, 0.,1.);
			q = reb_random_uniform(r, 0.,1.);
			g = q*q*pow(1.-q*q,7./2.);
		}while(0.1*x5>g);
		double ve = pow(2.,1./2.)*pow(1.+_r*_r,-1./4.);
		double v = q*ve;
		double x6 = reb_random_uniform(r, 0.,1.);
		double x7 = reb_random_uniform(r, 0.,2.*M_PI);
		star.vz = (1.-2.*x6)*v;
		star.vx = sqrt(v*v-star.vz*star.vz)*cos(x7);
		star.vy = sqrt(v*v-star.vz*star.vz)*sin(x7);
		
		star.x *= 3.*M_PI/64.*M*M/E;
		star.y *= 3.*M_PI/64.*M*M/E;
		star.z *= 3.*M_PI/64.*M*M/E;
		
		star.vx *= sqrt(E*64./3./M_PI/M);
		star.vy *= sqrt(E*64./3./M_PI/M);
		star.vz *= sqrt(E*64./3./M_PI/M);

		star.m = M/(double)_N;

		reb_add(r, star);
	}
}

double reb_tools_mod2pi(double f){
    const double pi2 = 2.*M_PI;
    return fmod(pi2 + fmod(f, pi2), pi2);
}

#define TINY 1.E-308 		///< Close to smallest representable floating point number, used for orbit calculation

struct reb_orbit reb_orbit_nan(void){
    struct reb_orbit o;
    o.d = nan("");
    o.v = nan("");
    o.h = nan("");
    o.P = nan("");
    o.n = nan("");
    o.a = nan("");
    o.e = nan("");
    o.inc = nan("");
    o.Omega = nan("");
    o.omega = nan("");
    o.pomega = nan("");
    o.f = nan("");
    o.M = nan("");
    o.l = nan("");
    o.theta = nan("");
    o.T = nan("");
    o.rhill = nan("");

    return o;
}

#define MIN_REL_ERROR 1.0e-12	///< Close to smallest relative floating point number, used for orbit calculation
#define MIN_INC 1.e-8		///< Below this inclination, the broken angles pomega and theta equal the corresponding 
							///< unbroken angles to within machine precision, so a practical boundary for planar orbits
							//
#define MIN_ECC 1.e-8       ///< Below this eccentricity, corrections at order e^2 are below machine precision, so we use
                            ///< stable expressions accurate to O(e) for the mean longitude below for near-circular orbits.
// returns acos(num/denom), using disambiguator to tell which quadrant to return.  
// will return 0 or pi appropriately if num is larger than denom by machine precision
// and will return 0 if denom is exactly 0.

static double acos2(double num, double denom, double disambiguator){
	double val;
	double cosine = num/denom;
	if(cosine > -1. && cosine < 1.){
		val = acos(cosine);
		if(disambiguator < 0.){
			val = - val;
		}
	}
	else{
		val = (cosine <= -1.) ? M_PI : 0.;
	}
	return val;
}

/***********************************
 * Variational Equations and Megno */

#define ROT32(x, y) ((x << y) | (x >> (32 - y))) // avoid effort
static uint32_t reb_murmur3_32(const char *key, uint32_t len, uint32_t seed) {
    // Source: Wikipedia
    static const uint32_t c1 = 0xcc9e2d51;
    static const uint32_t c2 = 0x1b873593;
    static const uint32_t r1 = 15;
    static const uint32_t r2 = 13;
    static const uint32_t m = 5;
    static const uint32_t n = 0xe6546b64;

    uint32_t hash = seed;

    const int nblocks = len / 4;
    const uint32_t *blocks = (const uint32_t *) key;
    int i;
    uint32_t k;
    for (i = 0; i < nblocks; i++) {
        k = blocks[i];
        k *= c1;
        k = ROT32(k, r1);
        k *= c2;

        hash ^= k;
        hash = ROT32(hash, r2) * m + n;
    }

    const uint8_t *tail = (const uint8_t *) (key + nblocks * 4);
    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
    case 2:
        k1 ^= tail[1] << 8;
    case 1:
        k1 ^= tail[0];

        k1 *= c1;
        k1 = ROT32(k1, r1);
        k1 *= c2;
        hash ^= k1;
    }

    hash ^= len;
    hash ^= (hash >> 16);
    hash *= 0x85ebca6b;
    hash ^= (hash >> 13);
    hash *= 0xc2b2ae35;
    hash ^= (hash >> 16);

    return hash;
}

uint32_t reb_hash(const char* str){
    const int reb_seed = 1983;
    return reb_murmur3_32(str,(uint32_t)strlen(str),reb_seed);
}

