/**
 * @file 	tools.h
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
#ifndef TOOLS_H
#define TOOLS_H

#include <stdint.h>

struct reb_simulation;
struct reb_particles;

/**
 * @brief Init random number generator based on time and process id.
 */
void reb_tools_init_srand(struct reb_simulation* r);

/**
 * @brief internal function to handle outputs for the Fast Simulation Restarter.
 */
void reb_fsr_heartbeat(struct reb_simulation* const r);

// Random sampling
double reb_random_uniform(struct reb_simulation* r, double min, double max);
double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope);
double reb_random_normal(struct reb_simulation* r, double variance);
double reb_random_rayleigh(struct reb_simulation* r, double sigma);

#endif 	// TOOLS_H
