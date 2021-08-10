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
 * 初始化随机种子
 */
void reb_tools_init_srand(struct reb_simulation* r);

// Random sampling
double reb_random_uniform(struct reb_simulation* r, double min, double max);
double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope);
double reb_random_normal(struct reb_simulation* r, double variance);
double reb_random_rayleigh(struct reb_simulation* r, double sigma);

/**
 * @brief internal function to handle outputs for the Fast Simulation Restarter.
 */
void reb_fsr_heartbeat(struct reb_simulation* const r);

// Diangnostic functions
// 返回系统总能量、总角动量
double reb_tools_energy(const struct reb_simulation* const r);
struct reb_vec3d reb_tools_angular_momentum(const struct reb_simulation* const r);

// Miscellaneous functions
// hash序列生成、mod、plummer生成、
uint32_t reb_hash(const char* str);
double reb_tools_mod2pi(double f);
void reb_tools_init_plummer(struct reb_simulation* r, int _N, double M, double R); // This function sets up a Plummer sphere, N=number of particles, M=total mass, R=characteristic radius
void reb_run_heartbeat(struct reb_simulation* const r);  // used internally

// Serialization functions.
// 主控离散数据转换为向量数据
void reb_serialize_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]); // NULL pointers will not be set.
void reb_set_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]); // Null pointers will be ignored.


#endif 	// TOOLS_H
