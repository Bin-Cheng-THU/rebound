/**
 * @file 	collision.h
 * @brief 	Collision search. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @section LICENSE
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
#ifndef _COLLISIONS_H
#define _COLLISIONS_H
/**
 * @brief Search for collisions and resolve them.
 * @param r REBOUND Simulation to consider
 */
void reb_collision_search(struct reb_simulation* const r);

/*
 * @brief halt the simulation.
 * @param r REBOUND Simulation to consider
 * @param c REBOUND collision event
int reb_collision_resolve_halt(struct reb_simulation* const r, struct reb_collision c);
 * @brief hard-sphere model for impact events.
 * @param r REBOUND Simulation to consider
 * @param c REBOUND collision event
int reb_collision_resolve_hardsphere(struct reb_simulation* const r, struct reb_collision c);
 * @brief Merge particles when impact.
 * @param r REBOUND Simulation to consider
 * @param c REBOUND collision event
int reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c);
*/

#endif // _COLLISIONS_H
