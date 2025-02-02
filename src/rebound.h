/**
 * @file    rebound.h
 * @brief   REBOUND API definition.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu
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

#ifndef _MAIN_H
#define _MAIN_H

#define REBOUND_RESTRICT restrict

#include <inttypes.h>
#include <stdint.h>
#include <sys/time.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include "mpi.h"
#endif // MPI
#ifndef GITHASH
#define GITHASH notavailable0000000000000000000000000001 
#endif // GITHASH

// 全局变量
extern const char* reb_build_str;   ///< Date and time build string.
extern const char* reb_version_str; ///< Version string.
extern const char* reb_githash_str; ///< Current git hash.
extern const char* reb_logo[26];    ///< Logo of rebound. 
extern volatile sig_atomic_t reb_sigint;  ///< Graceful global interrupt handler 

// Forward declarations
struct reb_simulation;
struct reb_display_data;
struct reb_treecell;

//rebound颗粒
struct reb_particle {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
    double m;                   // mass
    double r;                   // physical radius
    double lastcollision;       // Last time the particle had a physical collision.
    struct reb_treecell* c;     // Pointer to the cell the particle is currently in.
    uint32_t hash;              // Hash, can be used to identify particle.
    void* ap;                   // This pointer allows REBOUNDx to add additional properties to the particle.
    struct reb_simulation* sim; // Pointer to the parent simulation.
};

// Generic 3d vector
struct reb_vec3d {
    double x;
    double y;
    double z;
};

struct reb_ghostbox{
    double shiftx;
    double shifty;
    double shiftz;
    double shiftvx;
    double shiftvy;
    double shiftvz;
};

struct reb_collision{
    int p1;
    int p2;
    struct reb_ghostbox gb;
    int ri;
};

// Possible return values of of rebound_integrate
enum REB_STATUS {
    REB_RUNNING_PAUSED = -3,    // Simulation is paused by visualization.
    REB_RUNNING_LAST_STEP = -2, // Current timestep is the last one. Needed to ensure that t=tmax exactly.
    REB_RUNNING = -1,           // Simulation is current running, no error occurred.
    REB_EXIT_SUCCESS = 0,       // Integration finished successfully.
    REB_EXIT_ERROR = 1,         // A generic error occurred and the integration was not successful.
    REB_EXIT_NOPARTICLES = 2,   // The integration ends early because no particles are left in the simulation.
    REB_EXIT_ENCOUNTER = 3,     // The integration ends early because two particles had a close encounter (see exit_min_distance)
    REB_EXIT_ESCAPE = 4,        // The integration ends early because a particle escaped (see exit_max_distance)  
    REB_EXIT_USER = 5,          // User caused exit, simulation did not finish successfully.
    REB_EXIT_SIGINT = 6,        // SIGINT received. Simulation stopped.
    REB_EXIT_COLLISION = 7,     // The integration ends early because two particles collided. 
};

// IDs for content of a binary field. Used to read and write binary files.
enum REB_BINARY_FIELD_TYPE {
    REB_BINARY_FIELD_TYPE_T = 0,
    REB_BINARY_FIELD_TYPE_G = 1,
    REB_BINARY_FIELD_TYPE_SOFTENING = 2,
    REB_BINARY_FIELD_TYPE_DT = 3,
    REB_BINARY_FIELD_TYPE_N = 4,
    REB_BINARY_FIELD_TYPE_NVAR = 5,
    REB_BINARY_FIELD_TYPE_VARCONFIGN = 6,
    REB_BINARY_FIELD_TYPE_NACTIVE = 7,
    REB_BINARY_FIELD_TYPE_TESTPARTICLETYPE = 8,
    REB_BINARY_FIELD_TYPE_HASHCTR = 9, 
    REB_BINARY_FIELD_TYPE_OPENINGANGLE2 = 10,
    REB_BINARY_FIELD_TYPE_STATUS = 11,
    REB_BINARY_FIELD_TYPE_EXACTFINISHTIME = 12,
    REB_BINARY_FIELD_TYPE_FORCEISVELOCITYDEP = 13,
    REB_BINARY_FIELD_TYPE_GRAVITYIGNORETERMS = 14,
    REB_BINARY_FIELD_TYPE_OUTPUTTIMINGLAST = 15,
    REB_BINARY_FIELD_TYPE_SAVEMESSAGES = 16,
    REB_BINARY_FIELD_TYPE_EXITMAXDISTANCE = 17,
    REB_BINARY_FIELD_TYPE_EXITMINDISTANCE = 18,
    REB_BINARY_FIELD_TYPE_USLEEP = 19,
    REB_BINARY_FIELD_TYPE_TRACKENERGYOFFSET = 20,
    REB_BINARY_FIELD_TYPE_ENERGYOFFSET = 21,
    REB_BINARY_FIELD_TYPE_BOXSIZE = 22, 
    REB_BINARY_FIELD_TYPE_BOXSIZEMAX = 23, 
    REB_BINARY_FIELD_TYPE_ROOTSIZE = 24,
    REB_BINARY_FIELD_TYPE_ROOTN = 25,
    REB_BINARY_FIELD_TYPE_ROOTNX = 26, 
    REB_BINARY_FIELD_TYPE_ROOTNY = 27,
    REB_BINARY_FIELD_TYPE_ROOTNZ = 28,
    REB_BINARY_FIELD_TYPE_NGHOSTX = 29,
    REB_BINARY_FIELD_TYPE_NGHOSTY = 30,
    REB_BINARY_FIELD_TYPE_NGHOSTZ = 31,
    REB_BINARY_FIELD_TYPE_COLLISIONRESOLVEKEEPSORTED = 32,
    REB_BINARY_FIELD_TYPE_MINIMUMCOLLISIONVELOCITY = 33,
    REB_BINARY_FIELD_TYPE_COLLISIONSPLOG = 34, 
    REB_BINARY_FIELD_TYPE_MAXRADIUS = 35, 
    REB_BINARY_FIELD_TYPE_COLLISIONSNLOG = 36, 
    REB_BINARY_FIELD_TYPE_CALCULATEMEGNO = 37, 
    REB_BINARY_FIELD_TYPE_MEGNOYS = 38, 
    REB_BINARY_FIELD_TYPE_MEGNOYSS = 39, 
    REB_BINARY_FIELD_TYPE_MEGNOCOVYT = 40,
    REB_BINARY_FIELD_TYPE_MEGNOVART = 41, 
    REB_BINARY_FIELD_TYPE_MEGNOMEANT = 42, 
    REB_BINARY_FIELD_TYPE_MEGNOMEANY = 43, 
    REB_BINARY_FIELD_TYPE_MEGNON = 44,
    REB_BINARY_FIELD_TYPE_SASIZEFIRST = 45,
    REB_BINARY_FIELD_TYPE_SASIZESNAPSHOT = 46,
    REB_BINARY_FIELD_TYPE_SAAUTOINTERVAL = 47,
    REB_BINARY_FIELD_TYPE_SAAUTOWALLTIME = 102,
    REB_BINARY_FIELD_TYPE_SANEXT = 48,
    REB_BINARY_FIELD_TYPE_COLLISION = 50,
    REB_BINARY_FIELD_TYPE_INTEGRATOR = 51,
    REB_BINARY_FIELD_TYPE_BOUNDARY = 52,
    REB_BINARY_FIELD_TYPE_GRAVITY = 53,
    REB_BINARY_FIELD_TYPE_PARTICLES = 85,
    REB_BINARY_FIELD_TYPE_VARCONFIG = 86,
    REB_BINARY_FIELD_TYPE_FUNCTIONPOINTERS = 87,
    REB_BINARY_FIELD_TYPE_VISUALIZATION = 107,
    REB_BINARY_FIELD_TYPE_SAVERSION = 125,
    REB_BINARY_FIELD_TYPE_WALLTIME = 126,
    REB_BINARY_FIELD_TYPE_PYTHON_UNIT_L = 130,
    REB_BINARY_FIELD_TYPE_PYTHON_UNIT_M = 131,
    REB_BINARY_FIELD_TYPE_PYTHON_UNIT_T = 132,
    REB_BINARY_FIELD_TYPE_SAAUTOSTEP = 135,
    REB_BINARY_FIELD_TYPE_SANEXTSTEP = 136,
    REB_BINARY_FIELD_TYPE_STEPSDONE = 137,
    REB_BINARY_FIELD_TYPE_WHFAST_CORRECTOR2 = 143,
    REB_BINARY_FIELD_TYPE_WHFAST_KERNEL = 144,
    REB_BINARY_FIELD_TYPE_DTLASTDONE = 145,
    REB_BINARY_FIELD_TYPE_RAND_SEED = 154,
    REB_BINARY_FIELD_TYPE_TESTPARTICLEHIDEWARNINGS = 155,

    REB_BINARY_FIELD_TYPE_HEADER = 1329743186,  // Corresponds to REBO (first characters of header text)
    REB_BINARY_FIELD_TYPE_SABLOB = 9998,        // SA Blob
    REB_BINARY_FIELD_TYPE_END = 9999,
};

// This structure is used to save and load binary files.
struct reb_binary_field {
    uint32_t type;  // enum of REB_BINARY_FIELD_TYPE
    uint64_t size;  // Size in bytes of field (only counting what follows, not the binary field, itself).
};

// Holds a particle's hash and the particle's index in the particles array. Used for particle_lookup_table.
struct reb_hash_pointer_pair{
    uint32_t hash;
    int index;
};

// Main REBOUND Simulation structure
struct reb_simulation {
    double  t;
    double  G;
    double  softening;
    double  dt;
    double  dt_last_done;
    unsigned long long steps_done;
    int     N;
    int     N_active;
    int     testparticle_type;
    int     testparticle_hidewarnings;
    struct reb_hash_pointer_pair* particle_lookup_table; // Array of pairs that map particles' hashes to their index in the particles array.
    int     hash_ctr;               // Counter for number of assigned hashes to assign unique values.
    int     N_lookup;               // Number of entries in the particle lookup table.
    int     allocatedN_lookup;      // Number of lookup table entries allocated.
    int     allocatedN;             // Current maximum space allocated in the particles array on this node. 
    struct reb_particle* particles;
    struct reb_vec3d* gravity_cs;   // Containing the information for compensated gravity summation 
    int     gravity_cs_allocatedN;
    struct reb_treecell** tree_root;// Pointer to the roots of the trees. 
    int     tree_needs_update;      // Flag to force a tree update (after boundary check)
    double opening_angle2;
    enum REB_STATUS status;
    int     exact_finish_time;

    unsigned int force_is_velocity_dependent;
    unsigned int gravity_ignore_terms;
    double output_timing_last;      // Time when reb_output_timing() was called the last time. 
    unsigned long display_clock;    // Display clock, internal variable for timing refreshs.
    int save_messages;              // Set to 1 to ignore messages (used in python interface).
    char** messages;                // Array of strings containing last messages (only used if save_messages==1). 
    double exit_max_distance;
    double exit_min_distance;
    double usleep;
    struct reb_display_data* display_data; // Datastructure stores visualization related data. Does not have to be modified by the user. 
    int track_energy_offset;        // 是否记录能量变化
    double energy_offset;           // 能量变化
    double walltime;
    uint32_t python_unit_l;         // Only used for when working with units in python.
    uint32_t python_unit_m;         // Only used for when working with units in python.
    uint32_t python_unit_t;         // Only used for when working with units in python.
    
    // Ghost boxes 
    struct  reb_vec3d boxsize;      // Size of the entire box, root_x*boxsize. 
    double  boxsize_max;            // Maximum size of the entire box in any direction. Set in box_init().
    double  root_size;              // Size of a root box. 
    int     root_n;                 // Total number of root boxes in all directions, root_nx*root_ny*root_nz. Default: 1. Set in box_init().
    int     root_nx;                // Number of ghost boxes in x direction. Do not change manually.
    int     root_ny;
    int     root_nz;
    int     nghostx;
    int     nghosty;
    int     nghostz;

#ifdef MPI
    int    mpi_id;                              // Unique id of this node (starting at 0). Used for MPI only.
    int    mpi_num;                             // Number of MPI nodes. Used for MPI only.
    MPI_Datatype mpi_particle;                  // MPI datatype corresponding to the C struct reb_particle. 
    struct reb_particle** particles_send;       // Send buffer for particles. There is one buffer per node. 
    int*   particles_send_N;                    // Current length of particle send buffer. 
    int*   particles_send_Nmax;                 // Maximal length of particle send beffer before realloc() is needed. 
    struct reb_particle** particles_recv;       // Receive buffer for particles. There is one buffer per node. 
    int*   particles_recv_N;                    // Current length of particle receive buffer. 
    int*   particles_recv_Nmax;                 // Maximal length of particle receive beffer before realloc() is needed. */

    MPI_Datatype mpi_cell;                      // MPI datatype corresponding to the C struct reb_treecell. 
    struct reb_treecell** tree_essential_send;  // Send buffer for cells. There is one buffer per node. 
    int*   tree_essential_send_N;               // Current length of cell send buffer. 
    int*   tree_essential_send_Nmax;            // Maximal length of cell send beffer before realloc() is needed. 
    struct reb_treecell** tree_essential_recv;  // Receive buffer for cells. There is one buffer per node. 
    int*   tree_essential_recv_N;               // Current length of cell receive buffer. 
    int*   tree_essential_recv_Nmax;            // Maximal length of cell receive beffer before realloc() is needed. 
#endif // MPI

    int collision_resolve_keep_sorted;
    struct reb_collision* collisions;       ///< Array of all collisions. 
    int collisions_allocatedN;
    double minimum_collision_velocity;
    double collisions_plog;
    double max_radius[2];               // Two largest particle radii, set automatically, needed for collision search.
    long collisions_Nlog;

    unsigned int rand_seed; // seed for random number generator
    
     // SimulationArchive 
    int    simulationarchive_version;               // Version of the SA binary format (1=original/, 2=incremental)
    long   simulationarchive_size_first;            // (Deprecated SAV1) Size of the initial binary file in a SA
    long   simulationarchive_size_snapshot;         // (Deprecated SAV1) Size of a snapshot in a SA (other than 1st), in bytes
    double simulationarchive_auto_interval;         // Current sampling cadence, in code units
    double simulationarchive_auto_walltime;         // Current sampling cadence, in wall time
    unsigned long long simulationarchive_auto_step; // Current sampling cadence, in time steps
    double simulationarchive_next;                  // Next output time (simulation tim or wall time, depending on wether auto_interval or auto_walltime is set)
    unsigned long long simulationarchive_next_step; // Next output step (only used if auto_steps is set)
    char*  simulationarchive_filename;              // Name of output file

    // Modules
    enum {
        REB_VISUALIZATION_NONE = 0,     // No visualization (default if OPENGL compiler flag is turned off)
        REB_VISUALIZATION_OPENGL = 1,   // OpenGL visualization (default if OPENGL compiler flag is turned on)
        REB_VISUALIZATION_WEBGL = 2,    // WebGL visualization, only usable from Jupyter notebook widget
        } visualization;
    enum {
        REB_COLLISION_NONE = 0,     // Do not search for collisions (default)
        REB_COLLISION_DIRECT = 1,   // Direct collision search O(N^2)
        REB_COLLISION_TREE = 2,     // Tree based collision search O(N log(N))
        REB_COLLISION_LINE = 4,     // Direct collision search O(N^2), looks for collisions by assuming a linear path over the last timestep
        REB_COLLISION_LINETREE = 5, // Tree-based collision search O(N log(N)), looks for collisions by assuming a linear path over the last timestep
        } collision;
    enum {
        REB_INTEGRATOR_LEAPFROG = 4, // LEAPFROG integrator, simple, 2nd order, symplectic
        REB_INTEGRATOR_NONE = 7,     // Do not integrate anything
        } integrator;
    enum {
        REB_BOUNDARY_NONE = 0,      // Do not check for anything (default)
        REB_BOUNDARY_OPEN = 1,      // Open boundary conditions. Removes particles if they leave the box 
        REB_BOUNDARY_PERIODIC = 2,  // Periodic boundary conditions
        REB_BOUNDARY_SHEAR = 3,     // Shear periodic boundary conditions, needs OMEGA variable
        } boundary;
    enum {
        REB_GRAVITY_NONE = 0,       // Do not calculate graviational forces
        REB_GRAVITY_BASIC = 1,      // Basic O(N^2) direct summation algorithm, choose this for shearing sheet and periodic boundary conditions
        REB_GRAVITY_COMPENSATED = 2,// Direct summation algorithm O(N^2) but with compensated summation, slightly slower than BASIC but more accurate
        REB_GRAVITY_TREE = 3,       // Use the tree to calculate gravity, O(N log(N)), set opening_angle2 to adjust accuracy.
        } gravity;

    // Integrators

     // Callback functions
    void (*additional_forces) (struct reb_simulation* const r);
    void (*pre_timestep_modifications) (struct reb_simulation* const r);    // used by REBOUNDx
    void (*post_timestep_modifications) (struct reb_simulation* const r);   // used by REBOUNDx
    void (*heartbeat) (struct reb_simulation* r);
    void (*display_heartbeat) (struct reb_simulation* r);
    double (*coefficient_of_restitution) (const struct reb_simulation* const r, double v); 
    int (*collision_resolve) (struct reb_simulation* const r, struct reb_collision);
    void (*extras_cleanup) (struct reb_simulation* r);
    void* extras; // Pointer to connect additional (optional) libraries, e.g., reboundx
};


// Simulation life cycle
struct reb_simulation* reb_create_simulation(void);     // allocates memory, then calls reb_init_simulation
void reb_init_simulation(struct reb_simulation* r);    
void reb_free_simulation(struct reb_simulation* const r);
struct reb_simulation* reb_copy_simulation(struct reb_simulation* r);
void reb_free_pointers(struct reb_simulation* const r);
void reb_reset_temporary_pointers(struct reb_simulation* const r);
int reb_reset_function_pointers(struct reb_simulation* const r); // Returns 1 if one ore more function pointers were not NULL before.
// Configure the boundary/root box
void reb_configure_box(struct reb_simulation* const r, const double boxsize, const int root_nx, const int root_ny, const int root_nz);

// Messages and control functions
void reb_exit(const char* const msg); // Print out an error message, then exit in a semi-nice way.
void reb_warning(struct reb_simulation* const r, const char* const msg);   // Print or store a warning message, then continue.
void reb_error(struct reb_simulation* const r, const char* const msg);     // Print or store an error message, then continue.
int reb_get_next_message(struct reb_simulation* const r, char* const buf); // Get the next stored warning message. Used only if save_messages==1. Return value is 0 if no messages are present, 1 otherwise.

// Timestepping
void reb_step(struct reb_simulation* const r);
void reb_steps(struct reb_simulation* const r, unsigned int N_steps);
enum REB_STATUS reb_integrate(struct reb_simulation* const r, double tmax);
void reb_integrator_synchronize(struct reb_simulation* r);
void reb_integrator_reset(struct reb_simulation* r);
void reb_update_acceleration(struct reb_simulation* r);

// Compare simulations
// If r1 and r2 are exactly equal to each other then 0 is returned, otherwise 1. Walltime is ignored.
// If output_option=1, then output is printed on the screen. If 2, only return value os given. 

// Collision resolve functions
int reb_collision_resolve_halt(struct reb_simulation* const r, struct reb_collision c);
int reb_collision_resolve_halt(struct reb_simulation* const r, struct reb_collision c);
int reb_collision_resolve_hardsphere(struct reb_simulation* const r, struct reb_collision c);
int reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c);

// Output functions
int reb_output_check(struct reb_simulation* r, double interval);
void reb_output_timing(struct reb_simulation* r, const double tmax);
void reb_output_orbits(struct reb_simulation* r, char* filename);
void reb_output_binary(struct reb_simulation* r, const char* filename);
void reb_output_ascii(struct reb_simulation* r, char* filename);
void reb_output_binary_positions(struct reb_simulation* r, const char* filename);
void reb_output_velocity_dispersion(struct reb_simulation* r, char* filename);

// Input functions
struct reb_simulation* reb_create_simulation_from_binary(char* filename);

// Possible errors that might occur during binary file reading.
enum reb_input_binary_messages {
    REB_INPUT_BINARY_WARNING_NONE = 0,
    REB_INPUT_BINARY_ERROR_NOFILE = 1,
    REB_INPUT_BINARY_WARNING_VERSION = 2,
    REB_INPUT_BINARY_WARNING_POINTERS = 4,
    REB_INPUT_BINARY_WARNING_PARTICLES = 8,
    REB_INPUT_BINARY_ERROR_FILENOTOPEN = 16,
    REB_INPUT_BINARY_ERROR_OUTOFRANGE = 32,
    REB_INPUT_BINARY_ERROR_SEEK = 64,
    REB_INPUT_BINARY_WARNING_FIELD_UNKOWN = 128,
    REB_INPUT_BINARY_ERROR_INTEGRATOR = 256,
    REB_INPUT_BINARY_WARNING_CORRUPTFILE = 512,
};

// Functions to add and initialize particles
struct reb_particle reb_particle_nan(void); // Returns a reb_particle structure with fields/hash/ptrs initialized to nan/0/NULL. 
void reb_add(struct reb_simulation* const r, struct reb_particle pt);

// Functions to access and remove particles
void reb_remove_all(struct reb_simulation* const r);
int reb_remove(struct reb_simulation* const r, int index, int keepSorted);
int reb_remove_by_hash(struct reb_simulation* const r, uint32_t hash, int keepSorted);
struct reb_particle* reb_get_particle_by_hash(struct reb_simulation* const r, uint32_t hash);
int reb_get_particle_index(struct reb_particle* p); // Returns a particle's index in the simulation it's in. Needs to be in the simulation its sim pointer is pointing to. Otherwise -1 returned.
struct reb_particle reb_get_jacobi_com(struct reb_particle* p); // Returns the Jacobi center of mass for a given particle. Used by python. Particle needs to be in a simulation.

// Simulation Archive
struct reb_simulationarchive_blob {  // Used in the binary file to identify data blobs
    int32_t index;                   // Index of previous blob (binary file is 0, first blob is 1)
    int16_t offset_prev;             // Offset to beginning of previous blob (size of previous blob).
    int16_t offset_next;             // Offset to end of following blob (size of following blob).
};

struct reb_simulationarchive{
    FILE* inf;                   // File pointer (will be kept open)
    char* filename;              // Filename of open file
    int version;                 // SimulationArchive version
    long size_first;             // Size of first snapshot (only used for version 1)
    long size_snapshot;          // Size of snapshot (only used for version 1)
    double auto_interval;        // Interval setting used to create SA (if used)
    double auto_walltime;        // Walltime setting used to create SA (if used)
    unsigned long long auto_step;// Steps in-between SA snapshots (if used)
    long nblobs;                 // Total number of snapshots (including initial binary)
    uint32_t* offset;            // Index of offsets in file (length nblobs)
    double* t;                   // Index of simulation times in file (length nblobs)
};
struct reb_simulation* reb_create_simulation_from_simulationarchive(struct reb_simulationarchive* sa, long snapshot);
void reb_create_simulation_from_simulationarchive_with_messages(struct reb_simulation* r, struct reb_simulationarchive* sa, long snapshot, enum reb_input_binary_messages* warnings);
struct reb_simulationarchive* reb_open_simulationarchive(const char* filename);
void reb_close_simulationarchive(struct reb_simulationarchive* sa);
void reb_simulationarchive_snapshot(struct reb_simulation* r, const char* filename);
void reb_simulationarchive_automate_interval(struct reb_simulation* const r, const char* filename, double interval);
void reb_simulationarchive_automate_walltime(struct reb_simulation* const r, const char* filename, double walltime);
void reb_simulationarchive_automate_step(struct reb_simulation* const r, const char* filename, unsigned long long step);
void reb_free_simulationarchive_pointers(struct reb_simulationarchive* sa);


// Functions to between coordinate systems

#ifdef MPI
void reb_mpi_init(struct reb_simulation* const r);
void reb_mpi_finalize(struct reb_simulation* const r);
#endif // MPI

#ifdef OPENMP
// Wrapper method to set number of OpenMP threads from python.
void reb_omp_set_num_threads(int num_threads);
#endif // OPENMP

// The following stuctures are related to OpenGL/WebGL visualization. Nothing to be changed by the user.
struct reb_quaternion {
    double x, y, z, w;
};
struct reb_particle_opengl {
    float x,y,z;
    float vx,vy,vz;
    float r;
};
struct reb_orbit_opengl {
    float x,y,z;
    float a, e, f;
    float omega, Omega, inc;
};

struct reb_display_data {
    struct reb_simulation* r;
    struct reb_simulation* r_copy;
    struct reb_particle_opengl* particle_data;
    struct reb_particle* particles_copy;
    struct reb_particle* p_jh_copy;
    unsigned long allocated_N;
    unsigned long allocated_N_whfast;
    unsigned int opengl_enabled;
    double scale;
    double mouse_x;
    double mouse_y;
    double retina;
    pthread_mutex_t mutex;          // Mutex to guarantee non-flickering
    int spheres;                    // Switches between point sprite and real spheres.
    int pause;                      // Pauses visualization, but keep simulation running
    int wire;                       // Shows/hides orbit wires.
    int onscreentext;               // Shows/hides onscreen text.
    int onscreenhelp;               // Shows/hides onscreen help.
    int multisample;                // Turn off/on multisampling.
    int clear;                      // Toggles clearing the display on each draw.
    int ghostboxes;                 // Shows/hides ghost boxes.
    int reference;                  // reb_particle used as a reference for centering.
    unsigned int mouse_action;      
    unsigned int key_mods;      
    struct reb_quaternion view;
    unsigned int simplefont_tex;
    unsigned int simplefont_shader_program;
    unsigned int simplefont_shader_vao;
    unsigned int simplefont_shader_pos_location;
    unsigned int simplefont_shader_ypos_location;
    unsigned int simplefont_shader_scale_location;
    unsigned int simplefont_shader_aspect_location;
    unsigned int simplefont_shader_charval_buffer;
    unsigned int box_shader_program;
    unsigned int box_shader_box_vao;
    unsigned int box_shader_cross_vao;
    unsigned int box_shader_mvp_location;
    unsigned int box_shader_color_location;
    unsigned int point_shader_mvp_location;
    unsigned int point_shader_color_location;
    unsigned int point_shader_program;
    unsigned int point_shader_particle_vao;
    unsigned int sphere_shader_mvp_location;
    unsigned int sphere_shader_program;
    unsigned int sphere_shader_particle_vao;
    unsigned int sphere_shader_vertex_count;
    unsigned int orbit_shader_mvp_location;
    unsigned int orbit_shader_program;
    unsigned int orbit_shader_particle_vao;
    unsigned int orbit_shader_vertex_count;
};
#endif // _MAIN_H
