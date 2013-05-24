/*
    tsp.h -- header file for tsp, a graphical TSP example using fgen.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2013, Harm Hanemaaijer

    This file is part of fgen.

    fgen is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    fgen is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with fgen.  If not, see <http://www.gnu.org/licenses/>.

*/

#define MAX_CITIES 1000
#define AREA_SIZE 100.0

#define MUTATION_TYPE_SWAP 0
#define MUTATION_TYPE_SWAP_NEIGHBOURS 1
#define MUTATION_TYPE_INSERT 2
#define MUTATION_TYPE_INVERT 3
#define MUTATION_TYPE_NOOP 4
#define CROSSOVER_TYPE_ORDER1 0
#define CROSSOVER_TYPE_PBX 1
#define CROSSOVER_TYPE_NOOP 2

typedef struct {
	double x;
	double y;
} City;

extern City city[MAX_CITIES];
extern int nu_cities;
extern int best_route[MAX_CITIES];
extern int best_route_valid;
extern int ga_busy;
extern int ga_generation;
extern int population_size;
extern int mutation_type;
extern int crossover_type;
extern int reporting_frequency;
extern float mutation_probability;
extern float crossover_probability;
extern int steady_state;

void initialize_cities();

void gui_initialize(int *argc, char ***argv);
void gui_create_window_layout();
void gui_run();
void gui_handle_events();
void gui_draw_window();
void gui_draw_and_show_window();
double best_distance();

void start_ga();
void stop_ga();
void benchmark_current_settings();
void benchmark_operators();
void benchmark_rates();


