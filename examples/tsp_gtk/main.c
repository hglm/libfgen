/*
    main.c -- main file of tsp, a graphical TSP example using fgen.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2012, Harm Hanemaaijer

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "tsp.h"

#ifdef _WIN32
#define random() rand()
#define srandom(x) srand(x)
#endif

City city[MAX_CITIES];
int nu_cities;
int best_route[MAX_CITIES];
int best_route_valid = 0;
int ga_busy = 0;
int ga_generation;
int population_size = 256;
int mutation_type = MUTATION_TYPE_INVERT;
int crossover_type = CROSSOVER_TYPE_ORDER1;
int reporting_frequency = 1;
float mutation_probability = 0.2;
float crossover_probability = 0.2;
int steady_state = 0;

static double random_d(double range) {
    return (double)random() * range / ((double)RAND_MAX + 1);
}

void initialize_cities() {
    int i;
    for (i = 0; i < nu_cities; i++) {
        city[i].x = random_d(AREA_SIZE);
        city[i].y = random_d(AREA_SIZE);
    }
    best_route_valid = 0;
}

double best_distance() {
    double dist = 0;
    for (int i = 0; i < nu_cities - 1; i++) {
        dist += sqrt((city[best_route[i + 1]].x - city[best_route[i]].x) * (city[best_route[i + 1]].x -
            city[best_route[i]].x) + (city[best_route[i + 1]].y - city[best_route[i]].y) * (city[best_route[i + 1]].y
            - city[best_route[i]].y));
    }
    return dist;
}

int main(int argc, char **argv) {
    gui_initialize(&argc, &argv);
    nu_cities = 100;
    srandom(0);
    initialize_cities();
    gui_create_window_layout();
    gui_draw_window();
    gui_run();
}

