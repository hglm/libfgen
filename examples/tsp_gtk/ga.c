/*
    ga.c -- genetic algorithm interface of tsp, a graphical TSP example using fgen.

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
#include <string.h>
#include <math.h>
#include <limits.h>
#include <fgen.h>
#include "tsp.h"

static FgenPopulation *pop;

static void benchmark_generation_callback(FgenPopulation *pop, int generation);

static double ProblemCalculateFitness(const FgenPopulation *pop, const unsigned char *bitstring) {
    int *perm = (int *)bitstring;
    double dist;
    int i;
    dist = 0;
    for (i = 0; i < nu_cities - 1; i++)
        dist += sqrt((city[perm[i + 1]].x - city[perm[i]].x) * (city[perm[i + 1]].x - city[perm[i]].x) +
            (city[perm[i + 1]].y - city[perm[i]].y) * (city[perm[i + 1]].y - city[perm[i]].y));
    return AREA_SIZE * AREA_SIZE / dist;
}

static void ProblemGenerationCallback(FgenPopulation *pop, int generation) {
    if (generation % reporting_frequency != 0) {
        gui_handle_events();
        return;
    }
    FgenIndividual *best = fgen_best_individual_of_population(pop);
    int *perm = (int *)best->bitstring;
    memcpy(best_route, perm, sizeof(int) * nu_cities);
    best_route_valid = 1;
    ga_generation = generation;
    gui_draw_and_show_window();
    gui_handle_events();
}

// Custom mutation operator that swaps neighbouring cities in the permutation.

static void mutation_permutation_swap_neighbours(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	int index1, index2;
	int r = fgen_random_16(pop->rng);
	if (r >= pop->mutation_probability)
		return;
	index1 = fgen_random_n(pop->rng, pop->permutation_size - 1);
	index2 = index1 + 1;
	*(int *)&child[index1 * 4] = *(int *)&parent[index2 * 4];
	*(int *)&child[index2 * 4] = *(int *)&parent[index1 * 4];
}

FgenPopulation *setup_ga(int benchmark_mode) {
    FgenCrossoverFunc crossover_func;
    switch (crossover_type) {
    case CROSSOVER_TYPE_ORDER1 :
        crossover_func = fgen_crossover_permutation_order_based;
        break;
    case CROSSOVER_TYPE_PBX :
        crossover_func = fgen_crossover_permutation_position_based;
        break;
    case CROSSOVER_TYPE_NOOP :
        crossover_func = fgen_crossover_noop;
        break;
    }
    FgenMutationFunc mutation_func;
    mutation_func = fgen_mutation_permutation_swap;
    if (mutation_type == MUTATION_TYPE_SWAP_NEIGHBOURS)
        mutation_func = mutation_permutation_swap_neighbours;
    if (mutation_type == MUTATION_TYPE_INSERT)
        mutation_func = fgen_mutation_permutation_insert;
    if (mutation_type == MUTATION_TYPE_INVERT)
        mutation_func = fgen_mutation_permutation_invert;
    FgenGenerationCallbackFunc generation_callback_func = ProblemGenerationCallback;
    if (benchmark_mode)
        generation_callback_func = benchmark_generation_callback;
    pop = fgen_create(
        population_size,	/* Population size. */
        32 * nu_cities,		/* Individual size in bits. */
        32,			/* Data element size. */
        generation_callback_func,
        ProblemCalculateFitness,
        fgen_seed_permutation_random,
        mutation_func,
        crossover_func);
    float mut_prob = mutation_probability;
    if (mutation_type == MUTATION_TYPE_NOOP)
        mut_prob = 0;
    if (steady_state)
        fgen_set_parameters(
            pop,
            FGEN_KILL_TOURNAMENT,	/* Selection type. */
            FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
            crossover_probability,	/* Crossover probability. */
            mut_prob,			/* Per individual for permutation mutation. */
            0);				/* Macro-mutation probability. */
    else
        fgen_set_parameters(
            pop,
            FGEN_ELITIST_SUS,		/* Selection type. */
            FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
            crossover_probability,	/* Crossover probability. */
            mut_prob,			/* Per individual for permutation mutation. */
            0);				/* Macro-mutation probability. */
    fgen_set_permutation_size(pop, nu_cities);
    fgen_random_seed_with_timer(fgen_get_rng(pop));
    return pop;
}

void start_ga() {
    FgenPopulation *pop = setup_ga(0);
    ga_busy = 1;
    if (steady_state)
        fgen_run_steady_state(pop, - 1);
    else
        fgen_run(pop, - 1);
    fgen_destroy(pop);
    ga_busy = 0;
    gui_draw_and_show_window();
}

void stop_ga() {
    fgen_signal_stop(pop);
}

// Benchmark functionality.

static void benchmark_generation_callback(FgenPopulation *pop, int generation) {
    gui_handle_events();
    if (generation < 5000)
        return;
    FgenIndividual *best = fgen_best_individual_of_population(pop);
    int *perm = (int *)best->bitstring;
    memcpy(best_route, perm, sizeof(int) * nu_cities);
    best_route_valid = 1;
    ga_generation = generation;
    fgen_signal_stop(pop);
}

static char *mutation_type_str[] = { "SWAP", "SWPN", "INST", "INVT", "NOOP" };
static char *crossover_type_str[] = { "OX1 ", "PBX ", "NOOP" };

void benchmark_current_settings() {
    ga_busy = 1;
    printf("Running benchmarks (#cities = %d, population size = %d, 5000 generations):\n", nu_cities, population_size);
    printf("Mutation probability: %.3lf Crossover probability: %.3lf\n", mutation_probability, crossover_probability);
    printf("Mutation type: %s Crossover type: %s\n", mutation_type_str[mutation_type], crossover_type_str[crossover_type]);
    double total_dist = 0;
    for (int i = 0; i < 10; i++) {
        FgenPopulation *pop = setup_ga(1);
        fgen_run(pop, - 1);
        fgen_destroy(pop);
        total_dist += best_distance();
    }
    double average_dist = total_dist / 10;
    printf("Average distance: %8.3lf\n", average_dist);
    ga_busy = 0;
}

void benchmark_operators() {
    ga_busy = 1;
    printf("Running benchmarks (#cities = %d, population size = %d, 5000 generations):\n", nu_cities, population_size);
    printf("Mutation probability: %.3lf Crossover probability: %.3lf\n", mutation_probability, crossover_probability);
    int saved_mutation_type = mutation_type;
    int saved_crossover_type = crossover_type;  
    for (mutation_type = 0; mutation_type <= 4; mutation_type++)
        for (crossover_type = 0; crossover_type <= 2; crossover_type++) {
            double total_dist = 0;
            for (int i = 0; i < 10; i++) {
                FgenPopulation *pop = setup_ga(1);
                fgen_run(pop, - 1);
                fgen_destroy(pop);
                total_dist += best_distance();
            }
            double average_dist = total_dist / 10;
            printf("Mutation type: %s Crossover type: %s Average distance: %8.3lf\n", mutation_type_str[mutation_type],
                crossover_type_str[crossover_type], average_dist);
        }
    mutation_type = saved_mutation_type;
    crossover_type = saved_crossover_type;
    ga_busy = 0;
}

void benchmark_rates() {
    ga_busy = 1;
    printf("Running benchmarks (#cities = %d, population size = %d, 5000 generations):\n", nu_cities, population_size);
    printf("Mutation type: %s Crossover type: %s\n", mutation_type_str[mutation_type], crossover_type_str[crossover_type]);
    float saved_mutation_probability = mutation_probability;
    float saved_crossover_probability = crossover_probability;  
    for (mutation_probability = 0; mutation_probability < 1.01; mutation_probability += 0.10)
        for (crossover_probability = 0.2; crossover_probability < 1.01; crossover_probability += 0.2) {
            double total_dist = 0;
            for (int i = 0; i < 10; i++) {
                FgenPopulation *pop = setup_ga(1);
                fgen_run(pop, - 1);
                fgen_destroy(pop);
                total_dist += best_distance();
            }
            double average_dist = total_dist / 10;
            printf("Mutation prob: %.3lf Crossover prob: %.3lf Average distance: %8.3lf\n",
                mutation_probability, crossover_probability, average_dist);
        }
    mutation_probability = saved_mutation_probability;
    crossover_probability = saved_crossover_probability;
    ga_busy = 0;
}
