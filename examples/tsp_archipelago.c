/*
    sample_tsp_archipelago.c -- implementation of the Travelling Salesman Problem using an archipelago of
                                genetic algorithms

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


/*
 * Sample problem.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "fgen.h"

/*
 * Implement travelling salesman problem with GA archipelago.
 */


/*
 * Archipelago settings.
 */

#define NU_ISLANDS 8
/* One of the following three must be defined:
 * NO_CACHE		No fitness cache.
 * SEPERATE_CACHE	Seperate cache for each island.
 * SHARED_CACHE		One cache shared between all islands (works only unthreaded).
 *
 * The cache gets a respectable hit-rate, but because the fitness function is fast, the uncached version is
 * faster.
 */
#define NO_CACHE

/*
 * TSP definitions.
 */

#define NU_CITIES 32
#define AREA_SIZE 100.0

typedef struct {
	double x;
	double y;
} City;

static City city[NU_CITIES];
int steady_state;

static double random_d(double range) {
	return (double)rand() * range / ((double)RAND_MAX + 1);
}

static void InitializeCities() {
	int i;
	for (i = 0; i < NU_CITIES; i++) {
		city[i].x = random_d(AREA_SIZE);
		city[i].y = random_d(AREA_SIZE);
	}
}

static double ProblemCalculateFitness(const FgenPopulation *pop, const unsigned char *bitstring) {
	int *perm = (int *)bitstring;
	double dist;
	int i;
	dist = 0;
	for (i = 0; i < NU_CITIES - 1; i++)
		dist += sqrt((city[perm[i + 1]].x - city[perm[i]].x) * (city[perm[i + 1]].x - city[perm[i]].x) +
			(city[perm[i + 1]].y - city[perm[i]].y) * (city[perm[i + 1]].y - city[perm[i]].y));
	return AREA_SIZE * AREA_SIZE / dist;
}


void ProblemPrintIndividual(const FgenIndividual *ind) {
	int *perm = (int *)ind->bitstring;
	int i;
	for (i = 0; i < NU_CITIES; i++)
		printf("[%d]", perm[i]);
	printf(" dist = %lf\n", AREA_SIZE * AREA_SIZE / ind->fitness);
}

void PrintBestIndividual(const FgenIndividual *best, int generation) {
	printf("generation = %d, best fitness = %lf, solution:\n", generation,
		best->fitness);
	ProblemPrintIndividual(best);
}

void ProblemGenerationCallback(FgenPopulation *pop, int generation) {
	int interval = 100;
	if (steady_state)
		interval = 50000;
	if (generation % interval == 0) {
		printf("Island = %d, ", pop->island);
		FgenIndividual *best = fgen_best_individual_of_population(pop);
		PrintBestIndividual(best, generation);
		if (fgen_is_cached(pop)
#ifdef SHARED_CACHE
		&& fgen_get_island(pop) == NU_ISLANDS - 1
#endif
		)
			printf("Fitness cache hit-rate = %lf.\n", fgen_get_cache_hit_rate(pop));
	}
}

int main(int argc, char **argv) {
	FgenPopulation *pops[NU_ISLANDS];
	int i;

	steady_state = 0;
	if (argc > 1)
		if (strcmp(argv[1], "-s") == 0) {
			steady_state = 1;
			printf("Running with steady-state evolution.\n");
		}
	
	InitializeCities();
	for (i = 0; i < NU_ISLANDS; i++) {
		pops[i] = fgen_create(
			1024,			/* Population size. */
			32 * NU_CITIES,		/* Individual size in bits. */
			32,			/* Data element size. */
			ProblemGenerationCallback,
			ProblemCalculateFitness,
			fgen_seed_permutation_random,
			fgen_mutation_permutation_invert,
			fgen_crossover_permutation_order_based);
		if (steady_state)
			fgen_set_parameters(
				pops[i],
				FGEN_KILL_TOURNAMENT,		/* Selection type. */
				FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
				0.500,				/* Crossover probability. */
				0.500,				/* Per individual for permutation mutation. */
				0);				/* Macro-mutation probability. */
		else
			fgen_set_parameters(
				pops[i],
				FGEN_ELITIST_SUS,		/* Selection type. */
				FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
				0.200,				/* Crossover probability. */
				0.200,				/* Per individual for permutation mutation. */
				0);				/* Macro-mutation probability. */
		fgen_set_permutation_size(pops[i], NU_CITIES);
		fgen_set_migration_probability(pops[i], 0.001);
		if (steady_state) {
			fgen_set_generation_callback_interval(pops[i], 50000);
			fgen_set_migration_interval(pops[i], 10000);
		}
		else {
			fgen_set_generation_callback_interval(pops[i], 100);
			fgen_set_migration_interval(pops[i], 20);
		}
#ifdef SEPERATE_CACHE
		fgen_enable_cache(pops[i], 1024 * 128);	/* 128K entries. */
#endif
	}
#ifdef SHARED_CACHE
	fgen_enable_cache_on_archipelago(NU_ISLANDS, pops, 1024 * 1024);	/* 1MB entries. */
#endif
	/* Optionally seed the random number generator of the first population with system timer. */
	fgen_random_seed_with_timer(fgen_get_rng(pops[0]));

	int max_generations = 1000;
	if (steady_state) {
		max_generations = 500000;
		fgen_run_steady_state_archipelago_threaded(NU_ISLANDS, pops, max_generations);
	}
	else
		fgen_run_archipelago_threaded(NU_ISLANDS, pops, max_generations);
	printf("All islands: ");
	PrintBestIndividual(fgen_best_individual_of_archipelago(NU_ISLANDS, pops), max_generations);
        for (i = 0; i < NU_ISLANDS; i++)
		fgen_destroy(pops[i]);
    	exit(0);
}

