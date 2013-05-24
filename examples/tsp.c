/*
    sample_tsp.c -- implementation of the Travelling Salesman Problem

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
#include <math.h>
#include <limits.h>
#include "fgen.h"

/*
 * Implement travelling salesman problem.
 */

#define NU_CITIES 32
#define AREA_SIZE 100.0

typedef struct {
	double x;
	double y;
} City;

static City city[NU_CITIES];

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


void ProblemPrintIndividual(const FgenPopulation *pop, const FgenIndividual *ind) {
	int *perm = (int *)ind->bitstring;
	int i;
	for (i = 0; i < NU_CITIES; i++)
		printf("[%d]", perm[i]);
	printf(" dist = %lf\n", AREA_SIZE * AREA_SIZE / ind->fitness);
}

void ProblemGenerationCallback(FgenPopulation *pop, int generation) {
	if (generation % 1000 == 0) {
		FgenIndividual *best = fgen_best_individual_of_population(pop);
		printf("Generation = %d, best fitness = %lf, solution:\n", generation, best->fitness);
		ProblemPrintIndividual(pop, best);
		if (fgen_is_cached(pop))
			printf("Fitness cache hit-rate = %lf.\n", fgen_get_cache_hit_rate(pop));
	}
}

int main(int argc, char **argv) {
	FgenPopulation *pop;
	InitializeCities();
	pop = fgen_create(
		256,			/* Population size. */
		32 * NU_CITIES,		/* Individual size in bits. */
		32,			/* Data element size. */
		ProblemGenerationCallback,
		ProblemCalculateFitness,
		fgen_seed_permutation_random,
		fgen_mutation_permutation_invert,
		fgen_crossover_permutation_order_based);
	fgen_set_parameters(
		pop,
		FGEN_ELITIST_SUS,		/* Selection type. */
		FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
		0.200,				/* Crossover probability. */
		0.200,				/* Per individual for permutation mutation. */
		0);				/* Macro-mutation probability. */
	fgen_set_permutation_size(pop, NU_CITIES);
/*	fgen_enable_cache(pop, 1024 * 1024); */
	fgen_random_seed_with_timer(fgen_get_rng(pop));
	fgen_run(pop, 10000);
	fgen_destroy(pop);
    	exit(0);
}


