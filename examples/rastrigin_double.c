/*
    sampe_rastrigin_double.c -- fgen example with real-valued chromosomes for the Rastrigin function.

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
 * Rastrigin function example using fgen with custom seeding and mutation operator with real-valued
 * chromosomes.
 */

#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits.h>
#include "fgen.h"

#define NU_DIMENSIONS 20
#define DOMAIN_MIN - 5.12
#define DOMAIN_MAX 5.12

static void RastriginGenerationCallback(FgenPopulation *pop, int generation);
static double RastriginCalculateFitness(const FgenPopulation *pop, const unsigned char *bitstring);
static void RastriginSeed(FgenPopulation *pop, unsigned char *bitstring);
static void RastriginMutation(FgenPopulation *pop, const unsigned char *parent, unsigned char *child);
static void RastriginMutationSmall(FgenPopulation *pop, const unsigned char *parent, unsigned char *child);
static void RastriginMutationLarge(FgenPopulation *pop, const unsigned char *parent, unsigned char *child);
static void RastriginCrossover(FgenPopulation *pop, const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2);

int main(int argc, char **argv) {
	FgenPopulation *pop;
	pop = fgen_create(
		1024,			/* Population size. */
		64 * NU_DIMENSIONS,	/* Individual size in bits (NU_DIMENSIONS doubles)*/
		64,			/* Data element size. */
		RastriginGenerationCallback,
		RastriginCalculateFitness,
		RastriginSeed,				/* Custom seeding operator. */
		RastriginMutation,			/* Custom mutation operator. */
		RastriginCrossover);			/* Custom crossover operator. */
	/* Because we use a custom mutation function, the mutation probabilities are not used. */
	/* The custom crossover function is called with a frequency related to the crossover probability. */
	fgen_set_parameters(
		pop,
		FGEN_ELITIST_SUS,		/* Selection type. */
		FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
		0.800,				/* Crossover probability. */
		0.015,				/* Mutation probability (not used). */
		0);				/* Macro-mutation probability (not used). */
	fgen_run(pop, INT_MAX);
	fgen_destroy(pop);
    	exit(0);
}

/*
 * Example problem:
 *
 * Minimize the n-dimensional Rastrigin function:
 *
 * A * n + SIGMA[i = 1 to n](xi ^ 2 - A * cos(2 * PI * xi))
 *
 * where A = 10.
 *
 * The two-dimensional version is given by:
 *
 * f(x, y) = 10 * 2 + x ^ 2 + y ^ 2 - 10 * cos(2 * PI * x) - 10 * cos(2 * PI * y)
 */

#define sqr(x) ((x) * (x))

static double RastriginCalculateFitness(const FgenPopulation *pop, const unsigned char *bitstring) {
	double *param = (double *)bitstring;
	double result;
	int i;
	const double A = 10.0;
	const double n = NU_DIMENSIONS;
	result = A * n;
	for (i = 0; i < n; i++)
		result += sqr(param[i]) - A * cos(2 * M_PI * (param[i]));
	/* Scale to fitness so that lower function values have higher fitness */
	return 1.0 / result;
}

static void RastriginPrintIndividual(const double *param) {
	int i;
	printf("(");
	for (i = 0; i < NU_DIMENSIONS - 1; i++)
		printf("%lf, ", param[i]);
	printf("%lf)\n", param[i]);
}

static void RastriginGenerationCallback(FgenPopulation *pop, int generation) {
	FgenIndividual *best = fgen_best_individual_of_population(pop);
	if (generation % 200 == 0) {
		printf("Generation = %d, best error = %lf, solution ", generation, 1.0 / best->fitness);
		RastriginPrintIndividual((double *)best->bitstring);
	}
	if (generation == 1000)
		fgen_signal_stop(pop);
}

/* Custom seeding function that assignd parameters with random double floating point values within the range. */

static void RastriginSeed(FgenPopulation *pop, unsigned char *bitstring) {
	double *params = (double *)bitstring;
	int i;
	FgenRNG *rng;
	rng = fgen_get_rng(pop);
	for (i = 0; i < NU_DIMENSIONS; i++)
		params[i] = fgen_random_from_range_d(rng, DOMAIN_MIN, DOMAIN_MAX); 
}

/* Custom mutation function that assigns random values within the range to parameters with a mutation */
/* probability. When the function is called child already contains a copy of parent */

static void RastriginMutationLarge(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	double *params_parent = (double *)parent;
	double *params_child = (double *)child;
	FgenRNG *rng;
	int i;
	rng = fgen_get_rng(pop); 
	/* Mutate each parameter with a probability of 0.03. */
	for (i = 0; i < NU_DIMENSIONS; i++)
		if (fgen_random_f(rng, 1.0) < 0.03)
			params_child[i] = fgen_random_from_range_d(rng, DOMAIN_MIN, DOMAIN_MAX);
}

/* Alternative custom mutation function that tweaks the parameter by a maximum of 1/20th of the range. */

static void RastriginMutationSmall(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	double *params_parent = (double *)parent;
	double *params_child = (double *)child;
	FgenRNG *rng;
	int i;
	rng = fgen_get_rng(pop);
	/* Mutate each parameter with a probability of 0.03. */
	for (i = 0; i < NU_DIMENSIONS; i++)
		if (fgen_random_f(rng, 1.0) < 0.03) {
			double tweak = fgen_random_from_range_d(rng, - (DOMAIN_MAX - DOMAIN_MIN) / 20,
				(DOMAIN_MAX - DOMAIN_MIN) / 20);
			params_child[i] = params_parent[i] + tweak;
			/* Clamp to range. */
			if (params_child[i] < DOMAIN_MIN)
				params_child[i] = DOMAIN_MIN;
			if (params_child[i] > DOMAIN_MAX)
				params_child[i] = DOMAIN_MAX;
		}
}

/* Custom mutation function that is a combination of small and large mutations. */

static void RastriginMutation(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	RastriginMutationSmall(pop, parent, child);
	RastriginMutationLarge(pop, parent, child);
}

/* Custom real-valued crossover function. Intermediate recombination. */

static void RastriginCrossover(FgenPopulation *pop, const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	double *params_parent1 = (double *)parent1;
	double *params_parent2 = (double *)parent2;
	double *params_child1 = (double *)child1;
	double *params_child2 = (double *)child2;
	int i;
	FgenRNG *rng;
	rng = fgen_get_rng(pop);
	for (i = 0; i < NU_DIMENSIONS; i++) {
		double alpha;
		/* First offspring. */
		alpha = fgen_random_from_range_d(rng, - 0.25, 1.25);
		params_child1[i] = params_parent1[i] + alpha * (params_parent2[i] - params_parent1[i]);
		/* Clamp to range. */
		if (params_child1[i] < DOMAIN_MIN)
			params_child1[i] = DOMAIN_MIN;
		if (params_child1[i] > DOMAIN_MAX)
			params_child1[i] = DOMAIN_MAX;
		/* Second offspring. */
		alpha = fgen_random_from_range_d(rng, - 0.25, 1.25);
		params_child2[i] = params_parent2[i] + alpha * (params_parent1[i] - params_parent2[i]);
		/* Clamp to range. */
		if (params_child2[i] < DOMAIN_MIN)
			params_child2[i] = DOMAIN_MIN;
		if (params_child2[i] > DOMAIN_MAX)
			params_child2[i] = DOMAIN_MAX;
	}
}


