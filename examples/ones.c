/*
    sample_ones.c -- trivial example using fgen.

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
#include <fgen.h>

#define INDIVIDUAL_SIZE 1024

/* Prototypes of characteristic functions defined for the problem. */

static double ProblemCalculateFitness(const FgenPopulation *, const unsigned char *);

static void ProblemGenerationCallback(FgenPopulation *pop, int generation) {
	if (generation == 200) {
		printf("Generation = %d, best fitness = %lf.\n", generation,
			fgen_best_individual_of_population(pop)->fitness);
		if (fgen_is_cached(pop))
			printf("Fitness cache hit-rate = %lf.\n", fgen_get_cache_hit_rate(pop));
	}
}

int main(int argc, char **argv) {
	int i;
	FgenPopulation *pop;
	pop = fgen_create(
		1024,			/* Population size. */
		INDIVIDUAL_SIZE,	/* Individual size in bits. */
		1,			/* Data element size. */
		ProblemGenerationCallback,
		ProblemCalculateFitness,
		fgen_seed_random,
		fgen_mutation_per_bit_fast,
		fgen_crossover_one_point_per_bit);
	fgen_set_parameters(
		pop,
		FGEN_ELITIST_TOURNAMENT,	/* Selection type. */
		FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
		0.600,				/* Crossover probability. */
		0.001,				/* Mutation probability per bit. */
		0);				/* Macro-mutation probability. */
	fgen_random_seed_with_timer(fgen_get_rng(pop));
	fgen_set_initialization_type(pop, FGEN_INITIALIZATION_SEED);
	/* Do runs with several different parameters. */
	for (i = 2; i <= 4; i++) {
		fgen_set_tournament_size(pop, i);
		printf("Elitist tournament selection (tournament size = %d).\n", i);
		fgen_run(pop, 200);	/* Run for 200 generations. */
	}
	fgen_set_selection_type(pop, FGEN_TOURNAMENT);
	for (i = 2; i <= 4; i++) {
		fgen_set_tournament_size(pop, i);
		printf("Tournament selection (tournament size = %d).\n", i);
		fgen_run(pop, 200);	/* Run for 200 generations. */
	}
	fgen_set_selection_type(pop, FGEN_SUS);
	printf("Stochastic universal sampling selection.\n");
	fgen_run(pop, 200);
	fgen_set_selection_type(pop, FGEN_ELITIST_SUS);
	printf("Elitist stochastic universal sampling selection.\n");
	fgen_run(pop, 200);
	fgen_set_selection_type(pop, FGEN_RANK);
	printf("Rank-based selection.\n");
	fgen_run(pop, 200);
	fgen_set_selection_type(pop, FGEN_ELITIST_RANK);
	printf("Elitist rank-based selection.\n");
	fgen_run(pop, 200);
	fgen_destroy(pop);
    	exit(0);
}


/* Fitness is number of '1' bits. */

static double ProblemCalculateFitness(const FgenPopulation *pop, const unsigned char *bitstring) {
	int bits_to_go;
	int fitness;
	bits_to_go = INDIVIDUAL_SIZE;
	fitness = 0;
	while (bits_to_go >= 8) {
		unsigned int bits;
		bits = *bitstring;
		fitness += bits & 0x01;
		fitness += (bits & 0x02) >> 1;
		fitness += (bits & 0x04) >> 2;
		fitness += (bits & 0x08) >> 3;
		fitness += (bits & 0x10) >> 4;
		fitness += (bits & 0x20) >> 5;
		fitness += (bits & 0x40) >> 6;
		fitness += (bits & 0x80) >> 7;
		bitstring++;
		bits_to_go -= 8;
	}
	return fitness;
}


