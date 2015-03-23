/*
    population.c -- genetic algorithm population handling functions.

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
#include <limits.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include "fgen.h"
#include "parameters.h"
#include "bitstring.h"
#include "error.h"
#include "population.h"
#include "cache.h"
#include "win32_compat.h"

FgenIndividual **AllocatePopulation(FgenPopulation *pop) {
	return (FgenIndividual **)malloc(sizeof(FgenIndividual *) * pop->size);
}

void CreateInitialPopulation(FgenPopulation *pop) {
	int i;
	if (pop->ind == NULL)
		// Allocate a population.
		pop->ind = AllocatePopulation(pop);
	else {
		// A population already exists.
		if (pop->initialization_type == FGEN_INITIALIZATION_CONTINUE) {
			for (int i = 0; i < pop->size; i++)
				pop->ind[i]->fitness_is_valid = 0;
			return;
		}
		// FGEN_INITIALIZATION_SEED is set; free the individuals.
		FreePopulation(pop, pop->ind);
	}
	for (i = 0; i < pop->size; i++) {
		pop->ind[i] = NewIndividual(pop);
		pop->fgen_seed_func(pop, pop->ind[i]->bitstring);
	}
}

void CalculateIndividualFitness(FgenPopulation *pop, FgenIndividual *ind) {
	double cache_fitness;
	int hash;
	if (ind->fitness_is_valid) {
		return;
	}
	if (pop->cache != NULL)
		if (CacheHit(pop, ind, &cache_fitness, &hash)) {
			ind->fitness = cache_fitness;
			ind->fitness_is_valid = 1;
			return;
		}
	ind->fitness = pop->fgen_calculate_fitness_func(pop, ind->bitstring);
	ind->fitness_is_valid = 1;
	if (pop->cache != NULL)
		AddToCache(pop, ind, hash);
}

void CalculatePopulationFitness(FgenPopulation *pop, double *sum_out,
double *min_fitness_out ) {
	int i;
        double sum, min_fitness;
	sum = 0;
	min_fitness = POSITIVE_INFINITY_DOUBLE;
	for (i = 0; i < pop->size; i++) {
		CalculateIndividualFitness(pop, pop->ind[i]);
		// Be careful not to overflow on systems were infinity is finite. 
		if (pop->ind[i]->fitness == POSITIVE_INFINITY_DOUBLE)
			sum = POSITIVE_INFINITY_DOUBLE;
		else
			sum += pop->ind[i]->fitness;
		if (pop->ind[i]->fitness < min_fitness)
			min_fitness = pop->ind[i]->fitness;
	}
	if (sum_out != NULL)
		*sum_out = sum;
	if (min_fitness_out != NULL)
		*min_fitness_out = min_fitness;
}

/* Threaded version. Doesn't work with cache enabled. */

typedef struct {
	FgenPopulation *pop;
	int start_index;
	int end_index;
} ThreadData;

static void *CalculatePartPopulationFitness(void *threadarg) {
	ThreadData *thread_data = (ThreadData *)threadarg;
	FgenPopulation *pop = thread_data->pop;
	int i;
	for (i = thread_data->start_index; i < thread_data->end_index; i++)
		CalculateIndividualFitness(pop, pop->ind[i]);
	return NULL;
}

void CalculatePopulationFitnessThreaded(FgenPopulation *pop) {
	pthread_t *thread;
	ThreadData *thread_data;
	int i;
	thread = (pthread_t *)malloc(sizeof(pthread_t) * pop->max_threads);
	thread_data = (ThreadData *)malloc(sizeof(ThreadData) * pop->max_threads);
	/* Divide the population into n parts. */
	for (i = 0; i < pop->max_threads; i++) {
		thread_data[i].pop = pop;
		thread_data[i].start_index = pop->size * i / pop->max_threads;
		thread_data[i].end_index = pop->size * (i + 1) / pop->max_threads;
	}
	/* Create the first n - 1 threads. */
	for (i = 0; i < pop->max_threads - 1; i++)
		pthread_create(&thread[i], NULL, CalculatePartPopulationFitness, &thread_data[i]);
	/* The current thread is the n-th thread. */
	CalculatePartPopulationFitness(&thread_data[pop->max_threads - 1]);
	/* Wait for the threads to finish. */
	for (i = 0; i < pop->max_threads - 1; i++)
		pthread_join(thread[i], NULL);
	free(thread_data);
	free(thread);
}

/* Heavily threaded version that creates a new thread for each fitness calculation. */

typedef struct {
	FgenPopulation *pop;
	int index;
} ThreadData2;

void *CalculateIndividualFitnessThreaded(void *threadarg) {
	ThreadData2 *thread_data = (ThreadData2 *)threadarg;
	FgenPopulation *pop = thread_data->pop;
	int index = thread_data->index;
	CalculateIndividualFitness(pop, pop->ind[index]);
	return NULL;
}

void CalculatePopulationFitnessHeavilyThreaded(FgenPopulation *pop) {
	pthread_t *thread;
	int max_threads;
	int nu_threads;
	int i;
	ThreadData2 *thread_data;
	max_threads = 8;
	thread = (pthread_t *)malloc(sizeof(pthread_t) * max_threads);
	thread_data = (ThreadData2 *)malloc(sizeof(ThreadData) * max_threads);
	nu_threads = 0;
	for (i = 0; i < pop->size; i++) {
		if (pop->ind[i]->fitness_is_valid)
			continue;
		if (nu_threads == max_threads) {
			int j;
			for (j = 0; j < max_threads; j++)
				pthread_join(thread[j], NULL);
			nu_threads = 0;
		}
		thread_data[nu_threads].pop = pop;
		thread_data[nu_threads].index = i;
		pthread_create(&thread[nu_threads], NULL, CalculateIndividualFitnessThreaded,
			&thread_data[nu_threads]);
		nu_threads++;
	}
	for (i = 0; i < nu_threads; i++)
		pthread_join(thread[i], NULL);
	free(thread);
	free(thread_data);
}


void CopyIndividualBitstring(FgenPopulation *pop, FgenIndividual *src, FgenIndividual *dest ) {
	memcpy(dest->bitstring, src->bitstring, INDIVIDUAL_SIZE_IN_BYTES(pop));
}

void FreeIndividual(FgenIndividual *ind) {
	if (ind->refcount < 1)
		gen_error("FreeIndividual: Reference count already zero.");
	ind->refcount--;
	if (ind->refcount == 0) {
		free(ind->bitstring);
		free(ind);
	}
}

FgenIndividual *NewIndividual(FgenPopulation *pop) {
	FgenIndividual *ind;
	ind = (FgenIndividual *)malloc(sizeof(FgenIndividual));
	ind->bitstring = (unsigned char *)malloc(INDIVIDUAL_SIZE_IN_BYTES(pop));
	ind->fitness_is_valid = 0;
	ind->is_elite = 0;
	ind->refcount = 1;	/* Assumption. */
	return ind;
}

void FreePopulation(FgenPopulation *pop, FgenIndividual **ind) {
	int i;
	for (i = 0; i < pop->size; i++)
		FreeIndividual(ind[i]);
}

/**
 * Returns the best individual of the population.
 *
 * @param pop The population.
 * @return A pointer to the best individual.
 */

FgenIndividual *fgen_best_individual_of_population(FgenPopulation *pop) {
	int i, best;
        double bestfitness;
        CalculatePopulationFitness(pop, NULL, NULL); /* Make sure every fitness is updated. */
	bestfitness = NEGATIVE_INFINITY_DOUBLE;
	for (i = 0; i < pop->size; i++)
		if (pop->ind[i]->fitness > bestfitness) {
			bestfitness = pop->ind[i]->fitness;
			best = i;
		}
	return pop->ind[best];
}

/**
 * Returns the worst individual of the population.
 *
 * @param pop The population.
 * @return A pointer to the worst individual.
 */

FgenIndividual *fgen_worst_individual_of_population(FgenPopulation *pop) {
	int i, worst;
        double worst_fitness;
        CalculatePopulationFitness(pop, NULL, NULL); /* Make sure every fitness is updated. */
	worst_fitness = POSITIVE_INFINITY_DOUBLE;
	for (i = 0; i < pop->size; i++)
		if (pop->ind[i]->fitness < worst_fitness) {
			worst_fitness = pop->ind[i]->fitness;
			worst = i;
		}
	return pop->ind[worst];
}

/**
 * Returns the best individual of an archipelago of populations.
 *
 * @param nu_pops The number of islands in the archipelago.
 * @param pops The archipelago (an array of populations).
 * @return A pointer to the best individual.
 */

FgenIndividual *fgen_best_individual_of_archipelago(int nu_pops, FgenPopulation **pops) {
	int i;
	double best_fitness;
	FgenIndividual *best_ind;
	best_fitness = NEGATIVE_INFINITY_DOUBLE;
	for (i = 0; i < nu_pops; i++) {
		FgenIndividual *ind;
		ind = fgen_best_individual_of_population(pops[i]);
		if (ind->fitness > best_fitness) {
			best_fitness = ind->fitness;
			best_ind = ind;
		}
	}
	return best_ind;
}

/**
 * Returns the best individual of an archipelago of populations and the island to which it
 * belongs.
 *
 * @param nu_pops The number of islands in the archipelago.
 * @param pops The archipelago (an array of populations).
 * @param island A pointer to the integer where the island index will be stored.
 * @return A pointer to the best individual.
 */

FgenIndividual *fgen_best_individual_and_island_of_archipelago(int nu_pops,
FgenPopulation **pops, int *island) {
	int i;
	double best_fitness;
	FgenIndividual *best_ind;
	best_fitness = NEGATIVE_INFINITY_DOUBLE;
	int best_island = 0;
	for (i = 0; i < nu_pops; i++) {
		FgenIndividual *ind;
		ind = fgen_best_individual_of_population(pops[i]);
		if (ind->fitness > best_fitness) {
			best_fitness = ind->fitness;
			best_ind = ind;
			best_island = i;
		}
	}
	*island = best_island;
	return best_ind;
}

/**
 * Update the fitness of every individual in the population.
 */

void fgen_update_population_fitness(FgenPopulation *pop) {
	CalculatePopulationFitness(pop, NULL, NULL);
}

/**
 * Invalidate the fitness of every individual so that it will be recalculated (because the fitness
 * function has changed).
 */

void fgen_invalidate_population_fitness(FgenPopulation *pop) {
	int i;
	for (i = 0; i < pop->size; i++)
		pop->ind[i]->fitness_is_valid = 0;
	if (pop->cache != NULL)
		fgen_invalidate_cache(pop);
}

/**
 * Returns the size of the bitstring (genetic data) of an individual in bytes.
 */

int fgen_individual_size_in_bytes(const FgenPopulation *pop) {
    return INDIVIDUAL_SIZE_IN_BYTES(pop);
}

