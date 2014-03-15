/*
    steady_state.c -- implements a steady-state genetic algorithm.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2013, Harm Hanemaaijer <fgenfb at yahoo.com>

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
#ifdef __GNUC__
#include <unistd.h>
#endif
#include <math.h>
#include <limits.h>
#include <float.h>
#include <pthread.h>
#include "fgen.h"
#include "parameters.h"
#include "population.h"
#include "selection.h"
#include "crossover.h"
#include "mutation.h"
#include "random.h"
#include "cache.h"
#include "migration.h"
#include "error.h"
#include "win32_compat.h"

static int TournamentSelection(FgenPopulation *pop) {
	int winner;
	int j;
	double best_fitness = NEGATIVE_INFINITY_DOUBLE;
	FgenIndividual **ind = pop->ind;
	for (j = 0; j < pop->tournament_size; j++) {
		int k = fgen_random_n(pop->rng, pop->size);
		CalculateIndividualFitness(pop, ind[k]);
		if (ind[k]->fitness > best_fitness) {
			winner = k;
			best_fitness = ind[k]->fitness;
		}
	}
	return winner;
}

static int TournamentSelectionWorst(FgenPopulation *pop) {
	int winner;
	int j;
	double worst_fitness = POSITIVE_INFINITY_DOUBLE;
	FgenIndividual **ind = pop->ind;
	for (j = 0; j < pop->tournament_size; j++) {
		int k = fgen_random_n(pop->rng, pop->size);
		CalculateIndividualFitness(pop, ind[k]);
		if (ind[k]->fitness < worst_fitness) {
			winner = k;
			worst_fitness = ind[k]->fitness;
		}
	}
	return winner;
}

static void DoCrossoverSteadyState(FgenPopulation *pop, int i, int j, int k, int l) {
	if (pop->crossover_probability == 0 || pop->fgen_crossover_func == fgen_crossover_noop)
		return;
	int r;
	FgenIndividual **ind = pop->ind;
	r = RandomBits(pop->rng, 8);
	if (r < pop->crossover_probability) {
		/* Perform crossover operation. */
		FgenIndividual *new_ind1 = NewIndividual(pop);
		FgenIndividual *new_ind2 = NewIndividual(pop);
		pop->fgen_crossover_func(pop, ind[i]->bitstring, ind[j]->bitstring, new_ind1->bitstring,
			new_ind2->bitstring);
		/* Replace two random individuals. */
		FreeIndividual(ind[k]);
		ind[k] = new_ind1;
		FreeIndividual(ind[l]);
		ind[l] = new_ind2;
	}
	else {
		FreeIndividual(ind[k]);
		ind[k] = ind[i];
		ind[i]->refcount++;
		FreeIndividual(ind[l]);
		ind[l] = ind[j];
		ind[j]->refcount++;
	}
}

static void DoMutationSteadyState(FgenPopulation *pop, int i, int j) {
	FgenIndividual *new_ind = NewIndividual(pop);
	CopyIndividualBitstring(pop, pop->ind[i], new_ind);
	pop->fgen_mutation_func(pop, pop->ind[i]->bitstring, new_ind->bitstring);
	FreeIndividual(pop->ind[i]);
	pop->ind[i] = new_ind;
	new_ind = NewIndividual(pop);
	CopyIndividualBitstring(pop, pop->ind[j], new_ind);
	pop->fgen_mutation_func(pop, pop->ind[j]->bitstring, new_ind->bitstring);
	FreeIndividual(pop->ind[j]);
	pop->ind[j] = new_ind;
}

static void SteadyStateGeneration(FgenPopulation *pop) {
	/* Selection: select two random unique individuals using tournament selection. */
	int i = TournamentSelection(pop);
	int j;
	do {
		j = TournamentSelection(pop);
	} while (j == i);
	/* Select two random unique individuals that will be replaced using tournament selection */
	/* where the worst individual in the tournament is selected. */
	int k;
	do {
		if (pop->selection_type & FGEN_KILL_TOURNAMENT_ELEMENT)
			k = TournamentSelectionWorst(pop);
		else
			k = fgen_random_n(pop->rng, pop->size);
	} while (k == i || k == j);
	int l;
	do {
		if (pop->selection_type & FGEN_KILL_TOURNAMENT_ELEMENT)
			l = TournamentSelectionWorst(pop);
		else
			l = fgen_random_n(pop->rng, pop->size);
	} while (l == i || l == j || l == k);

	DoCrossoverSteadyState(pop, i, j, k, l);

	DoMutationSteadyState(pop, k, l);
}

/**
 * Run the steady state version of the genetic algorithm. In steady state evolution, only two individuals are selected
 * for possible crossover and mutation in each "generation" and replace two random individuals in the population. Only
 * tournament selection is supported, and elitism and extinction are not supported. The selection type must be
 * FGEN_TOURNAMENT or FGEN_KILL_TOURNAMENT; in the latter case, the random individuals selected for replacement are selected
 * with tournament selection of the worst individual. The generation callback interval should be set to a high value for
 * efficiency, as a rule of thumb generation intervals should multiplied by half the population size compared to the
 * regular genetic algorithm.
 *
 * @param pop The population.
 * @param max_generation The maximum generation. When this is reached (or when stop is signalled) the function exits.
 * A value of - 1 means the function only returns when stop is signalled.
 */

void fgen_run_steady_state(FgenPopulation *pop, int max_generation) {
	if ((pop->selection_type & FGEN_STOCHASTIC_TYPE_MASK) != FGEN_TOURNAMENT)
		gen_error("Error: fgen_run_steady_state: only tournament selection allowed.");
	if (pop->selection_type & FGEN_ELITIST_ELEMENT)
		gen_error("Error: fgen_run_steady_state: elitism not supported in steady-state evolution.");
	if (pop->selection_type & FGEN_EXTINCTION_ELEMENT)
		gen_error("Error: fgen_run_steady_state: extinction not supported in steady-state evolution.");

	CreateInitialPopulation(pop);

	pop->generation = 0;
	pop->stop_signalled = 0;
        pop->fgen_generation_callback_func(pop, pop->generation);
	if (pop->stop_signalled)
		goto end;

	for (;;) {
		SteadyStateGeneration(pop);

		pop->generation++;
		if (pop->generation % pop->generation_callback_interval == 0)
			pop->fgen_generation_callback_func(pop, pop->generation);
		if ((max_generation != - 1 && pop->generation >= max_generation) || pop->stop_signalled)
			break;
	}
end :
	CalculatePopulationFitness(pop, NULL, NULL);
}

/**
 * Run the steady state version of the genetic algorithm on an archipelago. In steady state evolution, only two individuals
 * are selected for possible crossover and mutation in each "generation" and replace two random individuals in the population.
 * Only tournament selection is supported, and elitism and extinction are not supported. The selection type must be
 * FGEN_TOURNAMENT or FGEN_KILL_TOURNAMENT; in the latter case, the random individuals selected for replacement are selected
 * with tournament selection of the worst individual. The generation callback interval should be set to a high value for
 * efficiency. The seeds of the random number generators of the islands after the first one are initialized with random seeds
 * derived from the first island's random number generator.
 *
 * @param nu_pops Number of islands.
 * @param pops An array of populations.
 * @param max_generation The maximum generation. When this is reached (or when stop is signalled) the function exits.
 * A value of - 1 means the function only returns when stop is signalled.
 */

void fgen_run_steady_state_archipelago(int nu_pops, FgenPopulation **pops, int max_generation) {
	int i;
	for (i = 0; i < nu_pops; i++)  {
		if ((pops[i]->selection_type & FGEN_STOCHASTIC_TYPE_MASK) != FGEN_TOURNAMENT)
			gen_error("Error: fgen_run_steady_state_archipelago: only tournament selection allowed.");
		if (pops[i]->selection_type & FGEN_ELITIST_ELEMENT)
			gen_error("Error: fgen_run_steady_state_archipelago: elitism not supported in steady-state evolution.");
		if (pops[i]->selection_type & FGEN_EXTINCTION_ELEMENT)
			gen_error("Error: fgen_run_steady_state_archipelago: extinction not supported in steady-state "
				"evolution.");
	}

	/* Randomize the seeds of the random number generators of the populations after */
	/* the first one based on random number from the first random number generator. */
	for (i = 1; i < nu_pops; i++)
		fgen_random_seed_rng(pops[i]->rng, fgen_random_32(pops[0]->rng));

	for (i = 0; i < nu_pops; i++) {
		CreateInitialPopulation(pops[i]);
		pops[i]->generation = 0;
		pops[i]->stop_signalled = 0;
		pops[i]->island = i;
	        pops[i]->fgen_generation_callback_func(pops[i], pops[i]->generation);
	}
	for (i = 0; i < nu_pops; i++)
		if (pops[i]->stop_signalled)
			goto end;

	for (;;) {
		for (i = 0; i < nu_pops; i++) {
			SteadyStateGeneration(pops[i]);

			pops[i]->generation++;
		}

		if (pops[0]->migration_interval != 0)
			if (pops[0]->generation % pops[0]->migration_interval == 0)
				DoMigration(nu_pops, pops);

		for (i = 0; i < nu_pops; i++)
			if (pops[i]->generation % pops[i]->generation_callback_interval == 0)
				pops[i]->fgen_generation_callback_func(pops[i], pops[i]->generation);
		for (i = 0; i < nu_pops; i++)
			if ((max_generation != - 1 && pops[i]->generation >= max_generation) || pops[i]->stop_signalled)
				goto end;
	}
end :
	for (i = 0; i < nu_pops; i++)
		CalculatePopulationFitness(pops[i], NULL, NULL);
}

/* Threaded steady-state archipelago. This is almost identical to the generational threaded archipelago in ga.c. */

static void *fgen_steady_state_create_initial_population_thread(void *threadarg) {
	FgenPopulation *pop = (FgenPopulation *)threadarg;
	CreateInitialPopulation(pop);
	CalculatePopulationFitness(pop, NULL, NULL);
	pthread_exit(NULL);
	return NULL;
}

typedef struct {
	FgenPopulation *pop;
	int nu_generations;
	pthread_cond_t *condition_ready;
	pthread_cond_t *condition_start;
	pthread_cond_t *condition_end;
	pthread_cond_t *condition_continue;
	pthread_mutex_t *mutex_ready;
	pthread_mutex_t *mutex_start;
	pthread_mutex_t *mutex_end;
	pthread_mutex_t *mutex_continue;
	int *ready_signalled_count;
	int *start_signalled;
	int *end_signalled_count;
	int *continue_signalled;
} ThreadData;

static void *fgen_steady_state_do_multiple_generations_cond_thread(void *threadarg) {
	ThreadData *thread_data = (ThreadData *)threadarg;
	FgenPopulation *pop = thread_data->pop;
	for (;;) {
		int nu_generations;
		int i;
		/* Give the ready signal. */
		pthread_mutex_lock(thread_data->mutex_ready);
		(*thread_data->ready_signalled_count)++;
		pthread_cond_signal(thread_data->condition_ready);
		pthread_mutex_unlock(thread_data->mutex_ready);
		/* Wait for the start signal. */
		pthread_mutex_lock(thread_data->mutex_start);
		while (!(*thread_data->start_signalled))
			pthread_cond_wait(thread_data->condition_start, thread_data->mutex_start);
		pthread_mutex_unlock(thread_data->mutex_start);
		/* Do the work. */
		nu_generations = thread_data->nu_generations;
		if (nu_generations == 0)
			pthread_exit(NULL);
		i = 0;
		for (i = 0; i < nu_generations; i++) {
			SteadyStateGeneration(pop);
			pop->generation++;
			/* A stop may be signalled from the fitness function. */
			if (pop->stop_signalled) {
				pop->generation += nu_generations - i - 1;
				break;
			}
		}
		/* Take advantage of the multiple threads to update the fitness now, instead of */
		/* later in single-threaded mode. */
		CalculatePopulationFitness(pop, NULL, NULL);
		/* Signal the end of the work. */
		pthread_mutex_lock(thread_data->mutex_end);
		(*thread_data->end_signalled_count)++;
		pthread_cond_signal(thread_data->condition_end);
		pthread_mutex_unlock(thread_data->mutex_end);
		/* Wait for the continue signal. */
		pthread_mutex_lock(thread_data->mutex_continue);
		while (!(*thread_data->continue_signalled))
			pthread_cond_wait(thread_data->condition_continue, thread_data->mutex_continue);
		pthread_mutex_unlock(thread_data->mutex_continue);
	}
}

/**
 * Run the steady state genetic algorithm on an archipelago of populations. Threaded version, uses one thread per island.
 * The seeds of the random number generators of the islands after the first one are initialized with random seeds
 * derived from the first island's random number generator.
 *
 * @param nu_pops Number of islands.
 * @param pops An array of populations.
 * @param max_generation The maximum generation. When this is reached (or when stop is signalled) the function exits.
 * A value of - 1 means the function only returns when stop is signalled.
 */

void fgen_run_steady_state_archipelago_threaded(int nu_pops, FgenPopulation **pops, int max_generation) {
	int i;
	pthread_t *thread;
	int max_threads;
	ThreadData *thread_data;
	pthread_cond_t *condition_ready;
	pthread_cond_t *condition_start;
	pthread_cond_t *condition_end;
	pthread_cond_t *condition_continue;
	pthread_mutex_t *mutex_ready;
	pthread_mutex_t *mutex_start;
	pthread_mutex_t *mutex_end;
	pthread_mutex_t *mutex_continue;
	int ready_signalled_count;
	int *start_signalled;
	int end_signalled_count;
	int continue_signalled;

	for (i = 0; i < nu_pops; i++)  {
		if ((pops[i]->selection_type & FGEN_STOCHASTIC_TYPE_MASK) != FGEN_TOURNAMENT)
			gen_error("Error: fgen_run_steady_state_archipelago_threaded_: only tournament selection allowed.");
		if (pops[i]->selection_type & FGEN_ELITIST_ELEMENT)
			gen_error("Error: fgen_run_steady_state_archipelago_threaded: elitism not supported in steady-state "
				"evolution.");
		if (pops[i]->selection_type & FGEN_EXTINCTION_ELEMENT)
			gen_error("Error: fgen_run_steady_state_archipelago_threaded: extinction not supported "
				"in steady-state evolution.");
	}

	/* Randomize the seeds of the random number generators of the populations after */
	/* the first one based on random number from the first random number generator. */
	for (i = 1; i < nu_pops; i++)
		fgen_random_seed_rng(pops[i]->rng, fgen_random_32(pops[0]->rng));

	/* In practice, just running a lot of threads seems to be at least as efficient as imposing a maximum */
	/* on the amount of concurrently active threads. */
//	max_threads = sysconf(_SC_NPROCESSORS_ONLN);
//	for (i = 0; i < nu_pops; i++)
//		pops[i]->max_threads = max_threads;
//	printf("fgen: number of CPU cores detected = %d.\n", max_threads);
	max_threads = nu_pops;

	thread = (pthread_t *)malloc(sizeof(pthread_t) * nu_pops);
	thread_data = (ThreadData *)malloc(sizeof(ThreadData) * nu_pops);
	condition_ready = (pthread_cond_t *)malloc(sizeof(pthread_cond_t));
	condition_start = (pthread_cond_t *)malloc(sizeof(pthread_cond_t) * nu_pops);
	condition_end = (pthread_cond_t *)malloc(sizeof(pthread_cond_t));
	condition_continue = (pthread_cond_t *)malloc(sizeof(pthread_cond_t));
	mutex_ready = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
	mutex_start = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
	mutex_end = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
	mutex_continue = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
	start_signalled = (int *)malloc(sizeof(int) * nu_pops);

	/* Create the initial populations and calculate fitness, threaded. */
	for (i = 0; i < nu_pops; i++) {
		pops[i]->generation = 0;
		pops[i]->stop_signalled = 0;
		pops[i]->island = i;
		pthread_create(&thread[i], NULL, fgen_steady_state_create_initial_population_thread, pops[i]);
	}
	for (i = 0; i < nu_pops; i++)
		pthread_join(thread[i], NULL);
	for (i = 0; i < nu_pops; i++)
	        pops[i]->fgen_generation_callback_func(pops[i], pops[i]->generation);
	for (i = 0; i < nu_pops; i++)
		if (pops[i]->stop_signalled)
			goto end2;

	/* Initialize the condition variables and mutex. */
	pthread_cond_init(condition_ready, NULL);
	for (i = 0; i < nu_pops; i++)
		pthread_cond_init(&condition_start[i], NULL);
	pthread_cond_init(condition_end, NULL);
	pthread_cond_init(condition_continue, NULL);
	pthread_mutex_init(mutex_ready, NULL);
	pthread_mutex_init(mutex_start, NULL);
	pthread_mutex_init(mutex_end, NULL);
	pthread_mutex_init(mutex_continue, NULL);

	/* Create a worker thread for each island. */
	ready_signalled_count = 0;
	for (i = 0; i < nu_pops; i++) {
		start_signalled[i] = 0;
		thread_data[i].pop = pops[i];
		thread_data[i].condition_ready = condition_ready;
		thread_data[i].condition_start = &condition_start[i];
		thread_data[i].condition_end = condition_end;
		thread_data[i].condition_continue = condition_continue;
		thread_data[i].mutex_ready = mutex_ready;
		thread_data[i].mutex_start = mutex_start;
		thread_data[i].mutex_end = mutex_end;
		thread_data[i].mutex_continue = mutex_continue;
		thread_data[i].ready_signalled_count = &ready_signalled_count;
		thread_data[i].start_signalled = &start_signalled[i];
		thread_data[i].end_signalled_count = &end_signalled_count;
		thread_data[i].continue_signalled = &continue_signalled;
		pthread_create(&thread[i], NULL, fgen_steady_state_do_multiple_generations_cond_thread,
			&thread_data[i]);
	}

	for (;;) {
		int nu_generations_to_run;
		int starting_thread;

		/* Figure out how many generations can be run without having to call the generation */
		/* callback function or do migration, and then run them concurrently. */
		if (max_generation == - 1)
			nu_generations_to_run = INT_MAX - pops[0]->generation;
		else
			nu_generations_to_run = max_generation - pops[0]->generation;
		for (i = pops[0]->generation + 1; i < pops[0]->generation + nu_generations_to_run; i++)
			if ((pops[0]->migration_interval != 0 && i % pops[0]->migration_interval == 0) ||
			i % pops[0]->generation_callback_interval == 0) {
				nu_generations_to_run = i - pops[0]->generation;
				break;
			}

		/* Check whether all the workers are ready. */
		pthread_mutex_lock(mutex_ready);
		while (ready_signalled_count < nu_pops) {
			pthread_cond_wait(condition_ready, mutex_ready);
		}
		pthread_mutex_unlock(mutex_ready);
		ready_signalled_count = 0;
		continue_signalled = 0;

		/* Run max_threads concurrent threads. Of the nu_pops threads, only allow max_threads to be */
		/* active at the same time. */
		for (starting_thread = 0; starting_thread < nu_pops; starting_thread += max_threads) {
			/* Calculate the number of active threads. */
			int nu_threads = max_threads;
			if (starting_thread + max_threads > nu_pops)
				nu_threads = nu_pops - starting_thread;
			/* Activate the worker threads for each selected island. */
			for (i = starting_thread; i < starting_thread + max_threads && i < nu_pops; i++)
				thread_data[i].nu_generations = nu_generations_to_run;
			end_signalled_count = 0;
			for (i = starting_thread; i < starting_thread + max_threads && i < nu_pops; i++) {
				pthread_mutex_lock(mutex_start);
				start_signalled[i] = 1;
				pthread_cond_signal(&condition_start[i]);
				pthread_mutex_unlock(mutex_start);
			}

			/* Wait for all of them to finish. */
			pthread_mutex_lock(mutex_end);
			while (end_signalled_count < nu_threads) {
				pthread_cond_wait(condition_end, mutex_end);
			}
			pthread_mutex_unlock(mutex_end);
	
			for (i = starting_thread; i < starting_thread + max_threads && i < nu_pops; i++)
				start_signalled[i] = 0;
		}
		/* Give the signal to continue. */
		pthread_mutex_lock(mutex_continue);
		continue_signalled = 1;
		pthread_cond_broadcast(condition_continue);
		pthread_mutex_unlock(mutex_continue);

		if (pops[0]->migration_interval != 0)
			if (pops[0]->generation % pops[0]->migration_interval == 0)
				DoMigration(nu_pops, pops);

		for (i = 0; i < nu_pops; i++) {
			if (pops[i]->generation % pops[i]->generation_callback_interval == 0)
				pops[i]->fgen_generation_callback_func(pops[i], pops[i]->generation);
		}
		for (i = 0; i < nu_pops; i++)
			if ((max_generation != - 1 && pops[i]->generation >= max_generation) || pops[i]->stop_signalled)
				goto end;
	}
end :
	/* Instruct the threads to exit. */
	for (i = 0; i < nu_pops; i++)
		thread_data[i].nu_generations = 0;
	for (i = 0; i < nu_pops; i++) {
		pthread_mutex_lock(mutex_start);
		start_signalled[i] = 1;
		pthread_cond_signal(&condition_start[i]);
		pthread_mutex_unlock(mutex_start);
	}
	for (i = 0; i < nu_pops; i++)
		pthread_join(thread[i], NULL);
	/* Free all data structures. */
	pthread_cond_destroy(condition_ready);
	for (i = 0; i < nu_pops; i++)
		pthread_cond_destroy(&condition_start[i]);
	pthread_cond_destroy(condition_end);
	pthread_cond_destroy(condition_continue);
	pthread_mutex_destroy(mutex_ready);
	pthread_mutex_destroy(mutex_start);
	pthread_mutex_destroy(mutex_end);
	pthread_mutex_destroy(mutex_continue);
end2 :
	free(thread);
	free(thread_data);
	free(condition_ready);
	free(condition_start);
	free(condition_end);
	free(condition_continue);
	free(mutex_ready);
	free(mutex_start);
	free(mutex_end);
	free(mutex_continue);
	free(start_signalled);
}

