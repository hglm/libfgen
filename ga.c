/*
    ga.c -- main module of the genetic algorithm implementation.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2012-13, Harm Hanemaaijer <fgenfb at yahoo.com>

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
 * This file contains the main module of the program; it controls the
 * high-level flow of the genetic algorithm.
 */

/** \mainpage 
 *
 * \section section_intro Introduction
 * 
 * libfgen is a library that implements an efficient and customizable genetic algorithm (GA). It also
 * provides particle swarm optimization (PSO) functionality and an interface for real-valued function minimization
 * or model fitting. It is written in C with a partial C++ wrapper.
 *
 * \section section_fgen fgen
 *
 * \link group_fgen fgen \endlink is the interface to the genetic algorithm implementation. In addition to bitstring
 * operators, permutation-based operators are also provided and custom seeding, mutation and
 * crossover operators can be defined by the application, also allowing use of real-valued parameters. Archipelagos
 * of multiple concurrent GAs are supported. Steady-state evolution is also supported. Internal threading allows
 * efficient parallel computation. The library is thread-safe.
 *
 * \section section_fgenpp fgenpp
 *
 * \link group_fgenpp fgenpp \endlink is the C++ wrapper interface to fgen (the genetic algorithm implementation).
 * It offers direct equivalents to the functions in the C version. A sub-class must be used to define the
 * application specific operators.
 *
 * \section section_fpso fpso
 *
 * \link group_fpso fpso \endlink is the interface to the particle swarm optimization implementation, which
 * provides a simple PSO with GBEST or LBEST topology and other customizable parameters. Although it can converge
 * quickly for simple problems, for more complex higher dimensional problems the genetic algorithm generally performs
 * much better.
 *
 * \section section_ffit ffit
 *
 * \link group_ffit ffit \endlink is the interface to straightforward multi-dimensional function minimization. It
 * can be used for model fitting where a set of floating point parameters is evaluated resulting in an error value.
 * It is also convenient for general optimization problems with real-valued parameters. A genetic algorithm or PSO
 * can be used to find the global minimum. Parameters can be distributed within a range with a profile other than linear
 * (for example logarithmic or binomial) to optimize the searching process. When using the genetic algorithm ffit uses
 * 16, 32 or 64-bit bitstring representations mapped to real values, or alternatively normalized real values.
 * 
 * \section section_samples Examples
 *
 * - ones.c is a trivial GA problem using bitstrings.
 * - func.c is a trivial function maximization problem using a GA.
 * - tsp.c implements the Travelling Salesman Problem using a genetic algorithm with permutation
 * representation and specialized mutation and crossover operators.
 * - tsp_archipelago.c does the same with an archipelago of concurrent genetic algorithms. It also supports
 * steady-state evolution.
 * - tsp_archipelago_cpp.cpp is a C++ version of the latter using the C++ wrapper API. 
 * - pi_pso.c is a trivial example of PSO.
 * - ffit_rastrigin.c uses the ffit API to find the minimum of a multi-dimensional Rastrigin function using a
 * genetic algorithm or PSO.
 * - ffit_rosenbrock.c does the same for the two-dimensional Rosenbrock-function.
 * - rastrigin_double.c is an example of real-valued chromosomes with custom operators applied to the
 *  Rastrigin function.
 * - pict is graphical example using ffit which decomposes an image into filled rectangles. It also supports texture
 * compression to ETC1 or DXT1 format. It requires GTK+ 3.
 * - tsp is a graphical example of the Travelling Salesman Problem. It requires GTK+ 3.
 * - gp is an example program that implements linear genetic programming.
 *
 * libfgen can be downloaded from the sourceforge.net page for the project:
 * http://sourceforge.net/projects/libfgen/.
 *
 * \section section_end License
 * fgen is copyright 2012-2013 Harm Hanemaaijer (fgenfb at yahoo.com) and released under the GNU Lesser General
 * Public License.
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
#include "win32_compat.h"

/**
 * Create an fgen population. The random number generator belonging to the population is seeded with
 * zero.
 * 
 * @param population_size Size of the population.
 * @param individual_size_in_bits Size of each individual in bits, must be a multiple of 8.
 * @param data_element_size Size of each data element in bits, preferably a multiple of 8 or 1.
 * @param fgen_generation_callback_func The callback function that is called by fgen after every n generations,
 * where n is defined by the generation callback interval.
 * @param fgen_calculate_fitness_func The fitness evaluation function that calculates the fitness of an individual.
 * It should return a double floating point value where higher values indicate higher fitness. In general fitness values
 * greater or equal to zero are advisable, although negative fitness values are allowed in some cases (when
 * FGEN_SUBTRACT_MIN_FITNESS is set with SUS selection or with any kind of tournament selection). Infinity is allowed
 * (INFINITY in Linux), negative infinity and NaN are not. In Windows positive infinity is defined as DBL_MAX. When a
 * fitness evaluation returns this value, the algorithm will continue but evolution may stop, and in the case of
 * threaded concurrent islands skip generations and generate an early generation_callback.
 * @param fgen_seed_func The seeding operator. This is called for every member of the population when the algorithm
 * starts (see fgen_run()).
 * @param fgen_mutation_func The mutation operator. When this is called the child is allocated and already contains
 * a copy of the parent. This function is called for every individual of the population; logic reflecting the mutation
 * probability should be used within the mutation operator function.
 * @param fgen_crossover_func The crossover operator. The children are allocated but not filled in. This function is
 * called depending on the crossover probability (so in general there is no need to use the crossover probability within
 * the operator).
 * @return The created population.
 */

FgenPopulation *fgen_create(int population_size, int individual_size_in_bits, int data_element_size,
FgenGenerationCallbackFunc fgen_generation_callback_func, FgenCalculateFitnessFunc fgen_calculate_fitness_func,
FgenSeedFunc fgen_seed_func, FgenMutationFunc fgen_mutation_func, FgenCrossoverFunc fgen_crossover_func) {
	FgenPopulation *pop = (FgenPopulation *)malloc(sizeof(FgenPopulation));
	fgen_initialize(pop, population_size, individual_size_in_bits, data_element_size,
		fgen_generation_callback_func, fgen_calculate_fitness_func, fgen_seed_func,
		fgen_mutation_func, fgen_crossover_func);
	return pop;
}

/**
 * Initialize an fgen population. pop must be allocated first. The random number generator belonging to the
 * population is seeded with zero.
 *
 * @param pop The population to initialize. 
 * @param population_size Size of the population.
 * @param individual_size_in_bits Size of each individual in bits, must be a multiple of 8.
 * @param data_element_size Size of each data element in bits, preferably a multiple of 8 or 1.
 * @param fgen_generation_callback_func The callback function that is called by fgen after each generation.
 * @param fgen_calculate_fitness_func The fitness evaluation function that calculates the fitness of an individual.
 * @param fgen_seed_func The seeding operator.
 * @param fgen_mutation_func The mutation operator.
 * @param fgen_crossover_func The crossover operator.
 * @return The created population.
 */

void fgen_initialize(FgenPopulation *pop, int population_size, int individual_size_in_bits,
int data_element_size, FgenGenerationCallbackFunc fgen_generation_callback_func, FgenCalculateFitnessFunc
fgen_calculate_fitness_func, FgenSeedFunc fgen_seed_func, FgenMutationFunc fgen_mutation_func, FgenCrossoverFunc
fgen_crossover_func) {
	pop->size = population_size;
	pop->population_size_shift = fgen_calculate_shift(population_size);
	pop->individual_size_in_bits = individual_size_in_bits;
        pop->individual_size_shift = fgen_calculate_shift(individual_size_in_bits);
	pop->data_element_size = data_element_size;
	pop->data_element_size_shift = fgen_calculate_shift(data_element_size);
        pop->nu_data_elements = pop->individual_size_in_bits / pop->data_element_size;
	pop->fgen_generation_callback_func = fgen_generation_callback_func;
	pop->fgen_calculate_fitness_func = fgen_calculate_fitness_func;
	pop->fgen_seed_func = fgen_seed_func;
	pop->fgen_mutation_func = fgen_mutation_func;
	pop->fgen_crossover_func = fgen_crossover_func;
	pop->cache = NULL;
	pop->rng = fgen_random_create_rng();
	pop->migration_probability = 0;
	pop->migration_probability_float = 0;
	pop->migration_interval = 1;
	fgen_set_tournament_size(pop, DEFAULT_TOURNAMENT_SIZE);
	if (pop->size < 64)
		fgen_set_number_of_elites(pop, 1);
	else
		fgen_set_number_of_elites(pop, pop->size / 64);
	pop->generation_callback_interval = 1;
	pop->fast_mutation_cumulative_chance = NULL;
	pop->ind = NULL;
	pop->initialization_type = FGEN_INITIALIZATION_SEED;
}

/**
 * Set the parameters for an fgen population.
 *
 * @param pop The population. It should be created or initialized before this function is called.
 * @param selection_type The selection type. See fgen_set_selection_type().
 * @param selection_fitness_type The fitness selection type. See fgen_set_selection_fitness_type().
 * @param crossover_probability_float The crossover probability per individual.
 * @param mutation_probability_float The mutation probability (per bit for bitstrings, per individual for permutations).
 * @param macro_mutation_probability_float The macro-mutation probability (per data element).
 */

void fgen_set_parameters(FgenPopulation *pop, int selection_type, int selection_fitness_type,
float crossover_probability_float, float mutation_probability_float, float macro_mutation_probability_float) {
	pop->selection_type = selection_type;
	pop->selection_fitness_type = selection_fitness_type;
	fgen_set_crossover_probability(pop, crossover_probability_float);
	fgen_set_mutation_probability(pop, mutation_probability_float);
	fgen_set_macro_mutation_probability(pop, macro_mutation_probability_float);
}

/**
 * Run the genetic algorithm. First the initial population is created using the seeding function. For every generation
 * selection, crossover and mutation are applied in that order. With a frequency defined by the generation callback
 * interval, the generation callback function is called. Use the generation callback function to keep track of evolution
 * (use fgen_best_individual_of_population() to get the individual with the best fitness) and print out results or draw
 * graphics reflecting the best solution, and use fgen_signal_stop() to let the algorithm terminate.
 *
 * @param pop The population.
 * @param max_generation The maximum generation. When this is reached (or when stop is signalled) the function exits.
 * If it is equal to - 1 the algorithm runs indefinitely and only stops when fgen_signal_stop() is called.
 */

void fgen_run(FgenPopulation *pop, int max_generation) {
	CreateInitialPopulation(pop);

	pop->generation = 0;
	pop->stop_signalled = 0;
        pop->fgen_generation_callback_func(pop, pop->generation);
	if (pop->stop_signalled)
		goto end;

	for (;;) {

		DoSelection(pop);

		DoCrossover(pop);

		DoMutation(pop);

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
 * Free all data structures associated with a population.
 *
 * @param pop The population.
 */

void fgen_destroy(FgenPopulation *pop) {
	if (pop->fast_mutation_cumulative_chance != NULL)
		free(pop->fast_mutation_cumulative_chance);
	if (pop->cache != NULL)
		DestroyCache(pop);
	fgen_random_destroy_rng(pop->rng);
	if (pop->ind != NULL) {
		FreePopulation(pop, pop->ind);
		free(pop->ind);
	}
	free(pop);
}

/**
 * Signal a stop to the algorithm from the generation callback function. After the callback functions returns,
 * the algorithm will terminate.
 *
 * @param pop The population.
 */

void fgen_signal_stop(FgenPopulation *pop) {
	pop->stop_signalled = 1;
}

/**
 * Run the genetic algorithm. Threaded version. Because the amount of parallellism is limited to concurrent fitness
 * calculations, the performance is generally better with the unthreaded version for simple problems.
 * However, with larger populations and more expensive fitness functions the speed-up can be respectable. In linux
 * the number of threads to use is obtained from a sysconf() system call returning the number of online processor cores
 * in the system. On windows the number of threads is set to 4. The fitness cache should not be enabled when using this
 * function.
 *
 * @param pop The population.
 * @param max_generation The maximum generation. When this is reached (or when stop is signalled) the function exits.
 * A value of - 1 means the function only returns when stop is signalled.
 */

void fgen_run_threaded(FgenPopulation *pop, int max_generation) {
//	RandomInformCommonRange(pop->individual_size_in_bits);
//	RandomInformCommonRange(pop->individual_size_in_bits / pop->data_element_size);
//	RandomInformCommonRange(pop->size);
#if defined(__GNUC__) && !defined(_WIN32)
	pop->max_threads = sysconf(_SC_NPROCESSORS_ONLN);
#else
	pop->max_threads = 4;
#endif
//	printf("fgen: number of CPU cores detected = %d.\n", pop->max_threads);

	CreateInitialPopulation(pop);
        CalculatePopulationFitnessThreaded(pop);

	pop->generation = 0;
	pop->stop_signalled = 0;
        pop->fgen_generation_callback_func(pop, pop->generation);
	if (pop->stop_signalled)
		goto end;

	for (;;) {

		DoSelection(pop);

		DoCrossover(pop);

		DoMutation(pop);

		/* The best we can do is to force the fitness of the population to be precalculated */
		/* with concurrency. */
		CalculatePopulationFitnessThreaded(pop);

		pop->generation++;
		if (pop->generation % pop->generation_callback_interval == 0)
			pop->fgen_generation_callback_func(pop, pop->generation);
		if ((max_generation != - 1 && pop->generation >= max_generation) || pop->stop_signalled)
			break;
	}
end : ;
}


/**
 * Run the genetic algorithm on an archipelago of populations. The seeds of the random number generators of the
 * islands after the first one are initialized with random seeds derived from the first island's random number
 * generator.
 *
 * @param nu_pops Number of islands.
 * @param pops An array of populations.
 * @param max_generation The maximum generation. When this is reached (or when stop is signalled) the function exits.
 * A value of - 1 means the function only returns when stop is signalled.
 */

void fgen_run_archipelago(int nu_pops, FgenPopulation **pops, int max_generation) {
	int i;
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
			DoSelection(pops[i]);

			DoCrossover(pops[i]);

			DoMutation(pops[i]);

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
	for (i = 0 ; i < nu_pops; i++)
		CalculatePopulationFitness(pops[i], NULL, NULL);
}

/* Threaded archipelago. */

static void *fgen_create_initial_population_thread(void *threadarg) {
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

static void *fgen_do_multiple_generations_cond_thread(void *threadarg) {
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
			DoSelection(pop);
			/* If there is an individual with infinite fitness, DoSelection will leave it as the */
			/* first individual. Make use of this to skip multiple generations if this is the case. */
			if (pop->ind[0]->fitness == POSITIVE_INFINITY_DOUBLE) {
				pop->generation += nu_generations - i;
				break;
			}
			DoCrossover(pop);
			DoMutation(pop);
			pop->generation++;
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
 * Run the genetic algorithm on an archipelago of populations. Threaded version, uses one thread per island.
 * The seeds of the random number generators of the islands after the first one are initialized with random seeds
 * derived from the first island's random number generator.
 *
 * @param nu_pops Number of islands.
 * @param pops An array of populations.
 * @param max_generation The maximum generation. When this is reached (or when stop is signalled) the function exits.
 * A value of - 1 means the function only returns when stop is signalled.
 */

void fgen_run_archipelago_threaded(int nu_pops, FgenPopulation **pops, int max_generation) {
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
		pthread_create(&thread[i], NULL, fgen_create_initial_population_thread, pops[i]);
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
		pthread_create(&thread[i], NULL, fgen_do_multiple_generations_cond_thread,
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

		if (pops[0]->migration_interval != 0 && pops[0]->generation != 0)
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


/**
 * Return the current island, to be called from the generation callback function.
 *
 * @return The current island.
 */

int fgen_get_island(const FgenPopulation *pop) {
	return pop->island;
}

/**
 * Indicates whether the given population is cached.
 *
 * @return 1 if cached, 0 otherwise.
 */

int fgen_is_cached(const FgenPopulation *pop) {
	return pop->cache != NULL;
}

/**
 * Returns the random number generator belonging to the population.
 */

FgenRNG *fgen_get_rng(const FgenPopulation *pop) {
	return pop->rng;
}

/**
 * Returns the current generation. May be useful in an operator.
 */

int fgen_get_generation(const FgenPopulation *pop) {
	return pop->generation;
}


