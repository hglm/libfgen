/*
    ffit.c -- ffit functionality.

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
 * Ffit functionality.
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#ifndef __GNUC__
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <limits.h>
#include <malloc.h>
#include "fgen.h"

static void ffit_fgen_generation_callback(FgenPopulation *pop, int generation);
static void ffit_fgen_archipelago_generation_callback(FgenPopulation *pop, int generation);
static double ffit_fgen_calculate_fitness(const FgenPopulation *pop, const unsigned char *bitstring);
static void ffit_fgen_extract_parameters(const Ffit *fit, const FgenPopulation *pop, const unsigned
char *bitstring, double *param);
// static void ffit_fgen_hill_climb(Ffit *fit, FgenPopulation *pop);
static double ffit_map_parameter(const Ffit *fit, int i, double f);
static void ffit_fgen_real_valued_seed(FgenPopulation *pop, unsigned char *bitstring);
static void ffit_fgen_real_valued_mutation(FgenPopulation *pop, const unsigned char *parent, unsigned char *child);
static void ffit_fgen_real_valued_crossover(FgenPopulation *pop, const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2);
static void ffit_fpso_generation_callback(FpsoPopulation *pop, int generation);
static double ffit_fpso_calculate_error(const FpsoPopulation *pop, const double *params);

/**
 * Create a new fit datastructure.
 *
 * @param nu_params The number of parameters in the solution.
 * @param generation_callback_func The generation callback function that is called after every n generations in the
 * optimization algorithm, where n is the generation callback interval (see ffit_set_generation_callback_interval()).
 * @param calculate_error_func The error evaluation function that returns the error value of the provided
 * parameter solution.
 * @return The created fit.
 */

Ffit *ffit_create(int nu_params, FfitGenerationCallbackFunc generation_callback_func,
FfitCalculateErrorFunc calculate_error_func) {
	Ffit *fit = (Ffit *)malloc(sizeof(Ffit));
	fit->nu_params = nu_params;
	fit->ffit_generation_callback_func = generation_callback_func;
	fit->ffit_calculate_error_func = calculate_error_func;
	fit->range_min = (double *)malloc(sizeof(double) * nu_params);
	fit->range_max = (double *)malloc(sizeof(double) * nu_params);
	fit->mapping = (int *)malloc(sizeof(int) * nu_params);
	fit->threading_level = FFIT_THREADING_DISABLED;
	fit->generation_callback_interval = 1;
	return fit;
}

/**
 * Set the range and mapping of a parameter.
 *
 * @param fit The fit to which the settings apply.
 * @param index The index (starting with 0) of the parameter.
 * @param range_min The lower bound of the range of the parameter.
 * @param range_max The upper bound (exclusive) of the range of the parameter.
 * @param mapping The mathematical mapping of the parameter based on the behaviour of standard functions in the
 * range [0, 1[.
 */

void ffit_set_parameter_range_and_mapping(Ffit *fit, int index, double range_min, double range_max, int mapping) {
	fit->range_min[index] = range_min;
	fit->range_max[index] = range_max;
	fit->mapping[index] = mapping;
}

/**
 * Run a fit using a genetic algorithm with default settings. Note that this function reseeds the random
 * number generator using the system timer so that every run is different. Default settings are:
 * - Population size 1024.
 * - 32 bits per parameter.
 * - Elitist universal stochastic sampling
 * - Uniform crossover per element.
 * - Crossover probability 0.9.
 * - Mutation probability per bit 0.015.
 * - Macro-mutation probability 0.050 per data element.
 *
 * @param fit The fit to run.
 */

void ffit_run_fgen_with_defaults(Ffit *fit) {
	ffit_run_fgen(fit, 1024, 32, FGEN_ELITIST_SUS, fgen_crossover_uniform_per_element, 0.9, 0.015, 0.050);
}

/**
 * Run a fit using a genetic algorithm with specified settings. Note that this function reseeds the random
 * number generator using the system timer so that every run is different.
 * @param fit The fit to run.
 * @param population_size The size of GA population.
 * @param nu_bits_per_param The number of bits used for the parameter representation. Must be 16, 32 or 64.
 * @param selection_type The fgen selection type.
 * @param crossover The crossover operator function.
 * @param crossover_probability_float The crossover rate.
 * @param mutation_per_bit_probability_float The mutation rate per bit.
 * @param macro_mutation_probability_float The macro-mutation probability per data element.
 */

void ffit_run_fgen(Ffit *fit, int population_size, int nu_bits_per_param, int selection_type,
FgenCrossoverFunc crossover, float crossover_probability_float, float mutation_per_bit_probability_float,
float macro_mutation_probability_float) {
	FgenPopulation *pop;
	pop = fgen_create(
		population_size,				/* Population size. */
		fit->nu_params * nu_bits_per_param,		/* Individual size in bits. */
		nu_bits_per_param,				/* Data element size. */
		ffit_fgen_generation_callback,
		ffit_fgen_calculate_fitness,
		fgen_seed_random,
		fgen_mutation_per_bit_plus_macro_mutation_fast,
		crossover);
	fgen_random_seed_with_timer(pop->rng);
	fgen_set_parameters(
		pop,
		selection_type,			/* Selection type. */
		FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
		crossover_probability_float,	/* Crossover probability. */
		mutation_per_bit_probability_float,	/* Mutation probability per bit. */
		macro_mutation_probability_float);	/* Macro-mutation probability. */
	fgen_set_generation_callback_interval(pop, fit->generation_callback_interval);
	fgen_set_user_data(pop, fit);
	fit->population_size = population_size;
	fit->nu_bits_per_param = nu_bits_per_param;
	fit->optimization_type = FFIT_OPTIMIZATION_FGEN;
        fit->stop_signalled = 0;
	fit->model_change_signalled = 0;
	fit->population = pop;
	if (fit->threading_level >= FFIT_THREADING_ENABLED)
		fgen_run_threaded(pop, - 1);
	else
		fgen_run(pop, - 1);
}

/**
 * Run a fit using a genetic algorithm with real-valued chromosomes with specified settings. Note that
 * ffit_run_fgen works very well with real valued parameters, despite using bitstring representation internally.
 * This function uses real-valued parameters (double) internally in the genetic algorithm. Note that this function
 * reseeds the random number generator using the system timer so that every run is different.
 * @param fit The fit to run.
 * @param population_size The size of GA population.
 * @param selection_type The fgen selection type.
 * @param crossover_probability_float The crossover rate.
 * @param mutation_per_bit_probability_float The mutation probability per data element.
 * @param macro_mutation_probability_float The macro-mutation probability per data element.
 */

void ffit_run_fgen_real_valued(Ffit *fit, int population_size, int selection_type,
float crossover_probability_float, float mutation_per_element_probability_float,
float macro_mutation_probability_float) {
	FgenPopulation *pop;
	pop = fgen_create(
		population_size,				/* Population size. */
		fit->nu_params * 64,				/* Individual size in bits. */
		64,						/* Data element size. */
		ffit_fgen_generation_callback,
		ffit_fgen_calculate_fitness,
		ffit_fgen_real_valued_seed,
		ffit_fgen_real_valued_mutation,
		ffit_fgen_real_valued_crossover);
	fgen_random_seed_with_timer(pop->rng);
	fgen_set_parameters(
		pop,
		selection_type,			/* Selection type. */
		FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
		crossover_probability_float,	/* Crossover probability. */
		mutation_per_element_probability_float,	/* Mutation probability per data element. */
		macro_mutation_probability_float);	/* Macro-mutation probability. */
	fgen_set_generation_callback_interval(pop, fit->generation_callback_interval);
	fgen_set_user_data(pop, fit);
	fit->population_size = population_size;
	fit->optimization_type = FFIT_OPTIMIZATION_FGEN_REAL_VALUED;
        fit->stop_signalled = 0;
	fit->model_change_signalled = 0;
	fit->population = pop;
	if (fit->threading_level >= FFIT_THREADING_ENABLED)
		fgen_run_threaded(pop, - 1);
	else
		fgen_run(pop, - 1);
}

#if 0

void ffit_run_fgen_with_hill_climb(Ffit *fit, int population_size, int nu_bits_per_param, int selection_type, FgenCrossoverFunc crossover, float crossover_probability_float, float mutation_per_bit_probability_float, float macro_mutation_probability_float) {
	FgenPopulation *pop;
	pop = fgen_create(
		population_size,				/* Population size. */
		fit->nu_params * nu_bits_per_param,		/* Individual size in bits. */
		nu_bits_per_param,				/* Data element size. */
		ffit_fgen_generation_callback,
		ffit_fgen_calculate_fitness,
		fgen_seed_random,
		fgen_mutation_per_bit_plus_macro_mutation,
		crossover);
	fgen_random_seed_with_timer(pop->rng);
	fgen_set_parameters(
		pop,
		selection_type,			/* Selection type. */
		FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
		crossover_probability_float,	/* Crossover probability. */
		mutation_per_bit_probability_float,	/* Mutation probability per bit. */
		macro_mutation_probability_float);	/* Macro-mutation probability. */
	fgen_set_generation_callback_interval(pop, fit->generation_callback_interval);
	fgen_set_user_data(pop, fit);
	fit->population_size = population_size;
	fit->nu_bits_per_param = nu_bits_per_param;
	fit->optimization_type = FFIT_OPTIMIZATION_FGEN_WITH_HILL_CLIMB;
        fit->stop_signalled = 0;
	fit->population = pop;
	if (fit->threading_level >= FFIT_THREADING_ENABLED)
		fgen_run_threaded(pop, INT_MAX);
	else
		fgen_run(pop, INT_MAX);
}

#endif

/**
 * Run a fit using an archipelago of genetic algorithms with migration. Note that this function reseeds the random
 * number generator using the system timer so that every run is different.
 * @param fit The fit to run.
 * @param nu_islands The number of islands in the archipelago.
 * @param population_size The size of each GA population.
 * @param nu_bits_per_param The number of bits used for the parameter representation. Must be 16, 32 or 64.
 * @param selection_type The fgen selection type.
 * @param crossover The crossover operator function.
 * @param crossover_probability_float The crossover rate.
 * @param mutation_probability_float The mutation rate per bit.
 * @param macro_mutation_probability_float The macro-mutation probability per data element.
 * @param migration_probability_float The migration probability per individual.
 * @param migration_interval The migration interval. Must be greater than 0. A value of 1 means migrations happens
 * every generation, a value of 50 means migrations only happens every 50th generation.
 */

void ffit_run_fgen_archipelago(Ffit *fit, int nu_islands, int population_size,
int nu_bits_per_param, int selection_type, FgenCrossoverFunc crossover, float
crossover_probability_float, float mutation_probability_float, float macro_mutation_probability_float,
float migration_probability_float, int migration_interval) {
	FgenPopulation **pops;
	int i;
	pops = (FgenPopulation **)malloc(sizeof(FgenPopulation *) * nu_islands);
	for (i = 0; i < nu_islands; i++) {
		pops[i] = fgen_create(
			population_size,				/* Population size. */
			fit->nu_params * nu_bits_per_param,		/* Individual size in bits. */
			nu_bits_per_param,				/* Data element size. */
			ffit_fgen_archipelago_generation_callback,
			ffit_fgen_calculate_fitness,
			fgen_seed_random,
			fgen_mutation_per_bit_plus_macro_mutation_fast,
			crossover);
		fgen_set_parameters(
			pops[i],
			selection_type,			/* Selection type. */
			FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
			crossover_probability_float,	/* Crossover probability. */
			mutation_probability_float,	/* Mutation probability per bit. */
			macro_mutation_probability_float);	/* Macro-mutation probability. */
		fgen_set_migration_probability(pops[i], migration_probability_float);
		fgen_set_migration_interval(pops[i], migration_interval);
		fgen_set_generation_callback_interval(pops[i], fit->generation_callback_interval);
		fgen_set_user_data(pops[i], fit);
	}
	fgen_random_seed_with_timer(pops[0]->rng);
	fit->population_size = population_size;
	fit->nu_bits_per_param = nu_bits_per_param;
	fit->optimization_type = FFIT_OPTIMIZATION_FGEN_ARCHIPELAGO;
        fit->stop_signalled = 0;
	fit->population = pops;
	fit->nu_islands = nu_islands;
	fit->best_island_params = (double *)malloc(sizeof(double) * fit->nu_params); 
	if (fit->threading_level >= FFIT_THREADING_ENABLED)
		fgen_run_archipelago_threaded(nu_islands, pops, - 1);
	else
		fgen_run_archipelago(nu_islands, pops, - 1);
}

/**
 * Run a fit using an archipelago of genetic algorithms with real-valued chromosomes internally, with migration.
 * Note that this function reseeds the random number generator using the system timer so that every run is different.
 * @param fit The fit to run.
 * @param nu_islands The number of islands in the archipelago.
 * @param population_size The size of each GA population.
 * @param selection_type The fgen selection type.
 * @param crossover_probability_float The crossover probability.
 * @param mutation_probability_float The mutation probability per data element.
 * @param macro_mutation_probability_float The macro-mutation probability per data element.
 * @param migration_probability_float The migration probability per individual.
 * @param migration_interval The migration interval. Must be greater than 0. A value of 1 means migrations happens
 * every generation, a value of 50 means migrations only happens every 50th generation.
 */

void ffit_run_fgen_real_valued_archipelago(Ffit *fit, int nu_islands, int population_size,
int selection_type, float crossover_probability_float, float mutation_probability_float,
float macro_mutation_probability_float, float migration_probability_float, int migration_interval) {
	FgenPopulation **pops;
	int i;
	pops = (FgenPopulation **)malloc(sizeof(FgenPopulation *) * nu_islands);
	for (i = 0; i < nu_islands; i++) {
		pops[i] = fgen_create(
			population_size,				/* Population size. */
			fit->nu_params * 64,				/* Individual size in bits. */
			64,						/* Data element size. */
			ffit_fgen_archipelago_generation_callback,
			ffit_fgen_calculate_fitness,
			ffit_fgen_real_valued_seed,
			ffit_fgen_real_valued_mutation,
			ffit_fgen_real_valued_crossover);
		fgen_set_parameters(
			pops[i],
			selection_type,			/* Selection type. */
			FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
			crossover_probability_float,	/* Crossover probability. */
			mutation_probability_float,	/* Mutation probability per element. */
			macro_mutation_probability_float);	/* Macro-mutation probability. */
		fgen_set_migration_probability(pops[i], migration_probability_float);
		fgen_set_migration_interval(pops[i], migration_interval);
		fgen_set_generation_callback_interval(pops[i], fit->generation_callback_interval);
		fgen_set_user_data(pops[i], fit);
	}
	fgen_random_seed_with_timer(pops[0]->rng);
	fit->population_size = population_size;
	fit->optimization_type = FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ARCHIPELAGO;
        fit->stop_signalled = 0;
	fit->population = pops;
	fit->nu_islands = nu_islands;
	fit->best_island_params = (double *)malloc(sizeof(double) * fit->nu_params); 
	if (fit->threading_level >= FFIT_THREADING_ENABLED)
		fgen_run_archipelago_threaded(nu_islands, pops, - 1);
	else
		fgen_run_archipelago(nu_islands, pops, - 1);
}

/**
 * Run a fit using particle swarm optimization. Note that this function reseeds the random
 * number generator using the system timer so that every run is different.
 * @param fit The fit to run.
 * @param population_size The size of the swarm population.
 * @param topology The fpso topology type, one of FPSO_TOPOLOGY_LBEST and FPSO_TOPOLOGY_GBEST.
 * @param bounding_strategy The fpso bounding strategy. One of the following:
 * - FPSO_BOUND_NOTHING
 * - FPSO_BOUND_POSITION
 * - FPSO_BOUND_VELOCITY
 * - FPSO_BOUND_POSITION_AND_VELOCITY
 * @param omega The omega value used in the PSO. A default value is defined as FPSO_DEFAULT_OMEGA.
 * @param phi1 The phi1 value used in the PSO. A default value is defined as FPSO_DEFAULT_PHI1.
 * @param phi2 the phi2 value used in the PSO. A default value is defined as FPSO_DEFAULT_PHI2.
 */ 

void ffit_run_fpso(Ffit *fit, int population_size, int topology, int bounding_strategy,
double omega, double phi1, double phi2) {
	FpsoPopulation *pop;
	int i;
	pop = fpso_create(
		population_size,
		fit->nu_params,
		ffit_fpso_generation_callback,
		ffit_fpso_calculate_error);
	fgen_random_seed_with_timer(pop->rng);
	fpso_set_parameters(
		pop,
		topology,
		bounding_strategy,
		omega,
		phi1,
		phi2);
	/* Use normalized parameter values. */
        for (i = 0; i < pop->nu_params; i++)
		fpso_set_parameter_bounds(pop, i, 0, 1.0);
	fpso_set_user_data(pop, fit);
	fit->population_size = population_size;
	fit->optimization_type = FFIT_OPTIMIZATION_FPSO;
        fit->stop_signalled = 0;
	fit->model_change_signalled = 0;
	fit->population = pop;
	fpso_run(pop, - 1);
}

/**
 * Free all data structures associated with a fit.
 * @param fit The fit to destroy.
 */

void ffit_destroy(Ffit *fit) {
	if ((fit->optimization_type & FFIT_OPTIMIZATION_FGEN_ELEMENT)
	&& !(fit->optimization_type & FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT))
		fgen_destroy((FgenPopulation *)fit->population);
	if (fit->optimization_type & FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT) {
                int i;
                for (i = 0; i < fit->nu_islands; i++) {
                    FgenPopulation **pops = (FgenPopulation **)fit->population;
                    fgen_destroy(pops[i]);
                }
		free(fit->best_island_params);
        }
	if (fit->optimization_type & FFIT_OPTIMIZATION_FPSO)
		fpso_destroy((FpsoPopulation *)fit->population);
	free(fit->range_min);
	free(fit->range_max);
	free(fit->mapping);
	free(fit);
}

/**
 * Return the size of the population used in a fit. Only valid after a fit has been started.
 * @param fit The fit.
 */

int ffit_get_population_size(const Ffit *fit) {
	if (!(fit->optimization_type & FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT))
		return fit->population_size;
	/* FFIT_OPTIMIZATION_FGEN_ARCHIPELAGO */
	return fit->population_size * fit->nu_islands;
}

/**
 * Return a pointer to the optimization algorithm population used in a fit.
 * @param fit The fit.
 * @return A pointer to the optimization algorithm population of the following type:
 * - FgenPopulation * if a genetic algorithm is running.
 * - FgenPopulation ** if an archipelago of genetic algorithms is running.
 * - FpsoPopulation * if a particle swarm optimization is running.
 */

void *ffit_get_population(const Ffit *fit) {
	return fit->population;
}

/**
 * Get the parameters of an individual from the optimization algorithm population.
 * @param fit The fit.
 * @param index The index (starting at 0) of the individual in the population.
 * @param params The location where the fetched parameter array of the individual is stored.
 */

void ffit_get_individual_params(const Ffit *fit, int index, double *params) {
	FgenPopulation **pops;
	int island;
	if ((fit->optimization_type & FFIT_OPTIMIZATION_FGEN_ELEMENT)
	&& !(fit->optimization_type & FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT)) {
		if (!(fit->optimization_type & FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ELEMENT)) {
			FgenPopulation *pop = (FgenPopulation *)fit->population;
			ffit_fgen_extract_parameters(fit, pop, pop->ind[index]->bitstring, params);
			return;
		}
		/* Real-valued genetic algorithm. */
		FgenPopulation *pop = (FgenPopulation *)fit->population;
		for (int i = 0; i < fit->nu_params; i++)
			params[i] = ffit_map_parameter(fit, i, ((double *)pop->ind[index]->bitstring)[i]);
		return;
	}
	if (fit->optimization_type & FFIT_OPTIMIZATION_FPSO_ELEMENT) {
		int i;
		FpsoPopulation *pop = (FpsoPopulation *)fit->population;
		for (i = 0; i < pop->nu_params; i++)
			params[i] = ffit_map_parameter(fit, i, pop->ind[index].position[i]);
		return;
	}
	/* FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT set. */
        pops = (FgenPopulation **)fit->population;
	island = index / fit->population_size;
	if (!(fit->optimization_type & FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ELEMENT)) {
		ffit_fgen_extract_parameters(fit, pops[island], pops[island]->ind[index % fit->population_size]->bitstring,
			params);
		return;
	}
	/* Real-valued genetic algorithm. */
	for (int i = 0; i < fit->nu_params; i++)
		params[i] = ffit_map_parameter(fit, i, ((double *)pops[island]->ind[
			index % fit->population_size]->bitstring)[i]); 
	return;
}

double ffit_map_parameter(const Ffit *fit, int i, double f) {
	double f2;
	switch (fit->mapping[i]) {
	case FFIT_MAPPING_LINEAR :
		return (fit->range_max[i] - fit->range_min[i]) * f + fit->range_min[i];
	case FFIT_MAPPING_SQUARE :
		return (fit->range_max[i] - fit->range_min[i]) * f * f + fit->range_min[i];
	case FFIT_MAPPING_CUBE :
		return (fit->range_max[i] - fit->range_min[i]) * f * f * f + fit->range_min[i];
	case FFIT_MAPPING_LOG :
		f2 = log(f * (M_E - 1) + 1);	/* Map to [0, 1[ with log. */
		return (fit->range_max[i] - fit->range_min[i]) * f + fit->range_min[i];
	case FFIT_MAPPING_BINOMIAL_1_TO_5 :
	        f2 = pow(2, (f * 4) + 1) - 2;	// Map to [0, 30[ with a bias to lower values.
		return f2 * (fit->range_max[i] - fit->range_min[i]) / 30.0 + fit->range_min[i];
	case FFIT_MAPPING_BINOMIAL_1_TO_7 :
	        f2 = pow(2, (f * 6) + 1) - 2;	// Map to [0, 126[ with a bias to lower values.
		return f2 * (fit->range_max[i] - fit->range_min[i]) / 126.0 + fit->range_min[i];
	}
	return 0;
}

/**
 * Signal stop to the fitting algorithm. When called from the generation callback function, the algorithm will
 * stop after the callback function returns.
 * @param fit The fit.
 */

void ffit_signal_stop(Ffit *fit) {
	fit->stop_signalled = 1;
}

/**
 * Signal that the calculate error function has changed, so that fitness values in the genetic algorithm are properly
 * invalidated and recalculated.
 */

void ffit_signal_model_change(Ffit *fit) {
	fit->model_change_signalled = 1;
}


/**
 * Set the threading level to either FFIT_THREADING_DISABLED or FFIT_THREADING_ENABLED. Threading allows threaded
 * executing of the genetic algorithm, which is useful with archipelagos or expensive fitness functions with large
 * populations.
 */

void ffit_set_threading(Ffit *fit, int level) {
	fit->threading_level = level;
}

/**
 * Set the interval at which the generation callback function is called. A value of 1 means it is called every
 * generation, a value of 50 means it is called only every 50th generation. The default is 1. Higher values are
 * beneficial for threaded execution when using a genetic algorithm archipelago.
 */

void ffit_set_generation_callback_interval(Ffit *fit, int interval) {
	fit->generation_callback_interval = interval;
}

/*
 * ffit fgen interface.
 */

static void ffit_fgen_extract_parameters(const Ffit *fit, const FgenPopulation *pop,
const unsigned char *bitstring, double *param) {
	int i;
	unsigned char *decoded;
	decoded = (unsigned char *)alloca(pop->individual_size_in_bits / 8);
	fgen_decode_from_gray(pop, bitstring, decoded);
	memcpy(decoded, bitstring, pop->individual_size_in_bits / 8);
	for (i = 0; i < fit->nu_params; i++) {
		double f;
		if (fit->nu_bits_per_param == 32) {
			unsigned int x;
			x = *(unsigned int *)&decoded[i * 4];
			f = (double)x / ((double)pow((double)2, 32));	/* Map to [0, 1[ */
		}
		else
		if (fit->nu_bits_per_param == 16) {
			unsigned short x;
			x = *(unsigned short *)&decoded[i * 2];
			f = (double)x / ((double)pow((double)2, 16));   /* Map to [0, 1[ */
                }
		else {	/* nu_bit_per_param = 64 */
			uint64_t x;
			x = *(uint64_t *)&decoded[i * 8];
			f = (double)x / ((double)pow((double)2, 64));	/* Map to [0, 1[. */
		}
		param[i] = ffit_map_parameter(fit, i, f);
	}
//	free(decoded);
}

/* The following three functions are used by both bitstring and real-valued fgen implementations. */

static void ffit_fgen_generation_callback(FgenPopulation *pop, int generation) { 
	FgenIndividual *best;
	Ffit *fit = (Ffit *)pop->user_data;
	double *best_param;
	double error;
//	if (fit->optimization_type & FFIT_OPTIMIZATION_HILL_CLIMB_ELEMENT)
//		ffit_fgen_hill_climb(fit, pop);
	best = fgen_best_individual_of_population(pop);
	best_param = (double *)alloca(sizeof(double) * fit->nu_params);
	if (!(fit->optimization_type & FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ELEMENT))
		ffit_fgen_extract_parameters(fit, pop, best->bitstring, best_param);
	else
		for (int i = 0; i < fit->nu_params; i++)
			best_param[i] = ffit_map_parameter(fit, i, ((double *)best->bitstring)[i]);
	error = fit->ffit_calculate_error_func(fit, best_param);
	fit->ffit_generation_callback_func(fit, generation, best_param, error);
//	free(best_param);
	if (fit->stop_signalled)
		fgen_signal_stop(pop);
	if (fit->model_change_signalled) {
		fgen_invalidate_population_fitness(pop);
		fit->model_change_signalled = 0;
	}
}

static void ffit_fgen_archipelago_generation_callback(FgenPopulation *pop, int generation) { 
	FgenIndividual *best;
	Ffit *fit = (Ffit *)pop->user_data;
	double *best_param;
	double error;
	best = fgen_best_individual_of_population(pop);
	best_param = (double *)alloca(sizeof(double) * fit->nu_params);
	if (!(fit->optimization_type & FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ELEMENT))
		ffit_fgen_extract_parameters(fit, pop, best->bitstring, best_param);
	else
		for (int i = 0; i < fit->nu_params; i++)
			best_param[i] = ffit_map_parameter(fit, i, ((double *)best->bitstring)[i]);
	error = fit->ffit_calculate_error_func(fit, best_param);
	if (pop->island == 0) {
		fit->best_island_error = error;
		memcpy(fit->best_island_params, best_param, sizeof(double) * fit->nu_params);
	}
	else
		if (error < fit->best_island_error) {
			fit->best_island_error = error;
			memcpy(fit->best_island_params, best_param, sizeof(double) * fit->nu_params);
		}
	if (pop->island == fit->nu_islands - 1) {
		fit->ffit_generation_callback_func(fit, generation, fit->best_island_params,
			fit->best_island_error);
		if (fit->stop_signalled)
			fgen_signal_stop(pop);
		if (fit->model_change_signalled) {
			FgenPopulation **pops = (FgenPopulation **)fit->population;
			for (int i = 0; i < fit->nu_islands; i++)
				fgen_invalidate_population_fitness(pops[i]);
			fit->model_change_signalled = 0;
		}
	}
//	free(best_param);
}

static double ffit_fgen_calculate_fitness(const FgenPopulation *pop, const unsigned char *bitstring) {
	Ffit *fit = (Ffit *)pop->user_data;;
	double *param = (double *)alloca(sizeof(double) * fit->nu_params);
        double error;
	if (!(fit->optimization_type & FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ELEMENT))
		ffit_fgen_extract_parameters(fit, pop, bitstring, param);
	else
		for (int i = 0; i < fit->nu_params; i++)
			param[i] = ffit_map_parameter(fit, i, ((double *)bitstring)[i]);
	error = fit->ffit_calculate_error_func(fit, param);
//	free(param);
	return 1000000 / error;
}

#if 0

/*
 * ffit hill climbing in normalized parameter space (every parameter has the range [0, 1[).
 */

static void ffit_fgen_extract_normalized_parameters(Ffit *fit, FgenPopulation *pop, unsigned char *bitstring,
double *param) {
	int i;
	unsigned char *decoded;
	decoded = (unsigned char *)malloc(pop->individual_size_in_bits / 8);
	fgen_decode_from_gray(pop, bitstring, decoded);
	for (i = 0; i < fit->nu_params; i++) {
		if (fit->nu_bits_per_param == 32) {
			unsigned int x;
			x = *(unsigned int *)&decoded[i * 4];
			param[i] = (double)x / ((double)pow(2, 32));	/* Map to [0, 1[ */
		}
		else
		if (fit->nu_bits_per_param == 16) {
			unsigned short x;
			x = *(unsigned short *)&decoded[i * 2];
			param[i] = (double)x / ((double)pow(2, 16));   /* Map to [0, 1[ */
                }
		else {	/* nu_bit_per_param = 64 */
			uint64_t x;
			x = *(uint64_t *)&decoded[i * 8];
			param[i] = (double)x / ((double)pow(2, 64));	/* Map to [0, 1[. */
		}
	}
	free(decoded);
}

static void ffit_fgen_encode_normalized_parameters(Ffit *fit, FgenPopulation *pop, double *param, unsigned char *bitstring) {
	int i;
	unsigned char *decoded;
	decoded = (unsigned char *)malloc(pop->individual_size_in_bits / 8);
	for (i = 0; i < fit->nu_params; i++) {
		if (fit->nu_bits_per_param == 32) {
			unsigned int x;
			x = param[i] * ((double)pow(2, 32));
			*(unsigned int *)&decoded[i * 4] = x;
		}
		else
		if (fit->nu_bits_per_param == 16) {
			unsigned short x;
			x = param[i] * ((double)pow(2, 16));
			*(unsigned short *)&decoded[i * 2] = x;
                }
		else {	/* nu_bit_per_param = 64 */
			uint64_t x;
			x = param[i] * ((double)pow(2, 64));
			*(uint64_t *)&decoded[i * 8] = x;
		}
	}
	fgen_encode_to_gray(pop, decoded, bitstring);
	free(decoded);
}

/* Calculate error based on given normalized parameters. */

static double ffit_hill_climb_calculate_error(Ffit *fit, double *normalized_param) {
	double *param;
	int i;
	double error;
	param = (double *)malloc(sizeof(double) * fit->nu_params);
	for (i = 0; i < fit->nu_params; i++)
		param[i] = ffit_map_parameter(fit, i, normalized_param[i]);
	error = fit->ffit_calculate_error_func(fit, param); 
	free(param);
	return error;
}

/* Do a hill-climb with given normalized parameters. */

static void ffit_do_hill_climb(Ffit *fit, FgenPopulation *pop, double *normalized_param) {
	int i;
	double max_step_size;
	double error;
	double *param;
	/* Try five random steps in the parameter space. */
	max_step_size = 0.0001;	/* Very arbitrary value. */
	error = ffit_hill_climb_calculate_error(fit, normalized_param);
	param = (double *)malloc(sizeof(double) * fit->nu_params);
	for (i = 0; i < 5; i++) {
		int param_index;
		double step;
		double new_error;
		/* Make a working copy of the old parameters. */
		memcpy(param, normalized_param, sizeof(double) * fit->nu_params);
		/* Choose a random parameter. */
		param_index = RandomInt(pop, fit->nu_params);
		/* Choose a random small step, positive or negative. */
		step = fgen_random_d(max_step_size * 2) - max_step_size;
		param[param_index] += step;
		if (param[param_index] < 0)
			param[param_index] = 0;
		if (param[param_index] > 1.0)
			param[param_index] = 0.99999999999999;
		new_error = ffit_hill_climb_calculate_error(fit, param);
		if (new_error < error) {
			/* The new parameters are better, replace the old parameters. */
			normalized_param[param_index] = param[param_index];
			error = new_error;
// printf("+");
		}
        }
	free(param);
}

static void print_params(Ffit *fit, double *param) {
	int i;
	printf("params: ");
	for (i = 0; i < fit->nu_params; i++)
		printf("%lf ", param[i]);
	printf("\n");
}

static void ffit_fgen_hill_climb(Ffit *fit, FgenPopulation *pop) {
	int i;
	double *normalized_param;
	normalized_param = (double *)malloc(sizeof(double) * fit->nu_params);
	for (i = 0; i < pop->size; i++) {
		ffit_fgen_extract_normalized_parameters(fit, pop, pop->ind[i]->bitstring, normalized_param);
// printf("Before hill-climb: ");
// print_params(fit, normalized_param);
		ffit_do_hill_climb(fit, pop, normalized_param);
// printf("After hill-climb: ");
// print_params(fit, normalized_param);
		ffit_fgen_encode_normalized_parameters(fit, pop, normalized_param, pop->ind[i]->bitstring);
	}
	free(normalized_param);
}

#endif

/*
 * ffit fgen real-valued interface.
 */

/* The parameters are normalized to the range [0, 1]. */

#define DOMAIN_MIN 0
#define DOMAIN_MAX 1.0

/* Custom seeding function that assignd parameters with random double floating point values within the range. */

static void ffit_fgen_real_valued_seed(FgenPopulation *pop, unsigned char *bitstring) {
	double *params = (double *)bitstring;
	int i;
	FgenRNG *rng;
	Ffit *fit = (Ffit *)pop->user_data;
	rng = fgen_get_rng(pop);
	for (i = 0; i < fit->nu_params; i++)
		params[i] = fgen_random_from_range_d(rng, DOMAIN_MIN, DOMAIN_MAX); 
}

/* Custom mutation function that assigns random values within the range to parameters with a mutation */
/* probability. When the function is called child already contains a copy of parent */

static void ffit_fgen_real_valued_mutation_large(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	double *params_child = (double *)child;
	int i;
	Ffit *fit = (Ffit *)pop->user_data; 
	/* Mutate each parameter with the given probability. */
	for (i = 0; i < fit->nu_params; i++)
		if (fgen_random_f(pop->rng, 1.0) < pop->macro_mutation_probability_float)
			params_child[i] = fgen_random_from_range_d(pop->rng, DOMAIN_MIN, DOMAIN_MAX);
}

/* Custom mutation function that tweaks the parameter by a maximum of 1/20th of the range. */

static void ffit_fgen_real_valued_mutation_small(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	double *params_parent = (double *)parent;
	double *params_child = (double *)child;
	int i;
	Ffit *fit = (Ffit *)pop->user_data;
	/* Mutate each parameter with the given probability. */
	for (i = 0; i < fit->nu_params; i++)
		if (fgen_random_f(pop->rng, 1.0) < pop->mutation_probability_float) {
			double tweak = fgen_random_from_range_d(pop->rng, - (DOMAIN_MAX - DOMAIN_MIN) / 20,
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

static void ffit_fgen_real_valued_mutation(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	ffit_fgen_real_valued_mutation_small(pop, parent, child);
	ffit_fgen_real_valued_mutation_large(pop, parent, child);
}

/* Custom real-valued crossover function. Intermediate recombination. */

static void ffit_fgen_real_valued_crossover(FgenPopulation *pop, const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	double *params_parent1 = (double *)parent1;
	double *params_parent2 = (double *)parent2;
	double *params_child1 = (double *)child1;
	double *params_child2 = (double *)child2;
	int i;
	Ffit *fit = (Ffit *)pop->user_data;
	for (i = 0; i < fit->nu_params; i++) {
		double alpha;
		/* First offspring. */
		alpha = fgen_random_from_range_d(pop->rng, - 0.25, 1.25);
		params_child1[i] = params_parent1[i] + alpha * (params_parent2[i] - params_parent1[i]);
		/* Clamp to range. */
		if (params_child1[i] < DOMAIN_MIN)
			params_child1[i] = DOMAIN_MIN;
		if (params_child1[i] > DOMAIN_MAX)
			params_child1[i] = DOMAIN_MAX;
		/* Second offspring. */
		alpha = fgen_random_from_range_d(pop->rng, - 0.25, 1.25);
		params_child2[i] = params_parent2[i] + alpha * (params_parent1[i] - params_parent2[i]);
		/* Clamp to range. */
		if (params_child2[i] < DOMAIN_MIN)
			params_child2[i] = DOMAIN_MIN;
		if (params_child2[i] > DOMAIN_MAX)
			params_child2[i] = DOMAIN_MAX;
	}
}

/*
 * ffit fpso interface.
 */

static void ffit_fpso_generation_callback(FpsoPopulation *pop, int generation) {
	Ffit *fit = (Ffit *)pop->user_data;
	if (generation % fit->generation_callback_interval != 0)
		return;
	double *best_param;
	double *best_normalized_param;
	int i;
	best_normalized_param = fpso_get_best_known_position(pop);
	best_param = (double *)alloca(sizeof(double) * pop->nu_params);
	/* Convert from normalized parameters to real parameters. */
	for (i = 0; i < pop->nu_params; i++)
		best_param[i] = ffit_map_parameter(fit, i, best_normalized_param[i]); 
	fit->ffit_generation_callback_func(fit, generation, best_param, fpso_get_best_known_error(pop));
//	free(best_param);
	if (fit->stop_signalled)
		fpso_signal_stop(pop);
}

static double ffit_fpso_calculate_error(const FpsoPopulation *pop, const double *params_in) {
	Ffit *fit = (Ffit *)pop->user_data;
	double *param = (double *)alloca(sizeof(double) * pop->nu_params);
        double error;
	int i;
	/* Convert from normalized parameters to real parameters. */
	for (i = 0; i < pop->nu_params; i++)
		param[i] = ffit_map_parameter(fit, i, params_in[i]); 
	error = fit->ffit_calculate_error_func(fit, param);
//	free(param);
	return error;
}

