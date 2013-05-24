/*
    parameters.c -- parameter setting functions.

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
 * This file contains the basic configurable parameters for the genetic
 * algorithm.
 */

#include <math.h>
#include "fgen.h"
#include "parameters.h"
#include "error.h"
#include "mutation.h"

/** 
 * Set the mutation probability to a floating value from 0 to 1.0. Note that the mutation operator function is called
 * for every individual regardless of the mutation probability. The operator itself uses the probability to apply
 * mutations. The precise meaning of the probability (per bit, or per individual for example) varies between operators.
 */

void fgen_set_mutation_probability(FgenPopulation *pop, float p) {
	pop->mutation_probability_float = p;
	pop->mutation_probability = floor(p * 65536);
	if (pop->fgen_mutation_func == fgen_mutation_per_bit_fast ||
	pop->fgen_mutation_func == fgen_mutation_per_bit_plus_macro_mutation_fast)
		SetupFastMutationCumulativeChanceArray(pop);
}

/**
 * Set the macro-mutation probability to a floating point value from 0 to 1.0. Only has meaning if
 * fgen_mutation_per_bit_plus_macro_mutation() or fgen_mutation_per_bit_plus_macro_mutation_fast() is the mutation
 * operator. May also be used as a parameter for a custom operator.
 */

void fgen_set_macro_mutation_probability(FgenPopulation *pop, float p) {
	pop->macro_mutation_probability_float = p;
	pop->macro_mutation_probability = floor(p * 65536);
}

/**
 * Set the crossover probability to a floating point value forom 0 to 1.0.
 */

void fgen_set_crossover_probability(FgenPopulation *pop, float p) {
	pop->crossover_probability_float = p;
	pop->crossover_probability = floor(p * 256);
}

/**
 * Set the selection fitness type.
 *
 * @param pop The population.
 * @param t One of FGEN_FITNESS_PROPORTIONAL, FGEN_SUBTRACT_MIN_FITNESS (the minimum fitness of the population
 * is subtracted from all fitness values), or FGEN_SUBTRACT_FITNESS_DIV_2 (the minimum fitness divided by two is
 * subtracted). In general FGEN_SUBTRACT_MIN_FITNESS is recommended.
 */

void fgen_set_selection_fitness_type(FgenPopulation *pop, int t) {
	pop->selection_fitness_type = t;
}

/**
 * Set the selection type. In general elitist SUS or elitist tournament selection is recommended.
 *
 * @param pop The population.
 * @param t One of the following:
 * - FGEN_STOCHASTIC: Use simple roulette-wheel selection.
 * - FGEN_SUS: Use stochastic universal sampling.
 * - FGEN_RANK: Use rank-based selection.
 * - FGEN_TOURNAMENT: Use tournament selection. The size is set with fgen_set_tournament_size().
 * - FGEN_ELITIST_SUS, FGEN_ELITIST_RANK, FGEN_ELITIST_TOURNAMENT: One or more of the best individuals
 * of the population are always selected. See fgen_set_number_of_elites().
 * - FGEN_ELITIST_SUS_WITH_EXTINCTION, FGEN_ELITIST_TOURNAMENT_WITH_EXTINCTION: Every 200 generations, the best individual
 * is preserved and the rest replaced with random new individuals.
 * - FGEN_KILL_TOURNAMENT: In steady state evolution using fgen_run_steady_state(), the random individuals selected for
 * replacement are selected with tournament selection of the worst individual.
 */

void fgen_set_selection_type(FgenPopulation *pop, int t) {
	pop->selection_type = t;
}

/**
 * Set the tournament size for tournament selection. The default value is 3.
 */

void fgen_set_tournament_size(FgenPopulation *pop, int n) {
	pop->tournament_size = n;
}

/**
 * Set the data element size in bits. Each individual can be divided into chromosomes (for example 32-bit values
 * that encode an integer or real number). The "per element" crossover operators as well as macro mutation operate
 * with data element boundaries.
 */

void fgen_set_data_element_size(FgenPopulation *pop, int n) {
	pop->data_element_size = n;
}

/**
 * Set the number of elite individuals (with the best fitness) that will always survive selection when elitism is enabled.
 * The default value is the population size divided by 64.
 */

void fgen_set_number_of_elites(FgenPopulation *pop, int n) {
	pop->nu_elites = n;
}

/**
 * Set the permutation size used by the permutation mutation and crossover operators.
 */

void fgen_set_permutation_size(FgenPopulation *pop, int n) {
	pop->permutation_size = n;
}

/**
 * Set the user_data pointer field of the population. This is used by ffit.
 */

void fgen_set_user_data(FgenPopulation *pop, void *data) {
	pop->user_data = data;
}

/**
 * Set the migration probability per individual during a generation when migration is allowed. This probability
 * must be set seperately for each island.
 */

void fgen_set_migration_probability(FgenPopulation *pop, float p) {
	pop->migration_probability_float = p;
	pop->migration_probability = floor(p * 65536);
}

/**
 * Set the migration interval. A value of 1 means migration happens every generation, a value of 50 means
 * migration happens only every 50th generation. The default is 1. A value of 0 means no migration. The
 * migration interval of the first population in an archipelago (pops[0]) is currently used for all islands.
 */

void fgen_set_migration_interval(FgenPopulation *pop, int n) {
	pop->migration_interval = n;
}

/**
 * Set the interval at which the generation callback function is called. A value of 1 means it is called every
 * generation, a value of 50 means it is called only every 50th generation. The default is 1. Higher values are
 * beneficial for threaded execution.
 */

void fgen_set_generation_callback_interval(FgenPopulation *pop, int n) {
	if (n < 1)
		gen_error("fgen_set_generation_callback_interval: interval should be 1 or greater.");
	pop->generation_callback_interval = n;
}

/**
 * Set the initialization type, which select whether to continue with an existing population at the start of fgen_run()
 * and variants, if it exists, or whether to seed a new population. FGEN_INITIALIZATION_SEED or
 * FGEN_INITIALIZATION_CONTINUE are valid values, the default is FGEN_INITIALIZATION_SEED.
 */

void fgen_set_initialization_type(FgenPopulation *pop, int v) {
	pop->initialization_type = v;
}

