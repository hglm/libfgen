/*
    fgenpp.cpp -- C++ wrapper API.

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
#include "fgenpp.h"


// Function hooks from the C library.

void fgenpp_generation_callback_func(FgenPopulation *pop, const int generation) {
	((FgenppPopulation *)pop->user_data)->GenerationCallback(generation);
}

double fgenpp_calculate_fitness_func(const FgenPopulation *pop, const unsigned char *bitstring) {
	return ((FgenppPopulation *)pop->user_data)->CalculateFitness(bitstring);
}

void fgenpp_seed_func(FgenPopulation *pop, unsigned char *bitstring) {
	((FgenppPopulation *)pop->user_data)->Seed(bitstring);
}

void fgenpp_mutation_func(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child) {
	((FgenppPopulation *)pop->user_data)->Mutate(parent, child);
}

void fgenpp_crossover_func(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	((FgenppPopulation *)pop->user_data)->Crossover(parent1, parent1, child1, child2);
}

// Default operator functions.

/**
 * The generation callback method. Usually overridden in a derived sub-class.
 */

void FgenppPopulation::GenerationCallback(int generation) {
}


/**
 * The calculate error method. The default one does nothing. It must be overridden by a derived
 * sub-class.
 */

double FgenppPopulation::CalculateFitness(const unsigned char *) {
}

/**
 * The seed method. The default is the random bitstring seeding operator SeedBitstringRandom. Can be overridden
 * by a derived sub-class.
 */

void FgenppPopulation::Seed(unsigned char *bitstring) {
	SeedBitstringRandom(bitstring);
}

/**
 * The mutate method. The default is mutation per bit plus macro mutation for bitstring
 * (MutatePerBitPlusMacroMutation). Can be overridden by a derived sub-class.
 */

void FgenppPopulation::Mutate(const unsigned char *parent, unsigned char *child) {
	MutatePerBitPlusMacroMutation(parent, child);
}

/**
 * The crossover method. The default is two-point bitwise crossover for bitstrings (CrossoverTwoPointPerBit).
 * Can be overridden in a derived sub-class. 
 */

void FgenppPopulation::Crossover(const unsigned char *parent1, const unsigned char *parent2, unsigned char *child1,
unsigned char *child2) {
	CrossoverTwoPointPerBit(parent1, parent2, child1, child2);
}

// Supplied operator functions.

void FgenppPopulation::SeedBitstringRandom(unsigned char *bitstring) {
	fgen_seed_random(pop, bitstring);
}

void FgenppPopulation::SeedPermutationRandom(unsigned char *bitstring) {
	fgen_seed_permutation_random(pop, bitstring);
}

void FgenppPopulation::MutatePerBitPlusMacroMutation(const unsigned char *parent, unsigned char *child) {
	fgen_mutation_per_bit_plus_macro_mutation(pop, parent, child);
}

void FgenppPopulation::MutatePerBit(const unsigned char *parent, unsigned char *child) {
	fgen_mutation_per_bit(pop, parent, child);
}

void FgenppPopulation::MutatePermutationSwap(const unsigned char *parent, unsigned char *child) {
	fgen_mutation_permutation_swap(pop, parent, child);
}

void FgenppPopulation::MutatePermutationInsert(const unsigned char *parent, unsigned char *child) {
	fgen_mutation_permutation_insert(pop, parent, child);
}

void FgenppPopulation::MutatePermutationInvert(const unsigned char *parent, unsigned char *child) {
	fgen_mutation_permutation_invert(pop, parent, child);
}

void FgenppPopulation::CrossoverOnePointPerBit(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_one_point_per_bit(pop, parent1, parent2, child1, child2);
}

void FgenppPopulation::CrossoverTwoPointPerBit(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_two_point_per_bit(pop, parent1, parent2, child1, child2);
}

void FgenppPopulation::CrossoverOnePointPerElement(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_one_point_per_element(pop, parent1, parent2, child1, child2);
}

void FgenppPopulation::CrossoverTwoPointPerElement(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_two_point_per_element(pop, parent1, parent2, child1, child2);
}

void FgenppPopulation::CrossoverUniformPerBit(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_uniform_per_bit(pop, parent1, parent2, child1, child2);
}

void FgenppPopulation::CrossoverUniformPerElement(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_uniform_per_element(pop, parent1, parent2, child1, child2);
}

void FgenppPopulation::CrossoverPermutationOrderBased(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_permutation_order_based(pop, parent1, parent2, child1, child2);
}

void FgenppPopulation::CrossoverPermutationPositionBased(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_permutation_position_based(pop, parent1, parent2, child1, child2);
}

void FgenppPopulation::CrossoverNoop(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	fgen_crossover_noop(pop, parent1, parent2, child1, child2);
}

// FgenppPopulation

/**
 * Initialize an FgenppPopulation. This does the same things as fgen_create(). The FgenppPopulation must be
 * instantiated.
 */

void FgenppPopulation::Initialize(int population_size, int individual_size_in_bits, int data_element_size) {
	pop = new FgenPopulation;
	fgen_initialize(
		pop,
		population_size,
		individual_size_in_bits,
		data_element_size,
		fgenpp_generation_callback_func,
		fgenpp_calculate_fitness_func,
		fgenpp_seed_func,
		fgenpp_mutation_func,
		fgenpp_crossover_func);
	// Trick to allow access to the class given an FgenPopulation.
	pop->user_data = this;
}

/**
 * Destroy data structures associated with a population. The FgenppPopulation class itself is not deallocated.
 */

void FgenppPopulation::Destroy() {
	fgen_destroy(pop);
}

void FgenppPopulation::SetParameters(
int selection_type,
int selection_fitness_type,
float crossover_probability_float,
float mutation_per_bit_probability_float,
float macro_mutation_probability_float) {
	fgen_set_parameters(pop,
		selection_type,
		selection_fitness_type,
		crossover_probability_float,
		mutation_per_bit_probability_float,
		macro_mutation_probability_float);

}

void FgenppPopulation::Run(int max_generation) {
	fgen_run(pop, max_generation);
}

void FgenppPopulation::RunThreaded(int max_generation) {
	fgen_run_threaded(pop, max_generation);
}

void FgenppPopulation::RunSteadyState(int max_generation) {
	fgen_run_steady_state(pop, max_generation);
}

void FgenppPopulation::SetMutationProbability(float prob) {
	fgen_set_mutation_probability(pop, prob);
}

void FgenppPopulation::SetMacroMutationProbability(float prob) {
	fgen_set_macro_mutation_probability(pop, prob);
}

void FgenppPopulation::SetCrossoverProbability(float prob) {
	fgen_set_crossover_probability(pop, prob);
}

void FgenppPopulation::SetSelectionFitnessType(int type) {
	fgen_set_selection_fitness_type(pop, type);
}

void FgenppPopulation::SetSelectionType(int type) {
	fgen_set_selection_type(pop, type);
}

void FgenppPopulation::SetTournamentSize(int size) {
	fgen_set_tournament_size(pop, size);
}

void FgenppPopulation::SetDataElementSize(int size) {
	fgen_set_data_element_size(pop, size);
}

void FgenppPopulation::SetNumberOfElites(int n) {
	fgen_set_number_of_elites(pop, n);
}

void FgenppPopulation::SetPermutationSize(int size) {
	fgen_set_permutation_size(pop, size);
}

void FgenppPopulation::SetMigrationInterval(int interval) {
	fgen_set_migration_interval(pop, interval);
}

void FgenppPopulation::SetMigrationProbability(float prob) {
	fgen_set_migration_probability(pop, prob);
}

void FgenppPopulation::SetGenerationCallbackInterval(int interval) {
	fgen_set_generation_callback_interval(pop, interval);
}

FgenppIndividual *FgenppPopulation::BestIndividual() {
	return (FgenppIndividual *)fgen_best_individual_of_population(pop);
}

FgenppIndividual *FgenppPopulation::WorstIndividual() {
	return (FgenppIndividual *)fgen_worst_individual_of_population(pop);
}

void FgenppPopulation::UpdateFitness() {
	fgen_update_population_fitness(pop);
}

void FgenppPopulation::SignalStop() {
	fgen_signal_stop(pop);
}

FgenRNG *FgenppPopulation::GetRNG() {
	return fgen_get_rng(pop);
}

int FgenppPopulation::GetIsland() {
	return fgen_get_island(pop);
}

int FgenppPopulation::GetIndividualSizeInBytes() {
	return fgen_individual_size_in_bytes(pop);
}

// FgenppArchipelago

/**
 * Constructor for an archipelago of genetic algorithms of the specified maximum size.
 */

FgenppArchipelago::FgenppArchipelago(int _max_size) {
	max_size = _max_size;
	size = 0;
	pops = new FgenPopulation *[max_size];
}

/**
 * Deconstructor for an archipelago that frees all data structures associated with the archipelago.
 */

FgenppArchipelago::~FgenppArchipelago() {
	for (int i = 0; i < size; i++)
		fgen_destroy(pops[i]);
	delete pops;
}

/**
 * Add an island of type FgenppPopulation * to an archipelago.
 */

void FgenppArchipelago::AddIsland(FgenppPopulation *pop) {
	pops[size] = pop->pop;
	size++;
}

/**
 * Run the genetic algorithm on the archipelago. 
 */

void FgenppArchipelago::Run(int max_generation) {
	fgen_run_archipelago(size, pops, max_generation);
}

void FgenppArchipelago::RunThreaded(int max_generation) {
	fgen_run_archipelago_threaded(size, pops, max_generation);
}

void FgenppArchipelago::RunSteadyState(int max_generation) {
	fgen_run_steady_state_archipelago(size, pops, max_generation);
}

void FgenppArchipelago::RunSteadyStateThreaded(int max_generation) {
	fgen_run_steady_state_archipelago_threaded(size, pops, max_generation);
}

/**
 * Returns the best individual of the archipelago (all islands).
 */

FgenppIndividual *FgenppArchipelago::BestIndividual() {
	return (FgenppIndividual *)fgen_best_individual_of_archipelago(size, pops);
}

/**
 * Returns the population (of type FgenppPopulation) corresponding to the given index.
 */

FgenppPopulation *FgenppArchipelago::GetPopulation(int index) {
	return (FgenppPopulation *)(pops[index]->user_data);
}

