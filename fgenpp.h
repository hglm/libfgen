/*
    fgenpp.h -- C++ wrapper API header file.

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
 * C++ wrapper API for the libfgen library.
 */

#ifndef __FGENPP_H__
#define __FGENPP_H__

#include <fgen.h>

/** \defgroup group_fgenpp fgenpp C++ wrapper API for the fgen (Genetic Algorithm) module
 * @{
 */

class FGEN_API FgenppIndividual : public FgenIndividual {
public :
	double GetFitness() const {
		return fitness;
	}
	unsigned char *GetBitstring() const {
		return bitstring;
	}
};

class FGEN_API FgenppPopulation {
private :
	FgenPopulation *pop;
public :
	FGEN_API void Initialize(int population_size, int individual_size_in_bits, int data_element_size);
	FGEN_API void Destroy();
	FGEN_API virtual void GenerationCallback(int generation);
	FGEN_API virtual double CalculateFitness(const unsigned char *bitstring);
	FGEN_API virtual void Seed(unsigned char *bitstring);
	FGEN_API virtual void Mutate(const unsigned char *parent, unsigned char *child);
	FGEN_API virtual void Crossover(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void SeedBitstringRandom(unsigned char *bitstring);
	FGEN_API void MutatePerBitPlusMacroMutation(const unsigned char *parent, unsigned char *child);
	FGEN_API void MutatePerBit(const unsigned char *parent, unsigned char *child);
	FGEN_API void CrossoverOnePointPerBit(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void CrossoverTwoPointPerBit(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void CrossoverOnePointPerElement(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void CrossoverTwoPointPerElement(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void CrossoverUniformPerBit(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void CrossoverUniformPerElement(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void CrossoverNoop(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void SeedPermutationRandom(unsigned char *bitstring);
	FGEN_API void MutatePermutationSwap(const unsigned char *parent, unsigned char *child);
	FGEN_API void MutatePermutationInsert(const unsigned char *parent, unsigned char *child);
	FGEN_API void MutatePermutationInvert(const unsigned char *parent, unsigned char *child);
	FGEN_API void CrossoverPermutationOrderBased(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void CrossoverPermutationPositionBased(const unsigned char *parent1, const unsigned char *parent2,
		unsigned char *child1, unsigned char *child2);
	FGEN_API void SetParameters(
		int selection_type,
		int selection_fitness_type,
		float crossover_probability_float,
		float mutation_per_bit_probability_float,
		float macro_mutation_probability_float);
	FGEN_API void Run(int max_generation);
	FGEN_API void RunThreaded(int max_generation);
	FGEN_API void RunSteadyState(int max_generation);
	FGEN_API void SetMutationProbability(float prob);
	FGEN_API void SetMacroMutationProbability(float prob);
	FGEN_API void SetCrossoverProbability(float prob);
	FGEN_API void SetSelectionFitnessType(int type);
	FGEN_API void SetSelectionType(int type);
	FGEN_API void SetTournamentSize(int size);
	FGEN_API void SetDataElementSize(int size);
	FGEN_API void SetNumberOfElites(int n);
	FGEN_API void SetPermutationSize(int size);
	FGEN_API void SetMigrationInterval(int interval);
	FGEN_API void SetMigrationProbability(float prob);
	FGEN_API void SetGenerationCallbackInterval(int interval);
	FGEN_API FgenppIndividual *BestIndividual();
	FGEN_API FgenppIndividual *WorstIndividual();
	FGEN_API void UpdateFitness();
	FGEN_API void SignalStop();
	FGEN_API FgenRNG *GetRNG();
	FGEN_API int GetIndividualSizeInBytes();
	FGEN_API int GetIsland();
	friend class FgenppArchipelago;
};

class FGEN_API FgenppArchipelago {
private :
	FgenPopulation **pops;
	int max_size;
	int size;
public :
	FGEN_API FgenppArchipelago(int max_size);
	FGEN_API ~FgenppArchipelago();
	FGEN_API void AddIsland(FgenppPopulation *pop);
	FGEN_API void Run(int max_generation);
	FGEN_API void RunThreaded(int max_generation);
	FGEN_API void RunSteadyState(int max_generation);
	FGEN_API void RunSteadyStateThreaded(int max_generation);
	FGEN_API FgenppIndividual *BestIndividual();
	FGEN_API FgenppPopulation *GetPopulation(int index);
};

/** @} */

#endif

