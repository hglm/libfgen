/*
    sample_tsp_archipelago_cpp.cpp -- implementation of the Travelling Salesman Problem using an archipelago
				      of genetic algorithms using the C++ wrapper API

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
#include <string.h>
#include <math.h>
#include <limits.h>
#include <fgenpp.h>

/*
 * Implement travelling salesman problem with GA archipelago using C++ wrapper API.
 */


/*
 * Archipelago settings.
 */

#define NU_ISLANDS 8

/*
 * TSP definitions.
 */

#define NU_CITIES 32
#define AREA_SIZE 100.0

// Define a subclass of FgenppPopulation.

class MyPopulation : public FgenppPopulation {
public :
	// Override the following functions:
	void GenerationCallback(int generation);
	double CalculateFitness(const unsigned char *bitstring);
	void Seed(unsigned char *bitstring);
	void Mutate(const unsigned char *parent, unsigned char *child);
	void Crossover(const unsigned char *parent1, const unsigned char *parent2, unsigned char *child1,
		unsigned char *child2);
};

typedef struct {
	double x;
	double y;
} City;

static City city[NU_CITIES];
bool steady_state;

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

double MyPopulation::CalculateFitness(const unsigned char *bitstring) {
	int *perm = (int *)bitstring;
	double dist;
	int i;
	dist = 0;
	for (i = 0; i < NU_CITIES - 1; i++)
		dist += sqrt((city[perm[i + 1]].x - city[perm[i]].x) * (city[perm[i + 1]].x - city[perm[i]].x) +
			(city[perm[i + 1]].y - city[perm[i]].y) * (city[perm[i + 1]].y - city[perm[i]].y));
	return AREA_SIZE * AREA_SIZE / dist;
}


static void PrintIndividual(const FgenppIndividual *ind) {
	int *perm = (int *)ind->GetBitstring();
	int i;
	for (i = 0; i < NU_CITIES; i++)
		printf("[%d]", perm[i]);
	printf(" dist = %lf\n", AREA_SIZE * AREA_SIZE / ind->GetFitness());
}

static void PrintBestIndividual(const FgenppIndividual *best, int generation) {
	printf("generation = %d, best fitness = %lf, solution:\n", generation,
		best->GetFitness());
	PrintIndividual(best);
}

void MyPopulation::GenerationCallback(int generation) {
	if (generation % 100 == 0) {
		printf("Island = %d, ", GetIsland());
		FgenppIndividual *best = BestIndividual();
		PrintBestIndividual(best, generation);
	}
}

void MyPopulation::Seed(unsigned char *bitstring) {
	SeedPermutationRandom(bitstring);
}

void MyPopulation::Mutate(const unsigned char *parent, unsigned char *child) {
	MutatePermutationInvert(parent, child);
}

void MyPopulation::Crossover(const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	CrossoverPermutationOrderBased(parent1, parent2, child1, child2);
}

int main(int argc, char **argv) {
	steady_state = false;
	if (argc > 1)
		if (strcmp(argv[1], "-s") == 0) {
			steady_state = true;
			printf("Running with steady-state evolution.\n");
		}
	InitializeCities();
	FgenppArchipelago *arch = new FgenppArchipelago(NU_ISLANDS);
	MyPopulation *pop;
	for (int i = 0; i < NU_ISLANDS; i++) {
		pop = new MyPopulation;
		pop->Initialize(
			1024,			/* Population size. */
			32 * NU_CITIES,		/* Individual size in bits. */
			32);			/* Data element size. */
		if (steady_state)
			pop->SetParameters(
				FGEN_KILL_TOURNAMENT,		/* Selection type. */
				FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
				0.500,				/* Crossover probability. */
				0.500,				/* Per individual for permutation mutation. */
				0);				/* Macro-mutation probability. */
		else
			pop->SetParameters(
				FGEN_ELITIST_SUS,		/* Selection type. */
				FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
				0.200,				/* Crossover probability. */
				0.200,				/* Per individual for permutation mutation. */
				0);				/* Macro-mutation probability. */
		pop->SetPermutationSize(NU_CITIES);
		pop->SetMigrationProbability(0.001);
		if (steady_state) {
			pop->SetMigrationInterval(10000);
			pop->SetGenerationCallbackInterval(50000);
		}
		else {
			pop->SetMigrationInterval(20);
			pop->SetGenerationCallbackInterval(100);
		}
		// Seed the random number generator of the first population with the timer.
		if (i == 0)
			fgen_random_seed_with_timer(pop->GetRNG());
		arch->AddIsland(pop);
	}

	int max_generations = 1000;
	if (steady_state) {
		max_generations = 500000;
		arch->RunSteadyStateThreaded(max_generations);
	}
	else
		arch->RunThreaded(max_generations);
	printf("All islands: ");
	PrintBestIndividual(arch->BestIndividual(), max_generations);
        delete arch;
    	exit(0);
}

