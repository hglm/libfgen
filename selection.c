/*
    selection.c -- Selection stage of the genetic algorithm.

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
#include <math.h>
#include <limits.h>
#include <float.h>
#include "fgen.h"
#include "parameters.h"
#include "population.h"
#include "random.h"
#include "error.h"
#include "win32_compat.h"

static int CompareFitness(const void *p1, const void *p2) {
	FgenIndividual *q1 = *(FgenIndividual **)p1;
	FgenIndividual *q2 = *(FgenIndividual **)p2;
	if (q1->fitness > q2->fitness)
            return -1;
        if (q1->fitness < q2->fitness)
            return 1;
        return 0;
}

/* In elitist selection, select the best individuals beforehand. They will also be allowed */
/* to be selected extra times based on fitness. */

static int DoElitism(FgenPopulation *pop, int *random_order, FgenIndividual **ind1, FgenIndividual **ind2) {
        if (!(pop->selection_type & FGEN_ELITIST_ELEMENT))
		return 0;
	for (int i = 0; i < pop->size; i++)
		ind1[i]->is_elite = 0;
	// During ranked-based selection, the population was already sorted.
	if (!((pop->selection_type & FGEN_STOCHASTIC_TYPE_MASK) == FGEN_RANK))
		// Sort the population on fitness.
        	qsort(ind1, pop->size, sizeof(ind1[0]), CompareFitness);
       	for (int i = 0; i < pop->nu_elites; i++) {
		// Have to mostly create new individuals to avoid problems with propagating
		// the is_elite == true characteristic to other individuals that share
		// the same data structure.
		// However, if the current elite is identical to the previous one, don't create a new
		// individual.
		if (i > 0)
			if (ind1[i] == ind1[i - 1]) {
				ind2[random_order[i]] = ind2[random_order[i - 1]];
				ind2[random_order[i - 1]]->refcount++;
				continue;
			}
		ind2[random_order[i]] = NewIndividual(pop);
		CopyIndividualBitstring(pop, ind1[i], ind2[random_order[i]]);
		ind2[random_order[i]]->fitness = ind1[i]->fitness;
		ind2[random_order[i]]->fitness_is_valid = 1;
		ind2[random_order[i]]->is_elite = 1;
	}
	return pop->nu_elites;
}


/* Stochastic selection with replacement. */

static void StochasticSelection(FgenPopulation *pop, FgenIndividual **ind1, FgenIndividual **ind2, double sum,
double fitness_cut) {
	int i;
	int *random_order;
	/*
	 * Inefficient roulette wheel selection.
	 * This traverses the whole roulette wheel for each individual
	 * selected. A more efficient way would be to pick individuals
	 * corresponding to equally distant markers (pop->size in
	 * number) along the roulette wheel (implemented in universal
	 * stochastic sampling).
	 */

	/* Calculate random order for target population. */
	random_order = (int *)malloc(sizeof(int) * pop->size);
	CalculateRandomOrder(pop->rng, random_order, pop->size);

	/* In elitist selection, select the best individuals beforehand. They will also be allowed */
	/* to be selected extra times based on fitness. */
	int nu_elites = DoElitism(pop, random_order, ind1, ind2);

        double sum2 = 0;
        for (i = 0; i < pop->size; i++)
            sum2 += ind1[i]->fitness - fitness_cut;
	for (i = nu_elites; i < pop->size; i++) {
		int j;
		double draw, csum;
		/* Pick a spot in the cumulative sum. */
		draw = fgen_random_d(pop->rng, sum2);
		/* Find the individual. */
		csum = 0;
		for (j = 0; j < pop->size; j++) {
			csum += ind1[j]->fitness - fitness_cut;
			if (draw < csum)
				break;
		}
		/* Avoid rounding errors. */
		if (j == pop->size) {
			j = pop->size - 1;
		}
		ind2[random_order[i]] = ind1[j];
		ind1[j]->refcount++;
	}

	free(random_order);
}

/* Stochastic Universal Sampling. */

static void StochasticUniversalSamplingSelection(FgenPopulation *pop, FgenIndividual **ind1, FgenIndividual **ind2,
double sum, double fitness_cut ) {
	int i, j;
        double csum, draw, step;
	double *pop_csum;
	int *random_order;
	/*
	 * The individuals corresponding to equally distant markers
	 * (pop->size in number) along the roulette wheel are picked.
	 * This ensures that the each individual gets at least the integer
	 * part of its expected number.
	 *
	 */

	/* Calculate random order for target population. */
	random_order = (int *)malloc(sizeof(int) * pop->size);
	CalculateRandomOrder(pop->rng, random_order, pop->size);

	/* In elitist selection, select the best individuals beforehand. They will also be allowed */
	/* to be selected extra times based on fitness. */
	int nu_elites = DoElitism(pop, random_order, ind1, ind2);

	/* Calculate array of cumulative sums. */
	pop_csum = (double *)malloc(sizeof(double) * pop->size);
	csum = 0;
	for (i = 0; i < pop->size; i++) {
		csum += ind1[i]->fitness - fitness_cut;
		pop_csum[i] = csum;
	}

	if (pop->selection_type & FGEN_ELITIST_ELEMENT)
		step = csum / (pop->size - nu_elites);
	else
		step = csum / pop->size;
	/* Pick a spot in the cumulative sum. */
	draw = fgen_random_d(pop->rng, csum);

	/* Traverse the roulette wheel with the calculated step. */
	/* Initially the starting spot will be looked for. */
	i = 0;
	for (j = nu_elites; j < pop->size; j++) {
		for (;;) {
			if (draw <= pop_csum[i])
				break;
			i++;
			if (i >= pop->size) {
				/* gen_error("StochasticUniversalSamplingSelection: Round problem?"); */
				/* This should almost never happen, but when it does, handle it gracefully. */
				i = pop->size - 1;
                                break;
                        }
		}

		ind2[random_order[j]] = ind1[i];
		ind1[i]->refcount++;

		draw += step;
		if (draw >= csum) {
			/* Wrap around. */
			draw -= csum;
			i = 0;
		}
	}

	free(pop_csum);
	free(random_order);
}

/*
 * Tournament selection.
 */

static void TournamentSelection(FgenPopulation *pop, FgenIndividual **ind1, FgenIndividual **ind2) {
	int *random_order;
        int i;
	/* Calculate random order for target population. */
	random_order = (int *)malloc(sizeof(int) * pop->size);
	CalculateRandomOrder(pop->rng, random_order, pop->size);

	/* In elitist selection, select the best individuals beforehand. They will also be allowed */
	/* to be selected extra times based on fitness. */
	int nu_elites = DoElitism(pop, random_order, ind1, ind2);

        if (pop->population_size_shift < 0)
		goto not_power_of_two_population_size;
        RandomIntPowerOfTwoPrepareForRepeat(pop->rng, pop->size);
	for (i = nu_elites; i < pop->size; i++) {
		FgenIndividual *winner;
		int j;
		double best_fitness = NEGATIVE_INFINITY_DOUBLE;
		for (j = 0; j < pop->tournament_size; j++) {
			int k = RandomIntPowerOfTwoRepeat(pop->rng);
			if (ind1[k]->fitness > best_fitness) {
				winner = ind1[k];
				best_fitness = ind1[k]->fitness;
			}
		}
		ind2[random_order[i]] = winner;
		winner->refcount++;
	}
	free(random_order);
	return;

not_power_of_two_population_size :
        RandomIntGeneralPrepareForRepeat(pop->rng, pop->size);
	for (i = nu_elites; i < pop->size; i++) {
		FgenIndividual *winner;
		int j;
		double best_fitness = NEGATIVE_INFINITY_DOUBLE;
		for (j = 0; j < pop->tournament_size; j++) {
			int k = RandomIntGeneralRepeat(pop->rng);
			if (ind1[k]->fitness > best_fitness) {
				winner = ind1[k];
				best_fitness = ind1[k]->fitness;
			}
		}
		ind2[random_order[i]] = winner;
		winner->refcount++;
	}
	free(random_order);
}

/*
 * Rank-based selection
 */

static void RankSelection(FgenPopulation *pop, FgenIndividual **ind1, FgenIndividual **ind2) {
	int i, j;
        double csum, draw, step;
	double *pop_csum;
	int *random_order;
	/*
	 * The individuals corresponding to equally distant markers
	 * (pop->size in number) along the roulette wheel are picked.
	 * This ensures that the each individual gets at least the integer
	 * part of its expected number.
	 */

	/* Calculate random order for target population. */
	random_order = (int *)malloc(sizeof(int) * pop->size);
	CalculateRandomOrder(pop->rng, random_order, pop->size);

	/* Sort the population on fitness. */
       	qsort(ind1, pop->size, sizeof(ind1[0]), CompareFitness);

	/* In elitist selection, select the best individuals beforehand. They will also be allowed */
	/* to be selected extra times based on fitness. */
	int nu_elites = DoElitism(pop, random_order, ind1, ind2);

	/* Calculate array of cumulative sums. */
	pop_csum = (double *)malloc(sizeof(double) * pop->size);
	csum = 0;
	for (i = 0; i < pop->size; i++) {
		// The fitness used for selection is proportional to the rank.
		csum += pop->size - i;
		pop_csum[i] = csum;
	}

	if (pop->selection_type & FGEN_ELITIST_ELEMENT)
		step = csum / (pop->size - nu_elites);
	else
		step = csum / pop->size;
	/* Pick a spot in the cumulative sum. */
	draw = fgen_random_d(pop->rng, csum);

	/* Traverse the roulette wheel with the calculated step. */
	/* Initially the starting spot will be looked for. */
	i = 0;
	for (j = nu_elites; j < pop->size; j++) {
		for (;;) {
			if (draw <= pop_csum[i])
				break;
			i++;
			if (i >= pop->size) {
				/* gen_error("RankSelection: Round problem?"); */
				/* This should almost never happen, but when it does, handle it gracefully. */
				i = pop->size - 1;
                                break;
			}
		}

		ind2[random_order[j]] = ind1[i];
		ind1[i]->refcount++;

		draw += step;
		if (draw >= csum) {
			/* Wrap around. */
			draw -= csum;
			i = 0;
		}
	}

	free(pop_csum);
	free(random_order);
}

/*
 * Extinction
 */

static void DoExtinction(FgenPopulation *pop, FgenIndividual **ind1, FgenIndividual **ind2) {
	int i;
	/* Create a population with a single copy of the best individual with the rest random new individuals. */
	FgenIndividual *best_ind = fgen_best_individual_of_population(pop);
	ind2[0] = NewIndividual(pop);
	CopyIndividualBitstring(pop, best_ind, ind2[0]);
	ind2[0]->fitness = best_ind->fitness;
	ind2[0]->fitness_is_valid = 1;
	ind2[0]->is_elite = 1;
	FreePopulation(pop, ind1);
	for (i = 1; i < pop->size; i++) {
		ind2[i] = NewIndividual(pop);
		pop->fgen_seed_func(pop, ind2[i]->bitstring);
	}
	free(ind1);
	pop->ind = ind2;
}


/*
 * Perform the selection step of the genetic algorithm.
 */

void DoSelection(FgenPopulation *pop) {
 	double sum, min_fitness, fitness_cut;
	FgenIndividual **ind1, **ind2;

	CalculatePopulationFitness(pop, &sum, &min_fitness);

	ind1 = pop->ind;
	ind2 = AllocatePopulation(pop);

	if (pop->selection_type & FGEN_EXTINCTION_ELEMENT) {
		if (pop->generation > 0 && pop->generation % 200 == 0) {
			DoExtinction(pop, ind1, ind2);
			return;
		}
        }

	/* If there is an individual with infinite fitness, copy it to the whole population. */
	if (sum == POSITIVE_INFINITY_DOUBLE) {
		int i, j;
		for (i = 0; i < pop->size; i++)
			if (ind1[i]->fitness == POSITIVE_INFINITY_DOUBLE)
				break;
		if (i == pop->size)
			gen_error("fgen: DoSelection: sum of fitness is infinite, but infinite fitness not found.");
		for (j = 0; j < pop->size; j++) {
			ind2[j] = ind1[i];
			ind1[i]->refcount++;
		}
		goto end;
	}
	if (isnan(sum) || sum == NEGATIVE_INFINITY_DOUBLE)
		gen_error("fgen: DoSelection: fitness with NaN or negative infinity detected.");
	/* If the fitness of every individual is identical, copy the population unchanged. */
	int i;
        for (i = 0; i < pop->size - 1; i++)
		if (ind1[i]->fitness != ind1[i + 1]->fitness)
			break;
        if (i == pop->size - 1) {
		int j;
		for (j = 0; j < pop->size; j++) {
			ind2[j] = ind1[j];
			ind1[j]->refcount++;
                }
		goto end;
        }

	if (pop->selection_fitness_type == FGEN_FITNESS_PROPORTIONAL)
		fitness_cut = 0;
	else
	if (pop->selection_fitness_type == FGEN_SUBTRACT_MIN_FITNESS)
		fitness_cut = min_fitness;
	else
	if (pop->selection_fitness_type == FGEN_SUBTRACT_MIN_FITNESS_DIV_2)
		fitness_cut = min_fitness / 2;
	else
		gen_error("DoSelection: Unsupported selection type.");

	if ((pop->selection_type & FGEN_STOCHASTIC_TYPE_MASK) == FGEN_STOCHASTIC)
		StochasticSelection(pop, ind1, ind2, sum, fitness_cut);
	else
	if ((pop->selection_type & FGEN_STOCHASTIC_TYPE_MASK) == FGEN_SUS)
		StochasticUniversalSamplingSelection(pop, ind1, ind2, sum, fitness_cut);
	else
	if ((pop->selection_type & FGEN_STOCHASTIC_TYPE_MASK) == FGEN_TOURNAMENT)
		TournamentSelection(pop, ind1, ind2);
	else
	if ((pop->selection_type & FGEN_STOCHASTIC_TYPE_MASK) == FGEN_RANK)
		RankSelection(pop, ind1, ind2);
	else
		gen_error("DoSelection: Unsupported selection stochastic type.");

end : 
	FreePopulation(pop, ind1);
 	free(ind1);
	pop->ind = ind2;
}

