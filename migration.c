/*
    migration.c -- genetic algorithm migration implementation.

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
#include "random.h"
#include "parameters.h"
#include "error.h"
#include "population.h"
#include "win32_compat.h"

/*
 * Each individual has a chance to migrate to the previous island.
 * The least fit individual on the previous island is replaced (ring topology).
 */

static int IndexOfWorstIndividualOfPopulation(FgenPopulation *pop) {
	int i, worst;
        double worst_fitness;
        CalculatePopulationFitness(pop, NULL, NULL); /* Make sure every fitness is updated. */
	worst_fitness = POSITIVE_INFINITY_DOUBLE;
	for (i = 0; i < pop->size; i++)
		if (pop->ind[i]->fitness < worst_fitness) {
			worst_fitness = pop->ind[i]->fitness;
			worst = i;
		}
	return worst;
}

static void Migrate(FgenIndividual *ind, FgenPopulation *dest_pop) {
	int worst_index;
	worst_index = IndexOfWorstIndividualOfPopulation(dest_pop);
	FreeIndividual(dest_pop->ind[worst_index]);
	dest_pop->ind[worst_index] = NewIndividual(dest_pop);
	CopyIndividualBitstring(dest_pop, ind, dest_pop->ind[worst_index]);
	// It is possible to copy the fitness value here, but that wouldn't work
	// in the academic case of the fitness function being different between islands
	// (which is unlikely).
}

void DoMigration(int nu_pops, FgenPopulation **pops) {
	int i;
	for (i = 0; i < nu_pops; i++) {
		int previous_i;
		int j;
		previous_i = i - 1;
		if (previous_i < 0)
			previous_i = nu_pops - 1;
		for (j = 0; j < pops[i]->size; j++)  {
			int r = RandomBits(pops[i]->rng, 16);
			if (r < pops[i]->migration_probability)
				Migrate(pops[i]->ind[j], pops[previous_i]);
		}
	}
}

