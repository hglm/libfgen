/*
    seed.c -- seeding functions for the genetic algorithm

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
#include "fgen.h"		/* General includes. */
#include "parameters.h"
#include "population.h"		/* Module includes. */
#include "bitstring.h"
#include "random.h"

/**
 * Seeding (initialization) operator that initializes a bitstring with random values.
 */

void fgen_seed_random(FgenPopulation *pop, unsigned char *bitstring) {
	CreateRandomBitstring(pop->rng, bitstring, INDIVIDUAL_SIZE_IN_BYTES(pop));
}


/**
 * Seeding operator that initializes an individual with a random permutation corresponding to the size
 * set by fgen_set_permutation_size().
 */

void fgen_seed_permutation_random(FgenPopulation *pop, unsigned char *bitstring) {
	CalculateRandomOrder(pop->rng, (int *)bitstring, pop->permutation_size);
}

