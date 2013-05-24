/*
    parameters.h -- default parameter values.

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
 * This file contains the definitions of the basic configurable parameters
 * for the genetic algorithm.
 */

/*
 * These are the parameter defaults; they can be overridden by the parameter setting functions.
 * The mutation probability is for one bit mutated in a single individual.
 */

#define DEFAULT_MAX_GENERATION		100

#define DEFAULT_POPULATION_SIZE		128

#define DEFAULT_INDIVIDUAL_SIZE_IN_BITS	64

#define DEFAULT_CROSSOVER_PROBABILITY	0.600

#define DEFAULT_MUTATION_PROBABILITY	0.005
#define DEFAULT_MACRO_MUTATION_PROBABILITY 0.050

#define DEFAULT_SELECTION_FITNESS_TYPE	FGEN_FITNESS_PROPORTIONAL

#define DEFAULT_SELECTION_TYPE		FGEN_SUS
#define DEFAULT_TOURNAMENT_SIZE		3

#define INDIVIDUAL_SIZE_IN_BYTES(pop) ((pop->individual_size_in_bits + 7) / 8)

