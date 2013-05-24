/*
    population.h -- prototypes of functions defined in population.c.

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

void CreateInitialPopulation(FgenPopulation *pop);

/* Fill in the fitness of the individual (calculate if necessary). */

void CalculateIndividualFitness(FgenPopulation *, FgenIndividual *);

/* Allocate a population (array of pointers to individuals). */

FgenIndividual **AllocatePopulation(FgenPopulation *pop);

/* Fill in the fitness of all individuals in the population, and */
/* calculate the sum and minimum fitness. */

void CalculatePopulationFitness(FgenPopulation *, double *sum_out, double *min_fitness_out);

/* Fill in the fitness of all individuals using threads. */

void CalculatePopulationFitnessThreaded(FgenPopulation *);
void CalculatePopulationFitnessHeavilyThreaded(FgenPopulation *);

/* Copy a bitstring from one individual to another (not shared). */

void CopyIndividualBitstring(FgenPopulation *pop, FgenIndividual *src, FgenIndividual *dest);

/* Reduce the reference count of the individual, and deallocate if it zero. */

void FreeIndividual(FgenIndividual *ind);

/* Allocate a new individual and initialize the fields. */

FgenIndividual *NewIndividual(FgenPopulation *pop);

/* Free all invididuals in the population. */

void FreePopulation(FgenPopulation *pop, FgenIndividual **ind);

/* Return the best individual of the population. */

FgenIndividual *BestIndividualOfPopulation(FgenPopulation *pop);

