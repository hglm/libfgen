/*
    cache.h -- prototypes of functions defined in cache.c.

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


int CacheLookup(FgenPopulation *pop, FgenIndividual *ind, int *hash_out, double *fitness_out);

int CacheHit(FgenPopulation *pop, FgenIndividual *ind, double *cache_fitness, int *hash_out);

void AddToCache(FgenPopulation *pop, FgenIndividual *ind, int hash);

void DestroyCache(FgenPopulation *pop);

void InvalidateCache(FgenPopulation *pop);

