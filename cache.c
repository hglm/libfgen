/*
    cache.c -- fitness cache implementation.

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
#include <string.h>
#include <limits.h>
#include <math.h>
#include "fgen.h"
#include "parameters.h"
#include "bitstring.h"
#include "error.h"
#include "population.h"
#include "cache.h"

#define CACHE_LOOK_AHEAD 4

/*
 * Test result, TSP problem:
 *
 * Cache size = 2MB entries:
 *
 * Look-ahead 8: 0.280483 (hit-rate)
 * Look-ahead 4: 0.279863
 * Look-ahead 2: 0.278402
 * Look-ahead 1: 0.275419
 *
 * Cache size = 16K entries:
 *
 * Look-ahead 32: 0.213105
 * Look-ahead 8: 0.213153
 * Look-ahead 4: 0.213110
 * Look-ahead 2: 0.212984
 * Look-ahead 1: 0.210340
 *
 * Cache size = 1K entries:
 *
 * Look-ahead 32: 0.182706
 * Look-ahead 4: 0.179864
 * Look-ahead 1: 0.155287
 */

/**
 * Enable fitness cache on population. Size is the number of entries. Memory used is roughly the number of entries
 * times the size of the bitstring of an individual. The cache imposes overhead, but for combinatorial problems
 * with an expensive fitness function it should be beneficial. For real-valued optimization problems, it is not
 * effective. When running an archipelago of GAs, it is safe to use a seperate cache for each population.
 */

void fgen_enable_cache(FgenPopulation *pop, int size) {
	int i;
	pop->cache = (FgenCache *)malloc(sizeof(FgenCache));
	pop->cache->size = size;
	pop->cache->nu_accesses = 0;
	pop->cache->nu_hits = 0;
	pop->cache->entry = (FgenCacheEntry **)malloc(sizeof(FgenCacheEntry *) * size);
	for (i = 0; i < size; i++) {
		pop->cache->entry[i] = (FgenCacheEntry *)malloc(sizeof(FgenCacheEntry));
		pop->cache->entry[i]->bitstring = (unsigned char *)malloc(INDIVIDUAL_SIZE_IN_BYTES(pop));
		pop->cache->entry[i]->date_mru = - 1;
	}
	pop->cache->refcount = 1;
}

/**
 * Enable one cache shared between all archipelagos. This only works with the unthreaded fgen_run_archipelago().
 */

void fgen_enable_cache_on_archipelago(int nu_pops, FgenPopulation **pops, int size) {
	int i;
	fgen_enable_cache(pops[0], size);
	for (i = 1; i < nu_pops; i++) {
		pops[i]->cache = pops[0]->cache;
		pops[0]->cache->refcount++;
	}
}

/**
 * Get the cache hit-rate.
 */

double fgen_get_cache_hit_rate(const FgenPopulation *pop) {
	return pop->cache->hit_rate;
}

/**
 * Invalidate the fitness cache. Should be run when the fitness functions changes in between calls to fgen_run with
 * the same cache enabled. When using a shared cache on an archipelago this function only needs to be called for one
 * island.
 */

void fgen_invalidate_cache(FgenPopulation *pop) {
	for (int i = 0; i < pop->cache->size; i++)
		pop->cache->entry[i]->date_mru = - 1;
}

void DestroyCache(FgenPopulation *pop) {
	int i;
	pop->cache->refcount--;
	if (pop->cache->refcount > 0)
		return;
	for (i = 0; i < pop->cache->size; i++) {
		free(pop->cache->entry[i]->bitstring);
		free(pop->cache->entry[i]);
	}
	free(pop->cache->entry);
	free(pop->cache);
}

static unsigned int jenkins_one_at_a_time_hash(unsigned char *key, size_t len)
{
    unsigned int hash, i;
    for (hash = i = 0; i < len; ++i)
    {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

static int Hash(FgenPopulation *pop, FgenIndividual *ind) {
	return jenkins_one_at_a_time_hash(ind->bitstring, INDIVIDUAL_SIZE_IN_BYTES(pop)) % pop->cache->size;
}

int CacheLookup(FgenPopulation *pop, FgenIndividual *ind, int *hash_out, double *fitness_out) {
	int j;
	int i = Hash(pop, ind);
	*hash_out = i;
	for (j = 0; j < CACHE_LOOK_AHEAD; j++) {
		int index = (i + j) % pop->cache->size;
		if (pop->cache->entry[index]->date_mru == -1)
			return 0;
		if (memcmp(ind->bitstring, pop->cache->entry[index]->bitstring, INDIVIDUAL_SIZE_IN_BYTES(pop)) == 0) {
			*fitness_out = pop->cache->entry[index]->fitness;
			return 1;
		}
	}
	return 0;
}

int CacheHit(FgenPopulation *pop, FgenIndividual *ind, double *cache_fitness, int *hash_out) {
	int hash;
	double fitness;
	pop->cache->nu_accesses++;
	if (CacheLookup(pop, ind, &hash, &fitness)) {
		*cache_fitness = fitness;
		*hash_out = hash;
		pop->cache->nu_hits++;
		return 1;
	}
	*hash_out = hash;
	return 0;
}

void AddToCache(FgenPopulation *pop, FgenIndividual *ind, int hash) {
	int i;
	int cache_index;
	int oldest_index;
	int oldest_mru = INT_MAX;
	for (i = 0; i < CACHE_LOOK_AHEAD; i++) {
		cache_index = (hash + i) % pop->cache->size;
		if (pop->cache->entry[cache_index]->date_mru == -1)
			goto skip;
		if (pop->cache->entry[cache_index]->date_mru < oldest_mru) {
			oldest_mru = pop->cache->entry[cache_index]->date_mru;
			oldest_index = cache_index;
		}
	}
	cache_index = oldest_index;
skip:
	memcpy(pop->cache->entry[cache_index]->bitstring, ind->bitstring,
		INDIVIDUAL_SIZE_IN_BYTES(pop));
	pop->cache->entry[cache_index]->fitness = ind->fitness;
	pop->cache->entry[cache_index]->date_mru = pop->generation;
	if (pop->cache->nu_accesses == 0)
		pop->cache->hit_rate = 0;
	else
		pop->cache->hit_rate = (double)pop->cache->nu_hits / pop->cache->nu_accesses;
}


